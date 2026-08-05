# A collections of functions called in other R scripts 


#######################
#####  Packages   #####
#######################


# Install missing packages and load required packages (if required)
UsePackages <- function( pkgs, load=TRUE, locn="http://cran.rstudio.com/" ) {
  # Identify missing (i.e., not yet installed) packages
  newPkgs <- pkgs[!(pkgs %in% installed.packages( )[, "Package"])]
  # Install missing packages if required
  if( length(newPkgs) )  install.packages( newPkgs, repos=locn )
  # Load packages
  if( load ){
    # Loop over all packages
    for( i in 1:length(pkgs) ) {
      # Load required packages using 'library'
      eval( parse(text=paste("library(", pkgs[i], ")", sep="")) )
    } # End i loop over package names
  }  # End if load
}  # End UsePackages function



#######################
#####  Data prep  #####
#######################
Mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Get land data for the study area
GetLandData <- function( doThin=TRUE, dir ) {
  # shapefile name
  coast <- list.files(path = dir, pattern = ".shp$")
  coast <- sub(".shp","",coast)
  # Load shapefile
  bcPoly <- readOGR( dsn=dir, layer=coast, verbose=FALSE )
  # Convert to spatial polygons (i.e., drop the data)
  bcPoly <- as( bcPoly, "SpatialPolygons" )
  # Set the coordinate reference system
  if( wkt( bcPoly ) != geoCRS ){
    # project to match geoCRS
    bcPoly <- spTransform( bcPoly, CRS(geoCRS) )
  }
  # Thin the polygons if requested
  if( doThin ){
    bcPoly <- maptools::thinnedSpatialPoly( SP=bcPoly, tolerance=100, minarea=100 )
  }
  # Convert the land layer to points
  bcDF <- fortify( model=bcPoly )
  # Set new column names, and drop others
  bcDF <- transmute( .data=bcDF, Longitude=long, Latitude=lat, group=group )
  # Return the data
  return( list(bcPoly=bcPoly, bcDF=bcDF) )
}  # End GetLandData function


# Calculate collinearity of environmental variables
CalcVIFs <- function( dat, VIFThresh = 10 ) {
  
  # Remove factors
  dat <- dat[ !names(dat) %in% facVars ] 
  
  # Calculate overall collinearity (VIF == variance inflation factor)
  # vif() selects variables for a linear model
  # returns a subset of variables for building a linear model
  VIFs <- usdm::vif( x=dat )
  
  # If there are collinear predictors
  if( any( VIFs$VIF > VIFThresh ) ) {
    # Message
    message( "\nCollinear predictors found: '", 
             paste(VIFs$Variables[VIFs$VIF > VIFThresh], collapse = "' '"), "'\n" )
    # Stepwise procedure to reduce number of predictors
    stepVIFs <- usdm::vifstep( x=dat, th=VIFThresh )
    # Predictors to remove
    exPred <<- stepVIFs@excluded
  } else {  # End if there are collinear predictors, otherwise
    exPred <<- NULL
    # Message
    message( "\nNo collinear predictors found\n" )
  }  # End if there are no collinear predictors
  
  # Return data
  return( VIFs )
}  # End CalcVIFs function



#######################
#####  R SQL Link #####
#######################

# Functions for establishing connections to DFO databases and editing sql code

sf_db_connection <- function(server = "WDC-SQL2016-P\\SIOSP01") {
  require(DBI)
  DBI::dbConnect(odbc::odbc(),
                 driver = "SQL Server",
                 server = server)
}

mdb_connection <- function(db_file_path)  {
  require(DBI)
  # Make sure that the file exists before attempting to connect
  if (!file.exists(db_file_path)) {
    stop("DB file does not exist at ", db_file_path)
  }
  # Connect
  DBI::dbConnect(odbc::odbc(),
                 .connection_string = 
                   paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};", 
                          paste0("DBQ=", db_file_path)))
}

##########################
#####  SDM Functions #####
##########################

tjur <- function(y, pred) {
  categories <- sort(unique(y))
  m1 <- mean(pred[which(y == categories[1])], na.rm = TRUE)
  m2 <- mean(pred[which(y == categories[2])], na.rm = TRUE)
  abs(m2 - m1)}

ll_binomial <- function(withheld_y, withheld_mu) {
  stats::dbinom(x = withheld_y, size = 1, prob = withheld_mu, log = TRUE) }

# cloglog inverse link from stats::family()
pclog <- function(x){
  pmax(pmin(-expm1(-exp(x)), 1 - .Machine$double.eps), .Machine$double.eps)
}

tjur_fun <- function(y, p) {
  mean(p[y == 1], na.rm = TRUE) - mean(p[y == 0], na.rm = TRUE)
}
rmse_fun <- function(y, p) {
  sqrt(mean((y - p)^2, na.rm = TRUE))
}
auc_fun <- function(y, p) {
  as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
}
tss_fun <- function(y, p, threshold) {
  pred_bin <- ifelse(p >= threshold, 1, 0)
  tp <- sum(pred_bin == 1 & y == 1)
  tn <- sum(pred_bin == 0 & y == 0)
  fp <- sum(pred_bin == 1 & y == 0)
  fn <- sum(pred_bin == 0 & y == 1)
  sens <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  spec <- ifelse(tn + fp == 0, 0, tn / (tn + fp))
  sens + spec - 1
}
get_tss_threshold <- function(obs, pred) {
  
  thresholds <- seq(0.01,0.99,0.01)
  
  tss <- sapply(thresholds,function(thr){
    
    pred_bin <- ifelse(pred >= thr,1,0)
    
    sens <- sensitivity_fun(obs,pred_bin)
    
    spec <- specificity_fun(obs,pred_bin)
    
    sens + spec - 1
    
  })
  
  thresholds[which.max(tss)]
  
}

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

logloss_fun <- function(obs, pred) {
  
  eps <- 1e-15
  
  pred <- pmin(pmax(pred, eps), 1 - eps)
  
  -mean(
    obs * log(pred) +
      (1 - obs) * log(1 - pred)
  )
}

brier_fun <- function(obs, pred) {
  
  mean((obs - pred)^2)
  
}

sensitivity_fun <- function(obs, pred_bin) {
  
  tp <- sum(obs == 1 & pred_bin == 1)
  
  fn <- sum(obs == 1 & pred_bin == 0)
  
  tp / (tp + fn)
  
}

specificity_fun <- function(obs, pred_bin) {
  
  tn <- sum(obs == 0 & pred_bin == 0)
  
  fp <- sum(obs == 0 & pred_bin == 1)
  
  tn / (tn + fp)
  
}

tss_metric <- function(obs, pred, threshold) {
  pred_bin <- as.numeric(pred >= threshold)
  cm <- table(factor(obs, levels = c(0,1)),
              factor(pred_bin, levels = c(0,1)))
  
  sens <- ifelse(sum(cm["1",]) > 0, cm["1","1"]/sum(cm["1",]), NA)
  spec <- ifelse(sum(cm["0",]) > 0, cm["0","0"]/sum(cm["0",]), NA)
  sens + spec - 1
}

confusion_stats <- function(obs,
                            pred,
                            threshold){
  
  pred_class <- ifelse(pred >= threshold,1,0)
  
  TP <- sum(pred_class == 1 & obs == 1)
  TN <- sum(pred_class == 0 & obs == 0)
  FP <- sum(pred_class == 1 & obs == 0)
  FN <- sum(pred_class == 0 & obs == 1)
  
  sensitivity <- ifelse(
    TP + FN == 0,
    NA,
    TP/(TP+FN)
  )
  specificity <- ifelse(
    TN + FP == 0,
    NA,
    TN/(TN+FP)
  )
  tss <- sensitivity + specificity - 1
  list(
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss
  )
}


# safe Tjur R2
tjur_r2 <- function(obs, pred) {
  if (length(unique(obs)) < 2) return(NA_real_)
  mean(pred[obs == 1], na.rm = TRUE) - mean(pred[obs == 0], na.rm = TRUE)
}
# safe AUC using ecospat
auc_safe <- function(obs, pred) {
  if (length(unique(obs)) < 2) return(NA_real_)
  ecospat::ecospat.auctest(pred[obs == 1], pred[obs == 0], nrep = 0)$auc
}
# threshold-based TSS
tss_from_threshold <- function(obs, pred, threshold) {
  pred_bin <- ifelse(pred >= threshold, 1, 0)
  
  tp <- sum(obs == 1 & pred_bin == 1, na.rm = TRUE)
  tn <- sum(obs == 0 & pred_bin == 0, na.rm = TRUE)
  fp <- sum(obs == 0 & pred_bin == 1, na.rm = TRUE)
  fn <- sum(obs == 1 & pred_bin == 0, na.rm = TRUE)
  
  sens <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
  spec <- ifelse((tn + fp) == 0, NA, tn / (tn + fp))
  
  sens + spec - 1
}

get_optimal_threshold <- function(obs, pred) {
  df <- data.frame(
    id = seq_along(obs),
    obs = obs,
    pred = pred
  )
  
  out <- tryCatch({
    PresenceAbsence::optimal.thresholds(
      df,
      opt.methods = "MaxSens+Spec"
    )$predicted[1]
  }, error = function(e) NA)
  
  if(is.null(out) || is.na(out)) {
    out <- mean(obs)  # prevalence fallback
  }
  
  return(out)
}

# forward feature selection glm
glm_ffs <- function(data, NumFolds){
  Fold_Outputs <- list()
  for(u in 1:numFolds){
    sp_data <- data %>% dplyr::filter(fold <=numFolds) %>% dplyr::filter(fold != u)
    Test_Data <- data %>% dplyr::filter(fold <=numFolds) %>% dplyr::filter(fold == u)
    #create a vector for fold membership
    foldmem <- seq(1:1:length(sp_data$fold))
    folds <- unique(sp_data$fold)
    new_folds <- seq(1:1:(numFolds-1))
    for (i in 1:length(sp_data$fold)){
      for (j in 1:(numFolds-1)){
        if (sp_data$fold[[i]] == folds[j]){
          foldmem[i] <- new_folds[j]}}}
    #Setting up an index list for the folds in the caret model training
    index_list <- list()
    for (i in 1:(numFolds-1)){
      index_list[[i]] <- which(foldmem == i)}
    #Setting up parameters for how my model is going to be fit
    fitControl <- caret::trainControl(method = "cv",
                                      number = (numFolds-1),
                                      index = index_list)
    set.seed(2024)
    #performing model selection by glmStepAIC 
    caret_model <- CAST::ffs(response = sp_data$presence, 
                             predictors = sp_data[,6:69], 
                             method = "glm", 
                             family = "binomial",
                             trControl = fitControl)
    #Create the final glm model using the above determined formula
    selectedVars <- caret_model$selectedvars
    final_formula <- paste0("presence~",selectedVars[1])
    for (i in 2: length(selectedVars)){
      final_formula <- paste0(final_formula,"+",selectedVars[i])}
    #fit the model
    caretmodel <- glm(as.formula(final_formula), data = sp_data, family = binomial(link = "logit"))
    #Calculate final model AUC on the testing fold
    pred.caretModel <- predict(caretmodel, newdata = Test_Data, type = "response")
    roc.caretModel <- pROC::roc(Test_Data$presence, pred.caretModel)
    auc.caretModel <- pROC::auc(roc.caretModel)
    Output <- list(caretmodel, auc.caretModel, caretmodel$formula)
    Fold_Outputs[[u]] <- Output
  }
  return(Fold_Outputs)  
}



#### nested spatial cross validation to select best smoothing terms 
# eelgrass
nested_sdmTMB_cv <- function(
    data,
    mesh,
    spatial_setting = FALSE,
    ocean_model = c("bccm", "nep"),
    outer_fold_col = "fold",
    inner_fold_col = "inner_fold",
    depth_k_values = list(
      "linear" = NULL,
      "2" = 2,
      "3" = 3,
      "4" = 4
    ),
    airtemp_k_values = list(
      "linear" = NULL,
      "2" = 2,
      "3" = 3,
      "4" = 4
    ),
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit")
) {
  
  library(dplyr)
  library(sdmTMB)
  library(pROC)
  
  ocean_model <- match.arg(ocean_model)
  
  #-------------------------------------------
  # choose oceanographic predictors
  #-------------------------------------------
  ocean_terms <- switch(
    ocean_model,
    bccm = "saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd",
    nep  = "saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd"
  )
  
  #-------------------------------------------
  # formula builder
  #-------------------------------------------
  build_formula <- function(depth_k, airtemp_k, ocean_terms) {
    depth_term <-
      if (is.null(depth_k))
        "depth_stnd"
    else
      paste0("s(depth_stnd, k=", depth_k, ")")
    
    airtemp_term <-
      if (is.null(airtemp_k))
        "airtempmin_stnd"
    else
      paste0("s(airtempmin_stnd, k=", airtemp_k, ")")
    
    as.formula(
      paste(
        "presence ~",
        depth_term,
        "+ substrate + slope_stnd + rei_stnd +",
        airtemp_term,
        "+ rsdsmin_stnd + prmin_stnd +",
        ocean_terms,
        "+ Survey"
      )
    )
  }
  
  
  #-------------------------------------------
  # helper metrics
  #-------------------------------------------
  calc_logloss <- function(actual, prob) {
    eps <- 1e-15
    prob_clip <- pmin(pmax(prob, eps), 1 - eps)
    -mean(actual * log(prob_clip) + (1 - actual) * log(1 - prob_clip), na.rm = TRUE)
  }
  
  calc_tjur <- function(actual, prob) {
    mean(prob[actual == 1], na.rm = TRUE) - mean(prob[actual == 0], na.rm = TRUE)
  }
  
  calc_brier <- function(actual, prob) {
    mean((prob - actual)^2, na.rm = TRUE)
  }
  
  calc_best_tss <- function(actual, prob) {
    thresholds <- sort(unique(prob))
    
    if (length(thresholds) == 0) {
      return(list(
        tss = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        threshold = NA_real_
      ))
    }
    
    tss_scores <- numeric(length(thresholds))
    sens_vals <- numeric(length(thresholds))
    spec_vals <- numeric(length(thresholds))
    
    for (j in seq_along(thresholds)) {
      th <- thresholds[j]
      pred_class <- as.integer(prob >= th)
      
      sens <- if (sum(actual == 1) > 0) {
        sum(pred_class == 1 & actual == 1) / sum(actual == 1)
      } else {
        NA_real_
      }
      
      spec <- if (sum(actual == 0) > 0) {
        sum(pred_class == 0 & actual == 0) / sum(actual == 0)
      } else {
        NA_real_
      }
      
      tss_scores[j] <- sens + spec - 1
      sens_vals[j] <- sens
      spec_vals[j] <- spec
    }
    
    best_j <- which.max(tss_scores)
    
    list(
      tss = tss_scores[best_j],
      sensitivity = sens_vals[best_j],
      specificity = spec_vals[best_j],
      threshold = thresholds[best_j]
    )
  }
  
  safe_auc <- function(actual, prob) {
    tryCatch(as.numeric(pROC::auc(actual, prob)), error = function(e) NA_real_)
  }
  
  #-------------------------------------------
  # checks
  #-------------------------------------------
  if (!outer_fold_col %in% names(data)) {
    stop("Outer fold column not found: ", outer_fold_col)
  }
  if (!inner_fold_col %in% names(data)) {
    stop("Inner fold column not found: ", inner_fold_col)
  }
  
  results <- list()
  
  for (outer_fold in sort(unique(data[[outer_fold_col]]))) {
    cat("Outer Fold:", outer_fold, "\n")
    
    test_data <- data[data[[outer_fold_col]] == outer_fold, ]
    train_data <- data[data[[outer_fold_col]] != outer_fold, ]
    
    param_grid <- expand.grid(
      depth_k = names(depth_k_values),
      airtemp_k = names(airtemp_k_values),
      stringsAsFactors = FALSE
    )
    
    inner_res <- list()
    
    for (i in seq_len(nrow(param_grid))) {
      params <- param_grid[i, , drop = FALSE]
      aucs <- c()
      loglosses <- c()
      
      for (inner_fold in sort(unique(train_data[[inner_fold_col]]))) {
        valid <- train_data[train_data[[inner_fold_col]] == inner_fold, ]
        train <- train_data[train_data[[inner_fold_col]] != inner_fold, ]
        
        depth_k <- depth_k_values[[params$depth_k]]
        airtemp_k <- airtemp_k_values[[params$airtemp_k]]
        
        formula <- build_formula(depth_k, airtemp_k, ocean_terms)
        
        fit <- try(
          sdmTMB(
            formula = formula,
            mesh = mesh,
            family = family,
            priors = priors,
            spatial = spatial_setting,
            data = train,
            silent = TRUE
          ),
          silent = TRUE
        )
        
        if (!inherits(fit, "try-error")) {
          pred <- predict(fit, newdata = valid)
          prob <- plogis(pred$est)
          
          auc_val <- try(pROC::auc(valid$presence, prob), silent = TRUE)
          if (!inherits(auc_val, "try-error")) {
            aucs <- c(aucs, as.numeric(auc_val))
          }
          
          ll <- calc_logloss(valid$presence, prob)
          loglosses <- c(loglosses, ll)
        }
      }
      
      inner_res[[i]] <- data.frame(
        depth_k = params$depth_k,
        airtemp_k = params$airtemp_k,
        mean_auc = mean(aucs, na.rm = TRUE),
        mean_logloss = mean(loglosses, na.rm = TRUE)
      )
    }
    
    inner_df <- dplyr::bind_rows(inner_res)
    
    if (all(is.na(inner_df$mean_logloss))) {
      best_params <- data.frame(
        depth_k = NA_character_,
        airtemp_k = NA_character_,
        mean_auc = NA_real_,
        mean_logloss = NA_real_
      )
    } else {
      best_i <- which.min(inner_df$mean_logloss)
      best_params <- inner_df[best_i, ]
    }
    
    cat("Best inner params for outer fold", outer_fold, ":\n")
    print(best_params)
    
    best_depth_k <- depth_k_values[[best_params$depth_k]]
    best_airtemp_k <- airtemp_k_values[[best_params$airtemp_k]]
    
    formula <- build_formula(
      best_depth_k,
      best_airtemp_k,
      ocean_terms
    )
    
    fit <- try(
      sdmTMB(
        formula = formula,
        mesh = mesh,
        family = family,
        priors = priors,
        spatial = spatial_setting,
        data = train_data,
        silent = TRUE
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      pred <- predict(fit, newdata = test_data)
      prob <- plogis(pred$est)
      actual <- test_data$presence
      
      tss_out <- calc_best_tss(actual, prob)
      auc <- safe_auc(actual, prob)
      tjur_r2 <- calc_tjur(actual, prob)
      brier_score <- calc_brier(actual, prob)
      logloss <- calc_logloss(actual, prob)
    } else {
      auc <- tjur_r2 <- brier_score <- logloss <- NA_real_
      tss_out <- list(
        tss = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        threshold = NA_real_
      )
    }
    
    results[[as.character(outer_fold)]] <- data.frame(
      outer_fold = outer_fold,
      ocean_model = ocean_model,
      spatial = spatial_setting,
      depth_k = best_params$depth_k,
      airtemp_k = best_params$airtemp_k,
      test_auc = auc,
      tjur_r2 = tjur_r2,
      brier_score = brier_score,
      sensitivity = tss_out$sensitivity,
      specificity = tss_out$specificity,
      tss = tss_out$tss,
      logloss = logloss,
      tss_threshold = tss_out$threshold
    )
  }
  
  final_results <- dplyr::bind_rows(results)
  
  summary_results <- final_results %>%
    summarise(
      mean_auc = mean(test_auc, na.rm = TRUE),
      mean_tjur = mean(tjur_r2, na.rm = TRUE),
      mean_brier = mean(brier_score, na.rm = TRUE),
      mean_tss = mean(tss, na.rm = TRUE),
      mean_logloss = mean(logloss, na.rm = TRUE)
    )
  
  list(
    fold_results = final_results,
    summary_results = summary_results
  )
}


nested_sdmTMB_cv_surf <- function(
    data,
    mesh,
    spatial_setting = FALSE,
    ocean_model = c("bccm", "nep"),
    outer_fold_col = "fold",
    inner_fold_col = "inner_fold",
    depth_k_values = list(
      "linear" = NULL,
      "2" = 2,
      "3" = 3,
      "4" = 4
    ),
    rei_sqrt_k_values = list(
      "linear" = NULL,
      "2" = 2,
      "3" = 3,
      "4" = 4
    ),
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit")
) {
  
  library(dplyr)
  library(sdmTMB)
  library(pROC)
  
  ocean_model <- match.arg(ocean_model)
  
  #-------------------------------------------
  # choose oceanographic predictors
  #-------------------------------------------
  ocean_terms <- switch(
    ocean_model,
    bccm = "saltmean_bccm_stnd + tempmean_bccm_stnd + surftempcv_bccm_stnd",
    nep  = "saltmean_nep_stnd + tempmean_nep_stnd + surftempcv_nep_stnd"
  )
  
  #-------------------------------------------
  # formula builder
  #-------------------------------------------
  build_formula <- function(depth_k, rei_k, ocean_terms) {
    depth_term <-
      if (is.null(depth_k))
        "depth_stnd"
    else
      paste0("s(depth_stnd, k=", depth_k, ")")
    
    rei_term <-
      if (is.null(rei_k))
        "rei_sqrt_stnd"
    else
      paste0("s(rei_sqrt_stnd, k=", rei_k, ")")
    
    as.formula(
      paste(
        "presence ~",
        depth_term,
        "+ substrate + tidal_sqrt_stnd +",
        rei_term,
        "+ airtempcv_stnd + prmean_stnd + rsdsmin_stnd +",
        ocean_terms,
        "+ Survey"
      )
    )
  }
  
  
  #-------------------------------------------
  # helper metrics
  #-------------------------------------------
  calc_logloss <- function(actual, prob) {
    eps <- 1e-15
    prob_clip <- pmin(pmax(prob, eps), 1 - eps)
    -mean(actual * log(prob_clip) + (1 - actual) * log(1 - prob_clip), na.rm = TRUE)
  }
  
  calc_tjur <- function(actual, prob) {
    mean(prob[actual == 1], na.rm = TRUE) - mean(prob[actual == 0], na.rm = TRUE)
  }
  
  calc_brier <- function(actual, prob) {
    mean((prob - actual)^2, na.rm = TRUE)
  }
  
  calc_best_tss <- function(actual, prob) {
    thresholds <- sort(unique(prob))
    
    if (length(thresholds) == 0) {
      return(list(
        tss = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        threshold = NA_real_
      ))
    }
    
    tss_scores <- numeric(length(thresholds))
    sens_vals <- numeric(length(thresholds))
    spec_vals <- numeric(length(thresholds))
    
    for (j in seq_along(thresholds)) {
      th <- thresholds[j]
      pred_class <- as.integer(prob >= th)
      
      sens <- if (sum(actual == 1) > 0) {
        sum(pred_class == 1 & actual == 1) / sum(actual == 1)
      } else {
        NA_real_
      }
      
      spec <- if (sum(actual == 0) > 0) {
        sum(pred_class == 0 & actual == 0) / sum(actual == 0)
      } else {
        NA_real_
      }
      
      tss_scores[j] <- sens + spec - 1
      sens_vals[j] <- sens
      spec_vals[j] <- spec
    }
    
    best_j <- which.max(tss_scores)
    
    list(
      tss = tss_scores[best_j],
      sensitivity = sens_vals[best_j],
      specificity = spec_vals[best_j],
      threshold = thresholds[best_j]
    )
  }
  
  safe_auc <- function(actual, prob) {
    tryCatch(as.numeric(pROC::auc(actual, prob)), error = function(e) NA_real_)
  }
  
  #-------------------------------------------
  # checks
  #-------------------------------------------
  if (!outer_fold_col %in% names(data)) {
    stop("Outer fold column not found: ", outer_fold_col)
  }
  if (!inner_fold_col %in% names(data)) {
    stop("Inner fold column not found: ", inner_fold_col)
  }
  
  results <- list()
  
  for (outer_fold in sort(unique(data[[outer_fold_col]]))) {
    cat("Outer Fold:", outer_fold, "\n")
    
    test_data <- data[data[[outer_fold_col]] == outer_fold, ]
    train_data <- data[data[[outer_fold_col]] != outer_fold, ]
    
    param_grid <- expand.grid(
      depth_k = names(depth_k_values),
      rei_k = names(rei_sqrt_k_values),
      stringsAsFactors = FALSE
    )
    
    inner_res <- list()
    
    for (i in seq_len(nrow(param_grid))) {
      params <- param_grid[i, , drop = FALSE]
      aucs <- c()
      loglosses <- c()
      
      for (inner_fold in sort(unique(train_data[[inner_fold_col]]))) {
        valid <- train_data[train_data[[inner_fold_col]] == inner_fold, ]
        train <- train_data[train_data[[inner_fold_col]] != inner_fold, ]
        
        depth_k <- depth_k_values[[params$depth_k]]
        rei_k <- rei_sqrt_k_values[[params$rei_k]]
        
        formula <- build_formula(depth_k, rei_k, ocean_terms)
        
        fit <- try(
          sdmTMB(
            formula = formula,
            mesh = mesh,
            family = family,
            priors = priors,
            spatial = spatial_setting,
            data = train,
            silent = TRUE
          ),
          silent = TRUE
        )
        
        if (!inherits(fit, "try-error")) {
          pred <- predict(fit, newdata = valid)
          prob <- plogis(pred$est)
          
          auc_val <- try(pROC::auc(valid$presence, prob), silent = TRUE)
          if (!inherits(auc_val, "try-error")) {
            aucs <- c(aucs, as.numeric(auc_val))
          }
          
          ll <- calc_logloss(valid$presence, prob)
          loglosses <- c(loglosses, ll)
        }
      }
      
      inner_res[[i]] <- data.frame(
        depth_k = params$depth_k,
        rei_k = params$rei_k,
        mean_auc = mean(aucs, na.rm = TRUE),
        mean_logloss = mean(loglosses, na.rm = TRUE)
      )
    }
    
    inner_df <- dplyr::bind_rows(inner_res)
    
    if (all(is.na(inner_df$mean_logloss))) {
      best_params <- data.frame(
        depth_k = NA_character_,
        rei_k = NA_character_,
        mean_auc = NA_real_,
        mean_logloss = NA_real_
      )
    } else {
      best_i <- which.min(inner_df$mean_logloss)
      best_params <- inner_df[best_i, ]
    }
    
    cat("Best inner params for outer fold", outer_fold, ":\n")
    print(best_params)
    
    best_depth_k <- depth_k_values[[best_params$depth_k]]
    best_rei_k <- rei_sqrt_k_values[[best_params$rei_k]]
    
    formula <- build_formula(
      best_depth_k,
      best_rei_k,
      ocean_terms
    )
    
    fit <- try(
      sdmTMB(
        formula = formula,
        mesh = mesh,
        family = family,
        priors = priors,
        spatial = spatial_setting,
        data = train_data,
        silent = TRUE
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      pred <- predict(fit, newdata = test_data)
      prob <- plogis(pred$est)
      actual <- test_data$presence
      
      tss_out <- calc_best_tss(actual, prob)
      auc <- safe_auc(actual, prob)
      tjur_r2 <- calc_tjur(actual, prob)
      brier_score <- calc_brier(actual, prob)
      logloss <- calc_logloss(actual, prob)
    } else {
      auc <- tjur_r2 <- brier_score <- logloss <- NA_real_
      tss_out <- list(
        tss = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        threshold = NA_real_
      )
    }
    
    results[[as.character(outer_fold)]] <- data.frame(
      outer_fold = outer_fold,
      ocean_model = ocean_model,
      spatial = spatial_setting,
      depth_k = best_params$depth_k,
      rei_k = best_params$rei_k,
      test_auc = auc,
      tjur_r2 = tjur_r2,
      brier_score = brier_score,
      sensitivity = tss_out$sensitivity,
      specificity = tss_out$specificity,
      tss = tss_out$tss,
      logloss = logloss,
      tss_threshold = tss_out$threshold
    )
  }
  
  final_results <- dplyr::bind_rows(results)
  
  summary_results <- final_results %>%
    summarise(
      mean_auc = mean(test_auc, na.rm = TRUE),
      mean_tjur = mean(tjur_r2, na.rm = TRUE),
      mean_brier = mean(brier_score, na.rm = TRUE),
      mean_tss = mean(tss, na.rm = TRUE),
      mean_logloss = mean(logloss, na.rm = TRUE)
    )
  
  list(
    fold_results = final_results,
    summary_results = summary_results
  )
}

# Variable importance
# Method from SDMtune R package: 
#  'The function randomly permutes one variable at time (using training and
#  absence/background datasets) and computes the decrease in training AUC. The
#  result is normalized to percentages. Same implementation of MaxEnt java
#  software but with the additional possibility of running several permutations
#  to obtain a better estimate of the permutation importance. In case of more
#  than one permutation (default is 10) the average of the decrease in training
#  AUC is computed.'
### We estimated the relative influence of covariates using a permutation method, based on the method implemented in MaxEnt software (Phillips et al. 2006). For each covariate we: (1) randomized the covariate with respect to the observations, (2) fit a model with the randomized covariate and all other non-randomized covariates, and (3) accessed model performance with the area under the receiver operating characteristic curve (AUC) metric. We completed steps (1) through (3) 10 times and returned a mean AUC value from the 10 permutations. We then calculated the AUC, the difference between the non-randomized model AUC and the permuted mean AUC from the randomized models. Finally, for each covariate, we divided the AUC by the sum of AUC values from all covariates to obtain the relative influence. A large AUC indicates that the randomized covariate has a large influence, while a small AUC indicates that the covariate has little influence on the model fit. For spatial random fields, we adjusted the procedure as randomization of sampling location (latitude and longitude) was not appropriate. We calculated the influence of the spatial random field by dropping it from the model, then measuring AUC between the model with random fields and the model without.


varImp_sdmTMB <- function(model, dat, preds, permute = 10) {
  library(sdmTMB)
  library(ModelMetrics)
  
  # ------------------------------
  # Checks
  # ------------------------------
  if (!all(preds %in% names(dat))) {
    stop("Some predictors are missing from the dataset.")
  }
  
  if (!"presence" %in% names(dat)) {
    stop("Column 'presence' is missing from the dataset.")
  }
  
  obs <- dat$presence
  
  # ------------------------------
  # Handle spatiotemporal models ONLY if needed
  # ------------------------------
  if (!is.null(model$spatiotemporal) && model$spatiotemporal != "off") {
    time_var <- model$time
    if (is.null(time_var) || !time_var %in% names(dat)) {
      stop("Model is spatiotemporal but time variable is missing from data.")
    }
    dat[[time_var]] <- as.integer(dat[[time_var]])
  }
  
  # ------------------------------
  # Base AUC
  # ------------------------------
  base_pred <- predict(model, newdata = dat, type = "response", re_form = NA)
  base_auc <- auc(obs, base_pred$est)
  rm(base_pred); gc()
  
  # ------------------------------
  # Individual variable importance
  # ------------------------------
  aucs_ind <- vector("list", length(preds))
  names(aucs_ind) <- preds
  
  for (v in preds) {
    perm_results <- numeric(permute)
    
    for (i in seq_len(permute)) {
      dat_perm <- dat
      dat_perm[[v]] <- sample(dat_perm[[v]])
      
      p <- try(
        predict(model, newdata = dat_perm, type = "response", re_form = NA),
        silent = TRUE
      )
      
      perm_results[i] <- ifelse(
        inherits(p, "try-error"),
        NA,
        auc(obs, p$est)
      )
      
      rm(p, dat_perm); gc()
    }
    
    aucs_ind[[v]] <- perm_results
  }
  
  perm_auc_ind <- sapply(aucs_ind, mean, na.rm = TRUE)
  imp_ind <- pmax(0, base_auc - perm_auc_ind)
  
  if (sum(imp_ind, na.rm = TRUE) == 0) {
    rel_ind <- rep(0, length(imp_ind))
  } else {
    rel_ind <- round(100 * imp_ind / sum(imp_ind, na.rm = TRUE), 1)
  }
  
  ind_df <- data.frame(
    term = preds,
    mean_permuted_auc = as.numeric(perm_auc_ind),
    importance = as.numeric(imp_ind),
    relimp = as.numeric(rel_ind),
    stringsAsFactors = FALSE
  )
  
  # ------------------------------
  # Return individual importance only
  # ------------------------------
  list(
    base_auc = base_auc,
    individual = ind_df
  )
}

get_oof_predictions <- function(cv_object, cv_list, response_col) {
  
  oof_pred <- rep(NA_real_, nrow(cv_object$data))
  
  for(i in seq_along(cv_object$models)) {
    
    test_idx <- cv_list[[i]]$test
    
    pred <- predict(
      cv_object$models[[i]],
      newdata = cv_object$data[test_idx, ]
    )
    
    oof_pred[test_idx] <- plogis(pred$est)
  }
  
  data.frame(
    obs = cv_object$data[[response_col]],
    pred = oof_pred
  )
}

get_optimal_threshold <- function(actual, prob) {
  
  keep <- !is.na(actual) & !is.na(prob)
  
  actual <- actual[keep]
  prob <- prob[keep]
  
  thresholds <- sort(unique(prob))
  
  best_tss <- -Inf
  best_threshold <- NA_real_
  
  for(th in thresholds) {
    
    pred_class <- as.integer(prob >= th)
    
    sens <- sum(pred_class == 1 & actual == 1) /
      sum(actual == 1)
    
    spec <- sum(pred_class == 0 & actual == 0) /
      sum(actual == 0)
    
    tss <- sens + spec - 1
    
    if(tss > best_tss) {
      best_tss <- tss
      best_threshold <- th
    }
  }
  
  best_threshold
}


evalStats <- function(folds, m, CV, response_col, fixed_threshold) {
  traintest.df <- list()
  
  # Helper: sensitivity, specificity, TSS from threshold
  class_metrics <- function(obs, pred, threshold) {
    pred_bin <- as.numeric(pred >= threshold)
    
    cm <- table(
      factor(obs, levels = c(0, 1)),
      factor(pred_bin, levels = c(0, 1))
    )
    
    sensitivity <- ifelse(sum(cm["1", ]) > 0, cm["1", "1"] / sum(cm["1", ]), NA)
    specificity <- ifelse(sum(cm["0", ]) > 0, cm["0", "0"] / sum(cm["0", ]), NA)
    tss <- sensitivity + specificity - 1
    
    list(
      sensitivity = sensitivity,
      specificity = specificity,
      tss = tss
    )
  }
  
  # Helper: log loss
  log_loss_fun <- function(obs, pred, eps = 1e-15) {
    pred <- pmin(pmax(pred, eps), 1 - eps)
    -mean(obs * log(pred) + (1 - obs) * log(1 - pred), na.rm = TRUE)
  }
  
  for (i in folds) {
    # Train/test indices
    train_idx <- CV[[i]][["train"]]
    test_idx  <- CV[[i]][["test"]]
    
    # Observed values
    obs_train <- m$data[[response_col]][train_idx]
    obs_test  <- m$data[[response_col]][test_idx]
    
    # Predictions
    pred_train <- plogis(predict(m$models[[i]])$est[train_idx])
    pred_test  <- plogis(predict(m$models[[i]])$est[test_idx])
    
    threshold <- fixed_threshold
    
    # Train classification metrics
    train_cls <- class_metrics(obs_train, pred_train, threshold)
    
    # Train metrics
    train_auc         <- ModelMetrics::auc(obs_train, pred_train)
    train_tjur        <- tjur(obs_train, pred_train)
    train_rmse        <- sqrt(mean((obs_train - pred_train)^2, na.rm = TRUE))
    train_brier       <- mean((obs_train - pred_train)^2, na.rm = TRUE)
    train_logloss     <- log_loss_fun(obs_train, pred_train)
    train_sensitivity <- train_cls$sensitivity
    train_specificity <- train_cls$specificity
    train_tss         <- train_cls$tss
    
    # Test classification metrics
    test_cls <- class_metrics(obs_test, pred_test, threshold)
    
    # Test metrics
    test_auc         <- ModelMetrics::auc(obs_test, pred_test)
    test_tjur        <- tjur(obs_test, pred_test)
    test_rmse        <- sqrt(mean((obs_test - pred_test)^2, na.rm = TRUE))
    test_brier       <- mean((obs_test - pred_test)^2, na.rm = TRUE)
    test_logloss     <- log_loss_fun(obs_test, pred_test)
    test_sensitivity <- test_cls$sensitivity
    test_specificity <- test_cls$specificity
    test_tss         <- test_cls$tss
    
    # Sum log-likelihood (from full model)
    sum_loglikelihood <- m$sum_loglik
    
    # Save fold metrics including threshold
    traintest.df[[i]] <- data.frame(
      fold = i,
      sum_loglikelihood = sum_loglikelihood,
      train_auc = train_auc,
      train_tjur = train_tjur,
      train_rmse = train_rmse,
      train_brier = train_brier,
      train_logloss = train_logloss,
      train_sensitivity = train_sensitivity,
      train_specificity = train_specificity,
      train_tss = train_tss,
      test_auc = test_auc,
      test_tjur = test_tjur,
      test_rmse = test_rmse,
      test_brier = test_brier,
      test_logloss = test_logloss,
      test_sensitivity = test_sensitivity,
      test_specificity = test_specificity,
      test_tss = test_tss,
      threshold = threshold
    )
  }
  
  # Combine folds and calculate mean metrics
  traintest.df <- do.call(rbind, traintest.df)
  
  mean_metrics <- traintest.df %>%
    dplyr::summarise(
      mean_sum_loglikelihood = mean(sum_loglikelihood, na.rm = TRUE),
      mean_train_auc = mean(train_auc, na.rm = TRUE),
      mean_train_tjur = mean(train_tjur, na.rm = TRUE),
      mean_train_rmse = mean(train_rmse, na.rm = TRUE),
      mean_train_brier = mean(train_brier, na.rm = TRUE),
      mean_train_logloss = mean(train_logloss, na.rm = TRUE),
      mean_train_sensitivity = mean(train_sensitivity, na.rm = TRUE),
      mean_train_specificity = mean(train_specificity, na.rm = TRUE),
      mean_train_tss = mean(train_tss, na.rm = TRUE),
      mean_test_auc = mean(test_auc, na.rm = TRUE),
      mean_test_tjur = mean(test_tjur, na.rm = TRUE),
      mean_test_rmse = mean(test_rmse, na.rm = TRUE),
      mean_test_brier = mean(test_brier, na.rm = TRUE),
      mean_test_logloss = mean(test_logloss, na.rm = TRUE),
      mean_test_sensitivity = mean(test_sensitivity, na.rm = TRUE),
      mean_test_specificity = mean(test_specificity, na.rm = TRUE),
      mean_test_tss = mean(test_tss, na.rm = TRUE),
      mean_threshold = mean(threshold, na.rm = TRUE)
    )
  
  return(list(per_fold = traintest.df, summary = mean_metrics))
}


# Calculate threshold
calcThresh <- function( x ){
  obspred <- data.frame( PlotID=1:nrow(x),
                         Observed=x[,'presence'],
                         Predicted=x[,'fitted_vals'] )
  # Select optimal threshold
  thresh <- PresenceAbsence::optimal.thresholds(
    obspred, na.rm=T)
  # Return threshold
  return( thresh )
  return(obspred)
}

find_best_threshold <- function(obs,
                                pred,
                                step = 0.01){
  
  thresholds <- seq(0,1,by = step)
  
  scores <- lapply(
    thresholds,
    function(th){
      
      confusion_stats(
        obs,
        pred,
        th
      )
      
    })
  
  tss <- sapply(scores,
                function(x) x$TSS)
  
  best <- which.max(tss)
  
  list(
    
    threshold = thresholds[best],
    
    sensitivity = scores[[best]]$sensitivity,
    
    specificity = scores[[best]]$specificity,
    
    TSS = scores[[best]]$TSS
    
  )
  
}

# evalulate forecast 
find_best_tss_threshold <- function(obs, pred) {
  thresholds <- sort(unique(pred))
  
  tss_scores <- sapply(thresholds, function(th) {
    pred_bin <- as.integer(pred >= th)
    
    sens <- if (sum(obs == 1) > 0) {
      sum(pred_bin == 1 & obs == 1) / sum(obs == 1)
    } else {
      NA_real_
    }
    
    spec <- if (sum(obs == 0) > 0) {
      sum(pred_bin == 0 & obs == 0) / sum(obs == 0)
    } else {
      NA_real_
    }
    
    sens + spec - 1
  })
  
  best_i <- which.max(tss_scores)
  
  list(
    threshold = thresholds[best_i],
    TSS = tss_scores[best_i]
  )
}

get_global_cv_threshold_sdmTMB <- function(data_train,
                                           formula,
                                           coastline,
                                           spatial = FALSE,
                                           family = binomial(link = "logit"),
                                           v = 10,
                                           seed = 123) {
  
  data_train <- data_train %>%
    dplyr::mutate(row_id_cv = dplyr::row_number())
  
  set.seed(seed)
  folds <- rsample::vfold_cv(data_train, v = v, strata = "presence")
  
  pred_oof <- rep(NA_real_, nrow(data_train))
  
  for (i in seq_along(folds$splits)) {
    split <- folds$splits[[i]]
    analysis_data <- rsample::analysis(split)
    assessment_data <- rsample::assessment(split)
    
    mesh <- make_mesh(data = analysis_data, xy_cols = c("X", "Y"), cutoff = 55)
    barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = FALSE)
    
    fit <- sdmTMB(
      formula = formula,
      mesh = barrier_mesh,
      priors = sdmTMBpriors(b = normal(0, 1)),
      family = family,
      spatial = spatial,
      data = analysis_data
    )
    
    pred <- plogis(predict(fit, newdata = assessment_data)$est)
    pred_oof[assessment_data$row_id_cv] <- pred
  }
  
  keep <- !is.na(pred_oof)
  
  best <- find_best_tss_threshold(
    obs = data_train$presence[keep],
    pred = pred_oof[keep]
  )
  
  list(
    threshold = best$threshold,
    TSS = best$TSS,
    oof_pred = pred_oof
  )
}

evaluate_forecast <- function(obs_train,
                              pred_train,
                              obs_test,
                              pred_test,
                              threshold) {
  
  obs_train  <- as.numeric(obs_train)
  obs_test   <- as.numeric(obs_test)
  pred_train <- as.numeric(pred_train)
  pred_test  <- as.numeric(pred_test)
  threshold  <- as.numeric(threshold)
  
  pred_train_bin <- as.numeric(pred_train >= threshold)
  pred_test_bin  <- as.numeric(pred_test >= threshold)
  
  cm_train <- table(
    factor(obs_train, levels = c(0, 1)),
    factor(pred_train_bin, levels = c(0, 1))
  )
  cm_test <- table(
    factor(obs_test, levels = c(0, 1)),
    factor(pred_test_bin, levels = c(0, 1))
  )
  
  sens_train <- ifelse(sum(cm_train["1", ]) > 0, cm_train["1", "1"] / sum(cm_train["1", ]), NA)
  spec_train <- ifelse(sum(cm_train["0", ]) > 0, cm_train["0", "0"] / sum(cm_train["0", ]), NA)
  sens_test  <- ifelse(sum(cm_test["1", ]) > 0, cm_test["1", "1"] / sum(cm_test["1", ]), NA)
  spec_test  <- ifelse(sum(cm_test["0", ]) > 0, cm_test["0", "0"] / sum(cm_test["0", ]), NA)
  
  TSS_train <- sens_train + spec_train - 1
  TSS_test  <- sens_test + spec_test - 1
  
  Tjur_train <- mean(pred_train[obs_train == 1], na.rm = TRUE) -
    mean(pred_train[obs_train == 0], na.rm = TRUE)
  Tjur_test <- mean(pred_test[obs_test == 1], na.rm = TRUE) -
    mean(pred_test[obs_test == 0], na.rm = TRUE)
  
  RMSE_train <- sqrt(mean((obs_train - pred_train)^2, na.rm = TRUE))
  RMSE_test  <- sqrt(mean((obs_test - pred_test)^2, na.rm = TRUE))
  
  Brier_train <- mean((obs_train - pred_train)^2, na.rm = TRUE)
  Brier_test  <- mean((obs_test - pred_test)^2, na.rm = TRUE)
  
  eps <- 1e-15
  pred_train_clip <- pmin(pmax(pred_train, eps), 1 - eps)
  pred_test_clip  <- pmin(pmax(pred_test, eps), 1 - eps)
  
  LogLoss_train <- -mean(
    obs_train * log(pred_train_clip) +
      (1 - obs_train) * log(1 - pred_train_clip),
    na.rm = TRUE
  )
  LogLoss_test <- -mean(
    obs_test * log(pred_test_clip) +
      (1 - obs_test) * log(1 - pred_test_clip),
    na.rm = TRUE
  )
  
  AUC_train <- ModelMetrics::auc(obs_train, pred_train)
  AUC_test  <- ModelMetrics::auc(obs_test, pred_test)
  
  data.frame(
    dataset = c("training", "forecast"),
    AUC = c(AUC_train, AUC_test),
    TjurR2 = c(Tjur_train, Tjur_test),
    RMSE = c(RMSE_train, RMSE_test),
    Brier = c(Brier_train, Brier_test),
    LogLoss = c(LogLoss_train, LogLoss_test),
    sensitivity = c(sens_train, sens_test),
    specificity = c(spec_train, spec_test),
    TSS = c(TSS_train, TSS_test),
    threshold = c(threshold, threshold)
  )
}

get_cv_threshold_sdmTMB <- function(data_train,
                                    formula,
                                    mesh_builder,
                                    coastline,
                                    spatial = FALSE,
                                    family = binomial(link = "logit"),
                                    v = 10,
                                    seed = 123) {
  
  set.seed(seed)
  folds <- rsample::vfold_cv(data_train, v = v, strata = "presence")
  
  pred_oof <- rep(NA_real_, nrow(data_train))
  
  for (i in seq_along(folds$splits)) {
    split <- folds$splits[[i]]
    analysis_data <- rsample::analysis(split)
    assessment_data <- rsample::assessment(split)
    
    mesh <- make_mesh(data = analysis_data, xy_cols = c("X", "Y"), cutoff = 55)
    barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = FALSE)
    
    fit <- sdmTMB(
      formula = formula,
      mesh = barrier_mesh,
      priors = sdmTMBpriors(b = normal(0, 1)),
      family = family,
      spatial = spatial,
      data = analysis_data
    )
    
    pred <- plogis(predict(fit, newdata = assessment_data)$est)
    
    pred_oof[assessment_data$row_id] <- pred
  }
  
  keep <- !is.na(pred_oof)
  best <- find_best_tss_threshold(
    obs = data_train$presence[keep],
    pred = pred_oof[keep]
  )
  
  list(
    threshold = best$threshold,
    TSS = best$TSS,
    oof_pred = pred_oof
  )
}

evalfmod <- function(x, thresh, sp = NA) {
  
  eval.df <- data.frame()
  
  obspred <- data.frame(
    PlotID = 1:nrow(x),
    Observed = x[, "presence"],
    Predicted = x[, "fitted_vals"]
  )
  
  # Threshold for max sensitivity + specificity
  tss_thresh <- thresh$Predicted[thresh$Method == "MaxSens+Spec"]
  
  # Get confusion matrix based on TSS threshold
  cmx_tss <- PresenceAbsence::cmx(
    DATA = obspred,
    which.model = 1,
    thresh = tss_thresh
  )
  
  true_neg <- cmx_tss[1, 1]
  false_neg <- cmx_tss[1, 2]
  false_pos <- cmx_tss[2, 1]
  true_pos <- cmx_tss[2, 2]
  
  sensitivity <- true_pos / (true_pos + false_neg)
  specificity <- true_neg / (true_neg + false_pos)
  TSS <- sensitivity + specificity - 1
  
  # Get confusion matrix based on Kappa
  cmx_kappa <- PresenceAbsence::cmx(
    DATA = obspred,
    which.model = 1,
    thresh = thresh$Predicted[thresh$Method == "MaxKappa"]
  )
  
  kappa <- PresenceAbsence::Kappa(CMX = cmx_kappa, st.dev = TRUE)
  
  # Brier score
  brier <- mean((x$fitted_vals - x$presence)^2, na.rm = TRUE)
  
  # Log loss
  eps <- 1e-15
  pred_clipped <- pmin(pmax(x$fitted_vals, eps), 1 - eps)
  log_loss <- -mean(
    x$presence * log(pred_clipped) +
      (1 - x$presence) * log(1 - pred_clipped),
    na.rm = TRUE
  )
  
  miller <- modEvA::MillerCalib(
    model = NULL,
    obs = x$presence,
    pred = x$fitted_vals
  )
  
  eer <- modEvA::errorMeasures(
    model = NULL,
    obs = x$presence,
    pred = x$fitted_vals
  )
  
  hlgof_quant <- modEvA::HLfit(
    model = NULL,
    obs = x$presence,
    pred = x$fitted_vals,
    bin.method = "quantiles",
    n.bins = 3000
  )
  
  hlgof_prob <- modEvA::HLfit(
    model = NULL,
    obs = x$presence,
    pred = x$fitted_vals,
    bin.method = "prob.bins"
  )
  
  hlgof_nbin <- modEvA::HLfit(
    model = NULL,
    obs = x$presence,
    pred = x$fitted_vals,
    bin.method = "n.bins",
    n.bins = 8
  )
  
  eval.df <- data.frame(
    kappa = if (is.list(kappa)) kappa$Kappa else kappa,
    TSS = TSS,
    sensitivity = sensitivity,
    specificity = specificity,
    brier = brier,
    log_loss = log_loss,
    miller = miller,
    eer = eer,
    hlgof_quant_chisq = hlgof_quant$chi.sq,
    hlgof_quant_pvalue = hlgof_quant$p.value,
    hlgof_quant_RMSE = hlgof_quant$RMSE,
    hlgof_prob_chisq = hlgof_prob$chi.sq,
    hlgof_prob_pvalue = hlgof_prob$p.value,
    hlgof_prob_RMSE = hlgof_prob$RMSE,
    hlgof_nbin_chisq = hlgof_nbin$chi.sq,
    hlgof_nbin_pvalue = hlgof_nbin$p.value,
    hlgof_nbin_RMSE = hlgof_nbin$RMSE,
    species = sp
  )
  
  return(eval.df)
}





# sdmtmb predictions

PredictSDM_bySurvey <- function(env, model, survey_type, species, model_name) {
  
  outdir <- file.path(
    "code/output_data/seagrass_predictions/survey",
    species,
    model_name
  )
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (s in survey_type) {
    message("Predicting survey: ", s)
    
    env$Survey <- as.factor(s)
    
    preds <- predict(model, newdata = env)
    sims  <- predict(model, newdata = env, nsim = 100)
    preds$SE <- apply(sims, 1, sd)
    
    pred <- env %>%
      dplyr::select(X_m, Y_m, X, Y, ID, region, Survey) %>%
      bind_cols(preds %>% dplyr::select(est, SE))
    
    # Save survey-level predictions with model name
    save(
      pred,
      file = file.path(
        outdir,
        paste0("Prediction_", species, "_", model_name, "_", s, ".RData")
      )
    )
    
    rm(pred, preds)
    gc()
  }
  
  invisible(TRUE)
}

AverageSurveyPredictions <- function(species, model_name) {
  
  indir <- file.path(
    "code/output_data/seagrass_predictions/survey",
    species,
    model_name
  )
  
  files <- list.files(
    indir,
    pattern = "^Prediction_.*\\.RData$",
    full.names = TRUE
  )
  
  message("Loading ", length(files), " survey prediction files...")
  
  all_preds <- lapply(files, function(f) {
    load(f)
    pred
  }) |> bind_rows()
  
  mean_pred <- all_preds %>%
    group_by(X_m, Y_m, ID, region) %>%
    summarise(
      est = mean(est, na.rm = TRUE),
      SE  = mean(SE,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(ID)
  
  save(
    mean_pred,
    file = file.path(
      indir,
      paste0("MeanSurveyPreds_", species, "_", model_name, ".RData")
    )
  )
  
  mean_pred
}


#Rasterizing function
rasterize_eelgrass_year <- function(year, poly_data, line_data, point_data, template_rast) {
  poly_year <- poly_data[poly_data$Year == year, ]
  line_year <- line_data[line_data$Year == year, ]
  point_year <- point_data[point_data$Year == year, ]
  
  poly_r <- rasterize(poly_year, template_rast, field = 1, background = NA, touches = TRUE)
  line_r <- rasterize(line_year, template_rast, field = 1, background = NA, touches = TRUE)
  point_r <- rasterize(point_year, template_rast, field = 1, background = NA, touches = TRUE)
  
  result_r <- ifel(!is.na(poly_r) | !is.na(line_r) | !is.na(point_r), 1, NA)
  return(result_r)
}


evaluate_occurrence <- function(df, observed_col, predicted_col) {
  
  # Convert column names to strings
  obs_name <- deparse(substitute(observed_col))
  pred_name <- deparse(substitute(predicted_col))
  
  # Filter out NAs
  df <- df %>%
    filter(!is.na(.data[[obs_name]]), !is.na(.data[[pred_name]])) %>%
    mutate(
      Presence = factor(.data[[obs_name]], levels = c(0,1), labels = c("Absent", "Present"))
    )
  
  # Extract vectors for statistical tests
  pred_vec <- df[[pred_name]]
  group_vec <- df$Presence
  
  # Mann-Whitney U test
  wilcox_test <- wilcox.test(pred_vec ~ group_vec, exact = FALSE)
  
  # Cliff's delta
  cliffs_delta <- cliff.delta(pred_vec ~ group_vec)
  
  # Summary stats per group
  summary_stats <- df %>%
    group_by(Presence) %>%
    summarise(
      n = n(),
      mean_pred = mean(.data[[pred_name]], na.rm = TRUE),
      median_pred = median(.data[[pred_name]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Return results as a list
  list(
    wilcox_test = wilcox_test,
    cliffs_delta = cliffs_delta,
    summary_stats = summary_stats
  )
}


#### ML functions 
###############################################################################
nested_gbm_spatial <- function(
    dat,
    predictors,
    response = "presence",
    outer_col = "outer_fold",
    inner_col = "inner_fold",
    gbm_grid){
  
  library(dplyr)
  library(gbm)
  library(pROC)
  
  calc_logloss <- function(actual, prob) {
    eps <- 1e-15
    prob_clip <- pmin(pmax(prob, eps), 1 - eps)
    -mean(actual * log(prob_clip) + (1 - actual) * log(1 - prob_clip), na.rm = TRUE)
  }
  
  calc_tjur <- function(actual, prob) {
    mean(prob[actual == 1], na.rm = TRUE) - mean(prob[actual == 0], na.rm = TRUE)
  }
  
  calc_brier <- function(actual, prob) {
    mean((prob - actual)^2, na.rm = TRUE)
  }
  
  class_metrics <- function(actual, prob, threshold) {
    pred_class <- as.integer(prob >= threshold)
    
    sensitivity <- if (sum(actual == 1) > 0) {
      sum(pred_class == 1 & actual == 1) / sum(actual == 1)
    } else {
      NA_real_
    }
    
    specificity <- if (sum(actual == 0) > 0) {
      sum(pred_class == 0 & actual == 0) / sum(actual == 0)
    } else {
      NA_real_
    }
    
    tss <- sensitivity + specificity - 1
    
    list(
      sensitivity = sensitivity,
      specificity = specificity,
      tss = tss
    )
  }
  
  get_optimal_threshold <- function(actual, prob) {
    keep <- !is.na(actual) & !is.na(prob)
    
    actual <- actual[keep]
    prob <- prob[keep]
    
    thresholds <- sort(unique(prob))
    
    if (length(thresholds) == 0) {
      return(NA_real_)
    }
    
    tss_scores <- numeric(length(thresholds))
    
    for (j in seq_along(thresholds)) {
      th <- thresholds[j]
      pred_class <- as.integer(prob >= th)
      
      sens <- if (sum(actual == 1) > 0) {
        sum(pred_class == 1 & actual == 1) / sum(actual == 1)
      } else {
        NA_real_
      }
      
      spec <- if (sum(actual == 0) > 0) {
        sum(pred_class == 0 & actual == 0) / sum(actual == 0)
      } else {
        NA_real_
      }
      
      tss_scores[j] <- sens + spec - 1
    }
    
    thresholds[which.max(tss_scores)]
  }
  
  safe_auc <- function(actual, prob) {
    tryCatch(as.numeric(pROC::auc(actual, prob)), error = function(e) NA_real_)
  }
  
  outer_results <- list()
  outer_folds <- sort(unique(dat[[outer_col]]))
  
  for (ofold in outer_folds) {
    
    message("Outer Fold: ", ofold)
    
    train_outer <- dat %>%
      filter(.data[[outer_col]] != ofold)
    
    test_outer <- dat %>%
      filter(.data[[outer_col]] == ofold)
    
    inner_scores <- list()
    
    for (g in seq_len(nrow(gbm_grid))) {
      oof_pred <- rep(NA_real_, nrow(train_outer))
      oof_obs <- train_outer[[response]]
      pars <- gbm_grid[g, ]
      inner_res <- c()
      inner_best_trees <- c()
      inner_folds <- sort(unique(train_outer[[inner_col]]))
      
      for (ifold in inner_folds) {
        
        train_inner <- train_outer %>%
          filter(.data[[inner_col]] != ifold)
        
        valid_inner <- train_outer %>%
          filter(.data[[inner_col]] == ifold)
        
        mod <- gbm(
          formula = as.formula(
            paste(response, "~", paste(predictors, collapse = "+"))
          ),
          data = train_inner,
          distribution = "bernoulli",
          n.trees = 5000,
          interaction.depth = pars$interaction.depth,
          shrinkage = pars$shrinkage,
          n.minobsinnode = pars$n.minobsinnode,
          bag.fraction = 0.7,
          train.fraction = 0.8,
          verbose = FALSE
        )
        
        best_trees <- gbm.perf(
          mod,
          method = "OOB",
          plot.it = FALSE
        )
        
        if (is.null(best_trees) || is.na(best_trees)) {
          best_trees <- 5000
        }
        
        pred <- predict(
          mod,
          valid_inner,
          n.trees = best_trees,
          type = "response"
        )
        
        valid_rows <- which(train_outer[[inner_col]] == ifold)
        
        oof_pred[valid_rows] <- pred
        
        inner_res <- c(
          inner_res,
          calc_logloss(valid_inner[[response]], pred)
        )
        
        inner_best_trees <- c(
          inner_best_trees,
          best_trees
        )
      }
      
      oof_threshold <- get_optimal_threshold(
        actual = oof_obs,
        prob = oof_pred
      )
      
      inner_scores[[g]] <- data.frame(
        interaction.depth = pars$interaction.depth,
        shrinkage = pars$shrinkage,
        n.minobsinnode = pars$n.minobsinnode,
        mean_logloss = mean(inner_res, na.rm = TRUE),
        sd_logloss = sd(inner_res, na.rm = TRUE),
        mean_trees = mean(inner_best_trees, na.rm = TRUE),
        threshold = oof_threshold
      )
    }
    
    tuning <- bind_rows(inner_scores) %>%
      arrange(mean_logloss)
    
    best <- tuning[1, ]
    
    final_mod <- gbm(
      formula = as.formula(
        paste(response, "~", paste(predictors, collapse = "+"))
      ),
      data = train_outer,
      distribution = "bernoulli",
      n.trees = 5000,
      interaction.depth = best$interaction.depth,
      shrinkage = best$shrinkage,
      n.minobsinnode = best$n.minobsinnode,
      bag.fraction = 0.7,
      train.fraction = 0.8,
      verbose = FALSE
    )
    
    final_trees <- round(best$mean_trees)
    
    pred_train <- predict(
      final_mod,
      train_outer,
      n.trees = final_trees,
      type = "response"
    )
    
    pred_test <- predict(
      final_mod,
      test_outer,
      n.trees = final_trees,
      type = "response"
    )
    
    obs_train <- train_outer[[response]]
    obs_test <- test_outer[[response]]
    
    threshold <- best$threshold
    
    train_cls <- class_metrics(obs_train, pred_train, threshold)
    test_cls <- class_metrics(obs_test, pred_test, threshold)
    
    outer_results[[length(outer_results) + 1]] <- data.frame(
      model = "GBM",
      outer_fold = ofold,
      interaction.depth = best$interaction.depth,
      shrinkage = best$shrinkage,
      n.minobsinnode = best$n.minobsinnode,
      n.trees = final_trees,
      threshold = threshold,
      train_auc = safe_auc(obs_train, pred_train),
      train_tjur = calc_tjur(obs_train, pred_train),
      train_brier = calc_brier(obs_train, pred_train),
      train_logloss = calc_logloss(obs_train, pred_train),
      train_sensitivity = train_cls$sensitivity,
      train_specificity = train_cls$specificity,
      train_tss = train_cls$tss,
      test_auc = safe_auc(obs_test, pred_test),
      test_tjur = calc_tjur(obs_test, pred_test),
      test_brier = calc_brier(obs_test, pred_test),
      test_logloss = calc_logloss(obs_test, pred_test),
      test_sensitivity = test_cls$sensitivity,
      test_specificity = test_cls$specificity,
      test_tss = test_cls$tss
    )
  }
  
  fold_metrics <- bind_rows(outer_results)
  
  list(
    fold_metrics = fold_metrics,
    summary = fold_metrics %>%
      summarise(
        mean_threshold = mean(threshold, na.rm = TRUE),
        mean_train_auc = mean(train_auc, na.rm = TRUE),
        mean_train_tjur = mean(train_tjur, na.rm = TRUE),
        mean_train_brier = mean(train_brier, na.rm = TRUE),
        mean_train_logloss = mean(train_logloss, na.rm = TRUE),
        mean_train_sensitivity = mean(train_sensitivity, na.rm = TRUE),
        mean_train_specificity = mean(train_specificity, na.rm = TRUE),
        mean_train_tss = mean(train_tss, na.rm = TRUE),
        mean_test_auc = mean(test_auc, na.rm = TRUE),
        mean_test_tjur = mean(test_tjur, na.rm = TRUE),
        mean_test_brier = mean(test_brier, na.rm = TRUE),
        mean_test_logloss = mean(test_logloss, na.rm = TRUE),
        mean_test_sensitivity = mean(test_sensitivity, na.rm = TRUE),
        mean_test_specificity = mean(test_specificity, na.rm = TRUE),
        mean_test_tss = mean(test_tss, na.rm = TRUE)
      )
  )
}


nested_xgb_spatial <- function(
    dat,
    predictors,
    response = "presence",
    outer_col = "outer_fold",
    inner_col = "inner_fold",
    xgb_grid
) {
  
  library(dplyr)
  library(xgboost)
  library(pROC)
  
  calc_logloss <- function(actual, prob) {
    eps <- 1e-15
    prob_clip <- pmin(pmax(prob, eps), 1 - eps)
    -mean(actual * log(prob_clip) + (1 - actual) * log(1 - prob_clip), na.rm = TRUE)
  }
  
  calc_tjur <- function(actual, prob) {
    mean(prob[actual == 1], na.rm = TRUE) - mean(prob[actual == 0], na.rm = TRUE)
  }
  
  calc_brier <- function(actual, prob) {
    mean((prob - actual)^2, na.rm = TRUE)
  }
  
  class_metrics <- function(actual, prob, threshold) {
    pred_class <- as.integer(prob >= threshold)
    
    sensitivity <- if (sum(actual == 1) > 0) {
      sum(pred_class == 1 & actual == 1) / sum(actual == 1)
    } else {
      NA_real_
    }
    
    specificity <- if (sum(actual == 0) > 0) {
      sum(pred_class == 0 & actual == 0) / sum(actual == 0)
    } else {
      NA_real_
    }
    
    tss <- sensitivity + specificity - 1
    
    list(
      sensitivity = sensitivity,
      specificity = specificity,
      tss = tss
    )
  }
  
  get_optimal_threshold <- function(actual, prob) {
    keep <- !is.na(actual) & !is.na(prob)
    
    actual <- actual[keep]
    prob <- prob[keep]
    
    thresholds <- sort(unique(prob))
    
    if (length(thresholds) == 0) {
      return(NA_real_)
    }
    
    tss_scores <- numeric(length(thresholds))
    
    for (j in seq_along(thresholds)) {
      th <- thresholds[j]
      pred_class <- as.integer(prob >= th)
      
      sens <- if (sum(actual == 1) > 0) {
        sum(pred_class == 1 & actual == 1) / sum(actual == 1)
      } else {
        NA_real_
      }
      
      spec <- if (sum(actual == 0) > 0) {
        sum(pred_class == 0 & actual == 0) / sum(actual == 0)
      } else {
        NA_real_
      }
      
      tss_scores[j] <- sens + spec - 1
    }
    
    thresholds[which.max(tss_scores)]
  }
  
  safe_auc <- function(actual, prob) {
    tryCatch(as.numeric(pROC::auc(actual, prob)), error = function(e) NA_real_)
  }
  
  outer_results <- list()
  outer_folds <- sort(unique(dat[[outer_col]]))
  
  for (ofold in outer_folds) {
    
    message("Outer Fold: ", ofold)
    
    train_outer <- dat %>%
      filter(.data[[outer_col]] != ofold)
    
    test_outer <- dat %>%
      filter(.data[[outer_col]] == ofold)
    
    inner_scores <- list()
    
    for (g in seq_len(nrow(xgb_grid))) {
      oof_pred <- rep(NA_real_, nrow(train_outer))
      oof_obs <- train_outer[[response]]
      pars <- xgb_grid[g, ]
      loglosses <- c()
      inner_folds <- sort(unique(train_outer[[inner_col]]))
      
      for (ifold in inner_folds) {
        
        train_inner <- train_outer %>%
          filter(.data[[inner_col]] != ifold)
        
        valid_inner <- train_outer %>%
          filter(.data[[inner_col]] == ifold)
        
        X_train <- model.matrix(
          ~ . - 1,
          train_inner[, predictors, drop = FALSE]
        )
        
        X_valid <- model.matrix(
          ~ . - 1,
          valid_inner[, predictors, drop = FALSE]
        )
        
        dtrain <- xgb.DMatrix(
          data = X_train,
          label = train_inner[[response]]
        )
        
        mod <- xgb.train(
          params = list(
            objective = "binary:logistic",
            eval_metric = "logloss",
            eta = pars$eta,
            max_depth = pars$max_depth,
            subsample = pars$subsample,
            colsample_bytree = pars$colsample_bytree
          ),
          data = dtrain,
          nrounds = pars$nrounds,
          verbose = 0
        )
        
        pred <- predict(mod, X_valid)
        
        valid_rows <- which(train_outer[[inner_col]] == ifold)
        
        oof_pred[valid_rows] <- pred
        
        loglosses <- c(
          loglosses,
          calc_logloss(valid_inner[[response]], pred)
        )
      }
      
      oof_threshold <- get_optimal_threshold(
        actual = oof_obs,
        prob = oof_pred
      )
      
      inner_scores[[g]] <- data.frame(
        eta = pars$eta,
        max_depth = pars$max_depth,
        subsample = pars$subsample,
        colsample_bytree = pars$colsample_bytree,
        nrounds = pars$nrounds,
        mean_logloss = mean(loglosses, na.rm = TRUE),
        sd_logloss = sd(loglosses, na.rm = TRUE),
        threshold = oof_threshold
      )
    }
    
    tuning <- bind_rows(inner_scores) %>%
      arrange(mean_logloss)
    
    best <- tuning[1, ]
    
    X_train_outer <- model.matrix(
      ~ . - 1,
      train_outer[, predictors, drop = FALSE]
    )
    
    X_test_outer <- model.matrix(
      ~ . - 1,
      test_outer[, predictors, drop = FALSE]
    )
    
    dtrain_outer <- xgb.DMatrix(
      data = X_train_outer,
      label = train_outer[[response]]
    )
    
    final_mod <- xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        eta = best$eta,
        max_depth = best$max_depth,
        subsample = best$subsample,
        colsample_bytree = best$colsample_bytree
      ),
      data = dtrain_outer,
      nrounds = best$nrounds,
      verbose = 0
    )
    
    pred_train <- predict(final_mod, X_train_outer)
    pred_test <- predict(final_mod, X_test_outer)
    
    obs_train <- train_outer[[response]]
    obs_test <- test_outer[[response]]
    
    threshold <- best$threshold
    
    train_cls <- class_metrics(obs_train, pred_train, threshold)
    test_cls <- class_metrics(obs_test, pred_test, threshold)
    
    outer_results[[length(outer_results) + 1]] <- data.frame(
      model = "XGBOOST",
      outer_fold = ofold,
      eta = best$eta,
      max_depth = best$max_depth,
      subsample = best$subsample,
      colsample_bytree = best$colsample_bytree,
      nrounds = best$nrounds,
      threshold = threshold,
      train_auc = safe_auc(obs_train, pred_train),
      train_tjur = calc_tjur(obs_train, pred_train),
      train_brier = calc_brier(obs_train, pred_train),
      train_logloss = calc_logloss(obs_train, pred_train),
      train_sensitivity = train_cls$sensitivity,
      train_specificity = train_cls$specificity,
      train_tss = train_cls$tss,
      test_auc = safe_auc(obs_test, pred_test),
      test_tjur = calc_tjur(obs_test, pred_test),
      test_brier = calc_brier(obs_test, pred_test),
      test_logloss = calc_logloss(obs_test, pred_test),
      test_sensitivity = test_cls$sensitivity,
      test_specificity = test_cls$specificity,
      test_tss = test_cls$tss
    )
  }
  
  fold_metrics <- bind_rows(outer_results)
  
  list(
    fold_metrics = fold_metrics,
    summary = fold_metrics %>%
      summarise(
        mean_threshold = mean(threshold, na.rm = TRUE),
        mean_train_auc = mean(train_auc, na.rm = TRUE),
        mean_train_tjur = mean(train_tjur, na.rm = TRUE),
        mean_train_brier = mean(train_brier, na.rm = TRUE),
        mean_train_logloss = mean(train_logloss, na.rm = TRUE),
        mean_train_sensitivity = mean(train_sensitivity, na.rm = TRUE),
        mean_train_specificity = mean(train_specificity, na.rm = TRUE),
        mean_train_tss = mean(train_tss, na.rm = TRUE),
        mean_test_auc = mean(test_auc, na.rm = TRUE),
        mean_test_tjur = mean(test_tjur, na.rm = TRUE),
        mean_test_brier = mean(test_brier, na.rm = TRUE),
        mean_test_logloss = mean(test_logloss, na.rm = TRUE),
        mean_test_sensitivity = mean(test_sensitivity, na.rm = TRUE),
        mean_test_specificity = mean(test_specificity, na.rm = TRUE),
        mean_test_tss = mean(test_tss, na.rm = TRUE)
      )
  )
}

get_global_cv_threshold_xgb <- function(
    data_train,
    predictors,
    response = "presence",
    eta,
    max_depth,
    subsample,
    colsample_bytree,
    nrounds,
    v = 10,
    seed = 123
) {
  
  library(dplyr)
  library(rsample)
  library(xgboost)
  
  data_train <- data_train %>%
    dplyr::mutate(row_id_cv = dplyr::row_number())
  
  set.seed(seed)
  folds <- rsample::vfold_cv(
    data_train,
    v = v,
    strata = response
  )
  
  pred_oof <- rep(NA_real_, nrow(data_train))
  
  for (i in seq_along(folds$splits)) {
    split <- folds$splits[[i]]
    analysis_data <- rsample::analysis(split)
    assessment_data <- rsample::assessment(split)
    
    X_analysis <- model.matrix(
      ~ . - 1,
      data = analysis_data[, predictors, drop = FALSE]
    )
    
    X_assessment <- model.matrix(
      ~ . - 1,
      data = assessment_data[, predictors, drop = FALSE]
    )
    
    dtrain <- xgb.DMatrix(
      data = X_analysis,
      label = analysis_data[[response]]
    )
    
    xgb_fit <- xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        eta = eta,
        max_depth = max_depth,
        subsample = subsample,
        colsample_bytree = colsample_bytree
      ),
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )
    
    pred <- predict(xgb_fit, X_assessment)
    pred_oof[assessment_data$row_id_cv] <- pred
  }
  
  keep <- !is.na(pred_oof)
  
  best <- find_best_tss_threshold(
    obs = data_train[[response]][keep],
    pred = pred_oof[keep]
  )
  
  list(
    threshold = best$threshold,
    TSS = best$TSS,
    oof_pred = pred_oof
  )
}

run_temporal_gbm_forecast <- function(
    data_pre2013,
    test_set,
    predictors,
    threshold,
    response = "presence",
    n.trees,
    interaction.depth,
    shrinkage,
    n.minobsinnode,
    v = 10,
    seed = 123,
    model_name = "GBM"
) {
  
  library(gbm)
  library(rsample)
  library(dplyr)
  
  set.seed(seed)
  
  folds <- rsample::vfold_cv(
    data_pre2013,
    v = v,
    strata = response
  )
  
  results <- list()
  
  for (i in seq_along(folds$splits)) {
    
    split <- folds$splits[[i]]
    train_data <- rsample::analysis(split)
    
    obs_train <- train_data[[response]]
    obs_test <- test_set[[response]]
    
    gbm_fit <- gbm(
      formula = as.formula(
        paste(response, "~", paste(predictors, collapse = "+"))
      ),
      distribution = "bernoulli",
      data = train_data,
      n.trees = n.trees,
      interaction.depth = interaction.depth,
      shrinkage = shrinkage,
      n.minobsinnode = n.minobsinnode,
      bag.fraction = 0.7,
      train.fraction = 1,
      verbose = FALSE
    )
    
    pred_train <- predict(
      gbm_fit,
      newdata = train_data,
      n.trees = n.trees,
      type = "response"
    )
    
    pred_test <- predict(
      gbm_fit,
      newdata = test_set,
      n.trees = n.trees,
      type = "response"
    )
    
    eval <- evaluate_forecast(
      obs_train = obs_train,
      pred_train = pred_train,
      obs_test = obs_test,
      pred_test = pred_test,
      threshold = threshold
    )
    
    eval$model <- model_name
    eval$fold <- i
    
    results[[i]] <- eval
  }
  
  bind_rows(results)
}

run_temporal_xgb_forecast <- function(
    data_pre2013,
    test_set,
    predictors,
    threshold,
    response = "presence",
    eta,
    max_depth,
    subsample,
    colsample_bytree,
    nrounds,
    v = 10,
    seed = 123,
    model_name = "XGBOOST"
) {
  
  library(xgboost)
  library(rsample)
  library(dplyr)
  
  set.seed(seed)
  
  folds <- rsample::vfold_cv(
    data_pre2013,
    v = v,
    strata = response
  )
  
  results <- list()
  
  for (i in seq_along(folds$splits)) {
    
    split <- folds$splits[[i]]
    train_data <- rsample::analysis(split)
    
    obs_train <- train_data[[response]]
    obs_test <- test_set[[response]]
    
    X_train <- model.matrix(
      ~ . - 1,
      data = train_data[, predictors, drop = FALSE]
    )
    
    X_test <- model.matrix(
      ~ . - 1,
      data = test_set[, predictors, drop = FALSE]
    )
    
    dtrain <- xgb.DMatrix(
      data = X_train,
      label = obs_train
    )
    
    xgb_fit <- xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        eta = eta,
        max_depth = max_depth,
        subsample = subsample,
        colsample_bytree = colsample_bytree
      ),
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )
    
    pred_train <- predict(xgb_fit, X_train)
    pred_test <- predict(xgb_fit, X_test)
    
    eval <- evaluate_forecast(
      obs_train = obs_train,
      pred_train = pred_train,
      obs_test = obs_test,
      pred_test = pred_test,
      threshold = threshold
    )
    
    eval$model <- model_name
    eval$fold <- i
    
    results[[i]] <- eval
  }
  
  bind_rows(results)
}

get_global_cv_threshold_gbm <- function(
    data_train,
    predictors,
    response = "presence",
    n.trees,
    interaction.depth,
    shrinkage,
    n.minobsinnode,
    v = 10,
    seed = 123
) {
  
  library(dplyr)
  library(rsample)
  library(gbm)
  
  data_train <- data_train %>%
    dplyr::mutate(row_id_cv = dplyr::row_number())
  
  set.seed(seed)
  folds <- rsample::vfold_cv(
    data_train,
    v = v,
    strata = response
  )
  
  pred_oof <- rep(NA_real_, nrow(data_train))
  
  for (i in seq_along(folds$splits)) {
    split <- folds$splits[[i]]
    analysis_data <- rsample::analysis(split)
    assessment_data <- rsample::assessment(split)
    
    gbm_fit <- gbm(
      formula = as.formula(
        paste(response, "~", paste(predictors, collapse = "+"))
      ),
      distribution = "bernoulli",
      data = analysis_data,
      n.trees = n.trees,
      interaction.depth = interaction.depth,
      shrinkage = shrinkage,
      n.minobsinnode = n.minobsinnode,
      bag.fraction = 0.7,
      train.fraction = 1,
      verbose = FALSE
    )
    
    pred <- predict(
      gbm_fit,
      newdata = assessment_data,
      n.trees = n.trees,
      type = "response"
    )
    
    pred_oof[assessment_data$row_id_cv] <- pred
  }
  
  keep <- !is.na(pred_oof)
  
  best <- find_best_tss_threshold(
    obs = data_train[[response]][keep],
    pred = pred_oof[keep]
  )
  
  list(
    threshold = best$threshold,
    TSS = best$TSS,
    oof_pred = pred_oof
  )
}






#relative importance for xgboost methods
perm_importance_xgb <- function(model, X, y, nrep = 10) {
  
  # baseline AUC
  baseline_pred <- predict(model, X)
  baseline_auc  <- as.numeric(pROC::auc(y, baseline_pred))
  
  p <- ncol(X)
  imp <- numeric(p)
  
  for (j in seq_len(p)) {
    aucs <- replicate(nrep, {
      X_perm <- X
      X_perm[, j] <- sample(X_perm[, j])  # permute column j
      pred <- predict(model, X_perm)
      as.numeric(pROC::auc(y, pred))
    })
    
    imp[j] <- baseline_auc - mean(aucs)
  }
  
  # clean up
  imp[imp < 0] <- 0
  rel_imp <- 100 * imp / sum(imp)
  
  data.frame(
    Feature = colnames(X),
    Importance = imp,
    RelImportance = rel_imp
  )
}



calc_independent_metrics <- function(data,
                                obs_col = "obs",
                                pa_set_col = "pa_set",
                                model_cols = NULL,
                                threshold = NULL,
                                threshold_method = "max_tss",
                                eps = 1e-15) {
  if (is.null(model_cols)) {
    model_cols <- data %>%
      st_drop_geometry() %>%
      select(where(is.numeric)) %>%
      select(-all_of(c(obs_col, pa_set_col))) %>%
      names()
  }
  dat <- data %>%
    sf::st_drop_geometry() %>%
    select(all_of(c(obs_col, pa_set_col, model_cols))) %>%
    mutate(
      "{obs_col}" := as.integer(.data[[obs_col]])
    )
  # presence rows have pa_set = NA, so repeat presences into each pseudoabsence set
  pa_sets <- sort(unique(dat[[pa_set_col]][!is.na(dat[[pa_set_col]])]))
  validation_by_set <- purrr::map(
    pa_sets,
    function(s) {
      dat %>%
        filter(.data[[obs_col]] == 1 | .data[[pa_set_col]] == s) %>%
        mutate("{pa_set_col}" := s)
    }
  ) %>%
    setNames(pa_sets)
  metric_results <- purrr::imap_dfr(
    validation_by_set,
    function(set_data, set_id) {
      purrr::map_dfr(
        model_cols,
        function(model) {
          df <- set_data %>%
            select(all_of(c(obs_col, model))) %>%
            filter(
              !is.na(.data[[obs_col]]),
              !is.na(.data[[model]])
            ) %>%
            mutate(
              obs_num = as.integer(.data[[obs_col]]),
              pred = pmin(pmax(as.numeric(.data[[model]]), eps), 1 - eps)
            )
          if (length(unique(df$obs_num)) < 2) {
            return(
              tibble(
                pa_set = as.integer(set_id),
                model = model,
                n = nrow(df),
                n_presence = sum(df$obs_num == 1),
                n_absence = sum(df$obs_num == 0),
                threshold = NA_real_,
                auc = NA_real_,
                tjur = NA_real_,
                brier = NA_real_,
                logloss = NA_real_,
                sensitivity = NA_real_,
                specificity = NA_real_,
                tss = NA_real_
              )
            )
          }
          auc_val <- yardstick::roc_auc_vec(
            truth = factor(df$obs_num, levels = c(0, 1)),
            estimate = df$pred,
            event_level = "second"
          )
          tjur_val <- mean(df$pred[df$obs_num == 1]) -
            mean(df$pred[df$obs_num == 0])
          brier_val <- mean((df$pred - df$obs_num)^2)
          logloss_val <- -mean(
            df$obs_num * log(df$pred) +
              (1 - df$obs_num) * log(1 - df$pred)
          )
          if (is.null(threshold)) {
            if (threshold_method == "max_tss") {
              threshold_grid <- sort(unique(df$pred))
              threshold_stats <- purrr::map_dfr(
                threshold_grid,
                function(th) {
                  pred_class <- ifelse(df$pred >= th, 1, 0)
                  tp <- sum(pred_class == 1 & df$obs_num == 1)
                  tn <- sum(pred_class == 0 & df$obs_num == 0)
                  fp <- sum(pred_class == 1 & df$obs_num == 0)
                  fn <- sum(pred_class == 0 & df$obs_num == 1)
                  sensitivity <- ifelse(
                    tp + fn > 0,
                    tp / (tp + fn),
                    NA_real_
                  )
                  specificity <- ifelse(
                    tn + fp > 0,
                    tn / (tn + fp),
                    NA_real_
                  )
                  tibble(
                    threshold = th,
                    sensitivity = sensitivity,
                    specificity = specificity,
                    tss = sensitivity + specificity - 1
                  )
                }
              )
              best_row <- threshold_stats %>%
                filter(tss == max(tss, na.rm = TRUE)) %>%
                arrange(desc(sensitivity), desc(specificity), threshold) %>%
                slice(1)
              threshold_use <- best_row$threshold
              sensitivity_val <- best_row$sensitivity
              specificity_val <- best_row$specificity
              tss_val <- best_row$tss
            } else {
              stop("Currently supported threshold_method: 'max_tss'.")
            }
          } else {
            if (length(threshold) == 1) {
              threshold_use <- threshold
            } else {
              threshold_use <- threshold[[model]]
            }
            pred_class <- ifelse(df$pred >= threshold_use, 1, 0)
            tp <- sum(pred_class == 1 & df$obs_num == 1)
            tn <- sum(pred_class == 0 & df$obs_num == 0)
            fp <- sum(pred_class == 1 & df$obs_num == 0)
            fn <- sum(pred_class == 0 & df$obs_num == 1)
            sensitivity_val <- ifelse(
              tp + fn > 0,
              tp / (tp + fn),
              NA_real_
            )
            specificity_val <- ifelse(
              tn + fp > 0,
              tn / (tn + fp),
              NA_real_
            )
            tss_val <- sensitivity_val + specificity_val - 1
          }
          tibble(
            pa_set = as.integer(set_id),
            model = model,
            n = nrow(df),
            n_presence = sum(df$obs_num == 1),
            n_absence = sum(df$obs_num == 0),
            threshold = threshold_use,
            auc = auc_val,
            tjur = tjur_val,
            brier = brier_val,
            logloss = logloss_val,
            sensitivity = sensitivity_val,
            specificity = specificity_val,
            tss = tss_val
          )
        }
      )
    }
  )
  metric_results
}




make_pseudoabsence_candidate_rast <- function(domain_rast,
                                              exclusion_rast,
                                              buffer_m = 100,
                                              filename = "",
                                              overwrite = TRUE) {
  # make sure rasters align
  if (!terra::compareGeom(domain_rast, exclusion_rast, stopOnError = FALSE)) {
    stop("domain_rast and exclusion_rast do not have the same geometry.")
  }
  # cells with mapped / observed / suspected eelgrass
  exclusion_binary <- terra::ifel(!is.na(exclusion_rast), 1, NA)
  # buffer exclusion cells by distance in map units
  # For EPSG:3005, map units are metres.
  if (buffer_m > 0) {
    exclusion_distance <- terra::distance(exclusion_binary)
    exclusion_buffer <- terra::ifel(
      exclusion_distance <= buffer_m,
      1,
      NA
    )
  } else {
    exclusion_buffer <- exclusion_binary
  }
  # candidate cells = valid eelgrass-domain cells not inside exclusion buffer
  candidate_rast <- terra::ifel(
    !is.na(domain_rast) & is.na(exclusion_buffer),
    1,
    NA,
    filename = filename,
    overwrite = overwrite
  )
  candidate_rast
}

sample_pseudoabsences_min_dist <- function(candidate_rast,
                                           n_pa,
                                           min_dist_m = 100,
                                           seed = 123,
                                           oversample_factor = 20,
                                           max_iter = 20) {
  set.seed(seed)
  
  candidate_cells <- which(!is.na(terra::values(candidate_rast)))
  candidate_n <- length(candidate_cells)
  
  if (candidate_n < n_pa) {
    stop(
      paste0(
        "Only ", candidate_n,
        " candidate cells available, but ", n_pa,
        " pseudoabsences requested."
      )
    )
  }
  
  selected_cells <- integer(0)
  selected_xy <- matrix(numeric(0), ncol = 2)
  
  remaining_cells <- candidate_cells
  
  for (iter in seq_len(max_iter)) {
    n_needed <- n_pa - length(selected_cells)
    
    if (n_needed <= 0) {
      break
    }
    
    sample_size <- min(
      length(remaining_cells),
      max(n_needed * oversample_factor, n_needed)
    )
    
    sampled_cells <- sample(remaining_cells, sample_size)
    sampled_xy <- terra::xyFromCell(candidate_rast, sampled_cells)
    
    sampled_order <- sample(seq_along(sampled_cells))
    
    for (j in sampled_order) {
      this_cell <- sampled_cells[j]
      this_xy <- sampled_xy[j, , drop = FALSE]
      
      if (length(selected_cells) == 0) {
        keep <- TRUE
      } else {
        d <- sqrt(
          (selected_xy[, 1] - this_xy[1, 1])^2 +
            (selected_xy[, 2] - this_xy[1, 2])^2
        )
        
        keep <- all(d >= min_dist_m)
      }
      
      if (keep) {
        selected_cells <- c(selected_cells, this_cell)
        selected_xy <- rbind(selected_xy, this_xy)
      }
      
      if (length(selected_cells) >= n_pa) {
        break
      }
    }
    
    remaining_cells <- setdiff(remaining_cells, selected_cells)
  }
  
  if (length(selected_cells) < n_pa) {
    stop(
      paste0(
        "Only selected ", length(selected_cells),
        " pseudoabsences with minimum distance of ", min_dist_m,
        " m, but ", n_pa,
        " were requested. Try reducing n_pa, reducing min_dist_m, ",
        "or increasing the candidate area."
      )
    )
  }
  
  selected_cells <- selected_cells[seq_len(n_pa)]
  selected_xy <- selected_xy[seq_len(n_pa), , drop = FALSE]
  
  pa_df <- data.frame(
    cell = selected_cells,
    x = selected_xy[, 1],
    y = selected_xy[, 2],
    obs = 0
  )
  
  pa_sf <- sf::st_as_sf(
    pa_df,
    coords = c("x", "y"),
    crs = terra::crs(candidate_rast)
  )
  
  list(
    pa_df = pa_df,
    pa_sf = pa_sf,
    candidate_n = candidate_n
  )
}

calc_field_metrics <- function(
    data,
    obs_col,
    model_cols = NULL,
    threshold = NULL,
    threshold_method = "max_tss",
    eps = 1e-15
){
  
  if (is.null(model_cols)) {
    model_cols <- names(data)[grepl("_pred$", names(data))]
  }
  
  dat <- data |>
    sf::st_drop_geometry() |>
    dplyr::mutate(
      obs = as.integer(.data[[obs_col]])
    )
  
  purrr::map_dfr(
    model_cols,
    function(model){
      
      df <- dat |>
        dplyr::select(obs, all_of(model)) |>
        dplyr::filter(
          !is.na(obs),
          !is.na(.data[[model]])
        ) |>
        dplyr::mutate(
          pred = pmin(
            pmax(as.numeric(.data[[model]]), eps),
            1 - eps
          )
        )
      
      if(length(unique(df$obs)) < 2){
        
        return(
          tibble::tibble(
            model = model,
            n = nrow(df),
            n_presence = sum(df$obs == 1),
            n_absence = sum(df$obs == 0),
            threshold = NA_real_,
            auc = NA_real_,
            tjur = NA_real_,
            brier = NA_real_,
            logloss = NA_real_,
            sensitivity = NA_real_,
            specificity = NA_real_,
            tss = NA_real_
          )
        )
        
      }
      
      auc_val <- yardstick::roc_auc_vec(
        truth = factor(df$obs, levels = c(0,1)),
        estimate = df$pred,
        event_level = "second"
      )
      
      tjur_val <-
        mean(df$pred[df$obs == 1]) -
        mean(df$pred[df$obs == 0])
      
      brier_val <-
        mean((df$pred - df$obs)^2)
      
      logloss_val <-
        -mean(
          df$obs * log(df$pred) +
            (1 - df$obs) * log(1 - df$pred)
        )
      
      if(is.null(threshold)){
        
        threshold_grid <- sort(unique(df$pred))
        
        threshold_stats <- purrr::map_dfr(
          threshold_grid,
          function(th){
            
            pred_class <- ifelse(df$pred >= th,1,0)
            
            tp <- sum(pred_class == 1 & df$obs == 1)
            tn <- sum(pred_class == 0 & df$obs == 0)
            fp <- sum(pred_class == 1 & df$obs == 0)
            fn <- sum(pred_class == 0 & df$obs == 1)
            
            sensitivity <-
              ifelse(tp + fn > 0,
                     tp/(tp+fn),
                     NA_real_)
            
            specificity <-
              ifelse(tn + fp > 0,
                     tn/(tn+fp),
                     NA_real_)
            
            tibble::tibble(
              threshold = th,
              sensitivity = sensitivity,
              specificity = specificity,
              tss = sensitivity + specificity - 1
            )
            
          }
        )
        
        best <- threshold_stats |>
          dplyr::filter(
            tss == max(tss, na.rm = TRUE)
          ) |>
          dplyr::arrange(
            desc(sensitivity),
            desc(specificity),
            threshold
          ) |>
          dplyr::slice(1)
        
      } else {
        
        th <- if(length(threshold) == 1)
          threshold
        else
          threshold[[model]]
        
        pred_class <- ifelse(df$pred >= th,1,0)
        
        tp <- sum(pred_class == 1 & df$obs == 1)
        tn <- sum(pred_class == 0 & df$obs == 0)
        fp <- sum(pred_class == 1 & df$obs == 0)
        fn <- sum(pred_class == 0 & df$obs == 1)
        
        best <- tibble::tibble(
          threshold = th,
          sensitivity = tp/(tp+fn),
          specificity = tn/(tn+fp),
          tss = tp/(tp+fn) + tn/(tn+fp) - 1
        )
        
      }
      
      tibble::tibble(
        model = model,
        n = nrow(df),
        n_presence = sum(df$obs == 1),
        n_absence = sum(df$obs == 0),
        threshold = best$threshold,
        auc = auc_val,
        tjur = tjur_val,
        brier = brier_val,
        logloss = logloss_val,
        sensitivity = best$sensitivity,
        specificity = best$specificity,
        tss = best$tss
      )
      
    }
  )
  
}

make_atm_plot <- function(var_name,
                          ylab,
                          ylim = NULL,
                          convert_kelvin = FALSE) {
  
  plot_df <- bind_rows(
    env_long %>%
      filter(variable == var_name) %>%
      mutate(source = "Hindcast"),
    
    obs_long %>%
      filter(variable == var_name) %>%
      mutate(source = "Survey")
  )
  
  if (convert_kelvin) {
    plot_df <- plot_df %>%
      mutate(value = value - 273.15)
  }
  
  p <- ggplot(
    plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Survey" = "#E69F00",
        "Hindcast" = "#56B4E9"
      )
    ) +
    labs(
      x = NULL,
      y = ylab,
      fill = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines")
    )
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  p
}

make_bccm_nep_plot <- function(bccm_var,
                               nep_var,
                               ylab,
                               ylim = NULL) {
  
  plot_df <- bind_rows(
    
    env_long %>%
      filter(variable %in% c(bccm_var, nep_var)) %>%
      mutate(
        source = case_when(
          variable == bccm_var ~ "Hindcast BCCM",
          variable == nep_var  ~ "Hindcast NEP36"
        )
      ),
    
    obs_long %>%
      filter(variable %in% c(bccm_var, nep_var)) %>%
      mutate(
        source = case_when(
          variable == bccm_var ~ "Survey BCCM",
          variable == nep_var  ~ "Survey NEP36"
        )
      )
  )
  
  p <- ggplot(
    plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Hindcast BCCM" = "#56B4E9",
        "Hindcast NEP36" = "#009E73",
        "Survey BCCM" = "#E69F00",
        "Survey NEP36" = "#CC79A7"
      )
    ) +
    labs(
      x = NULL,
      y = ylab,
      fill = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines")
    )
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  p
}