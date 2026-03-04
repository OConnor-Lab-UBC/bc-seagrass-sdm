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

# tjur <- function(y, pred){
#   mean(pred[y == 1], na.rm = TRUE) - mean(pred[y == 0], na.rm = TRUE)
# }


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
  thrs <- seq(0, 1, by = 0.01)
  tss_vals <- sapply(thrs, function(t) tss_fun(obs, pred, t))
  thrs[which.max(tss_vals)]
}

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

tss_metric <- function(obs, pred, threshold) {
  pred_bin <- as.numeric(pred >= threshold)
  cm <- table(factor(obs, levels = c(0,1)),
              factor(pred_bin, levels = c(0,1)))
  
  sens <- ifelse(sum(cm["1",]) > 0, cm["1","1"]/sum(cm["1",]), NA)
  spec <- ifelse(sum(cm["0",]) > 0, cm["0","0"]/sum(cm["0",]), NA)
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


varImp_sdmTMB <- function(model, dat, preds, groups = NULL, permute = 10) {
  library(sdmTMB)
  library(ModelMetrics)
  # ------------------------------
  # Checks
  # ------------------------------
  if(!all(preds %in% names(dat))) {
    stop("Some predictors are missing from the dataset.")
  }
  
  obs <- dat$presence
  
  # ------------------------------
  # Handle spatiotemporal models ONLY if needed
  # ------------------------------
  if(!is.null(model$spatiotemporal) && model$spatiotemporal != "off") {
    time_var <- model$time
    if(is.null(time_var) || !time_var %in% names(dat)) {
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
  
  # ======================================================
  # 1️⃣ INDIVIDUAL VARIABLE IMPORTANCE
  # ======================================================
  
  aucs_ind <- vector("list", length(preds))
  names(aucs_ind) <- preds
  
  for(v in preds) {
    perm_results <- numeric(permute)
    for(i in seq_len(permute)) {
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
  rel_ind <- round(100 * imp_ind / sum(imp_ind), 1)
  
  ind_df <- data.frame(
    term = preds,
    relimp = rel_ind,
    stringsAsFactors = FALSE
  )
  
  # ======================================================
  # 2️⃣ GROUPED VARIABLE IMPORTANCE (optional)
  # ======================================================
  
  grp_df <- NULL
  
  if(!is.null(groups)) {
    
    # groups must be a named list
    if(is.null(names(groups))) {
      stop("Groups must be a *named* list.")
    }
    
    aucs_grp <- vector("list", length(groups))
    names(aucs_grp) <- names(groups)
    
    for(g in names(groups)) {
      vars <- groups[[g]]
      
      if(!all(vars %in% preds)) {
        stop(paste("Group", g, "contains variables not in preds."))
      }
      
      perm_results <- numeric(permute)
      
      for(i in seq_len(permute)) {
        dat_perm <- dat
        
        # Permute ALL variables in the group
        for(v in vars) {
          dat_perm[[v]] <- sample(dat_perm[[v]])
        }
        
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
      
      aucs_grp[[g]] <- perm_results
    }
    
    perm_auc_grp <- sapply(aucs_grp, mean, na.rm = TRUE)
    imp_grp <- pmax(0, base_auc - perm_auc_grp)
    rel_grp <- round(100 * imp_grp / sum(imp_grp), 1)
    
    grp_df <- data.frame(
      group = names(groups),
      relimp = rel_grp,
      stringsAsFactors = FALSE
    )
  }
  
  # ------------------------------
  # Return both
  # ------------------------------
  list(
    base_auc = base_auc,
    individual = ind_df,
    grouped = grp_df
  )
}



evalStats <- function(folds, m, CV, response_col) {
  traintest.df <- list()
  
  for(i in folds){
    # Train/test indices
    train_idx <- CV[[i]][["train"]]
    test_idx  <- CV[[i]][["test"]]
    
    # Observed values
    obs_train <- m$data[[response_col]][train_idx]
    obs_test  <- m$data[[response_col]][test_idx]
    
    # Predictions
    pred_train <- plogis(predict(m$models[[i]])$est[train_idx])
    pred_test  <- plogis(predict(m$models[[i]])$est[test_idx])
    
    # Calculate TSS threshold from training data only
    threshold <- get_optimal_threshold(obs_train, pred_train)
    
    # Train metrics
    train_auc  <- ModelMetrics::auc(obs_train, pred_train)
    train_tjur <- tjur(obs_train, pred_train)
    train_rmse <- sqrt(mean((obs_train - pred_train)^2))
    train_tss  <- tss_metric(obs_train, pred_train, threshold)
    
    # Test metrics
    test_auc  <- ModelMetrics::auc(obs_test, pred_test)
    test_tjur <- tjur(obs_test, pred_test)
    test_rmse <- sqrt(mean((obs_test - pred_test)^2))
    test_tss  <- tss_metric(obs_test, pred_test, threshold)
    
    # Sum log-likelihood (from full model)
    sum_loglikelihood <- m$sum_loglik
    
    # Save fold metrics including threshold
    traintest.df[[i]] <- data.frame(
      fold = i,
      sum_loglikelihood = sum_loglikelihood,
      train_auc  = train_auc,
      train_tjur = train_tjur,
      train_rmse = train_rmse,
      train_tss  = train_tss,
      test_auc   = test_auc,
      test_tjur  = test_tjur,
      test_rmse  = test_rmse,
      test_tss   = test_tss,
      threshold  = threshold
    )
  }
  
  # Combine folds and calculate mean metrics
  traintest.df <- do.call(rbind, traintest.df)
  mean_metrics <- traintest.df %>%
    dplyr::summarise(
      mean_sum_loglikelihood = mean(sum_loglikelihood, na.rm = TRUE),
      mean_train_auc  = mean(train_auc, na.rm = TRUE),
      mean_train_tjur = mean(train_tjur, na.rm = TRUE),
      mean_train_rmse = mean(train_rmse, na.rm = TRUE),
      mean_train_tss  = mean(train_tss, na.rm = TRUE),
      mean_test_auc   = mean(test_auc, na.rm = TRUE),
      mean_test_tjur  = mean(test_tjur, na.rm = TRUE),
      mean_test_rmse  = mean(test_rmse, na.rm = TRUE),
      mean_test_tss   = mean(test_tss, na.rm = TRUE),
      mean_threshold  = mean(threshold, na.rm = TRUE)
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

# evalulate forecast 
evaluate_forecast <- function(obs_train,
                              pred_train,
                              obs_test,
                              pred_test,
                              threshold = NULL) {
  
  # Coerce to numeric
  obs_train <- as.numeric(obs_train)
  obs_test  <- as.numeric(obs_test)
  pred_train <- as.numeric(pred_train)
  pred_test  <- as.numeric(pred_test)
  
  # Use provided threshold or fallback
  if(is.null(threshold)) {
    threshold <- max(0.01, mean(obs_train))  # fallback for rare species
  }
  
  # Binary predictions
  pred_train_bin <- as.numeric(pred_train >= threshold)
  pred_test_bin  <- as.numeric(pred_test  >= threshold)
  
  # Confusion matrices
  cm_train <- table(factor(obs_train, levels = c(0,1)),
                    factor(pred_train_bin, levels = c(0,1)))
  cm_test  <- table(factor(obs_test, levels = c(0,1)),
                    factor(pred_test_bin, levels = c(0,1)))
  
  # TSS
  sens_train <- ifelse(sum(cm_train["1",]) > 0, cm_train["1","1"]/sum(cm_train["1",]), NA)
  spec_train <- ifelse(sum(cm_train["0",]) > 0, cm_train["0","0"]/sum(cm_train["0",]), NA)
  TSS_train  <- sens_train + spec_train - 1
  
  sens_test  <- ifelse(sum(cm_test["1",]) > 0, cm_test["1","1"]/sum(cm_test["1",]), NA)
  spec_test  <- ifelse(sum(cm_test["0",]) > 0, cm_test["0","0"]/sum(cm_test["0",]), NA)
  TSS_test   <- sens_test + spec_test - 1
  
  # Other metrics
  Tjur_train <- mean(pred_train[obs_train==1]) - mean(pred_train[obs_train==0])
  Tjur_test  <- mean(pred_test[obs_test==1]) - mean(pred_test[obs_test==0])
  RMSE_train <- sqrt(mean((obs_train - pred_train)^2))
  RMSE_test  <- sqrt(mean((obs_test - pred_test)^2))
  AUC_train  <- ModelMetrics::auc(obs_train, pred_train)
  AUC_test   <- ModelMetrics::auc(obs_test, pred_test)
  
  # Return
  data.frame(
    dataset = c("training","forecast"),
    AUC = c(AUC_train, AUC_test),
    TjurR2 = c(Tjur_train, Tjur_test),
    RMSE = c(RMSE_train, RMSE_test),
    TSS = c(TSS_train, TSS_test),
    threshold_used = threshold
  )
}

evalfmod <- function( x, thresh ){
  eval.df <- data.frame()
  obspred <- data.frame( PlotID=1:nrow(x),
                         Observed=x[,'presence'],
                         Predicted=x[,'fitted_vals'] )
  # Get confusion matrix based on TSS
  cmx_tss<- PresenceAbsence::cmx(DATA = obspred, which.model = 1, thresh = thresh$Predicted[thresh$Method == "MaxSens+Spec"])
  true_neg <- cmx_tss[1, 1]
  false_neg <- cmx_tss[1, 2]
  false_pos <- cmx_tss[2, 1]
  true_pos <- cmx_tss[2, 2]
  true_pos_rate <- true_pos / (true_pos + false_neg)
  true_neg_rate <- true_neg / (true_neg + false_pos)
  TSS <- true_pos_rate + true_neg_rate - 1
  # Get confusion matrix based on Kappa
  cmx_kappa<- PresenceAbsence::cmx(DATA = obspred, which.model = 1, thresh = thresh$Predicted[thresh$Method == "MaxKappa"])
  kappa <- PresenceAbsence::Kappa(CMX = cmx_kappa, st.dev = TRUE)
  miller<- modEvA::MillerCalib(model = NULL, obs = x$presence, pred = x$fitted_vals)
  eer<- modEvA::errorMeasures(model = NULL, obs = x$presence, pred = x$fitted_vals)
  hlgof_quant<- modEvA::HLfit(model = NULL, obs = x$presence, pred = x$fitted_vals, bin.method = "quantiles", n.bins = 3000) # these values are fine
  hlgof_prob<- modEvA::HLfit(model = NULL, obs = x$presence, pred = x$fitted_vals, bin.method = "prob.bins") # these values are not great but the highest probs don't have many values
  hlgof_nbin<- modEvA::HLfit(model = NULL, obs = x$presence, pred = x$fitted_vals, bin.method = "n.bins", n.bins = 8) # these values are fine
  eval.df <- data.frame(kappa=kappa, TSS=TSS, miller = miller, eer=eer, hlgof_quant$chi.sq, hlgof_quant$p.value, hlgof_quant$RMSE,  hlgof_prob$chi.sq, hlgof_prob$p.value, hlgof_prob$RMSE, hlgof_nbin$chi.sq, hlgof_nbin$p.value, hlgof_nbin$RMSE, species = sp)
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


#### biomod 2 functions 
###############################################################################
# BIOMOD2 Cross-Validation Wrapper
###############################################################################
run_biomod_cv <- function(data, config, ml_models) {
  
  dat <- data %>%
    dplyr::select(presence, X_m, Y_m, fold, Year, all_of(config$pred_vars)) %>%
    na.omit() %>%
    mutate(presence = as.numeric(presence))
  
  myBiomodData <- BIOMOD_FormatingData(
    resp.var = dat$presence,
    expl.var = dat %>% dplyr::select(all_of(config$pred_vars)),
    resp.xy = dat[, c("X_m", "Y_m")],
    resp.name = config$resp_name,
    PA.nb.rep = 1,
    PA.nb.absences = 0,
    PA.strategy = "none"
  )
  
  folds <- sort(unique(dat$fold))
  K <- length(folds)
  
  cv_table <- matrix(TRUE, nrow = nrow(dat), ncol = K + 1)
  for (k in seq_along(folds)) {
    cv_table[dat$fold == folds[k], k] <- FALSE
  }
  colnames(cv_table) <- c(paste0("_allData_RUN", seq_len(K)), "_allData_allRun")
  
  myOptions <- bm_ModelingOptions(
    data.type = "binary",
    models = ml_models,
    strategy = "bigboss"
  )
  
  myBiomodModelOut <- BIOMOD_Modeling(
    bm.format = myBiomodData,
    models = ml_models,
    CV.strategy = "user.defined",
    CV.user.table = cv_table,
    CV.do.full.models = TRUE,
    metric.eval = c("AUCroc", "TSS", "KAPPA"),
    var.import = 3,
    seed.val = 123
  )
  
  evals <- get_evaluations(myBiomodModelOut)
  
  tss_thresholds <- evals %>%
    filter(metric.eval == "TSS") %>%
    dplyr::select(algo, run, cutoff, validation) %>%
    group_by(algo) %>%
    summarise(
      threshold_mean = mean(cutoff, na.rm = TRUE),
      threshold_sd = sd(cutoff, na.rm = TRUE),
      TSS_validation = mean(validation, na.rm = TRUE)
    )
  
  auc_metrics <- bm_PlotEvalMean(
    bm.out = myBiomodModelOut,
    metric.eval = c("AUCroc", "TSS"),
    dataset = "validation"
  )$tab %>%
    as.data.frame() %>%
    setNames(c("algo", "AUC_mean", "TSS_mean", "AUC_sd", "TSS_sd"))
  
  var_imp <- get_variables_importance(myBiomodModelOut)
  
  return(list(
    data = dat,
    biomod_out = myBiomodModelOut,
    tss_thresholds = tss_thresholds,
    auc_metrics = auc_metrics,
    var_importance = var_imp
  ))
}
###############################################################################
# Independent GBM Cross-Validation (correct threshold)
###############################################################################
run_gbm_cv <- function(dat, predictors, response = "presence", fold_col = "fold") {
  
  folds <- sort(unique(dat[[fold_col]]))
  out <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    f <- folds[i]
    train <- dat[dat[[fold_col]] != f, ]
    test <- dat[dat[[fold_col]] == f, ]
    
    gbm_fit <- gbm(
      formula = as.formula(paste(response, "~", paste(predictors, collapse = " + "))),
      data = train[, c(response, predictors)],
      distribution = "bernoulli",
      n.trees = 5000,
      interaction.depth = 3,
      shrinkage = 0.005,
      bag.fraction = 0.5,
      n.minobsinnode = 10,
      train.fraction = 1,
      verbose = FALSE
    )
    
    p <- predict(gbm_fit, newdata = test[, predictors], n.trees = 5000, type = "response")
    y <- test[[response]]
    
    thr_fold <- get_tss_threshold(y, p)
    tss_fold <- tss_fun(y, p, thr_fold)
    
    out[[i]] <- data.frame(
      fold = f,
      threshold = thr_fold,
      RMSE = rmse_fun(y, p),
      Tjur = tjur_fun(y, p),
      AUC = auc_fun(y, p),
      TSS = tss_fold
    )
  }
  
  result <- bind_rows(out) %>%
    summarise(
      mean_threshold = mean(threshold),
      sd_threshold = sd(threshold),
      mean_RMSE = mean(RMSE),
      sd_RMSE = sd(RMSE),
      mean_Tjur = mean(Tjur),
      sd_Tjur = sd(Tjur),
      mean_AUC = mean(AUC),
      sd_AUC = sd(AUC),
      mean_TSS = mean(TSS),
      sd_TSS = sd(TSS)
    )
  
  result$algo <- "GBM_robust"
  return(result)
}
###############################################################################
# Temporal Forecasting (GBM uses the independent threshold)
###############################################################################
run_temporal_forecast <- function(dat, pred_vars, config,
                                  biomod_thresholds,
                                  gbm_threshold) {
  
  data_pre2010 <- dat %>% filter(Year < 2010)
  data_post2012 <- dat %>% filter(Year > 2012)
  
  if (nrow(data_pre2010) == 0 || nrow(data_post2012) == 0)
    return(NULL)
  
  models <- list()
  
  # RF
  models$RF <- randomForest(
    x = data_pre2010[, pred_vars],
    y = as.factor(data_pre2010$presence),
    ntree = 500
  )
  
  # GBM
  models$GBM <- gbm(
    formula = presence ~ .,
    distribution = "bernoulli",
    data = data_pre2010[, c("presence", pred_vars)],
    n.trees = 5000,
    interaction.depth = 3,
    shrinkage = 0.005,
    verbose = FALSE
  )
  
  # XGBOOST
  xgb_matrix_pre <- model.matrix(~ . - 1, data = data_pre2010[, pred_vars])
  xgb_matrix_post <- model.matrix(~ . - 1, data = data_post2012[, pred_vars])
  dtrain <- xgb.DMatrix(data = xgb_matrix_pre, label = data_pre2010$presence)
  models$XGBOOST <- xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 6,
      eta = 0.1,
      nthread = 4
    ),
    data = dtrain,
    nrounds = 200,
    verbose = 0
  )
  
  # ANN
  data_pre_num <- data_pre2010[, pred_vars] %>%
    mutate(across(everything(), as.numeric)) %>%
    replace(is.na(.), 0)
  
  data_post_num <- data_post2012[, pred_vars] %>%
    mutate(across(everything(), as.numeric)) %>%
    replace(is.na(.), 0)
  
  y_pre_1hot <- class.ind(as.factor(data_pre2010$presence))
  models$ANN <- nnet(
    x = as.matrix(data_pre_num),
    y = y_pre_1hot,
    size = 5,
    linout = FALSE,
    entropy = TRUE,
    maxit = 500,
    trace = FALSE
  )
  
  # CTA
  models$CTA <- rpart(
    formula = presence ~ .,
    data = cbind(presence = as.factor(data_pre2010$presence), data_pre_num),
    method = "class",
    control = rpart.control(minsplit = 5, cp = 0.001)
  )
  
  # Predictions
  predictions <- list(
    RF = predict(models$RF, data_post2012[, pred_vars], type = "prob")[, 2],
    GBM = predict(models$GBM, data_post2012[, pred_vars], n.trees = 5000, type = "response"),
    XGBOOST = predict(models$XGBOOST, newdata = xgb_matrix_post),
    ANN = predict(models$ANN, newdata = as.matrix(data_post_num))[, 2],
    CTA = predict(models$CTA, newdata = data_post_num, type = "prob")[, 2]
  )
  
  obs_fore <- data_post2012$presence
  
  out <- list()
  
  for (mod in names(predictions)) {
    pred <- predictions[[mod]]
    
    # BIOMOD for non‑GBM models
    if (mod != "GBM") {
      threshold <- biomod_thresholds$threshold_mean[biomod_thresholds$algo == mod]
      
      # fallback threshold
      if (length(threshold) == 0 || is.na(threshold)) {
        threshold <- 0.5
      }
      
    } else {
      # GBM gets the correct independent threshold
      threshold <- gbm_threshold
    }
    
    out[[mod]] <- list(
      AUC = ModelMetrics::auc(obs_fore, pred),
      Tjur = tjur_fun(obs_fore, pred),
      RMSE = rmse_fun(obs_fore, pred),
      TSS = tss_fun(obs_fore, pred, threshold),
      threshold = threshold,
      pred = pred,
      obs = obs_fore
    )
  }
  
  return(out)
}