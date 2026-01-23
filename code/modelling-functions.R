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

# Calculate threshold-independent statistics using train and test data
evalStats <- function( folds, m, CV ){
  traintest.df <- data.frame()
  for(i in folds){
    # Get train obs
    train <- CV[[i]][["train"]]
    sp_data_cv <- filter(seagrass_data_long, species == sp)
    trainobs <- sp_data_cv[ train, 76]
    # Get train preds
    trainpred <- plogis(predict(m$models[[i]])$est[train])
    # Calculate area under the receiver-operator curve (AUC))
    train.AUC <-  ModelMetrics::auc(trainobs, trainpred)
    # Calculate tjur R2
    train.tjur <- tjur(trainobs, trainpred )
    # Get test indices
    test <- CV[[i]][["test"]]
    testobs <- sp_data_cv[ test, 76]
    #testobs <- tobs$presence
    # Get test preds
    testpred <- plogis(predict(m$models[[i]])$est[test])
    # Calculate area under the receiver-operator curve (AUC))
    test.AUC <- ModelMetrics::auc(testobs, testpred)
    # Calculate tjur R2
    test.tjur <- tjur(testobs, testpred)
    # sum log likelihood
    ll<- m$sum_loglik
    mae <- Metrics::mae(testobs, testpred)
    bias <- Metrics::bias(testobs, testpred)
        traintest.df <- data.frame(ll=ll, mae=mae, bias=bias, train.AUC = train.AUC, test.AUC=test.AUC, train.tjur = train.tjur, test.tjur = test.tjur, species = sp, fold = i)
    traintest.df <- traintest.df %>% dplyr::summarise(mean_AUC_train = mean(train.AUC, na.rm = TRUE), mean_AUC_test = mean(test.AUC, na.rm = TRUE), mean_Tjur_train = mean(train.tjur, na.rm = TRUE), mean_Tjur_test = mean(test.tjur, na.rm = TRUE), sum_loglike = mean(ll, na.rm = TRUE), mean_mae = mean(mae, na.rm =TRUE), mean_bias = mean(bias, na.rm = TRUE))
  }
  return(traintest.df)
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

# PredictSDM <- function(env, model, survey_type, species, years) {
#   message("Predicting with environmental layers...")
#   
#   outdir <- file.path("code/output_data/seagrass_predictions/survey", species)
#   if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#   
#   #predictions for each survey
#   predbysurvey <- function(s, ...) {
#     env$Survey <- as.factor(s)
#     pname <- paste0("Prediction_", s)
#     
#     #average across years
#     env_all <- data.frame()
#     for (y in years) {
#       env$Year_factor <- as.factor(y)
#       env_all <- rbind(env_all, env)
#     }
#     
#     hold_all <- predict(model, newdata = env_all)
#     sims <- predict(model, newdata = env_all, nsim = 100)
#     hold_all$SE <- apply(sims, 1, sd)
#     
#     hold <- hold_all %>%
#       group_by(X_m, Y_m, ID, Survey) %>%
#       summarise(across(everything(), mean, na.rm = TRUE)) %>%
#       arrange(ID) %>%
#       as.data.frame()
#     
#     epreds <- env %>%
#       select(X_m, Y_m, X, Y, ID, Survey) %>%
#       left_join(hold %>% select(ID, est:SE, Survey), by = c("ID", "Survey"))
#     
#     save(epreds, file = file.path(outdir, paste0(pname, "_survey_preds.RData")))
#     return(epreds)
#   }
#   
#   predlist <- lapply(survey_type, FUN = predbysurvey)
#   message("Combining and averaging predictions across survey types...")
#   
#   all_preds <- do.call(rbind, predlist)
#   
#   mean_preds <- all_preds %>%
#     group_by(X_m, Y_m, ID) %>%
#     summarise(across(est:SE, mean, na.rm = TRUE)) %>%
#     arrange(ID) %>%
#     as.data.frame()
#   
#   save(mean_preds, file = file.path(outdir, paste0("MeanSurveyPreds_", species, ".RData")))
#   
#   return(mean_preds)
# }

# PredictSDM <- function(env, model, survey_type, species, model_name) {
#   message("Predicting with environmental layers...")
#   
#   outdir <- file.path(
#     "code/output_data/seagrass_predictions/survey",
#     species,
#     model_name
#   )
#   if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#   
#   pred_list <- list()
#   
#   for (s in survey_type) {
#     message("  Survey: ", s)
#     
#     env$Survey <- as.factor(s)
#     
#     preds <- predict(model, newdata = env)
#     sims  <- predict(model, newdata = env) # model, newdata = env, nsim = 100)
#     preds$SE <- apply(sims, 1, sd)
#     
#     pred <- env %>%
#       select(X_m, Y_m, X, Y, ID, Survey) %>%
#       bind_cols(preds %>% select(est, SE))
#     
#     # Save survey-level predictions with model name
#     save(
#       pred,
#       file = file.path(
#         outdir,
#         paste0("Prediction_", species, "_", model_name, "_", s, ".RData")
#       )
#     )
#     
#     pred_list[[s]] <- pred
#   }
#   
#   message("Averaging predictions across survey types...")
#   
#   all_preds <- bind_rows(pred_list)
#   
#   mean_pred <- all_preds %>%
#     group_by(X_m, Y_m, ID) %>%
#     summarise(
#       across(c(est, SE), \(x) mean(x, na.rm = TRUE)),
#       .groups = "drop"
#     ) %>%
#     arrange(ID)
#   
#   # Save mean prediction with model name
#   save(
#     mean_pred,
#     file = file.path(
#       outdir,
#       paste0("MeanSurveyPreds_", species, "_", model_name, ".RData")
#     )
#   )
#   
#   return(mean_pred)
# }

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
      select(X_m, Y_m, X, Y, ID, region, Survey) %>%
      bind_cols(preds %>% select(est, SE))
    
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
    group_by(X_m, Y_m, ID) %>%
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