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
    #Setting up parameters for how my model is going to be fitted
    fitControl <- caret::trainControl(method = "cv",
                                      number = (numFolds-1),
                                      index = index_list)
    set.seed(2024)
    #performing model selection by glmStepAIC 
    caret_model <- CAST::ffs(response = sp_data$presence, 
                             predictors = sp_data[,6:37], 
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

varImp <- function( model, nosp_model, dat, preds, permute ) {
  # Get inputs
  obs <- dat[, "presence"]
  env <- dat[, c("X","Y", preds)]
  preds <- predict(model, type="response")
  AUC <- ModelMetrics::auc( obs, preds$est )
  vars <- names(env)[!names(env) %in% c("X","Y")]
  var_reps <- rep(vars, each=permute)
  var_reps <- c(var_reps, "spatial")
  # Plan
  plan(multisession)
  # # Set seed
  # set.seed(42)
  # Randomize the order of variable v and make predictions
  rand_mod <- function(v, rand_j, obs){
    # require
    require(sdmTMB)
    require(ModelMetrics)
    # predict
    # randomize var
    rand_j[v] <- sample( rand_j[[v]] )
    # Predict with rand_j
    p <- predict(model, newdata = rand_j, type="response")
    # Calculate AUC
    AUC <- ModelMetrics::auc(obs, p$est)
    # Return AUC
    return(AUC)
  }
  # Get AUC values from preds form randomized env dataset
  aucs <- unlist(future_lapply( var_reps, rand_mod, rand_j=env, obs=obs, future.seed=42 ))
  names(aucs) <- var_reps
  # Mean auc values across permutations
  permuted_auc <- aggregate(aucs, by=list(names(aucs)), FUN=mean)
  # Difference between original model AUC and randomized AUC
  # Return 0 if negative (pmax == 0)
  perm_imp <- pmax(0, (AUC - permuted_auc$x))
  rel_imp <- 100 * perm_imp / sum(perm_imp) # as percentage
  rel_imp <- round(rel_imp, 1)
  # Return dataframe
  relauc <- data.frame(terms = permuted_auc[,1], relimp = rel_imp,
                       stringsAsFactors = FALSE)
  # Return
  return( relauc )
}

# Calculate threshold-independent statistics using train and test data
evalStats <- function( folds, m, CV ){
  traintest.df <- data.frame()
  for(i in folds){
    # Get train obs
    train <- CV[[i]][["train"]]
    sp_data_cv <- filter(seagrass_data_long, species == sp)
    trainobs <- sp_data_cv[ train, 42]
    # Get train preds
    trainpred <- plogis(predict(m$models[[i]])$est[train])
    # Calculate area under the receiver-operator curve (AUC))
    train.AUC <-  ModelMetrics::auc(trainobs, trainpred)
    # Calculate tjur R2
    train.tjur <- tjur(trainobs, trainpred )
    # Get test indices
    test <- CV[[i]][["test"]]
    testobs <- sp_data_cv[ test, 42]
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
