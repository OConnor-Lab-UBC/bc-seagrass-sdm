###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# fit sdm models with biomod2 
###############################################################################


library(biomod2)
library(dplyr)
library(terra)

load("code/output_data/seagrass_model_inputs.RData")
source("code/modelling-functions.R")

seagrass_data_long <- seagrass_data_long %>% dplyr::select(-saltmean_bccm_sq_stnd, -saltmean_nep_sq_stnd, -slope_sqrt_stnd, -saltmin_bccm_sq_stnd, -saltmin_nep_sq_stnd)
seagrass_data_long <- seagrass_data_long %>%
  mutate(Survey = as.factor(substr(HKey, 1, 3)),
         HKey = as.factor(HKey),
         Year_factor = as.factor(Year))

####Eelgrass model with BCCM####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
data <- data %>% mutate(Survey = as.factor(Survey)) # was still not recognizing i changed survye to factor so did it again

pred_vars <- c("depth_stnd", "slope_stnd", "rei_stnd", "substrate", "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd", "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempmin_bccm_stnd")

dat2 <- data %>%
  dplyr::select(presence, X_m, Y_m, fold, Year, all_of(pred_vars)) %>%
  na.omit()

dat2$presence <- as.numeric(dat2$presence)

myBiomodData <- BIOMOD_FormatingData(
  resp.var  = dat2$presence,
  expl.var  = dat2 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = dat2[, c("X_m", "Y_m")],
  resp.name = "eelgrass",
  #PA.nb.rep = 0,
  PA.nb.rep = 1, # <-- REQUIRED for CV to work 
  PA.nb.absences = 0, # <-- keeps your true absences 
  PA.strategy = "none") # <-- prevents BIOMOD from generating fake absences


folds <- sort(unique(dat2$fold))
K <- length(folds)

cv_table <- matrix(TRUE, nrow = nrow(dat2), ncol = K + 1)  # extra column for full model
for (k in seq_along(folds)) {
  test_fold <- folds[k]
  cv_table[dat2$fold == test_fold, k] <- FALSE  # use your fold column
}
colnames(cv_table) <- c(paste0("_allData_RUN", seq_len(K)), "_allData_allRun")

ml_models <- c("GBM", "RF", "XGBOOST", "ANN", "CTA")

myOptions <- bm_ModelingOptions(
  data.type = "binary",
  models = ml_models,
  strategy = "bigboss"  # stronger defaults than the "default" set
)

# need to change var import to 10-20 for final runs
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

# Extract validation thresholds that maximized TSS
tss_thresholds_df <- evals %>%
  filter(metric.eval == "TSS") %>%
  select(algo, run, cutoff, validation)

print(tss_thresholds_df)

biomod_threshold_summary <- tss_thresholds_df %>%
  group_by(algo) %>%
  summarise(
    TSS_threshold_mean = mean(cutoff, na.rm = TRUE),
    TSS_threshold_sd   = sd(cutoff, na.rm = TRUE),
    mean_TSS_validation = mean(validation, na.rm = TRUE)
  )

print(biomod_threshold_summary)


plot_output <- bm_PlotEvalMean(bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = 'validation') 
ggdat <- plot_output$tab
names(ggdat) <- c("algo", "AUC_mean","TSS_mean","AUC_sd","TSS_sd")
ggdat$algo <- as.factor(ggdat$algo)



get_variables_importance(myBiomodModelOut)



bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# random forest here can get cut out, TSS is horrible
# none of the ML models do better than the sdmTMB models, BUT becasue the GBM bccm model forecast is comparable want to add without biomod
#the RMSE and Tjur cross validation metrics, which is not possible to do through biomod after much testing
auc_tab <- bm_PlotEvalMean( bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = "validation" )$tab
cv_tbl <- auc_tab %>% rename( algo = name, mean_AUC_validation = mean1, mean_TSS_validation = mean2, sd_AUC_validation = sd1, sd_TSS_validation = sd2 )

save(cv_tbl, file = "code/output_data/model_results/eelgrass_bccm_cv_metrics_biomod2.RData")

# GBM CV without biomod to get Tjur and RMSE and also becasue default values are dumn
library(gbm)
library(pROC)

tjur_fun <- function(y, p) {
  mean(p[y == 1], na.rm = TRUE) - mean(p[y == 0], na.rm = TRUE)
}

rmse_fun <- function(y, p) {
  sqrt(mean((y - p)^2, na.rm = TRUE))
}

auc_fun <- function(y, p) {
  as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
}

tss_fun <- function(y, p) {
  roc_obj <- pROC::roc(y, p, quiet = TRUE)
  # threshold that maximizes sensitivity + specificity
  coords_obj <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  tss <- coords_obj["sensitivity"] + coords_obj["specificity"] - 1
  as.numeric(tss)
}

gbm_cv_metrics <- function(dat, predictors, response = "presence", fold_col = "fold") {
  
  folds <- sort(unique(dat[[fold_col]]))
  out <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    f <- folds[i]
    
    train <- dat[dat[[fold_col]] != f, ]
    test  <- dat[dat[[fold_col]] == f, ]
    
    gbm_fit <- gbm(
      formula = as.formula(
        paste(response, "~", paste(predictors, collapse = " + "))
      ),
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
    
    out[[i]] <- data.frame(
      algo = "GBM_robust",
      fold = f,
      RMSE = rmse_fun(y, p),
      Tjur = tjur_fun(y, p),
      AUC  = auc_fun(y, p),
      TSS  = tss_fun(y, p)
    )
  }
  
  bind_rows(out)
}

gbm_cv_res <- gbm_cv_metrics(dat2, predictors = pred_vars)

gbm_cv_res          # per-fold metrics
gbm_summary <- gbm_cv_res %>%
  summarise(
    algo = "GBM_robust",
    mean_RMSE = mean(RMSE),
    mean_Tjur = mean(Tjur),
    mean_AUC  = mean(AUC),
    mean_TSS  = mean(TSS)
  )

gbm_summary

save(gbm_summary, file = "code/output_data/model_results/eelgrass_bccm_cv_metricsextraGBM_biomod2.RData")

# test forecasting
data_pre2013 <- dat2 %>% filter(Year < 2010)
data_post2012 <- dat2 %>% filter(Year > 2012)

# Format data for BIOMOD2
myBiomodData_pre2013 <- BIOMOD_FormatingData(
  resp.var  = data_pre2013$presence,
  expl.var  = data_pre2013 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = data_pre2013[, c("X_m", "Y_m")],
  resp.name = "eelgrass",
  PA.nb.rep = 0
)

n <- nrow(data_pre2013)
cv_table <- matrix(TRUE, nrow = n, ncol = 1)  # 1 column = 1 run
colnames(cv_table) <- "_allData_allRun"       # required by BIOMOD2


myBiomodModelOut_pre2013 <- BIOMOD_Modeling(
  bm.format       = myBiomodData_pre2013,
  models          = ml_models,
  CV.strategy     = "user.defined",
  CV.user.table   = cv_table,
  CV.do.full.models = TRUE,  # ensures BIOMOD2 builds the “full model”
  metric.eval     = c("AUCroc","TSS","KAPPA"),
  var.import      = 3,
  seed.val        = 123,
  modeling.id     = "pre2013_ML"
)


library(randomForest)
library(xgboost)
library(gbm)
library(ModelMetrics)  
library(nnet)       
library(MASS) 
library(rpart) 


# Random Forest
rf_mod <- randomForest(
  x = data_pre2013[, pred_vars],
  y = as.factor(data_pre2013$presence),
  ntree = 500
)

# GBM
gbm_mod <- gbm(
  formula = presence ~ .,
  distribution = "bernoulli",
  data = data_pre2013[, c("presence", pred_vars)],
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.005,
  verbose = FALSE
)

# XGBoost
xgb_matrix_pre <- model.matrix(~ . -1, data = data_pre2013[, pred_vars])
xgb_matrix_post <- model.matrix(~ . -1, data = data_post2012[, pred_vars])

dtrain <- xgb.DMatrix(
  data = xgb_matrix_pre,
  label = data_pre2013$presence
)

xgb_params <- list(
  objective = "binary:logistic",  # presence/absence
  eval_metric = "auc",
  max_depth = 6,
  eta = 0.1,
  nthread = 4
)

xgb_mod <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 200,
  verbose = 1
)

# --- ANN ---

data_pre2013_num <- data_pre2013 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.)))) %>%
  replace(is.na(.), 0)  
data_post2012_num <- data_post2012 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.))))%>%
  replace(is.na(.), 0)  
#data_pre2013_num$substrate <- as.factor(data_pre2013$substrate)

y_pre_1hot <- class.ind(as.factor(data_pre2013$presence))  # one-hot encoding: 0 -> [1,0], 1 -> [0,1]

# -------------------------------
# 3️⃣ Train ANN
# -------------------------------
set.seed(123)
ann_mod <- nnet(
  x = as.matrix(data_pre2013_num),
  y = y_pre_1hot,
  size = 5,       # hidden neurons
  linout = FALSE,
  entropy = TRUE,
  maxit = 500,
  trace = FALSE
)



# --- CTA (Classification Tree Analysis) ---
cta_mod <- rpart(
  formula = presence ~ .,
  data = cbind(presence = as.factor(data_pre2013$presence), data_pre2013_num),
  method = "class",
  control = rpart.control(minsplit = 5, cp = 0.001)
)




predictions <- list(
  RF = list(
    training = predict(rf_mod, data_pre2013[, pred_vars], type = "prob")[,2],
    forecast = predict(rf_mod, data_post2012[, pred_vars], type = "prob")[,2]
  ),
  GBM = list(
    training = predict(gbm_mod, data_pre2013[, pred_vars], n.trees = 2000, type = "response"),
    forecast = predict(gbm_mod, data_post2012[, pred_vars], n.trees = 2000, type = "response")
  ),
  XGBOOST = list(
    training = predict(xgb_mod, newdata = xgb_matrix_pre),
    forecast = predict(xgb_mod, newdata = xgb_matrix_post)
  ),
  ANN = list(
    training = predict(ann_mod, newdata = as.matrix(data_pre2013_num))[,2],
    forecast = predict(ann_mod, newdata = as.matrix(data_post2012_num))[,2]
  ),
  CTA = list(
    training = predict(cta_mod, newdata = data_pre2013_num, type = "prob")[,2],
    forecast = predict(cta_mod, newdata = data_post2012_num, type = "prob")[,2]
  )
)


forecast_results <- lapply(names(predictions), function(mod) {
  
  obs_train <- data_pre2013$presence
  pred_train <- predictions[[mod]]$training
  
  obs_fore  <- data_post2012$presence
  pred_fore <- predictions[[mod]]$forecast
  
  # Model-specific threshold from pre-2013 data
  threshold_mod <- get_optimal_threshold(obs_train, pred_train)
  
  data.frame(
    model = mod,
    AUC_forecast  = ModelMetrics::auc(obs_fore, pred_fore),
    Tjur_forecast = tjur(obs_fore, pred_fore),
    RMSE_forecast = rmse(obs_fore, pred_fore),
    TSS_forecast  = tss_metric(obs_fore, pred_fore, threshold_mod),
    threshold_used = threshold_mod
  )
}) %>% bind_rows()

print(forecast_results)
save(forecast_results, file = "code/output_data/model_results/eelgrass_bccm_forecast_metrics_biomod2.RData")
# only the GBM model ML models is comparable to the sdm tmb model 


#### make projection of GBM model only
#this was not working in biomod!
load("code/output_data/prediction_model_inputs.RData")
env<- env_20m_all %>% dplyr::select(depth_stnd, slope_stnd, rei_stnd, substrate, airtempmin_stnd, rsdsmin_stnd, prmin_stnd, 
                             saltcv_bccm_stnd, NH4_bccm_stnd, tempmin_bccm_stnd)

set.seed(123)

# make a full gbm model without cv
gbm_fit <- gbm(
  presence ~ .,
  data = dat2[, c("presence", pred_vars)],
  distribution = "bernoulli",
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.005,
  bag.fraction = 0.5,
  n.minobsinnode = 10,
  train.fraction = 1,
  verbose = FALSE
)

gbm_prob <- predict(
  gbm_fit,
  newdata = env,
  n.trees = 5000,
  type = "response"   # returns probabilities for bernoulli 
)

range(gbm_prob)


env_20m_all$gbm_prob <- gbm_prob

outdir <- file.path("./raster/eelgrass/gbm_bccm")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_gbm_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_gbm_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_gbm_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_gbm_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_gbm_bccm.tif")), overwrite = TRUE)





####Eelgrass model with NEP####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
data <- data %>% mutate(Survey = as.factor(Survey)) # was still not recognizing i changed survye to factor so did it again

pred_vars <- c("depth_stnd", "slope_stnd", "rei_stnd", "substrate", "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd", "saltcv_nep_stnd", "tempcv_nep_stnd")

dat2 <- data %>%
  dplyr::select(presence, X_m, Y_m, fold, Year, all_of(pred_vars)) %>%
  na.omit()

dat2$presence <- as.numeric(dat2$presence)

myBiomodData <- BIOMOD_FormatingData(
  resp.var  = dat2$presence,
  expl.var  = dat2 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = dat2[, c("X_m", "Y_m")],
  resp.name = "eelgrass",
  #PA.nb.rep = 0,
  PA.nb.rep = 1, # <-- REQUIRED for CV to work 
  PA.nb.absences = 0, # <-- keeps your true absences 
  PA.strategy = "none") # <-- prevents BIOMOD from generating fake absences


folds <- sort(unique(dat2$fold))
K <- length(folds)

cv_table <- matrix(TRUE, nrow = nrow(dat2), ncol = K + 1)  # extra column for full model
for (k in seq_along(folds)) {
  test_fold <- folds[k]
  cv_table[dat2$fold == test_fold, k] <- FALSE  # use your fold column
}
colnames(cv_table) <- c(paste0("_allData_RUN", seq_len(K)), "_allData_allRun")

ml_models <- c("GBM", "RF", "XGBOOST", "ANN", "CTA")

myOptions <- bm_ModelingOptions(
  data.type = "binary",
  models = ml_models,
  strategy = "bigboss"  # stronger defaults than the "default" set
)

# need to change var import to 10-20 for final runs
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

# Extract validation thresholds that maximized TSS
tss_thresholds_df <- evals %>%
  filter(metric.eval == "TSS") %>%
  dplyr::select(algo, run, cutoff, validation)

print(tss_thresholds_df)

biomod_threshold_summary <- tss_thresholds_df %>%
  group_by(algo) %>%
  summarise(
    TSS_threshold_mean = mean(cutoff, na.rm = TRUE),
    TSS_threshold_sd   = sd(cutoff, na.rm = TRUE),
    mean_TSS_validation = mean(validation, na.rm = TRUE)
  )

print(biomod_threshold_summary)


plot_output <- bm_PlotEvalMean(bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = 'validation') 
ggdat <- plot_output$tab
names(ggdat) <- c("algo", "AUC_mean","TSS_mean","AUC_sd","TSS_sd")
ggdat$algo <- as.factor(ggdat$algo)



get_variables_importance(myBiomodModelOut)



bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# random forest here can get cut out, TSS is horrible
# none of the ML models do better than the sdmTMB models, 
#the RMSE and Tjur cross validation metrics, which is not possible to do through biomod after much testing
auc_tab <- bm_PlotEvalMean( bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = "validation" )$tab
cv_tbl <- auc_tab %>% rename( algo = name, mean_AUC_validation = mean1, mean_TSS_validation = mean2, sd_AUC_validation = sd1, sd_TSS_validation = sd2 )

save(cv_tbl, file = "code/output_data/model_results/eelgrass_nep_cv_metrics_biomod2.RData")

# GBM CV without biomod to get Tjur and RMSE and also becasue default values are dumn
library(gbm)
library(pROC)

tjur_fun <- function(y, p) {
  mean(p[y == 1], na.rm = TRUE) - mean(p[y == 0], na.rm = TRUE)
}

rmse_fun <- function(y, p) {
  sqrt(mean((y - p)^2, na.rm = TRUE))
}

auc_fun <- function(y, p) {
  as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
}

tss_fun <- function(y, p) {
  roc_obj <- pROC::roc(y, p, quiet = TRUE)
  # threshold that maximizes sensitivity + specificity
  coords_obj <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  tss <- coords_obj["sensitivity"] + coords_obj["specificity"] - 1
  as.numeric(tss)
}

gbm_cv_metrics <- function(dat, predictors, response = "presence", fold_col = "fold") {
  
  folds <- sort(unique(dat[[fold_col]]))
  out <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    f <- folds[i]
    
    train <- dat[dat[[fold_col]] != f, ]
    test  <- dat[dat[[fold_col]] == f, ]
    
    gbm_fit <- gbm(
      formula = as.formula(
        paste(response, "~", paste(predictors, collapse = " + "))
      ),
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
    
    out[[i]] <- data.frame(
      algo = "GBM_robust",
      fold = f,
      RMSE = rmse_fun(y, p),
      Tjur = tjur_fun(y, p),
      AUC  = auc_fun(y, p),
      TSS  = tss_fun(y, p)
    )
  }
  
  bind_rows(out)
}

gbm_cv_res <- gbm_cv_metrics(dat2, predictors = pred_vars)

gbm_cv_res          # per-fold metrics
gbm_summary <- gbm_cv_res %>%
  summarise(
    algo = "GBM_robust",
    mean_RMSE = mean(RMSE),
    mean_Tjur = mean(Tjur),
    mean_AUC  = mean(AUC),
    mean_TSS  = mean(TSS)
  )

gbm_summary

save(gbm_summary, file = "code/output_data/model_results/eelgrass_nep_cv_metricsextraGBM_biomod2.RData")

# test forecasting
data_pre2013 <- dat2 %>% filter(Year < 2010)
data_post2012 <- dat2 %>% filter(Year > 2012)

# Format data for BIOMOD2
myBiomodData_pre2013 <- BIOMOD_FormatingData(
  resp.var  = data_pre2013$presence,
  expl.var  = data_pre2013 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = data_pre2013[, c("X_m", "Y_m")],
  resp.name = "eelgrass",
  PA.nb.rep = 0
)

n <- nrow(data_pre2013)
cv_table <- matrix(TRUE, nrow = n, ncol = 1)  # 1 column = 1 run
colnames(cv_table) <- "_allData_allRun"       # required by BIOMOD2


myBiomodModelOut_pre2013 <- BIOMOD_Modeling(
  bm.format       = myBiomodData_pre2013,
  models          = ml_models,
  CV.strategy     = "user.defined",
  CV.user.table   = cv_table,
  CV.do.full.models = TRUE,  # ensures BIOMOD2 builds the “full model”
  metric.eval     = c("AUCroc","TSS","KAPPA"),
  var.import      = 3,
  seed.val        = 123,
  modeling.id     = "pre2013_ML"
)


library(randomForest)
library(xgboost)
library(gbm)
library(ModelMetrics)  
library(nnet)       
library(MASS) 
library(rpart) 


# Random Forest
rf_mod <- randomForest(
  x = data_pre2013[, pred_vars],
  y = as.factor(data_pre2013$presence),
  ntree = 500
)

# GBM
gbm_mod <- gbm(
  formula = presence ~ .,
  distribution = "bernoulli",
  data = data_pre2013[, c("presence", pred_vars)],
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.005,
  verbose = FALSE
)

# XGBoost
xgb_matrix_pre <- model.matrix(~ . -1, data = data_pre2013[, pred_vars])
xgb_matrix_post <- model.matrix(~ . -1, data = data_post2012[, pred_vars])

dtrain <- xgb.DMatrix(
  data = xgb_matrix_pre,
  label = data_pre2013$presence
)

xgb_params <- list(
  objective = "binary:logistic",  # presence/absence
  eval_metric = "auc",
  max_depth = 6,
  eta = 0.1,
  nthread = 4
)

xgb_mod <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 200,
  verbose = 1
)

# --- ANN ---

data_pre2013_num <- data_pre2013 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.)))) %>%
  replace(is.na(.), 0)  
data_post2012_num <- data_post2012 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.))))%>%
  replace(is.na(.), 0)  
#data_pre2013_num$substrate <- as.factor(data_pre2013$substrate)

y_pre_1hot <- class.ind(as.factor(data_pre2013$presence))  # one-hot encoding: 0 -> [1,0], 1 -> [0,1]

# -------------------------------
# 3️⃣ Train ANN
# -------------------------------
set.seed(123)
ann_mod <- nnet(
  x = as.matrix(data_pre2013_num),
  y = y_pre_1hot,
  size = 5,       # hidden neurons
  linout = FALSE,
  entropy = TRUE,
  maxit = 500,
  trace = FALSE
)



# --- CTA (Classification Tree Analysis) ---
cta_mod <- rpart(
  formula = presence ~ .,
  data = cbind(presence = as.factor(data_pre2013$presence), data_pre2013_num),
  method = "class",
  control = rpart.control(minsplit = 5, cp = 0.001)
)




predictions <- list(
  RF = list(
    training = predict(rf_mod, data_pre2013[, pred_vars], type = "prob")[,2],
    forecast = predict(rf_mod, data_post2012[, pred_vars], type = "prob")[,2]
  ),
  GBM = list(
    training = predict(gbm_mod, data_pre2013[, pred_vars], n.trees = 2000, type = "response"),
    forecast = predict(gbm_mod, data_post2012[, pred_vars], n.trees = 2000, type = "response")
  ),
  XGBOOST = list(
    training = predict(xgb_mod, newdata = xgb_matrix_pre),
    forecast = predict(xgb_mod, newdata = xgb_matrix_post)
  ),
  ANN = list(
    training = predict(ann_mod, newdata = as.matrix(data_pre2013_num))[,2],
    forecast = predict(ann_mod, newdata = as.matrix(data_post2012_num))[,2]
  ),
  CTA = list(
    training = predict(cta_mod, newdata = data_pre2013_num, type = "prob")[,2],
    forecast = predict(cta_mod, newdata = data_post2012_num, type = "prob")[,2]
  )
)


forecast_results <- lapply(names(predictions), function(mod) {
  
  obs_train <- data_pre2013$presence
  pred_train <- predictions[[mod]]$training
  
  obs_fore  <- data_post2012$presence
  pred_fore <- predictions[[mod]]$forecast
  
  # Model-specific threshold from pre-2013 data
  threshold_mod <- get_optimal_threshold(obs_train, pred_train)
  
  data.frame(
    model = mod,
    AUC_forecast  = ModelMetrics::auc(obs_fore, pred_fore),
    Tjur_forecast = tjur(obs_fore, pred_fore),
    RMSE_forecast = rmse(obs_fore, pred_fore),
    TSS_forecast  = tss_metric(obs_fore, pred_fore, threshold_mod),
    threshold_used = threshold_mod
  )
}) %>% bind_rows()

print(forecast_results)
save(forecast_results, file = "code/output_data/model_results/eelgrass_nep_forecast_metrics_biomod2.RData")
# only the GBM model ML models is comparable to the sdm tmb model 

#### make projection of GBM model only
#this was not working in biomod!
load("code/output_data/prediction_model_inputs.RData")
env<- env_20m_all %>% dplyr::select(depth_stnd, slope_stnd, rei_stnd, substrate, airtempmin_stnd, rsdsmin_stnd, prmin_stnd, 
                                    saltcv_nep_stnd, tempcv_nep_stnd)

set.seed(123)

# make a full gbm model without cv
gbm_fit <- gbm(
  presence ~ .,
  data = dat2[, c("presence", pred_vars)],
  distribution = "bernoulli",
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.005,
  bag.fraction = 0.5,
  n.minobsinnode = 10,
  train.fraction = 1,
  verbose = FALSE
)

gbm_prob <- predict(
  gbm_fit,
  newdata = env,
  n.trees = 5000,
  type = "response"   # returns probabilities for bernoulli 
)

range(gbm_prob)


env_20m_all$gbm_prob <- gbm_prob

outdir <- file.path("./raster/eelgrass/gbm_nep")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_gbm_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_gbm_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_gbm_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_gbm_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_gbm_nep.tif")), overwrite = TRUE)





####Surfgrass bccm model####
sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)
data <- data %>% mutate(Survey = as.factor(Survey))


pred_vars <- c("depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd", "substrate", "cul_eff_stnd", 
               "airtempcv_stnd", "rsdsmin_stnd", "prmean_stnd", 
               "saltcv_bccm_stnd", "NO3_bccm_stnd", "tempmin_bccm_stnd", "surftempcv_bccm_stnd", "Survey")

dat2 <- data %>%
  dplyr::select(presence, X_m, Y_m, fold, Year, all_of(pred_vars)) %>%
  na.omit()

dat2$presence <- as.numeric(dat2$presence)

myBiomodData <- BIOMOD_FormatingData(
  resp.var  = dat2$presence,
  expl.var  = dat2 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = dat2[, c("X_m", "Y_m")],
  resp.name = "surfgrass",
  #PA.nb.rep = 0,
  PA.nb.rep = 1, # <-- REQUIRED for CV to work 
  PA.nb.absences = 0, # <-- keeps your true absences 
  PA.strategy = "none") # <-- prevents BIOMOD from generating fake absences


folds <- sort(unique(dat2$fold))
K <- length(folds)

cv_table <- matrix(TRUE, nrow = nrow(dat2), ncol = K + 1)  # extra column for full model
for (k in seq_along(folds)) {
  test_fold <- folds[k]
  cv_table[dat2$fold == test_fold, k] <- FALSE  # use your fold column
}
colnames(cv_table) <- c(paste0("_allData_RUN", seq_len(K)), "_allData_allRun")

ml_models <- c("GBM", "RF", "XGBOOST", "ANN", "CTA")

myOptions <- bm_ModelingOptions(
  data.type = "binary",
  models = ml_models,
  strategy = "bigboss"  # stronger defaults than the "default" set
)

# need to change var import to 10-20 for final runs
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

# Extract validation thresholds that maximized TSS
tss_thresholds_df <- evals %>%
  filter(metric.eval == "TSS") %>%
  dplyr::select(algo, run, cutoff, validation)

print(tss_thresholds_df)

biomod_threshold_summary <- tss_thresholds_df %>%
  group_by(algo) %>%
  summarise(
    TSS_threshold_mean = mean(cutoff, na.rm = TRUE),
    TSS_threshold_sd   = sd(cutoff, na.rm = TRUE),
    mean_TSS_validation = mean(validation, na.rm = TRUE)
  )

print(biomod_threshold_summary)


plot_output <- bm_PlotEvalMean(bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = 'validation') 
ggdat <- plot_output$tab
names(ggdat) <- c("algo", "AUC_mean","TSS_mean","AUC_sd","TSS_sd")
ggdat$algo <- as.factor(ggdat$algo)



get_variables_importance(myBiomodModelOut)



bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# random forest here can get cut out, TSS is horrible
# none of the ML models do better than the sdmTMB models, BUT becasue the GBM bccm model forecast is comparable want to add without biomod
#the RMSE and Tjur cross validation metrics, which is not possible to do through biomod after much testing
auc_tab <- bm_PlotEvalMean( bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = "validation" )$tab
cv_tbl <- auc_tab %>% rename( algo = name, mean_AUC_validation = mean1, mean_TSS_validation = mean2, sd_AUC_validation = sd1, sd_TSS_validation = sd2 )

save(cv_tbl, file = "code/output_data/model_results/surfgrass_bccm_cv_metrics_biomod2.RData")

# GBM CV without biomod to get Tjur and RMSE and also becasue default values are dumn
library(gbm)
library(pROC)

tjur_fun <- function(y, p) {
  mean(p[y == 1], na.rm = TRUE) - mean(p[y == 0], na.rm = TRUE)
}

rmse_fun <- function(y, p) {
  sqrt(mean((y - p)^2, na.rm = TRUE))
}

auc_fun <- function(y, p) {
  as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
}

tss_fun <- function(y, p) {
  roc_obj <- pROC::roc(y, p, quiet = TRUE)
  # threshold that maximizes sensitivity + specificity
  coords_obj <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  tss <- coords_obj["sensitivity"] + coords_obj["specificity"] - 1
  as.numeric(tss)
}

gbm_cv_metrics <- function(dat, predictors, response = "presence", fold_col = "fold") {
  
  folds <- sort(unique(dat[[fold_col]]))
  out <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    f <- folds[i]
    
    train <- dat[dat[[fold_col]] != f, ]
    test  <- dat[dat[[fold_col]] == f, ]
    
    gbm_fit <- gbm(
      formula = as.formula(
        paste(response, "~", paste(predictors, collapse = " + "))
      ),
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
    
    out[[i]] <- data.frame(
      algo = "GBM_robust",
      fold = f,
      RMSE = rmse_fun(y, p),
      Tjur = tjur_fun(y, p),
      AUC  = auc_fun(y, p),
      TSS  = tss_fun(y, p)
    )
  }
  
  bind_rows(out)
}

gbm_cv_res <- gbm_cv_metrics(dat2, predictors = pred_vars)

gbm_cv_res          # per-fold metrics
gbm_summary <- gbm_cv_res %>%
  summarise(
    algo = "GBM_robust",
    mean_RMSE = mean(RMSE),
    mean_Tjur = mean(Tjur),
    mean_AUC  = mean(AUC),
    mean_TSS  = mean(TSS)
  )

gbm_summary

save(gbm_summary, file = "code/output_data/model_results/surfgrass_bccm_cv_metricsextraGBM_biomod2.RData")

# test forecasting
data_pre2013 <- dat2 %>% filter(Year < 2010)
data_post2012 <- dat2 %>% filter(Year > 2012)

# Format data for BIOMOD2
myBiomodData_pre2013 <- BIOMOD_FormatingData(
  resp.var  = data_pre2013$presence,
  expl.var  = data_pre2013 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = data_pre2013[, c("X_m", "Y_m")],
  resp.name = "surfgrass",
  PA.nb.rep = 0
)

n <- nrow(data_pre2013)
cv_table <- matrix(TRUE, nrow = n, ncol = 1)  # 1 column = 1 run
colnames(cv_table) <- "_allData_allRun"       # required by BIOMOD2


myBiomodModelOut_pre2013 <- BIOMOD_Modeling(
  bm.format       = myBiomodData_pre2013,
  models          = ml_models,
  CV.strategy     = "user.defined",
  CV.user.table   = cv_table,
  CV.do.full.models = TRUE,  # ensures BIOMOD2 builds the “full model”
  metric.eval     = c("AUCroc","TSS","KAPPA"),
  var.import      = 3,
  seed.val        = 123,
  modeling.id     = "pre2013_ML"
)


library(randomForest)
library(xgboost)
library(gbm)
library(ModelMetrics)  
library(nnet)       
library(MASS) 
library(rpart) 


# Random Forest
rf_mod <- randomForest(
  x = data_pre2013[, pred_vars],
  y = as.factor(data_pre2013$presence),
  ntree = 500
)

# GBM
gbm_mod <- gbm(
  formula = presence ~ .,
  distribution = "bernoulli",
  data = data_pre2013[, c("presence", pred_vars)],
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.005,
  verbose = FALSE
)

# XGBoost
xgb_matrix_pre <- model.matrix(~ . -1, data = data_pre2013[, pred_vars])
xgb_matrix_post <- model.matrix(~ . -1, data = data_post2012[, pred_vars])

dtrain <- xgb.DMatrix(
  data = xgb_matrix_pre,
  label = data_pre2013$presence
)

xgb_params <- list(
  objective = "binary:logistic",  # presence/absence
  eval_metric = "auc",
  max_depth = 6,
  eta = 0.1,
  nthread = 4
)

xgb_mod <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 200,
  verbose = 1
)

# --- ANN ---

data_pre2013_num <- data_pre2013 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.)))) %>%
  replace(is.na(.), 0)  
data_post2012_num <- data_post2012 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.))))%>%
  replace(is.na(.), 0)  
#data_pre2013_num$substrate <- as.factor(data_pre2013$substrate)

y_pre_1hot <- class.ind(as.factor(data_pre2013$presence))  # one-hot encoding: 0 -> [1,0], 1 -> [0,1]

# -------------------------------
# 3️⃣ Train ANN
# -------------------------------
set.seed(123)
ann_mod <- nnet(
  x = as.matrix(data_pre2013_num),
  y = y_pre_1hot,
  size = 5,       # hidden neurons
  linout = FALSE,
  entropy = TRUE,
  maxit = 500,
  trace = FALSE
)



# --- CTA (Classification Tree Analysis) ---
cta_mod <- rpart(
  formula = presence ~ .,
  data = cbind(presence = as.factor(data_pre2013$presence), data_pre2013_num),
  method = "class",
  control = rpart.control(minsplit = 5, cp = 0.001)
)




predictions <- list(
  RF = list(
    training = predict(rf_mod, data_pre2013[, pred_vars], type = "prob")[,2],
    forecast = predict(rf_mod, data_post2012[, pred_vars], type = "prob")[,2]
  ),
  GBM = list(
    training = predict(gbm_mod, data_pre2013[, pred_vars], n.trees = 2000, type = "response"),
    forecast = predict(gbm_mod, data_post2012[, pred_vars], n.trees = 2000, type = "response")
  ),
  XGBOOST = list(
    training = predict(xgb_mod, newdata = xgb_matrix_pre),
    forecast = predict(xgb_mod, newdata = xgb_matrix_post)
  ),
  ANN = list(
    training = predict(ann_mod, newdata = as.matrix(data_pre2013_num))[,2],
    forecast = predict(ann_mod, newdata = as.matrix(data_post2012_num))[,2]
  ),
  CTA = list(
    training = predict(cta_mod, newdata = data_pre2013_num, type = "prob")[,2],
    forecast = predict(cta_mod, newdata = data_post2012_num, type = "prob")[,2]
  )
)


forecast_results <- lapply(names(predictions), function(mod) {
  
  obs_train <- data_pre2013$presence
  pred_train <- predictions[[mod]]$training
  
  obs_fore  <- data_post2012$presence
  pred_fore <- predictions[[mod]]$forecast
  
  # Model-specific threshold from pre-2013 data
  threshold_mod <- get_optimal_threshold(obs_train, pred_train)
  
  data.frame(
    model = mod,
    AUC_forecast  = ModelMetrics::auc(obs_fore, pred_fore),
    Tjur_forecast = tjur(obs_fore, pred_fore),
    RMSE_forecast = rmse(obs_fore, pred_fore),
    TSS_forecast  = tss_metric(obs_fore, pred_fore, threshold_mod),
    threshold_used = threshold_mod
  )
}) %>% bind_rows()

print(forecast_results)
save(forecast_results, file = "code/output_data/model_results/surfgrass_bccm_forecast_metrics_biomod2.RData")
# none of the ML models are comparable to the sdm tmb model 





####surfgrass model with NEP####
sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)
data <- data %>% mutate(Survey = as.factor(Survey))


pred_vars <- c("depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd", "substrate", "cul_eff_stnd", 
               "airtempcv_stnd", "rsdsmin_stnd", "prmean_stnd", 
               "saltmean_nep_stnd", "tempmean_nep_stnd", "surftempcv_nep_stnd", "Survey")

dat2 <- data %>%
  dplyr::select(presence, X_m, Y_m, fold, Year, all_of(pred_vars)) %>%
  na.omit()

dat2$presence <- as.numeric(dat2$presence)

myBiomodData <- BIOMOD_FormatingData(
  resp.var  = dat2$presence,
  expl.var  = dat2 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = dat2[, c("X_m", "Y_m")],
  resp.name = "surfgrass",
  #PA.nb.rep = 0,
  PA.nb.rep = 1, # <-- REQUIRED for CV to work 
  PA.nb.absences = 0, # <-- keeps your true absences 
  PA.strategy = "none") # <-- prevents BIOMOD from generating fake absences


folds <- sort(unique(dat2$fold))
K <- length(folds)

cv_table <- matrix(TRUE, nrow = nrow(dat2), ncol = K + 1)  # extra column for full model
for (k in seq_along(folds)) {
  test_fold <- folds[k]
  cv_table[dat2$fold == test_fold, k] <- FALSE  # use your fold column
}
colnames(cv_table) <- c(paste0("_allData_RUN", seq_len(K)), "_allData_allRun")

ml_models <- c("GBM", "RF", "XGBOOST", "ANN", "CTA")

myOptions <- bm_ModelingOptions(
  data.type = "binary",
  models = ml_models,
  strategy = "bigboss"  # stronger defaults than the "default" set
)

# need to change var import to 10-20 for final runs
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

# Extract validation thresholds that maximized TSS
tss_thresholds_df <- evals %>%
  filter(metric.eval == "TSS") %>%
  dplyr::select(algo, run, cutoff, validation)

print(tss_thresholds_df)

biomod_threshold_summary <- tss_thresholds_df %>%
  group_by(algo) %>%
  summarise(
    TSS_threshold_mean = mean(cutoff, na.rm = TRUE),
    TSS_threshold_sd   = sd(cutoff, na.rm = TRUE),
    mean_TSS_validation = mean(validation, na.rm = TRUE)
  )

print(biomod_threshold_summary)


plot_output <- bm_PlotEvalMean(bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = 'validation') 
ggdat <- plot_output$tab
names(ggdat) <- c("algo", "AUC_mean","TSS_mean","AUC_sd","TSS_sd")
ggdat$algo <- as.factor(ggdat$algo)



get_variables_importance(myBiomodModelOut)



bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# random forest here can get cut out, TSS is horrible
# none of the ML models do better than the sdmTMB models, 
#the RMSE and Tjur cross validation metrics, which is not possible to do through biomod after much testing
auc_tab <- bm_PlotEvalMean( bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = "validation" )$tab
cv_tbl <- auc_tab %>% rename( algo = name, mean_AUC_validation = mean1, mean_TSS_validation = mean2, sd_AUC_validation = sd1, sd_TSS_validation = sd2 )

save(cv_tbl, file = "code/output_data/model_results/surfgrass_nep_cv_metrics_biomod2.RData")

# GBM CV without biomod to get Tjur and RMSE and also becasue default values are dumn
library(gbm)
library(pROC)

tjur_fun <- function(y, p) {
  mean(p[y == 1], na.rm = TRUE) - mean(p[y == 0], na.rm = TRUE)
}

rmse_fun <- function(y, p) {
  sqrt(mean((y - p)^2, na.rm = TRUE))
}

auc_fun <- function(y, p) {
  as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
}

tss_fun <- function(y, p) {
  roc_obj <- pROC::roc(y, p, quiet = TRUE)
  # threshold that maximizes sensitivity + specificity
  coords_obj <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  tss <- coords_obj["sensitivity"] + coords_obj["specificity"] - 1
  as.numeric(tss)
}

gbm_cv_metrics <- function(dat, predictors, response = "presence", fold_col = "fold") {
  
  folds <- sort(unique(dat[[fold_col]]))
  out <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    f <- folds[i]
    
    train <- dat[dat[[fold_col]] != f, ]
    test  <- dat[dat[[fold_col]] == f, ]
    
    gbm_fit <- gbm(
      formula = as.formula(
        paste(response, "~", paste(predictors, collapse = " + "))
      ),
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
    
    out[[i]] <- data.frame(
      algo = "GBM_robust",
      fold = f,
      RMSE = rmse_fun(y, p),
      Tjur = tjur_fun(y, p),
      AUC  = auc_fun(y, p),
      TSS  = tss_fun(y, p)
    )
  }
  
  bind_rows(out)
}

gbm_cv_res <- gbm_cv_metrics(dat2, predictors = pred_vars)

gbm_cv_res          # per-fold metrics
gbm_summary <- gbm_cv_res %>%
  summarise(
    algo = "GBM_robust",
    mean_RMSE = mean(RMSE),
    mean_Tjur = mean(Tjur),
    mean_AUC  = mean(AUC),
    mean_TSS  = mean(TSS)
  )

gbm_summary

save(gbm_summary, file = "code/output_data/model_results/surfgrass_nep_cv_metricsextraGBM_biomod2.RData")

# test forecasting
data_pre2013 <- dat2 %>% filter(Year < 2010)
data_post2012 <- dat2 %>% filter(Year > 2012)

# Format data for BIOMOD2
myBiomodData_pre2013 <- BIOMOD_FormatingData(
  resp.var  = data_pre2013$presence,
  expl.var  = data_pre2013 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = data_pre2013[, c("X_m", "Y_m")],
  resp.name = "surfgrass",
  PA.nb.rep = 0
)

n <- nrow(data_pre2013)
cv_table <- matrix(TRUE, nrow = n, ncol = 1)  # 1 column = 1 run
colnames(cv_table) <- "_allData_allRun"       # required by BIOMOD2


myBiomodModelOut_pre2013 <- BIOMOD_Modeling(
  bm.format       = myBiomodData_pre2013,
  models          = ml_models,
  CV.strategy     = "user.defined",
  CV.user.table   = cv_table,
  CV.do.full.models = TRUE,  # ensures BIOMOD2 builds the “full model”
  metric.eval     = c("AUCroc","TSS","KAPPA"),
  var.import      = 3,
  seed.val        = 123,
  modeling.id     = "pre2013_ML"
)


library(randomForest)
library(xgboost)
library(gbm)
library(ModelMetrics)  
library(nnet)       
library(MASS) 
library(rpart) 


# Random Forest
rf_mod <- randomForest(
  x = data_pre2013[, pred_vars],
  y = as.factor(data_pre2013$presence),
  ntree = 500
)

# GBM
gbm_mod <- gbm(
  formula = presence ~ .,
  distribution = "bernoulli",
  data = data_pre2013[, c("presence", pred_vars)],
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.005,
  verbose = FALSE
)

# XGBoost
xgb_matrix_pre <- model.matrix(~ . -1, data = data_pre2013[, pred_vars])
xgb_matrix_post <- model.matrix(~ . -1, data = data_post2012[, pred_vars])

dtrain <- xgb.DMatrix(
  data = xgb_matrix_pre,
  label = data_pre2013$presence
)

xgb_params <- list(
  objective = "binary:logistic",  # presence/absence
  eval_metric = "auc",
  max_depth = 6,
  eta = 0.1,
  nthread = 4
)

xgb_mod <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 200,
  verbose = 1
)

# --- ANN ---

data_pre2013_num <- data_pre2013 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.)))) %>%
  replace(is.na(.), 0)  
data_post2012_num <- data_post2012 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.))))%>%
  replace(is.na(.), 0)  
#data_pre2013_num$substrate <- as.factor(data_pre2013$substrate)

y_pre_1hot <- class.ind(as.factor(data_pre2013$presence))  # one-hot encoding: 0 -> [1,0], 1 -> [0,1]

# -------------------------------
# 3️⃣ Train ANN
# -------------------------------
set.seed(123)
ann_mod <- nnet(
  x = as.matrix(data_pre2013_num),
  y = y_pre_1hot,
  size = 5,       # hidden neurons
  linout = FALSE,
  entropy = TRUE,
  maxit = 500,
  trace = FALSE
)



# --- CTA (Classification Tree Analysis) ---
cta_mod <- rpart(
  formula = presence ~ .,
  data = cbind(presence = as.factor(data_pre2013$presence), data_pre2013_num),
  method = "class",
  control = rpart.control(minsplit = 5, cp = 0.001)
)




predictions <- list(
  RF = list(
    training = predict(rf_mod, data_pre2013[, pred_vars], type = "prob")[,2],
    forecast = predict(rf_mod, data_post2012[, pred_vars], type = "prob")[,2]
  ),
  GBM = list(
    training = predict(gbm_mod, data_pre2013[, pred_vars], n.trees = 2000, type = "response"),
    forecast = predict(gbm_mod, data_post2012[, pred_vars], n.trees = 2000, type = "response")
  ),
  XGBOOST = list(
    training = predict(xgb_mod, newdata = xgb_matrix_pre),
    forecast = predict(xgb_mod, newdata = xgb_matrix_post)
  ),
  ANN = list(
    training = predict(ann_mod, newdata = as.matrix(data_pre2013_num))[,2],
    forecast = predict(ann_mod, newdata = as.matrix(data_post2012_num))[,2]
  ),
  CTA = list(
    training = predict(cta_mod, newdata = data_pre2013_num, type = "prob")[,2],
    forecast = predict(cta_mod, newdata = data_post2012_num, type = "prob")[,2]
  )
)


forecast_results <- lapply(names(predictions), function(mod) {
  
  obs_train <- data_pre2013$presence
  pred_train <- predictions[[mod]]$training
  
  obs_fore  <- data_post2012$presence
  pred_fore <- predictions[[mod]]$forecast
  
  # Model-specific threshold from pre-2013 data
  threshold_mod <- get_optimal_threshold(obs_train, pred_train)
  
  data.frame(
    model = mod,
    AUC_forecast  = ModelMetrics::auc(obs_fore, pred_fore),
    Tjur_forecast = tjur(obs_fore, pred_fore),
    RMSE_forecast = rmse(obs_fore, pred_fore),
    TSS_forecast  = tss_metric(obs_fore, pred_fore, threshold_mod),
    threshold_used = threshold_mod
  )
}) %>% bind_rows()

print(forecast_results)
save(forecast_results, file = "code/output_data/model_results/surfgrass_nep_forecast_metrics_biomod2.RData")
# none of the models are comparable to sdmtmb
