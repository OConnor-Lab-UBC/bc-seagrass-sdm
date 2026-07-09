###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# fit sdm models with GBM and XGBOOST 
###############################################################################


library(dplyr)
library(terra)
library(gbm)
library(pROC)
library(randomForest)
library(xgboost)
library(ModelMetrics)  
library(nnet)       
library(MASS) 
library(rpart)

# Load functions
source("code/modelling-functions.R")
# Load data
load("code/output_data/seagrass_model_inputs.RData")
# Preprocess
seagrass_data_long <- seagrass_data_long %>%
  dplyr::select(
    -saltmean_bccm_sq_stnd, -saltmean_nep_sq_stnd,
    -slope_sqrt_stnd, -saltmin_bccm_sq_stnd, -saltmin_nep_sq_stnd
  ) %>%
  mutate(
    Survey = as.factor(substr(HKey, 1, 3)),
    HKey = as.factor(HKey),
    Year_factor = as.factor(Year)
  )

#having trouble with MSE survey causing b_j standard error warning as it was only absences for sea grass so combined it with BHm because people doing data would be similar and there would be care in iding correctly
seagrass_data_long <- seagrass_data_long %>%
  mutate(Survey = recode(Survey, MSE = "BHM")) %>%
  mutate(Survey = recode(Survey, GSU = "RSU")) %>%
  mutate(Survey = recode(Survey, ABL = "RSU")) %>%
  mutate(Survey = factor(Survey))



#GBM and GXBOOST are the best models so need to tune to find best hyper parameters
## eelgrass
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
dat_eel <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass, inner_fold = innerfold_eelgrass)

# variables
pred_vars_bccm <- c(
  "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
  "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
  "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempcv_bccm_stnd", "Survey"
)

pred_vars_nep <- c(
  "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
  "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
  "saltcv_nep_stnd", "NH4_nep_stnd", "tempcv_nep_stnd", "Survey"
)
 
# Hyperparameter Grids
gbm_grid <- expand.grid(
  interaction.depth = c(2,3,4),
  shrinkage = c(0.001,0.005,0.01),
  n.minobsinnode = c(5,10)
)

xgb_grid <- expand.grid(
  eta = c(0.01, 0.05),
  max_depth = c(3, 5, 7),
  subsample = c(0.7, 0.9),
  colsample_bytree = c(0.7, 0.9),
  nrounds = c(500, 1000, 2000)
)

### nested spatial cross validation to choose hyper parameters and report on metrics
gbm_eel_bccm <- nested_gbm_spatial(
  dat = dat_eel,
  predictors = pred_vars_bccm,
  outer_col = "fold",
  inner_col = "inner_fold",
  gbm_grid = gbm_grid
)
save(gbm_eel_bccm, file = "code/output_data/model_results/eelgrass_eval_cv_gbm_bccm.RData")
gbm_eel_bccm$fold_metrics %>%
  dplyr::count(interaction.depth, shrinkage, n.minobsinnode, n.trees, sort = TRUE)
gbm_eel_bccm$fold_metrics %>%
  dplyr::count(interaction.depth, shrinkage, n.minobsinnode, sort = TRUE)
gbm_eel_bccm$fold_metrics %>%
  dplyr::group_by(interaction.depth, shrinkage, n.minobsinnode) %>%
  dplyr::summarise(
    mean_test_auc = mean(test_auc, na.rm = TRUE),
    mean_test_logloss = mean(test_logloss, na.rm = TRUE),
    mean_test_tss = mean(test_tss, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_test_logloss)

xgb_eel_bccm <- nested_xgb_spatial(
  dat = dat_eel,
  predictors = pred_vars_bccm,
  outer_col = "fold",
  inner_col = "inner_fold",
  xgb_grid = xgb_grid
)

save(xgb_eel_bccm, file = "code/output_data/model_results/eelgrass_eval_cv_xgb_bccm.RData")
xgb_eel_bccm$fold_metrics %>%
  dplyr::count(eta, max_depth, subsample, colsample_bytree, nrounds, sort = TRUE)
xgb_eel_bccm$fold_metrics %>%
  dplyr::group_by(eta, max_depth, subsample, colsample_bytree, nrounds) %>%
  dplyr::summarise(
    mean_test_auc = mean(test_auc, na.rm = TRUE),
    mean_test_logloss = mean(test_logloss, na.rm = TRUE),
    mean_test_tss = mean(test_tss, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_test_logloss)


gbm_eel_nep <- nested_gbm_spatial(
  dat = dat_eel,
  predictors = pred_vars_nep,
  outer_col = "fold",
  inner_col = "inner_fold",
  gbm_grid = gbm_grid
)

save(gbm_eel_nep, file = "code/output_data/model_results/eelgrass_eval_cv_gbm_nep.RData")

gbm_eel_nep$fold_metrics %>%
  dplyr::count(interaction.depth, shrinkage, n.minobsinnode, n.trees, sort = TRUE)
gbm_eel_nep$fold_metrics %>%
  dplyr::count(interaction.depth, shrinkage, n.minobsinnode, sort = TRUE)
gbm_eel_nep$fold_metrics %>%
  dplyr::group_by(interaction.depth, shrinkage, n.minobsinnode) %>%
  dplyr::summarise(
    mean_test_auc = mean(test_auc, na.rm = TRUE),
    mean_test_logloss = mean(test_logloss, na.rm = TRUE),
    mean_test_tss = mean(test_tss, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_test_logloss)


xgb_eel_nep <- nested_xgb_spatial(
  dat = dat_eel,
  predictors = pred_vars_nep,
  outer_col = "fold",
  inner_col = "inner_fold",
  xgb_grid = xgb_grid
)

save(xgb_eel_nep, file = "code/output_data/model_results/eelgrass_eval_cv_xgb_nep.RData")
xgb_eel_nep$fold_metrics %>%
  dplyr::count(eta, max_depth, subsample, colsample_bytree, nrounds, sort = TRUE)
xgb_eel_nep$fold_metrics %>%
  dplyr::group_by(eta, max_depth, subsample, colsample_bytree, nrounds) %>%
  dplyr::summarise(
    mean_test_auc = mean(test_auc, na.rm = TRUE),
    mean_test_logloss = mean(test_logloss, na.rm = TRUE),
    mean_test_tss = mean(test_tss, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_test_logloss)




# temporal validation
#update variables to remove survey as cannot be included in temporal validation
pred_vars_bccm <- c(
  "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
  "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
  "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempcv_bccm_stnd"
)

pred_vars_nep <- c(
  "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
  "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
  "saltcv_nep_stnd", "NH4_nep_stnd", "tempcv_nep_stnd"
)

data_pre2013_eel <- dat_eel %>%
  filter(Year < 2010) %>%
  mutate(row_id = row_number())
names(data_pre2013_eel)

unique(data_pre2013_eel$Survey) # can't include survey becasue only RSU, Cuke and GDK in early years
test_set <- dat_eel %>% filter(Year > 2012)
obs_test <- test_set$presence


# based off most chosen parameters
best_gbm_bccm <- list(
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.001,
  n.minobsinnode = 10
)

best_gbm_nep <- list(
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.001,
  n.minobsinnode = 10
)

best_xgb_bccm <- list(
  eta = 0.01, 
  max_depth = 3, 
  subsample = 0.9, 
  colsample_bytree = 0.7, 
  nrounds = 1000
)

best_xgb_nep <- list(
  eta = 0.01, 
  max_depth = 5, 
  subsample = 0.7, 
  colsample_bytree = 0.9, 
  nrounds = 500
)

cv_thr_gbm_bccm <- get_global_cv_threshold_gbm(
  data_train = data_pre2013_eel,
  predictors = pred_vars_bccm,
  response = "presence",
  n.trees = best_gbm_bccm$n.trees,
  interaction.depth = best_gbm_bccm$interaction.depth,
  shrinkage = best_gbm_bccm$shrinkage,
  n.minobsinnode = best_gbm_bccm$n.minobsinnode
)
forecast_gbm_bccm_eelgrass <- run_temporal_gbm_forecast(
  data_pre2013 = data_pre2013_eel,
  test_set = test_set,
  predictors = pred_vars_bccm,
  threshold = cv_thr_gbm_bccm$threshold,
  response = "presence",
  n.trees = best_gbm_bccm$n.trees,
  interaction.depth = best_gbm_bccm$interaction.depth,
  shrinkage = best_gbm_bccm$shrinkage,
  n.minobsinnode = best_gbm_bccm$n.minobsinnode,
  model_name = "GBM_BCCM"
)

cv_thr_gbm_nep <- get_global_cv_threshold_gbm(
  data_train = data_pre2013_eel,
  predictors = pred_vars_nep,
  response = "presence",
  n.trees = best_gbm_nep$n.trees,
  interaction.depth = best_gbm_nep$interaction.depth,
  shrinkage = best_gbm_nep$shrinkage,
  n.minobsinnode = best_gbm_nep$n.minobsinnode
)

forecast_gbm_nep_eelgrass <- run_temporal_gbm_forecast(
  data_pre2013 = data_pre2013_eel,
  test_set = test_set,
  predictors = pred_vars_nep,
  threshold = cv_thr_gbm_nep$threshold,
  response = "presence",
  n.trees = best_gbm_nep$n.trees,
  interaction.depth = best_gbm_nep$interaction.depth,
  shrinkage = best_gbm_nep$shrinkage,
  n.minobsinnode = best_gbm_nep$n.minobsinnode,
  model_name = "GBM_NEP"
)

cv_thr_xgb_bccm <- get_global_cv_threshold_xgb(
  data_train = data_pre2013_eel,
  predictors = pred_vars_bccm,
  response = "presence",
  eta = best_xgb_bccm$eta,
  max_depth = best_xgb_bccm$max_depth,
  subsample = best_xgb_bccm$subsample,
  colsample_bytree = best_xgb_bccm$colsample_bytree,
  nrounds = best_xgb_bccm$nrounds
)
forecast_xgb_bccm_eelgrass <- run_temporal_xgb_forecast(
  data_pre2013 = data_pre2013_eel,
  test_set = test_set,
  predictors = pred_vars_bccm,
  threshold = cv_thr_xgb_bccm$threshold,
  response = "presence",
  eta = best_xgb_bccm$eta,
  max_depth = best_xgb_bccm$max_depth,
  subsample = best_xgb_bccm$subsample,
  colsample_bytree = best_xgb_bccm$colsample_bytree,
  nrounds = best_xgb_bccm$nrounds,
  model_name = "XGBOOST_BCCM"
)

cv_thr_xgb_nep <- get_global_cv_threshold_xgb(
  data_train = data_pre2013_eel,
  predictors = pred_vars_nep,
  response = "presence",
  eta = best_xgb_nep$eta,
  max_depth = best_xgb_nep$max_depth,
  subsample = best_xgb_nep$subsample,
  colsample_bytree = best_xgb_nep$colsample_bytree,
  nrounds = best_xgb_nep$nrounds
)
forecast_xgb_nep_eelgrass <- run_temporal_xgb_forecast(
  data_pre2013 = data_pre2013_eel,
  test_set = test_set,
  predictors = pred_vars_nep,
  threshold = cv_thr_xgb_nep$threshold,
  response = "presence",
  eta = best_xgb_nep$eta,
  max_depth = best_xgb_nep$max_depth,
  subsample = best_xgb_nep$subsample,
  colsample_bytree = best_xgb_nep$colsample_bytree,
  nrounds = best_xgb_nep$nrounds,
  model_name = "XGBOOST_NEP"
)

forecast_ml_eelgrass <- bind_rows(
  forecast_gbm_bccm_eelgrass,
  forecast_gbm_nep_eelgrass,
  forecast_xgb_bccm_eelgrass,
  forecast_xgb_nep_eelgrass
)

forecast_summary_ml_eelgrass <- forecast_ml_eelgrass %>%
  dplyr::filter(dataset == "forecast") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    AUC = mean(AUC, na.rm = TRUE),
    TjurR2 = mean(TjurR2, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    Brier = mean(Brier, na.rm = TRUE),
    LogLoss = mean(LogLoss, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    TSS = mean(TSS, na.rm = TRUE),
    Threshold = mean(threshold, na.rm = TRUE),
    .groups = "drop"
  )


save(forecast_ml_eelgrass, forecast_summary_ml_eelgrass, file = "code/output_data/model_results/forecast_ml_eelgrass_models.RData")










#surfgrass
sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
dat2 <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)

pred_vars <- c("depth_stnd", "tidal_sqrt_stnd", "rei_sqrt_stnd", "substrate", "cul_eff_stnd",
               "airtempcv_stnd", "prmean_stnd", "rsdsmin_stnd", "saltcv_bccm_stnd",
               "tempmin_bccm_stnd", "surftempcv_bccm_stnd")







