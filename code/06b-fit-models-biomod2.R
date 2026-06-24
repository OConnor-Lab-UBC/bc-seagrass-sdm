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

# TO DO 
### need to update to remove biomod form the picture. Just use GBM and XGBOOST 
## need nested spatial cross validation to choose hyperparameters for GBM and XGBOOST, need to add new metrics Brier, specificity, sensitivity as well as AUC, Tjur R2, and TSS
# for temporal validation need to add random folds and do 10 k folds
# make sure survey is included as fixed effect

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
 


###############################################################################
# Hyperparameter Grids
###############################################################################

gbm_grid <- expand.grid(
  
  n.trees = c(1000,2000,4000,6000),
  
  interaction.depth = c(2,3,4,5),
  
  shrinkage = c(0.001,0.005,0.01),
  
  n.minobsinnode = c(5,10,20)
  
)


xgb_grid <- expand.grid(
  
  eta = c(0.01,0.05),
  
  max_depth = c(3,5,7),
  
  subsample = c(0.7,0.9),
  
  colsample_bytree = c(0.7,0.9),
  
  nrounds = c(500,1000,2000)
  
)



gbm_eel_bccm <- nested_gbm_spatial(
  dat = dat_eel,
  predictors = pred_vars_bccm,
  outer_col = "fold",
  inner_col = "inner_fold",
  gbm_grid = gbm_grid
)

xgb_eel <- nested_xgb_spatial(
  dat = dat_eel,
  predictors = pred_vars,
  outer_col = "outer_fold",
  inner_col = "inner_fold",
  xgb_grid = xgb_grid
)



data_pre2013_eel <- dat_eel %>%
  filter(Year < 2010) %>%
  mutate(row_id = row_number())
names(data_pre2013)

unique(data_pre2013_eel$Survey)
test_set <- dat_eel %>% filter(Year > 2012)
obs_test <- test_set$presence



best_gbm <- list(
  n.trees = 4000,
  interaction.depth = 3,
  shrinkage = 0.005,
  n.minobsinnode = 10
)

best_xgb <- list(
  eta = 0.05,
  max_depth = 5,
  subsample = 0.7,
  colsample_bytree = 0.7,
  nrounds = 1500
)


forecast_gbm_bccm_eelgrass <- run_temporal_gbm(
  data_pre2013 = data_pre2013_eel,
  test_set = test_set,
  predictors = pred_vars_bccm,
  
  n.trees = best_gbm$n.trees,
  interaction.depth = best_gbm$interaction.depth,
  shrinkage = best_gbm$shrinkage,
  n.minobsinnode = best_gbm$n.minobsinnode
)

forecast_xgb_bccm_eelgrass <- run_temporal_xgb(
  data_pre2013 = data_pre2013_eel,
  test_set = test_set,
  predictors = pred_vars_bccm,
  
  eta = best_xgb$eta,
  max_depth = best_xgb$max_depth,
  subsample = best_xgb$subsample,
  colsample_bytree = best_xgb$colsample_bytree,
  nrounds = best_xgb$nrounds
)


forecast_ml <- bind_rows(
  forecast_gbm,
  forecast_xgb
)

forecast_summary <- forecast_ml %>%
  filter(dataset == "forecast") %>%
  group_by(model) %>%
  summarise(
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

gbm_temporal_threshold <-
  forecast_gbm %>%
  filter(dataset == "training") %>%
  summarise(
    threshold = mean(threshold)
  ) %>%
  pull()
``
xgb_temporal_threshold <-
  forecast_xgb %>%
  filter(dataset == "training") %>%
  summarise(
    threshold = mean(threshold)
  ) %>%
  pull()









#surfgrass
sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
dat2 <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)

pred_vars <- c("depth_stnd", "tidal_sqrt_stnd", "rei_sqrt_stnd", "substrate", "cul_eff_stnd",
               "airtempcv_stnd", "prmean_stnd", "rsdsmin_stnd", "saltcv_bccm_stnd",
               "tempmin_bccm_stnd", "surftempcv_bccm_stnd")







