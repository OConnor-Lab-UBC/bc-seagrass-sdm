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



#GBM and GXBOOST are the best models so need to tune to find best hyper parameters
## eelgrass
# sp = "ZO"
# numFolds <- length(unique(seagrass_data$fold_eelgrass))
# dat2 <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
# 
# pred_vars <- c(
#   "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
#   "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
#   "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempmin_bccm_stnd"
# )
# 
# #surfgrass
# sp = "PH"
# numFolds <- length(unique(seagrass_data$fold_seagrass))
# dat2 <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)
# 
# pred_vars <- c("depth_stnd", "tidal_sqrt_stnd", "rei_sqrt_stnd", "substrate", "cul_eff_stnd", 
#                  "airtempcv_stnd", "prmean_stnd", "rsdsmin_stnd", "saltcv_bccm_stnd", 
#                  "NO3_bccm_stnd", "tempmin_bccm_stnd", "surftempcv_bccm_stnd") 
# 
# gbm_grid <- expand.grid(
#   interaction.depth = c(2,3,4),
#   shrinkage = c(0.01, 0.005),
#   n.minobsinnode = c(5,10)
# )
# 
# gbm_tuning <- tune_gbm(
#   dat = dat2,
#   predictors = pred_vars,
#   gbm_grid = gbm_grid
# )
# 
# 
# head(gbm_tuning)
# 
# best_gbm <- gbm_tuning[1,]
# print (best_gbm)
# #eelgrass
# #    depth    lr minobs  trees       AUC     AUC_sd
# #    2 0.005     10 4696.2 0.9211391 0.01825039
# 
# #surfgrass
# #  depth    lr minobs  trees       AUC     AUC_sd
# #     4 0.005      5 2418.9 0.9446731 0.03416803
# 
# xgb_results <- tune_xgboost_spatial_parallel(
#   dat = dat2,                   # your full data frame with presence/absence
#   pred_vars = c(
#     "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
#     "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
#     "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempmin_bccm_stnd"),  # predictor columns
#   folds = dat2$fold,            # pre-defined fold column
#   eta_vals = c(0.01, 0.05, 0.1),   # learning rates to try
#   max_depth_vals = c(3, 5, 7),     # tree depths to try
#   subsample_vals = c(0.7, 1),      # row subsampling
#   colsample_bytree_vals = c(0.7, 1), # column subsampling
#   nrounds = 5000,                  # max trees
#   early_stop = 50,                 # early stopping rounds
#   n_cores = 4                      # number of parallel cores to use
# )
# 
# 
# #eelgrass and surgrass
# # Because multiple parameter combinations produced nearly identical AUC values, we selected a slightly less complex model (max_depth = 5) to reduce overfitting and improve model generalization.
# #Chose these parameters
# #eta = 0.05
# #max_depth = 5
# #subsample = 0.7
# #colsample_bytree = 0.7
# 



###############################################################################
# MODEL CONFIGURATION
###############################################################################
model_configs <- list(
  eelgrass_bccm = list(
    species = "ZO",
    fold_col = "fold_eelgrass",
    pred_vars = c(
      "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
      "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
      "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempmin_bccm_stnd"
    ),
    resp_name = "eelgrass",
    ocean_model = "bccm"
  ),
  eelgrass_nep = list(
    species = "ZO",
    fold_col = "fold_eelgrass",
    pred_vars = c(
      "depth_stnd", "slope_stnd", "rei_stnd", "substrate",
      "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
      "saltcv_nep_stnd", "tempcv_nep_stnd"
    ),
    resp_name = "eelgrass",
    ocean_model = "nep"
  ),
  surfgrass_bccm = list(
    species = "PH",
    fold_col = "fold_seagrass",
    pred_vars = c(
      "depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd",
      "substrate", "cul_eff_stnd", "airtempcv_stnd",
      "rsdsmin_stnd", "prmean_stnd", "saltcv_bccm_stnd",
      "NO3_bccm_stnd", "tempmin_bccm_stnd", "surftempcv_bccm_stnd"
    ),
    resp_name = "surfgrass",
    ocean_model = "bccm"
  ),
  surfgrass_nep = list(
    species = "PH",
    fold_col = "fold_seagrass",
    pred_vars = c(
      "depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd",
      "substrate", "cul_eff_stnd", "airtempcv_stnd",
      "rsdsmin_stnd", "prmean_stnd", "saltmean_nep_stnd",
      "tempmean_nep_stnd", "surftempcv_nep_stnd"
    ),
    resp_name = "surfgrass",
    ocean_model = "nep"
  )
)

ml_models <- c("GBM", "RF", "XGBOOST", "ANN", "CTA")

###############################################################################
# MAIN LOOP
###############################################################################
all_results <- list()

for (config_name in names(model_configs)) {
  
  message("\n========== Running: ", config_name, " ==========\n")
  config <- model_configs[[config_name]]
  
  # Filter and Prepare Data
  data <- seagrass_data_long %>%
    filter(species == config$species) %>%
    rename(fold = !!sym(config$fold_col))
  
  # BIOMOD2 CV
  message("Running BIOMOD2 CV...")
  biomod_results <- run_biomod_cv(data, config, ml_models)
  
  # Independent GBM CV
  message("Running independent GBM CV...")
  gbm_cv_results <- run_gbm_cv(
    dat = biomod_results$data,
    predictors = config$pred_vars,
    fold_col = "fold"
  )
  
  gbm_threshold <- gbm_cv_results$mean_threshold
  
  # Independent XGBoost CV
  message("Running independent XGBoost CV...")
  xgb_cv_results <- run_xgb_cv(
    dat = biomod_results$data,
    predictors = config$pred_vars,
    response = "presence",
    fold_col = "fold"
  )
  xgb_threshold <- xgb_cv_results$mean_threshold
  
  # Temporal Forecast
  forecast_results <- run_temporal_forecast(
    dat = biomod_results$data,
    pred_vars = config$pred_vars,
    config = config,
    biomod_thresholds = biomod_results$tss_thresholds,
    gbm_threshold = gbm_threshold,
    xgb_threshold = xgb_threshold
  )
  
  # Aggregate Metrics
  cv_metrics_combined <- biomod_results$auc_metrics %>%
    left_join(biomod_results$tss_thresholds, by = "algo") %>%
    bind_rows(gbm_cv_results) %>%
    bind_rows(xgb_cv_results)
  
  # Final Results Storage
  all_results[[config_name]] <- list(
    config = config,
    cv_metrics = cv_metrics_combined,
    forecast_results = forecast_results,
    var_importance = biomod_results$var_importance,
    biomod_model = biomod_results$biomod_out,
    data = biomod_results$data
  )
  
  # Save Metrics
  save(cv_metrics_combined, forecast_results,
       file = paste0("code/output_data/model_results/", config_name, "_all_metrics.RData"))
  
  message("Completed: ", config_name)
}

# Save all results combined
save(all_results, file = "code/output_data/model_results/all_biomod_results_combined.RData")


###############################################################################
# SUMMARY TABLES
###############################################################################
message("Building summary tables...")
# CV Summary
cv_summary <- bind_rows(
  lapply(names(all_results), function(nm) {
    all_results[[nm]]$cv_metrics %>%
      mutate(
        model_config = nm,
        species = all_results[[nm]]$config$resp_name,
        ocean_model = all_results[[nm]]$config$ocean_model
      )
  })
)
# Forecast Summary
forecast_summary <- bind_rows(
  lapply(names(all_results), function(nm) {
    fr <- all_results[[nm]]$forecast_results
    
    data.frame(
      model_config = nm,
      species = all_results[[nm]]$config$resp_name,
      ocean_model = all_results[[nm]]$config$ocean_model,
      
      model = names(fr),
      AUC_forecast = sapply(fr, function(x) x$AUC),
      Tjur_forecast = sapply(fr, function(x) x$Tjur),
      RMSE_forecast = sapply(fr, function(x) x$RMSE),
      TSS_forecast = sapply(fr, function(x) x$TSS),
      threshold_used = sapply(fr, function(x) x$threshold)
    )
  })
)


save(cv_summary, file = "code/output_data/model_results/seagrass_cv_metrics_biomod2.RData")
save(forecast_summary, file = "code/output_data/model_results/seagrass_forecast_metrics_biomod2.RData")





