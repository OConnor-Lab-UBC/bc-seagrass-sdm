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
  
  # Filter and prepare data
  data <- seagrass_data_long %>%
    filter(species == config$species) %>%
    rename(fold = !!sym(config$fold_col))
  
  # 1. Run BIOMOD2 CV
  message("Running BIOMOD2 CV...")
  biomod_results <- run_biomod_cv(data, config, ml_models)
  
  # 2. Independent GBM CV
  message("Running independent GBM CV...")
  gbm_cv_results <- run_gbm_cv(
    dat = biomod_results$data,
    predictors = config$pred_vars,
    fold_col = "fold"
  )
  gbm_threshold <- gbm_cv_results$mean_threshold
  
  # 3. Temporal forecast
  message("Running temporal forecast...")
  forecast_results <- run_temporal_forecast(
    dat = biomod_results$data,
    pred_vars = config$pred_vars,
    config = config,
    biomod_thresholds = biomod_results$tss_thresholds,
    gbm_threshold = gbm_threshold
  )
  
  # Combine metrics
  cv_metrics_combined <- biomod_results$auc_metrics %>%
    left_join(biomod_results$tss_thresholds, by = "algo") %>%
    bind_rows(gbm_cv_results)
  
  all_results[[config_name]] <- list(
    config = config,
    cv_metrics = cv_metrics_combined,
    forecast_results = forecast_results,
    var_importance = biomod_results$var_importance,
    biomod_model = biomod_results$biomod_out,
    data = biomod_results$data
  )
  
  save(
    cv_metrics_combined,
    forecast_results,
    file = paste0("code/output_data/model_results/", config_name, "_all_metrics.RData")
  )
  
  message("Completed: ", config_name)
}
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



#### make projection of GBM eelgrass models only
#this was not working in biomod!
load("code/output_data/prediction_model_inputs.RData")

#BCCM model
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



#NEP model
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



