###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# internal validation metrics table

###############################################################################
# packages 
library(terra)
library(ecospat)
library(Metrics)
library(dplyr)
library(purrr)
library(tidyr)


# create metric tables

# eelgrass ml cv metrics
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/eelgrass_eval_cv_gbm_nep.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/eelgrass_eval_cv_gbm_bccm.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/eelgrass_eval_cv_xgb_nep.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/eelgrass_eval_cv_xgb_bccm.RData")

#eelgrass forecast metrics
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/forecast_ml_eelgrass_models.RData")

# eelgrass regression cv metrics
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/eelgrass_eval_cv.RData")
names(eval_cv_list) <- c("GLM_bccm", "GLMM_bccm", "GLM_nep", "GLMM_nep")

# eelgrass regression forecast metrics
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/forecast_eelgrass_models.RData")

spatial_validation_eelgrass <- map_dfr(names(eval_cv_list), function(n) {
  s <- eval_cv_list[[n]]$summary
  data.frame(
    model = n,
    validation = "spatial",
    auc  = s$mean_test_auc,
    tjur = s$mean_test_tjur,
    brier = s$mean_test_brier,
    logloss = s$mean_test_logloss,
    sensitivity = s$mean_test_sensitivity,
    specificity = s$mean_test_specificity,
    tss  = s$mean_test_tss,
    threshold = s$mean_threshold)
})

gbm_bccm <- gbm_eel_bccm$summary %>%
  mutate(model = "GBM_bccm") %>%
  rename(auc  = mean_test_auc,
         tjur = mean_test_tjur,
         brier = mean_test_brier,
         logloss = mean_test_logloss,
         sensitivity = mean_test_sensitivity,
         specificity = mean_test_specificity,
         tss  = mean_test_tss,
         threshold = mean_threshold) %>%
  select(-mean_train_auc, -mean_train_tjur, -mean_train_brier, -mean_train_logloss, -mean_train_sensitivity, -mean_train_specificity, -mean_train_tss)

gbm_nep <- gbm_eel_nep$summary %>%
  mutate(model = "GBM_nep") %>%
  rename(auc  = mean_test_auc,
         tjur = mean_test_tjur,
         brier = mean_test_brier,
         logloss = mean_test_logloss,
         sensitivity = mean_test_sensitivity,
         specificity = mean_test_specificity,
         tss  = mean_test_tss,
         threshold = mean_threshold) %>%
  select(-mean_train_auc, -mean_train_tjur, -mean_train_brier, -mean_train_logloss, -mean_train_sensitivity, -mean_train_specificity, -mean_train_tss)

xgb_bccm <- xgb_eel_bccm$summary %>%
  mutate(model = "XGBoost_bccm") %>%
  rename(auc  = mean_test_auc,
         tjur = mean_test_tjur,
         brier = mean_test_brier,
         logloss = mean_test_logloss,
         sensitivity = mean_test_sensitivity,
         specificity = mean_test_specificity,
         tss  = mean_test_tss,
         threshold = mean_threshold) %>%
  select(-mean_train_auc, -mean_train_tjur, -mean_train_brier, -mean_train_logloss, -mean_train_sensitivity, -mean_train_specificity, -mean_train_tss)

xgb_nep <- xgb_eel_nep$summary %>%
  mutate(model = "XGBoost_nep") %>%
  rename(auc  = mean_test_auc,
         tjur = mean_test_tjur,
         brier = mean_test_brier,
         logloss = mean_test_logloss,
         sensitivity = mean_test_sensitivity,
         specificity = mean_test_specificity,
         tss  = mean_test_tss,
         threshold = mean_threshold) %>%
  select(-mean_train_auc, -mean_train_tjur, -mean_train_brier, -mean_train_logloss, -mean_train_sensitivity, -mean_train_specificity, -mean_train_tss)


spatial_validation_eelgrass_biomod <- xgb_nep %>% rbind(xgb_bccm, gbm_nep, gbm_bccm) %>% mutate(validation = "spatial") 

spatial_validation_eelgrass_biomod <- spatial_validation_eelgrass_biomod %>%
  select(model, validation, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold)

spatial_validation_eelgrass <- spatial_validation_eelgrass %>% rbind(spatial_validation_eelgrass_biomod)

forecast_table_biomod <- forecast_summary_ml_eelgrass %>%
  mutate(
    model = recode(model,
                   "XGBOOST_NEP" = "XGBoost_nep",
                   "XGBOOST_BCCM" = "XGBoost_bccm",
                   "GBM_NEP" = "GBM_nep",
                   "GBM_BCCM" = "GBM_bccm"),
    validation = "temporal",
    auc = AUC,
    tjur = TjurR2,
    brier = Brier,
    logloss = LogLoss,
    tss = TSS, 
    threshold = Threshold
  ) %>%
  select(model, validation, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold)



temporal_validation_eelgrass <- forecast_by_model %>%
  mutate(model = recode(model,
                        "BCCM_no_spatial" = "GLM_bccm",
                        "BCCM_spatial" = "GLMM_bccm",
                        "NEP_no_spatial" = "GLM_nep",
                        "NEP_spatial" = "GLMM_nep"))%>%
  rename(
    auc = AUC,
    tjur = TjurR2,
    tss = TSS,
    brier = Brier,
    logloss = LogLoss,
    threshold = Threshold
  ) %>%
  mutate (validation = "temporal") %>%
  select(model, validation, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold)

temporal_validation_eelgrass <- temporal_validation_eelgrass %>% rbind(forecast_table_biomod)

# Add suffixes to make clear which validation type
spatial_metrics_eelgrass <- spatial_validation_eelgrass %>%
  select(model, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold) %>%
  rename_with(~ paste0(., "_spatial"), -model)

temporal_metrics_eelgrass <- temporal_validation_eelgrass %>%
  select(model, validation, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold) %>%
  rename_with(~ paste0(., "_temporal"), -model)

# Join by model
combined_metrics_eelgrass <- full_join(spatial_metrics_eelgrass, temporal_metrics_eelgrass, by = "model")

save(combined_metrics_eelgrass, file = "code/output_data/model_results/internal_metrics_eelgrass.RData")




# HAVE NOT UPDATED YET FOR SURFGRASS!!!! Need to update all code similar to eelgrass
# surfgrass
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/eval_cv_surfgrass.RData")
names(eval_cv_list_surfgrass) <- c("bccm_nospatial", "bccm_spatial", "nep_nospatial", "nep_spatial")

load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/forecast_surfgrass_models.RData")

spatial_validation_surfgrass <- map_dfr(names(eval_cv_list_surfgrass), function(n) {
  s <- eval_cv_list_surfgrass[[n]]$summary
  data.frame(
    model = n,
    validation = "spatial",
    auc  = s$mean_test_auc,
    tjur = s$mean_test_tjur,
    rmse = s$mean_test_rmse,
    tss  = s$mean_test_tss,
    threshold = s$mean_threshold)
})

spatial_validation_surfgrass_biomod <- cv_summary %>%
  filter(species == "surfgrass") %>%
  mutate(
    algo_clean = recode(algo,
                        "GBM_robust" = "GBM",
                        "XGBOOST_robust" = "XGBOOST"
    ),
    
    model = paste0(algo_clean, "_", ocean_model),
    validation = "spatial",
    auc = coalesce(AUC_mean, mean_AUC),
    tjur = mean_Tjur,
    tss = coalesce(TSS_mean, mean_TSS),
    rmse = mean_RMSE,
    threshold = coalesce(mean_threshold, threshold_mean)
  ) %>%
  select(model, validation, auc, tjur, rmse, tss, threshold)

spatial_validation_surfgrass <- spatial_validation_surfgrass %>% rbind(spatial_validation_surfgrass_biomod)

forecast_table_biomod <- forecast_summary %>%
  filter(species == "surfgrass") %>%
  mutate(
    model = paste0(model, "_", ocean_model),
    validation = "temporal",
    auc = AUC_forecast,
    tjur = Tjur_forecast,
    rmse = RMSE_forecast,
    tss = TSS_forecast,
    threshold = threshold_used
  ) %>%
  select(model, validation, auc, tjur, rmse, tss, threshold)
rownames(forecast_table_biomod) <- NULL

temporal_validation_surfgrass <- forecast_results_all %>%
  filter(dataset == "forecast") %>%
  mutate(model = recode(model,
                        "BCCM_no_spatial" = "bccm_nospatial",
                        "BCCM_spatial" = "bccm_spatial",
                        "NEP_no_spatial" = "nep_nospatial",
                        "NEP_spatial" = "nep_spatial"))%>%
  rename(
    auc = AUC,
    tjur = TjurR2,
    rmse = RMSE,
    tss = TSS,
    threshold = threshold_used
  ) %>%
  mutate (validation = "temporal") %>%
  select(model, validation, auc, tjur, rmse, tss, threshold)

temporal_validation_surfgrass <- temporal_validation_surfgrass %>% rbind(forecast_table_biomod)

# Add suffixes to make clear which validation type
spatial_metrics_surfgrass <- spatial_validation_surfgrass %>%
  select(model, auc, tjur, rmse, tss, threshold) %>%
  rename_with(~ paste0(., "_spatial"), -model)

temporal_metrics_surfgrass <- temporal_validation_surfgrass %>%
  select(model, auc, tjur, rmse, tss) %>%
  rename_with(~ paste0(., "_temporal"), -model)

# Join by model
combined_metrics_surfgrass <- full_join(spatial_metrics_surfgrass, temporal_metrics_surfgrass, by = "model")






