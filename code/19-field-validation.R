###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# compare sdm predictions to sdm validations sites to choose best model
#
###############################################################################

#functions
source("code/modelling-functions.R")

# Load packages
library(sf)
library(tidyverse)
library(terra)
library(patchwork)
library(effsize)
library(ggpubr)
library(pROC)
library(forcats)
library(caret)
library(purrr)
library(tibble)

#themes for figures
boxed_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title.y = element_text(margin = margin(r = 8))
  )

load("code/output_data/model_results/metrics_eelgrass.RData")
#load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/combined_metrics_surfgrass.RData")

# load sdm validation dataset
load("code/output_data/field_validation/validation_dataset.RData")
validation_sf <- validation_sf %>%
  relocate(PH_freq, PH, ZM_freq, ZM, .after = PC_PH)

# need to also drop sites that after reviewing the database should not have been included as they were aborted or vis was too bad to see species
validation_sf <- validation_sf %>%
  filter(!HKey %in% c(74, 84, 85, 161))

print(paste("zm present in ", round((sum(validation_sf$ZM)/nrow(validation_sf))*100,2), "% of observations", sep = ""))
print(paste("ph present in ", round((sum(validation_sf$PH)/nrow(validation_sf))*100,2), "% of observations", sep = ""))


# determine difference between observed and modelled environmental layers
# if either observed or modelled layers were mixed these were always retained as we are just looking at predominant substrate and most sites are actually mixed in reality
validation_sf <- validation_sf %>%
  mutate(Depth_diff = abs(bathy_mod - avgCorDepth_obs),
         Slope_diff = abs(slope_mod - slope_obs),
         Substrate_diff =
           if_else(
             substrate_mod == "Mixed" | substrate_mode_obs == "Mixed", TRUE,
             substrate_mod == substrate_mode_obs))

validation_sf <- validation_sf %>% select (HKey, Region, Survey, Visibility, High.vegetation, avgCorDepth_obs, bathy_mod,
                                           Depth_diff, slope_obs, slope_mod, Slope_diff, substrate_mode_obs, 
                                           substrate_mod, Substrate_diff, 
                                           ZM, ZM_freq, PC_ZM, zo_bccm_nospatial_pred, zo_bccm_spatial_pred, zo_nep_nospatial_pred, zo_nep_spatial_pred, zo_XGBOOST_bccm_pred, zo_XGBOOST_nep_pred, zo_GBM_bccm_pred, zo_GBM_nep_pred,
                                           PH, PH_freq, PC_PH, ph_bccm_nospatial_pred, ph_bccm_spatial_pred, ph_nep_nospatial_pred, ph_nep_spatial_pred, ph_XGBOOST_bccm_pred, ph_XGBOOST_nep_pred, ph_GBM_bccm_pred, ph_GBM_nep_pred, geometry) # need to add in surfgrass predictions
summary(validation_sf)
# substrate update diff false is 130 and substrate diff is 175 false (shows that update with shorezone substrate data makes substrate more accurate)
# at least 26% of sites had substrate mismatch between modelled and obs during field validation study

#not going to remove sites at this time as this issue would persist across the whole province. good to know that we are having a substrate mismatch >25% of the time
# remove sites that we know will have wrong predictions becasue the environmental layers do not match reality
# validation_sf <- validation_sf %>%
#   filter(
#     (Depth_diff < 5 | is.na(Depth_diff)),
#     (Slope_diff < 20 | is.na(Slope_diff)), # this was just a guess, need to add something more thought out
#     #(Substrate_diff == TRUE | is.na(Substrate_diff)), #leaving this out for now as that means more get eliminated for zm and ph then necessary
#     (Substrate_diff == TRUE | is.na(Substrate_diff))
#   )
# 334 passed threshold 
# 3 sites >20 degrees slope difference, 40 sites has >5m depth difference, 112 sites had substrate mismatch in updated substrate



##ZM Eelgrass
zm_metrics <- calc_field_metrics(
  data = validation_sf,
  obs_col = "ZM",
  model_cols = grep("^zo_.*_pred$", names(validation_sf), value = TRUE)
)

# test with no high intertida, models do much better indicating we have overprediction into the intertidal in the sdms
validation_sf_nointertidal <- validation_sf %>%
  dplyr::filter(
    Survey != "Intertidal" | bathy_mod > -1
  )

zm_metrics_nointertidal <- calc_field_metrics(
  data = validation_sf_nointertidal,
  obs_col = "ZM",
  model_cols = grep("^zo_.*_pred$", names(validation_sf), value = TRUE)
)

eelgrass_field_metrics_summary <- zm_metrics %>%
  mutate(model = recode(model,
                        "zo_bccm_nospatial_pred" = "GLM_bccm",
                        "zo_bccm_spatial_pred" = "GLMM_bccm",
                        "zo_nep_nospatial_pred" = "GLM_nep",
                        "zo_nep_spatial_pred" = "GLMM_nep",
                        "zo_XGBOOST_nep_pred" = "XGBoost_nep",
                        "zo_XGBOOST_bccm_pred" = "XGBoost_bccm",
                        "zo_GBM_nep_pred" = "GBM_nep",
                        "zo_GBM_bccm_pred" = "GBM_bccm"))%>%
  select(model, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold) %>%
  rename_with(~ paste0(., "_field"), -model)

final_metrics_eelgrass <- full_join(combined_metrics_eelgrass, eelgrass_field_metrics_summary, by = "model")
save(final_metrics_eelgrass, file = "code/output_data/model_results/final_metrics_eelgrass.RData")


eelgrass_field_metrics_summary_nointertidal <- zm_metrics_nointertidal %>%
  mutate(model = recode(model,
                        "zo_bccm_nospatial_pred" = "GLM_bccm",
                        "zo_bccm_spatial_pred" = "GLMM_bccm",
                        "zo_nep_nospatial_pred" = "GLM_nep",
                        "zo_nep_spatial_pred" = "GLMM_nep",
                        "zo_XGBOOST_nep_pred" = "XGBoost_nep",
                        "zo_XGBOOST_bccm_pred" = "XGBoost_bccm",
                        "zo_GBM_nep_pred" = "GBM_nep",
                        "zo_GBM_bccm_pred" = "GBM_bccm"))%>%
  select(model, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold) %>%
  rename_with(~ paste0(., "_field"), -model)

final_metrics_eelgrass_nointertidal <- full_join(combined_metrics_eelgrass, eelgrass_field_metrics_summary_nointertidal, by = "model")
save(final_metrics_eelgrass_nointertidal, file = "code/output_data/model_results/final_metrics_eelgrass_fieldnointertidal.RData")










#still need to update surfgrass
#### PH surfgrass

ph_metrics <- calc_field_metrics(
  data = field_validation,
  obs_col = "PH",
  model_cols = grep("^ph_.*_pred$", names(field_validation), value = TRUE)
)



