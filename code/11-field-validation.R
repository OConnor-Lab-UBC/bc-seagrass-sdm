###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# compare sdm predictions to sdm validations to choose best model
#
###############################################################################


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

#themes for figures
boxed_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title.y = element_text(margin = margin(r = 8))
  )

# load sdm validation dataset
load("code/output_data/field_validation/validation_dataset.RData")


# need to also drop sites that after reviewing the database should not have been included as they were aborted or vis was too bad to see species
validation_sf <- validation_sf %>%
  filter(!HKey %in% c(74, 84, 85, 161))

# noticed error in database that have alerted to Sandie so remove this once it is changed
validation_sf$PC_ZM[validation_sf$HKey == "125"] <- "51-75"
validation_sf$PC_PH[validation_sf$HKey == "125"] <- "0"

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
                                           ZM, ZM_freq, PC_ZM, eelgrass_bccm_nospatial_pred, eelgrass_bccm_nospatial_se, eelgrass_bccm_spatial_pred, eelgrass_bccm_spatial_se, eelgrass_nep_nospatial_pred, eelgrass_nep_nospatial_pred, eelgrass_nep_spatial_pred, eelgrass_nep_spatial_se,
                                           PH, PH_freq, PC_PH, geometry) # need to add in surfgrass predictions
summary(validation_sf)
# substrate update diff false is 130 and substrate diff is 175 false (shows that update with shorezone substrate data makes substrate more accurate)
# at least 26% of sites had substrate mismatch between modelled and obs during field validation study

# remove sites that we know will have wrong predictions becasue the environmental layers do not match reality
validation_sf <- validation_sf %>%
  filter(
    (Depth_diff < 5 | is.na(Depth_diff)),
    (Slope_diff < 20 | is.na(Slope_diff)), # this was just a guess, need to add something more thought out
    #(Substrate_diff == TRUE | is.na(Substrate_diff)), #leaving this out for now as that means more get eliminated for zm and ph then necessary
    (Substrate_diff == TRUE | is.na(Substrate_diff))
  )
# 334 passed threshold 
# 3 sites >20 degrees slope difference, 40 sites has >5m depth difference, 112 sites had substrate mismatch in updated substrate
