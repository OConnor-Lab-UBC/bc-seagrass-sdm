###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# Generate psuedo absences based on eelgrass presences in netforce data all years and dfo dive data
#
###############################################################################

library(sf)
library(terra)
library(dplyr)
library(purrr)
library(tidyr)
library(Metrics)


source("code/modelling-functions.R")
template_rast <- rast("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif")

eelgrass_domain <- terra::vrt(c("raster/eelgrass/xgb_nep/eelgrass_predictions_hg_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ncc_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_qcs_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ss_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_wcvi_xgb_nep.tif"), "eelgrass_xgb_nep.vrt", overwrite=T)   # values 0–1


# load independent eelgrass rasters
eelgrass_indep <- rast("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif")
values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1
eelgrass_indep2 <- rast("code/output_data/independent_validation/BCeelgrass_netforce_1974_2012.tif")
values(eelgrass_indep2)[values(eelgrass_indep2) >= 1] <- 1
# load model inputs
load("code/output_data/seagrass_model_inputs.RData")
data <- seagrass_data_long %>%
  filter(species == "ZO", presence == 1) %>%
  select(X_m, Y_m) %>%
  distinct()
data_sf <- st_as_sf(data, coords = c("X_m", "Y_m"), crs = 3005)
data_sf$pres <- 1
data_vect <- terra::vect(data_sf)
training_pres_rast <- terra::rasterize(
  x = data_vect,
  y = template_rast,
  field = "pres",
  touches = TRUE
)
training_pres_rast <- ifel(!is.na(training_pres_rast), 1, NA)
freq(training_pres_rast)

combined_exclusion <- terra::ifel(
  !is.na(training_pres_rast) | !is.na(eelgrass_indep) | !is.na(eelgrass_indep2),
  1,
  NA
)

freq(combined_exclusion)



terra::values(eelgrass_domain) <- ifelse(!is.na(terra::values(eelgrass_domain)), 1, NA)

combined_exclusion_resample<- terra::crop(
  terra::resample(combined_exclusion, eelgrass_domain, method = "bilinear"),
  eelgrass_domain)
  
n_pres <- sum(terra::values(combined_exclusion) > 0, na.rm = TRUE)
eelgrass_pa <- generate_pseudoabsences(
  domain_rast = eelgrass_domain,
  exclusion_rast = combined_exclusion,
  n_pa = n_pres,
  buffer_cells = 5,
  seed = 101
)

