###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# Compare independent data sets to sdms

###############################################################################
# packages 
library(terra)
library(ecospat)
library(Metrics)
library(dplyr)
library(purrr)
library(tidyr)
library(yardstick)

#### Load modelling functions ####
source("code/modelling-functions.R")



model_names <- c("GLM_bccm", "GLMM_bccm", "GLM_nep", "GLMM_nep", "GBM_bccm", "GBM_nep", "XGBoost_bccm", "XGBoost_nep")
#eelgrass_threshold_df <- combined_metrics_eelgrass %>% select(model, threshold_spatial) %>% rename(mean_threshold = threshold_spatial)
#surfgrass_threshold_df <- combined_metrics_surfgrass %>% select(model, threshold_spatial) %>% rename(mean_threshold = threshold_spatial)

#### Eelgrass ####
#load netforce eelgrass  data 2013-2023 that has been thinned
load("code/output_data/model_results/internal_metrics_eelgrass.RData")
eelgrass_obs <- vect("code/output_data/independent_validation/thinned_netforce_2013_2023_eelgrassobs.shp")
psuedoabsences <- vect("code/output_data/independent_validation/pseudoabsences_all.shp")

#### get sdm predictions
# get predictions
parent_dir <- "raster/eelgrass"
folders <- list.dirs(parent_dir, recursive = FALSE)

# load and mosaic eelgrass SDM rasters from 4 different SDM predictions
mosaics <- list()

for (f in folders) {
  tif_files <- list.files(
    f,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  pred_files <- tif_files[!grepl("se", basename(tif_files), ignore.case = TRUE)]
  folder_name <- basename(f)
  mosaics[[folder_name]] <- list()
  # mosaic prediction files (no se)
  if (length(pred_files) > 0) {
    r_pred <- lapply(pred_files, rast)
    mosaics[[folder_name]]$pred <- do.call(mosaic, r_pred)
  }
}

eelgrass_stack <- c(
  mosaics[["bccm_nospatial"]]$pred,
  mosaics[["bccm_spatial"]]$pred,
  mosaics[["gbm_bccm"]]$pred,
  mosaics[["gbm_nep"]]$pred,
  mosaics[["nep_nospatial"]]$pred,
  mosaics[["nep_spatial"]]$pred,
  mosaics[["xgb_bccm"]]$pred,
  mosaics[["xgb_nep"]]$pred
)
names(eelgrass_stack) <- c(
  "GLM_bccm",
  "GLMM_bccm",
  "GBM_bccm",
  "GBM_nep",
  "GLM_nep",
  "GLMM_nep",
  "XGBoost_bccm",
  "XGBoost_nep"
)

#### prepare selected validation points
eelgrass_obs_sf <- sf::st_as_sf(eelgrass_obs) %>%
  mutate(
    obs = 1,
    pa_set = NA_integer_
  )
pseudoabsences_sf <- sf::st_as_sf(psuedoabsences) %>%
  mutate(obs = 0)
if (!"pa_set" %in% names(pseudoabsences_sf)) {
  pseudoabsences_sf$pa_set <- NA_integer_
}
validation_points <- bind_rows(
  eelgrass_obs_sf %>%
    select(obs, pa_set, geometry),
  pseudoabsences_sf %>%
    select(obs, pa_set, geometry)
)


#### project points to prediction CRS if needed
validation_points <- sf::st_transform(
  validation_points,
  crs = terra::crs(eelgrass_stack)
)

#### extract predictions at selected points only
prediction_extract <- terra::extract(
  eelgrass_stack,
  terra::vect(validation_points)
) %>%
  select(-ID)
eelgrass_sf <- validation_points %>%
  bind_cols(prediction_extract) %>%
  filter(if_all(all_of(names(eelgrass_stack)), ~ !is.na(.)))
summary(eelgrass_sf)
#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"

#extract substrate data to see if that is causing mismatch
substrate_extract <- terra::extract(substrate_all, eelgrass_sf)
eelgrass_sf$substrate <- substrate_extract$substrate
eelgrass_sf$substrate <- c("Rock", "Mixed", "Sand", "Mud")[eelgrass_sf$substrate]

eelgrass_sf$rock_group <- ifelse(eelgrass_sf$substrate == "Rock", "Rock", "Not Rock")

#load depth layers, 
bathy_all <- terra::vrt(c("raw_data/envlayers-20m-hg//bathymetry.tif", "raw_data/envlayers-20m-ncc/bathymetry.tif", "raw_data/envlayers-20m-qcs/bathymetry.tif", "raw_data/envlayers-20m-wcvi/bathymetry.tif", "raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif"), "bathy.vrt", overwrite=T)
bathy_extract <- terra::extract(bathy_all, eelgrass_sf)
eelgrass_sf$bathy <- bathy_extract$bathy

# check to see overlap between psudoabsence data sets.
#Very little overlap, only about 0.4% of locations reused across two sets.
summary_table <- bind_rows(
  eelgrass_sf %>%
    filter(obs == 1) %>%
    summarise(
      Group = "Observations",
      pa_set = NA_integer_,
      n = n()
    ),
  eelgrass_sf %>%
    filter(obs == 0) %>%
    count(pa_set, name = "n") %>%
    mutate(Group = "Pseudoabsences")
)

summary_table

pa <- eelgrass_sf %>%
  filter(obs == 0)

n_total <- nrow(pa)
n_unique <- nrow(distinct(pa, geometry))

c(
  total_points = n_total,
  unique_locations = n_unique
)

location_overlap <- eelgrass_sf %>%
  filter(obs == 0) %>%
  count(geometry, name = "n_sets") %>%
  count(n_sets)

location_overlap

eelgrass_df <- sf::st_drop_geometry(eelgrass_sf)

model_cols <- c(
  "GLM_bccm",
  "GLMM_bccm",
  "GBM_bccm",
  "GBM_nep",
  "GLM_nep",
  "GLMM_nep",
  "XGBoost_bccm",
  "XGBoost_nep"
)
validation_metrics <- calc_independent_metrics(
  data = eelgrass_df,
  obs_col = "obs",
  pa_set_col = "pa_set",
  model_cols = model_cols,
  threshold_method = "max_tss"
)
validation_metrics

eelgrass_indep_metrics_summary <- validation_metrics %>%
  group_by(model) %>%
  summarise(
    auc_sd = sd(auc, na.rm = TRUE),
    auc = mean(auc, na.rm = TRUE),
    tjur_sd = sd(tjur, na.rm = TRUE),
    tjur = mean(tjur, na.rm = TRUE),
    brier_sd = sd(brier, na.rm = TRUE),
    brier = mean(brier, na.rm = TRUE),
    logloss_sd = sd(logloss, na.rm = TRUE),
    logloss = mean(logloss, na.rm = TRUE),
    sensitivity_sd = sd(sensitivity, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    specificity_sd = sd(specificity, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    tss_sd = sd(tss, na.rm = TRUE),
    tss = mean(tss, na.rm = TRUE),
    threshold = mean(threshold, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  select (model, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold, auc_sd, tjur_sd, brier_sd, logloss_sd, sensitivity_sd, specificity_sd, tss_sd)
eelgrass_indep_metrics_summary
# log loss values are high becasue there are many psuedoabsences in areas predicted suitable for eelgrass

independent_metrics_eelgrass <- eelgrass_indep_metrics_summary %>%
  rename_with(~ paste0(., "_independent"), -model)


combined_metrics_eelgrass <- full_join(combined_metrics_eelgrass, independent_metrics_eelgrass, by = "model")

save(combined_metrics_eelgrass, file = "code/output_data/model_results/metrics_eelgrass.RData")


save(eelgrass_indep_metrics_summary, validation_metrics, eelgrass_df, eelgrass_sf, file = "code/output_data/model_results/eelgrass_independent_eval.RData")







#### Surfgrass
#load surfgrass independent data
# these use the shorezone line data extrapolated to the depth ribbon, so likely some overestimation of surfgrass distirbution

load("code/output_data/model_results/internal_metrics_surfgrass.RData")


surfgrass_obs <- vect("code/output_data/independent_validation/thinned_2013_2024_surfgrassobs.shp")
psuedoabsences <- vect("code/output_data/independent_validation/pseudoabsences_all_surfgrass.shp")

#### get sdm predictions
# get predictions
parent_dir <- "raster/surfgrass"
folders <- list.dirs(parent_dir, recursive = FALSE)

# load and mosaic surfgrass SDM rasters from 4 different SDM predictions
mosaics <- list()

for (f in folders) {
  tif_files <- list.files(
    f,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  pred_files <- tif_files[!grepl("se", basename(tif_files), ignore.case = TRUE)]
  folder_name <- basename(f)
  mosaics[[folder_name]] <- list()
  # mosaic prediction files (no se)
  if (length(pred_files) > 0) {
    r_pred <- lapply(pred_files, rast)
    mosaics[[folder_name]]$pred <- do.call(mosaic, r_pred)
  }
}

surfgrass_stack <- c(
  mosaics[["bccm_nospatial"]]$pred,
  mosaics[["bccm_spatial"]]$pred,
  mosaics[["gbm_bccm"]]$pred,
  mosaics[["gbm_nep"]]$pred,
  mosaics[["nep_nospatial"]]$pred,
  mosaics[["nep_spatial"]]$pred,
  mosaics[["xgb_bccm"]]$pred,
  mosaics[["xgb_nep"]]$pred
)
names(surfgrass_stack) <- c(
  "GLM_bccm",
  "GLMM_bccm",
  "GBM_bccm",
  "GBM_nep",
  "GLM_nep",
  "GLMM_nep",
  "XGBoost_bccm",
  "XGBoost_nep"
)

#### prepare selected validation points
surfgrass_obs_sf <- sf::st_as_sf(surfgrass_obs) %>%
  mutate(
    obs = 1,
    pa_set = NA_integer_
  )
pseudoabsences_sf <- sf::st_as_sf(psuedoabsences) %>%
  mutate(obs = 0)
if (!"pa_set" %in% names(pseudoabsences_sf)) {
  pseudoabsences_sf$pa_set <- NA_integer_
}
validation_points <- bind_rows(
  surfgrass_obs_sf %>%
    select(obs, pa_set, geometry),
  pseudoabsences_sf %>%
    select(obs, pa_set, geometry)
)


#### project points to prediction CRS if needed
validation_points <- sf::st_transform(
  validation_points,
  crs = terra::crs(surfgrass_stack)
)

#### extract predictions at selected points only
prediction_extract <- terra::extract(
  surfgrass_stack,
  terra::vect(validation_points)
) %>%
  select(-ID)
surfgrass_sf <- validation_points %>%
  bind_cols(prediction_extract) %>%
  filter(if_all(all_of(names(surfgrass_stack)), ~ !is.na(.)))
summary(surfgrass_sf)
#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"

#extract substrate data to see if that is causing mismatch
substrate_extract <- terra::extract(substrate_all, surfgrass_sf)
surfgrass_sf$substrate <- substrate_extract$substrate
surfgrass_sf$substrate <- c("Rock", "Mixed", "Sand", "Mud")[surfgrass_sf$substrate]

surfgrass_sf$rock_group <- ifelse(surfgrass_sf$substrate == "Rock", "Rock", "Not Rock")

#load depth layers, 
bathy_all <- terra::vrt(c("raw_data/envlayers-20m-hg//bathymetry.tif", "raw_data/envlayers-20m-ncc/bathymetry.tif", "raw_data/envlayers-20m-qcs/bathymetry.tif", "raw_data/envlayers-20m-wcvi/bathymetry.tif", "raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif"), "bathy.vrt", overwrite=T)
bathy_extract <- terra::extract(bathy_all, surfgrass_sf)
surfgrass_sf$bathy <- bathy_extract$bathy

# check to see overlap between psudoabsence data sets.
#Very little overlap, only about 0.4% of locations reused across two sets.
summary_table <- bind_rows(
  surfgrass_sf %>%
    filter(obs == 1) %>%
    summarise(
      Group = "Observations",
      pa_set = NA_integer_,
      n = n()
    ),
  surfgrass_sf %>%
    filter(obs == 0) %>%
    count(pa_set, name = "n") %>%
    mutate(Group = "Pseudoabsences")
)

summary_table

pa <- surfgrass_sf %>%
  filter(obs == 0)

n_total <- nrow(pa)
n_unique <- nrow(distinct(pa, geometry))

c(
  total_points = n_total,
  unique_locations = n_unique
)

location_overlap <- surfgrass_sf %>%
  filter(obs == 0) %>%
  count(geometry, name = "n_sets") %>%
  count(n_sets)

location_overlap

surfgrass_df <- sf::st_drop_geometry(surfgrass_sf)

model_cols <- c(
  "GLM_bccm",
  "GLMM_bccm",
  "GBM_bccm",
  "GBM_nep",
  "GLM_nep",
  "GLMM_nep",
  "XGBoost_bccm",
  "XGBoost_nep"
)
validation_metrics <- calc_independent_metrics(
  data = surfgrass_df,
  obs_col = "obs",
  pa_set_col = "pa_set",
  model_cols = model_cols,
  threshold_method = "max_tss"
)
validation_metrics

surfgrass_indep_metrics_summary <- validation_metrics %>%
  group_by(model) %>%
  summarise(
    auc_sd = sd(auc, na.rm = TRUE),
    auc = mean(auc, na.rm = TRUE),
    tjur_sd = sd(tjur, na.rm = TRUE),
    tjur = mean(tjur, na.rm = TRUE),
    brier_sd = sd(brier, na.rm = TRUE),
    brier = mean(brier, na.rm = TRUE),
    logloss_sd = sd(logloss, na.rm = TRUE),
    logloss = mean(logloss, na.rm = TRUE),
    sensitivity_sd = sd(sensitivity, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    specificity_sd = sd(specificity, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    tss_sd = sd(tss, na.rm = TRUE),
    tss = mean(tss, na.rm = TRUE),
    threshold = mean(threshold, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  select (model, auc, tjur, brier, logloss, sensitivity, specificity, tss, threshold, auc_sd, tjur_sd, brier_sd, logloss_sd, sensitivity_sd, specificity_sd, tss_sd)
surfgrass_indep_metrics_summary
# log loss values are high becasue there are many psuedoabsences in areas predicted suitable for surfgrass

independent_metrics_surfgrass <- surfgrass_indep_metrics_summary %>%
  rename_with(~ paste0(., "_independent"), -model)


combined_metrics_surfgrass <- full_join(combined_metrics_surfgrass, independent_metrics_surfgrass, by = "model")

save(combined_metrics_surfgrass, file = "code/output_data/model_results/metrics_surfgrass.RData")


save(surfgrass_indep_metrics_summary, validation_metrics, surfgrass_df, surfgrass_sf, file = "code/output_data/model_results/surfgrass_independent_eval.RData")




