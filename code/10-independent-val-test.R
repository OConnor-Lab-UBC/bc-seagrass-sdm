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

#### Load modelling functions ####
source("code/modelling-functions.R")

# load model thresholds created
load("code/output_data/model_results/eelgrass_eval_cv.RData")
model_names <- c("bccm_nospatial", "bccm_spatial", 
                 "nep_nospatial", "nep_spatial")
names(eval_cv_list) <- model_names

# gbm threshold
load("code/output_data/model_results/seagrass_cv_metrics_biomod2.RData")

spatial_thresholds <- sapply(model_names, function(x) {
  eval_cv_list[[x]]$summary$mean_threshold
})

spatial_threshold_df <- data.frame(
  model = model_names,
  mean_threshold = as.numeric(spatial_thresholds)
)

gbm_robust_df <- cv_summary[cv_summary$species == "eelgrass" & cv_summary$algo == "GBM_robust", ]
gbm_robust_df$model <- ifelse(gbm_robust_df$ocean_model == "bccm", "bccm_gbm", "nep_gbm")
gbm_robust_df <- gbm_robust_df[, c("model", "mean_threshold")]
eelgrass_threshold_df <- rbind(spatial_threshold_df, gbm_robust_df)

load("code/output_data/model_results/eval_cv_surfgrass.RData")
names(eval_cv_list_surfgrass) <- model_names
spatial_thresholds <- sapply(model_names, function(x) {
  eval_cv_list_surfgrass[[x]]$summary$mean_threshold
})

spatial_threshold_df <- data.frame(
  model = model_names,
  mean_threshold = as.numeric(spatial_thresholds)
)

gbm_robust_df <- cv_summary[cv_summary$species == "surfgrass" & cv_summary$algo == "GBM_robust", ]
gbm_robust_df$model <- ifelse(gbm_robust_df$ocean_model == "bccm", "bccm_gbm", "nep_gbm")
gbm_robust_df <- gbm_robust_df[, c("model", "mean_threshold")]
surfgrass_threshold_df <- rbind(spatial_threshold_df, gbm_robust_df)

model_names <- c("bccm_nospatial", "bccm_spatial", 
                 "nep_nospatial", "nep_spatial", "bccm_gbm", "nep_gbm")



#### Eelgrass ####
#load netforce eelgrass  data 2013-2023
eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))

# from looking previously there is no difference in the predictions if eelgrass has been observed more times so just change all to 1. A Spearman correlation of –0.0303 suggests that higher modelled probabilities do not correspond in any meaningful way to more years of observed presence.
values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1
names(eelgrass_indep)<-"obs"

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

# Ensure rasters align, just take model predictions where there is observed data
# we can only compare in areas where eelgrass has been observed as there is no absence data from independent dataset
mosaics_flat <- unlist(mosaics, recursive = FALSE)

eelgrass_sdm_resampled <- lapply(mosaics_flat, function(r) {
  if (!compareGeom(r, eelgrass_indep, crs = TRUE, stopOnError = FALSE)) {
    r <- project(r, eelgrass_indep)  }
  crop(
    resample(r, eelgrass_indep, method = "bilinear"),
    eelgrass_indep)})
#plot(eelgrass_sdm_resampled[[1]])
observed_cells <- which(values(eelgrass_indep) > 0)
coords <- xyFromCell(eelgrass_indep, observed_cells)

eelgrass_stack <- c(eelgrass_sdm_resampled[[1]], eelgrass_sdm_resampled[[2]], eelgrass_sdm_resampled[[3]], eelgrass_sdm_resampled[[4]], eelgrass_sdm_resampled[[5]], eelgrass_sdm_resampled[[6]])
names(eelgrass_stack) <- c("bccm_nospatial", "bccm_spatial", "bccm_gbm", "nep_gbm", "nep_nospatial", "nep_spatial")

r_points <- as.points(eelgrass_indep, na.rm = TRUE)
eelgrass_sf<- r_points %>% st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

prediction_extract <- terra::extract(eelgrass_stack, eelgrass_sf)
eelgrass_sf <- eelgrass_sf %>% bind_cols(prediction_extract)

eelgrass_sf <- eelgrass_sf %>%
  filter(if_all(everything(), ~ !is.na(.)))

#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"

#extract substrate data to see if that is causing mismatch
substrate_extract <- terra::extract(substrate_all, eelgrass_sf)
eelgrass_sf$substrate <- substrate_extract$substrate
eelgrass_sf$substrate <- c("Rock", "Mixed", "Sand", "Mud")[eelgrass_sf$substrate]

eelgrass_sf$rock_group <- ifelse(eelgrass_sf$substrate == "Rock", "Rock", "Not Rock")
# should remove areas with modelled rock substrate from comparison as the substrate model is likely wrong in these areas

eelgrass_df <- sf::st_drop_geometry(eelgrass_sf)

eelgrass_independent_results <- evaluate_independent_seagrass(
  independent = eelgrass_df,
  model_names = model_names,
  cv_thresholds = eelgrass_threshold_df,
  raster_stack = eelgrass_stack
)

#The best eelgrass model is bccm_gbm, followed closely by nep_gbm, with both substantially outperforming the sdm models based on Boyce index and suitability predictions at observed eelgrass locations.

save(eelgrass_independent_results, file = "code/output_data/model_results/eelgrass_independent_eval.RData")



#### Surfgrass
#load surfgrass independent data
# these use the shorezone line data extrapolated to the depth ribbon, so likely some overestimation of surfgrass distirbution
surfgrass_indep <- rast(c("code/output_data/independent_validation/surfgrass_validation_raster_2013_2024.tif"))

values(surfgrass_indep)[values(surfgrass_indep) >= 1] <- 1
names(surfgrass_indep)<-"obs"


#### get sdm predictions

# get predictions
parent_dir <- "raster/surfgrass"
folders <- list.dirs(parent_dir, recursive = FALSE)

# load and mosaic SDM rasters from 4 different SDM predictions
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

# Ensure rasters align, just take model predictions where there is observed data
# we can only compare in areas where surfgrass has been observed as there is no absence data from independent dataset
mosaics_flat <- unlist(mosaics, recursive = FALSE)

surfgrass_sdm_resampled <- lapply(mosaics_flat, function(r) {
  if (!compareGeom(r, surfgrass_indep, crs = TRUE, stopOnError = FALSE)) {
    r <- project(r, surfgrass_indep)  }
  crop(
    resample(r, surfgrass_indep, method = "bilinear"),
    surfgrass_indep)})

observed_cells <- which(values(surfgrass_indep) > 0)
coords <- xyFromCell(surfgrass_indep, observed_cells)

surfgrass_stack <- c(surfgrass_sdm_resampled[[1]], surfgrass_sdm_resampled[[2]], surfgrass_sdm_resampled[[3]], surfgrass_sdm_resampled[[4]])
names(surfgrass_stack) <- c("bccm_nospatial", "bccm_spatial", "nep_nospatial", "nep_spatial")

r_points <- as.points(surfgrass_indep, na.rm = TRUE)
surfgrass_sf<- r_points %>% st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

prediction_extract <- terra::extract(surfgrass_stack, surfgrass_sf)
surfgrass_sf <- surfgrass_sf %>% bind_cols(prediction_extract)

surfgrass_sf <- surfgrass_sf %>%
  filter(if_all(everything(), ~ !is.na(.)))

surfgrass_df <- sf::st_drop_geometry(surfgrass_sf)

# this is not working!!!!
surfgrass_independent_results <- evaluate_independent_seagrass(
  independent = surfgrass_df,
  model_names = model_names,
  cv_thresholds = surfgrass_threshold_df,
  raster_stack = surfgrass_stack
)

#The best surfgrass model is 

save(surfgrass_independent_results, file = "code/output_data/model_results/surfgrass_independent_eval.RData")






