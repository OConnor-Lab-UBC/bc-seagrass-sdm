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

#### Load modelling functions ####
source("code/modelling-functions.R")

load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/seagrass_cv_metrics_biomod2.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/seagrass_forecast_metrics_biomod2.RData")

# eelgrass
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/eelgrass_eval_cv.RData")
names(eval_cv_list) <- c("bccm_nospatial", "bccm_spatial", "nep_nospatial", "nep_spatial")

load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/forecast_eelgrass_models.RData")

spatial_validation_eelgrass <- map_dfr(names(eval_cv_list), function(n) {
  s <- eval_cv_list[[n]]$summary
  data.frame(
    model = n,
    validation = "spatial",
    auc  = s$mean_test_auc,
    tjur = s$mean_test_tjur,
    rmse = s$mean_test_rmse,
    tss  = s$mean_test_tss,
    threshold = s$mean_threshold)
})

spatial_validation_eelgrass_biomod <- cv_summary %>%
  filter(species == "eelgrass") %>%
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

spatial_validation_eelgrass <- spatial_validation_eelgrass %>% rbind(spatial_validation_eelgrass_biomod)

forecast_table_biomod <- forecast_summary %>%
  filter(species == "eelgrass") %>%
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

temporal_validation_eelgrass <- forecast_results_all %>%
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

temporal_validation_eelgrass <- temporal_validation_eelgrass %>% rbind(forecast_table_biomod)

# Add suffixes to make clear which validation type
spatial_metrics_eelgrass <- spatial_validation_eelgrass %>%
  select(model, auc, tjur, rmse, tss, threshold) %>%
  rename_with(~ paste0(., "_spatial"), -model)

temporal_metrics_eelgrass <- temporal_validation_eelgrass %>%
  select(model, auc, tjur, rmse, tss) %>%
  rename_with(~ paste0(., "_temporal"), -model)

# Join by model
combined_metrics_eelgrass <- full_join(spatial_metrics_eelgrass, temporal_metrics_eelgrass, by = "model")

# drop CTA, ANN, RF they all rpeform horribly and not carried forward to independent and field
combined_metrics_eelgrass <- combined_metrics_eelgrass %>%
  filter(!model %in% c("RF_nep", "RF_bccm", "ANN_bccm", "ANN_nep", "CTA_bccm", "CTA_nep"))



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

# drop CTA, ANN, RF they all rpeform horribly and not carried forward to independent and field
combined_metrics_surfgrass <- combined_metrics_surfgrass %>%
  filter(!model %in% c("RF_nep", "RF_bccm", "ANN_bccm", "ANN_nep", "CTA_bccm", "CTA_nep"))

eelgrass_threshold_df <- combined_metrics_eelgrass %>% select(model, threshold_spatial) %>% rename(mean_threshold = threshold_spatial)
surfgrass_threshold_df <- combined_metrics_surfgrass %>% select(model, threshold_spatial) %>% rename(mean_threshold = threshold_spatial)

model_names <- c("bccm_nospatial", "bccm_spatial", 
                 "nep_nospatial", "nep_spatial", "GBM_bccm", "GBM_nep", "XGBOOST_bccm", "XGBOOST_nep")


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

eelgrass_stack <- c(eelgrass_sdm_resampled[[1]], eelgrass_sdm_resampled[[2]], eelgrass_sdm_resampled[[3]], eelgrass_sdm_resampled[[4]], eelgrass_sdm_resampled[[5]], eelgrass_sdm_resampled[[6]], eelgrass_sdm_resampled[[7]], eelgrass_sdm_resampled[[8]])
names(eelgrass_stack) <- c("bccm_nospatial", "bccm_spatial", "GBM_bccm", "GBM_nep", "nep_nospatial", "nep_spatial", "XGBOOST_bccm", "XGBOOST_nep")

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


save(eelgrass_independent_results, file = "code/output_data/model_results/eelgrass_independent_eval.RData")

eelgrass_independent_clean <- eelgrass_independent_results %>%
  rename(model = Model) %>%
  rename_with(~ paste0(tolower(.x), "_independent"),
              c(MPS, FPPS, FNR, CBI)) %>%
  select(-Threshold)

combined_metrics_eelgrass <- combined_metrics_eelgrass %>%
  left_join(eelgrass_independent_clean, by = "model")

save(combined_metrics_eelgrass, file = "code/output_data/model_results/combined_metrics_eelgrass.RData")

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

surfgrass_stack <- c(surfgrass_sdm_resampled[[1]], surfgrass_sdm_resampled[[2]], surfgrass_sdm_resampled[[3]], surfgrass_sdm_resampled[[4]], surfgrass_sdm_resampled[[5]], surfgrass_sdm_resampled[[6]], surfgrass_sdm_resampled[[7]], surfgrass_sdm_resampled[[8]])
names(surfgrass_stack) <- c("bccm_nospatial", "bccm_spatial", "GBM_bccm", "GBM_nep", "nep_nospatial", "nep_spatial", "XGBOOST_bccm", "XGBOOST_nep")

r_points <- as.points(surfgrass_indep, na.rm = TRUE)
surfgrass_sf<- r_points %>% sf::st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

prediction_extract <- terra::extract(surfgrass_stack, surfgrass_sf)
surfgrass_sf <- surfgrass_sf %>% bind_cols(prediction_extract)

surfgrass_sf <- surfgrass_sf %>%
  filter(if_all(everything(), ~ !is.na(.)))

#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"

#extract substrate data to see if that is causing mismatch
substrate_extract <- terra::extract(substrate_all, surfgrass_sf)
surfgrass_sf$substrate <- substrate_extract$substrate
surfgrass_sf$substrate <- c("Rock", "Mixed", "Sand", "Mud")[surfgrass_sf$substrate]

surfgrass_sf$rock_group <- ifelse(surfgrass_sf$substrate == "Rock", "Rock", "Not Rock")

# load bathy
bathy_all <- terra::vrt(c("raw_data/envlayers-20m-hg/bathymetry.tif", "raw_data/envlayers-20m-ncc/bathymetry.tif", "raw_data/envlayers-20m-qcs/bathymetry.tif", "raw_data/envlayers-20m-wcvi/bathymetry.tif", "raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif"), "bathy.vrt", overwrite=T)
bathy_extract <- terra::extract(bathy_all, surfgrass_sf)
surfgrass_sf$bathy <- bathy_extract$bathy

# Looksed at distirbution of dataset of observations from dive datasets and majoriity of presences of PH were above 5 m CD.
# to better align the "mapped" PH we will only retain cells that have bathy above 5 m. 

surfgrass_sf <- surfgrass_sf %>% filter(bathy <= 5)

surfgrass_df <- sf::st_drop_geometry(surfgrass_sf)

# this is not working!!!!
surfgrass_independent_results <- evaluate_independent_seagrass(
  independent = surfgrass_df,
  model_names = model_names,
  cv_thresholds = surfgrass_threshold_df,
  raster_stack = surfgrass_stack
)


save(surfgrass_independent_results, file = "code/output_data/model_results/surfgrass_independent_eval.RData")


surfgrass_independent_clean <- surfgrass_independent_results %>%
  rename(model = Model) %>%
  rename_with(~ paste0(tolower(.x), "_independent"),
              c(MPS, FPPS, FNR, CBI)) %>%
  select(-Threshold)

combined_metrics_surfgrass <- combined_metrics_surfgrass %>%
  left_join(surfgrass_independent_clean, by = "model")

save(combined_metrics_surfgrass, file = "code/output_data/model_results/combined_metrics_surfgrass.RData")




