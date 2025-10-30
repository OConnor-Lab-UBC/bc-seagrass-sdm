###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
#
#
# Objective:
# ---------
# create prediction dataset, need to decide if to restrict it by range of observations 
#
###############################################################################

#load packages####
library(sf)
library(tidyverse)
library(terra)
library(GGally)
library(reproducible)
library(factoextra)
library(viridis)
library(gridExtra)
library(grid)

#load seagrass data
load("code/output_data/seagrass_model_inputs.RData")

# coastline
coastline_full <- st_read("raw_data/CHS_HWL2015_Coastline.gdb", layer = "Line_CHS_Pacific_HWL_2015_5028437")
coastline <- coastline_full %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_crop(st_bbox(coastline_full) + c(0, 100000, -10000, -80000)) %>%
  st_transform(crs = "EPSG:3005")

#read in 20m rasters####
bathy_hg <- rast("raw_data/envlayers-20m-hg//bathymetry.tif")
bathy_ncc <- rast("raw_data/envlayers-20m-ncc/bathymetry.tif")
bathy_qcs <- rast("raw_data/envlayers-20m-qcs/bathymetry.tif")
bathy_wcvi <- rast("raw_data/envlayers-20m-wcvi/bathymetry.tif")
bathy_ss <- rast("raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif")

rei_hg <- rast("raw_data/REI/rei_hg.tif")
rei_ncc <- rast("raw_data/REI/rei_ncc.tif")
rei_qcs <- rast("raw_data/REI/rei_qcs.tif")
rei_wcvi <- rast("raw_data/REI/rei_wcvi.tif")
rei_ss <- rast("raw_data/REI/rei_sog.tif")

slope_hg <- rast("raw_data/envlayers-20m-hg/slope.tif")
slope_ncc <- rast("raw_data/envlayers-20m-ncc/slope.tif")
slope_qcs <- rast("raw_data/envlayers-20m-qcs/slope.tif")
slope_wcvi <- rast("raw_data/envlayers-20m-wcvi/slope.tif")
slope_ss <- rast("raw_data/envlayers-20m-shelfsalishsea/slope.tif")

freshwater_hg <- rast("raw_data/freshwater-index/hg_freshwater_index.tif")
names(freshwater_hg) <- "freshwater"
freshwater_ncc <- rast("raw_data/freshwater-index/ncc_freshwater_index.tif")
names(freshwater_ncc) <- "freshwater"
freshwater_qcs <- rast("raw_data/freshwater-index/qcs_freshwater_index.tif")
names(freshwater_qcs) <- "freshwater"
freshwater_wcvi <- rast("raw_data/freshwater-index/wcvi_freshwater_index.tif")
names(freshwater_wcvi) <- "freshwater"
freshwater_ss <- rast("raw_data/freshwater-index/salish_sea_freshwater_index.tif")
names(freshwater_ss) <- "freshwater"

tidal_all <- vrt(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))
names(tidal_all)<-"tidal"
#change to index 0-1 scale
tidal_index_all <- tidal_all/(maxFn(tidal_all))
crs(tidal_index_all) <- "EPSG:3005"

substrate_hg <- rast("raw_data/substrate_20m/updated/hg_20m.tif")
substrate_ncc <- rast("raw_data/substrate_20m/updated/ncc_20m.tif")
substrate_qcs <- rast("raw_data/substrate_20m/updated/qcs_20m.tif")
substrate_wcvi <- rast("raw_data/substrate_20m/updated/wcvi_20m.tif")
substrate_ss <- rast("raw_data/substrate_20m/updated/sog_20m.tif")
names(substrate_hg) <- names(substrate_ncc) <- names(substrate_qcs) <- names(substrate_wcvi) <- names(substrate_ss) <- "substrate"
crs(substrate_hg) <- crs(substrate_ncc) <- crs(substrate_qcs) <- crs(substrate_wcvi) <- crs(substrate_ss) <- "EPSG:3005"

folder_2013_2023 <- "code/output_data/processed_ocean_variables/years_2013-2023"
files_2013_2023 <- list.files(folder_2013_2023, pattern = "\\.tif$", full.names = TRUE)
hindcast2013_2023 <- terra::rast(files_2013_2023)

culmulative_effects <- terra::rast("code/output_data/culmulative_effects_all_20m.tif")
names(culmulative_effects) <- "cul_eff"

####make prediction data####
max_depth <- quantile(seagrass_data$depth, probs = c(0.99))
min_depth <- quantile(seagrass_data$depth, probs = c(0.001))

#haida gwaii
env_20m_hg <- as.data.frame(bathy_hg, xy=TRUE)
names(env_20m_hg) <- c("X_m", "Y_m", "depth")
env_20m_hg <- env_20m_hg %>% filter(depth <= max_depth, depth >= min_depth)
env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m","Y_m"), crs = "EPSG:3005")
HG_layers <- c(rei_hg, slope_hg, substrate_hg, freshwater_hg)
hold <- terra::extract(x = HG_layers, y = env_20m_hg_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_hg$rei <- hold$rei
env_20m_hg$substrate <- hold$substrate
env_20m_hg$slope <- hold$slope
env_20m_hg$freshwater <- hold$freshwater
hold <- extract(x = tidal_index_all, y = env_20m_hg_sf)
env_20m_hg$tidal <- hold$tidal
hold <- extract(x = culmulative_effects, y = env_20m_hg_sf)
env_20m_hg$cul_eff <- hold$cul_eff
env_20m_hg <- filter(env_20m_hg, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_hg_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_hg$NH4_5m_mean_bccm <- hold$NH4_5m_mean_bccm
env_20m_hg$NH4_5m_mean_nep36 <- hold$NH4_5m_mean_nep36
env_20m_hg$NO3_5m_mean_bccm <- hold$NO3_5m_mean_bccm
env_20m_hg$NO3_5m_mean_nep36 <- hold$NO3_5m_mean_nep36
env_20m_hg$salt_5m_mean_bccm <- hold$salt_5m_mean_bccm
env_20m_hg$salt_5m_mean_nep36 <- hold$salt_5m_mean_nep36
env_20m_hg$salt_5m_min_bccm <- hold$salt_5m_min_bccm
env_20m_hg$salt_5m_min_nep36 <- hold$salt_5m_min_nep36
env_20m_hg$salt_5m_cv_bccm <- hold$salt_5m_cv_bccm
env_20m_hg$salt_5m_cv_nep36 <- hold$salt_5m_cv_nep36
env_20m_hg$PAR_5m_mean_bccm <- hold$PAR_5m_mean_bccm
env_20m_hg$PAR_5m_mean_nep36 <- hold$PAR_5m_mean_nep36
env_20m_hg$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_hg$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_hg$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_hg$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_hg$temp_s_mean_bccm <- hold$temp_s_mean_bccm
env_20m_hg$temp_s_mean_nep36 <- hold$temp_s_mean_nep36
env_20m_hg$temp_s_max_bccm <- hold$temp_s_max_bccm
env_20m_hg$temp_s_max_nep36 <- hold$temp_s_max_nep36
env_20m_hg$temp_s_min_bccm <- hold$temp_s_min_bccm
env_20m_hg$temp_s_min_nep36 <- hold$temp_s_min_nep36
env_20m_hg$temp_s_cv_bccm <- hold$temp_s_cv_bccm
env_20m_hg$temp_s_cv_nep36 <- hold$temp_s_cv_nep36
env_20m_hg$temp_s_diff_bccm <- hold$temp_s_diff_bccm
env_20m_hg$temp_s_diff_nep36 <- hold$temp_s_diff_nep36
env_20m_hg$temp_5m_mean_bccm <- hold$temp_5m_mean_bccm
env_20m_hg$temp_5m_mean_nep36 <- hold$temp_5m_mean_nep36
env_20m_hg$temp_5m_max_bccm <- hold$temp_5m_max_bccm
env_20m_hg$temp_5m_max_nep36 <- hold$temp_5m_max_nep36
env_20m_hg$temp_5m_min_bccm <- hold$temp_5m_min_bccm
env_20m_hg$temp_5m_min_nep36 <- hold$temp_5m_min_nep36
env_20m_hg$temp_5m_cv_bccm <- hold$temp_5m_cv_bccm
env_20m_hg$temp_5m_cv_nep36 <- hold$temp_5m_cv_nep36
env_20m_hg$temp_5m_diff_bccm <- hold$temp_5m_diff_bccm
env_20m_hg$temp_5m_diff_nep36 <- hold$temp_5m_diff_nep36
env_20m_hg$do_5m_mean_bccm <- hold$do_5m_mean_bccm
env_20m_hg$do_5m_mean_nep36 <- hold$do_5m_mean_nep36
env_20m_hg$do_5m_min_bccm <- hold$do_5m_min_bccm
env_20m_hg$do_5m_min_nep36 <- hold$do_5m_min_nep36
env_20m_hg$precip_cv <- hold$precip_cv
env_20m_hg$precip_max <- hold$precip_max
env_20m_hg$precip_mean <- hold$precip_mean
env_20m_hg$precip_min <- hold$precip_min
env_20m_hg$rsds_cv <- hold$rsds_cv
env_20m_hg$rsds_max <- hold$rsds_max
env_20m_hg$rsds_mean <- hold$rsds_mean
env_20m_hg$rsds_min <- hold$rsds_min
env_20m_hg$temp_air_cv <- hold$temp_air_cv
env_20m_hg$temp_air_max <- hold$temp_air_max
env_20m_hg$temp_air_mean <- hold$temp_air_mean
env_20m_hg$temp_air_min <- hold$temp_air_min

env_20m_hg <- env_20m_hg %>% filter(do_5m_min_bccm >= 100, temp_s_min_bccm >= 0, temp_s_min_nep36 >= 0, !is.na(NH4_5m_mean_bccm)) # areas where BCCM doesn't make good predictions

env_20m_hg$region <- "Haida Gwaii"
env_20m_hg$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_hg$substrate]
summary(env_20m_hg)

save(env_20m_hg, file = "code/output_data/regional_prediction_datasets/prediction_HG.RData")

ggplot(env_20m_hg, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()

#north central coast
env_20m_ncc <- as.data.frame(bathy_ncc, xy=TRUE)
names(env_20m_ncc) <- c("X_m", "Y_m", "depth")
env_20m_ncc <- env_20m_ncc %>% filter(depth <= max_depth, depth >= min_depth)
env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
NCC_layers <- c(rei_ncc, slope_ncc, substrate_ncc, freshwater_ncc)
hold <- terra::extract(x = NCC_layers, y = env_20m_ncc_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_ncc$rei <- hold$rei
env_20m_ncc$substrate <- hold$substrate
env_20m_ncc$slope <- hold$slope
env_20m_ncc$freshwater <- hold$freshwater
hold <- extract(x = tidal_index_all, y = env_20m_ncc_sf)
env_20m_ncc$tidal <- hold$tidal
hold <- extract(x = culmulative_effects, y = env_20m_ncc_sf)
env_20m_ncc$cul_eff <- hold$cul_eff
env_20m_ncc <- filter(env_20m_ncc, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_ncc_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_ncc$NH4_5m_mean_bccm <- hold$NH4_5m_mean_bccm
env_20m_ncc$NH4_5m_mean_nep36 <- hold$NH4_5m_mean_nep36
env_20m_ncc$NO3_5m_mean_bccm <- hold$NO3_5m_mean_bccm
env_20m_ncc$NO3_5m_mean_nep36 <- hold$NO3_5m_mean_nep36
env_20m_ncc$salt_5m_mean_bccm <- hold$salt_5m_mean_bccm
env_20m_ncc$salt_5m_mean_nep36 <- hold$salt_5m_mean_nep36
env_20m_ncc$salt_5m_min_bccm <- hold$salt_5m_min_bccm
env_20m_ncc$salt_5m_min_nep36 <- hold$salt_5m_min_nep36
env_20m_ncc$salt_5m_cv_bccm <- hold$salt_5m_cv_bccm
env_20m_ncc$salt_5m_cv_nep36 <- hold$salt_5m_cv_nep36
env_20m_ncc$PAR_5m_mean_bccm <- hold$PAR_5m_mean_bccm
env_20m_ncc$PAR_5m_mean_nep36 <- hold$PAR_5m_mean_nep36
env_20m_ncc$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_ncc$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_ncc$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_ncc$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_ncc$temp_s_mean_bccm <- hold$temp_s_mean_bccm
env_20m_ncc$temp_s_mean_nep36 <- hold$temp_s_mean_nep36
env_20m_ncc$temp_s_max_bccm <- hold$temp_s_max_bccm
env_20m_ncc$temp_s_max_nep36 <- hold$temp_s_max_nep36
env_20m_ncc$temp_s_min_bccm <- hold$temp_s_min_bccm
env_20m_ncc$temp_s_min_nep36 <- hold$temp_s_min_nep36
env_20m_ncc$temp_s_cv_bccm <- hold$temp_s_cv_bccm
env_20m_ncc$temp_s_cv_nep36 <- hold$temp_s_cv_nep36
env_20m_ncc$temp_s_diff_bccm <- hold$temp_s_diff_bccm
env_20m_ncc$temp_s_diff_nep36 <- hold$temp_s_diff_nep36
env_20m_ncc$temp_5m_mean_bccm <- hold$temp_5m_mean_bccm
env_20m_ncc$temp_5m_mean_nep36 <- hold$temp_5m_mean_nep36
env_20m_ncc$temp_5m_max_bccm <- hold$temp_5m_max_bccm
env_20m_ncc$temp_5m_max_nep36 <- hold$temp_5m_max_nep36
env_20m_ncc$temp_5m_min_bccm <- hold$temp_5m_min_bccm
env_20m_ncc$temp_5m_min_nep36 <- hold$temp_5m_min_nep36
env_20m_ncc$temp_5m_cv_bccm <- hold$temp_5m_cv_bccm
env_20m_ncc$temp_5m_cv_nep36 <- hold$temp_5m_cv_nep36
env_20m_ncc$temp_5m_diff_bccm <- hold$temp_5m_diff_bccm
env_20m_ncc$temp_5m_diff_nep36 <- hold$temp_5m_diff_nep36
env_20m_ncc$do_5m_mean_bccm <- hold$do_5m_mean_bccm
env_20m_ncc$do_5m_mean_nep36 <- hold$do_5m_mean_nep36
env_20m_ncc$do_5m_min_bccm <- hold$do_5m_min_bccm
env_20m_ncc$do_5m_min_nep36 <- hold$do_5m_min_nep36
env_20m_ncc$precip_cv <- hold$precip_cv
env_20m_ncc$precip_max <- hold$precip_max
env_20m_ncc$precip_mean <- hold$precip_mean
env_20m_ncc$precip_min <- hold$precip_min
env_20m_ncc$rsds_cv <- hold$rsds_cv
env_20m_ncc$rsds_max <- hold$rsds_max
env_20m_ncc$rsds_mean <- hold$rsds_mean
env_20m_ncc$rsds_min <- hold$rsds_min
env_20m_ncc$temp_air_cv <- hold$temp_air_cv
env_20m_ncc$temp_air_max <- hold$temp_air_max
env_20m_ncc$temp_air_mean <- hold$temp_air_mean
env_20m_ncc$temp_air_min <- hold$temp_air_min

env_20m_ncc <- env_20m_ncc %>% filter(do_5m_min_bccm >= 100, temp_s_min_bccm >= 0, temp_s_min_nep36 >= 0, !is.na(NH4_5m_mean_bccm)) # areas where BCCM doesn't make good predictions

env_20m_ncc$region <- "North Central Coast"
env_20m_ncc$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_ncc$substrate]

summary(env_20m_ncc)

save(env_20m_ncc, file = "code/output_data/regional_prediction_datasets/prediction_NCC.RData")
ggplot(env_20m_ncc, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()

#queen charlotte strait
env_20m_qcs <- as.data.frame(bathy_qcs, xy=TRUE)
names(env_20m_qcs) <- c("X_m", "Y_m", "depth")
env_20m_qcs <- env_20m_qcs %>% filter(depth <= max_depth, depth >= min_depth)
env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
QCS_layers <- c(rei_qcs, slope_qcs, substrate_qcs, freshwater_qcs)
hold <- terra::extract(x = QCS_layers, y = env_20m_qcs_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_qcs$rei <- hold$rei
env_20m_qcs$substrate <- hold$substrate
env_20m_qcs$slope <- hold$slope
env_20m_qcs$freshwater <- hold$freshwater
hold <- extract(x = tidal_index_all, y = env_20m_qcs_sf)
env_20m_qcs$tidal <- hold$tidal
hold <- extract(x = culmulative_effects, y = env_20m_qcs_sf)
env_20m_qcs$cul_eff <- hold$cul_eff
env_20m_qcs <- filter(env_20m_qcs, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_qcs_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_qcs$NH4_5m_mean_bccm <- hold$NH4_5m_mean_bccm
env_20m_qcs$NH4_5m_mean_nep36 <- hold$NH4_5m_mean_nep36
env_20m_qcs$NO3_5m_mean_bccm <- hold$NO3_5m_mean_bccm
env_20m_qcs$NO3_5m_mean_nep36 <- hold$NO3_5m_mean_nep36
env_20m_qcs$salt_5m_mean_bccm <- hold$salt_5m_mean_bccm
env_20m_qcs$salt_5m_mean_nep36 <- hold$salt_5m_mean_nep36
env_20m_qcs$salt_5m_min_bccm <- hold$salt_5m_min_bccm
env_20m_qcs$salt_5m_min_nep36 <- hold$salt_5m_min_nep36
env_20m_qcs$salt_5m_cv_bccm <- hold$salt_5m_cv_bccm
env_20m_qcs$salt_5m_cv_nep36 <- hold$salt_5m_cv_nep36
env_20m_qcs$PAR_5m_mean_bccm <- hold$PAR_5m_mean_bccm
env_20m_qcs$PAR_5m_mean_nep36 <- hold$PAR_5m_mean_nep36
env_20m_qcs$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_qcs$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_qcs$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_qcs$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_qcs$temp_s_mean_bccm <- hold$temp_s_mean_bccm
env_20m_qcs$temp_s_mean_nep36 <- hold$temp_s_mean_nep36
env_20m_qcs$temp_s_max_bccm <- hold$temp_s_max_bccm
env_20m_qcs$temp_s_max_nep36 <- hold$temp_s_max_nep36
env_20m_qcs$temp_s_min_bccm <- hold$temp_s_min_bccm
env_20m_qcs$temp_s_min_nep36 <- hold$temp_s_min_nep36
env_20m_qcs$temp_s_cv_bccm <- hold$temp_s_cv_bccm
env_20m_qcs$temp_s_cv_nep36 <- hold$temp_s_cv_nep36
env_20m_qcs$temp_s_diff_bccm <- hold$temp_s_diff_bccm
env_20m_qcs$temp_s_diff_nep36 <- hold$temp_s_diff_nep36
env_20m_qcs$temp_5m_mean_bccm <- hold$temp_5m_mean_bccm
env_20m_qcs$temp_5m_mean_nep36 <- hold$temp_5m_mean_nep36
env_20m_qcs$temp_5m_max_bccm <- hold$temp_5m_max_bccm
env_20m_qcs$temp_5m_max_nep36 <- hold$temp_5m_max_nep36
env_20m_qcs$temp_5m_min_bccm <- hold$temp_5m_min_bccm
env_20m_qcs$temp_5m_min_nep36 <- hold$temp_5m_min_nep36
env_20m_qcs$temp_5m_cv_bccm <- hold$temp_5m_cv_bccm
env_20m_qcs$temp_5m_cv_nep36 <- hold$temp_5m_cv_nep36
env_20m_qcs$temp_5m_diff_bccm <- hold$temp_5m_diff_bccm
env_20m_qcs$temp_5m_diff_nep36 <- hold$temp_5m_diff_nep36
env_20m_qcs$do_5m_mean_bccm <- hold$do_5m_mean_bccm
env_20m_qcs$do_5m_mean_nep36 <- hold$do_5m_mean_nep36
env_20m_qcs$do_5m_min_bccm <- hold$do_5m_min_bccm
env_20m_qcs$do_5m_min_nep36 <- hold$do_5m_min_nep36
env_20m_qcs$precip_cv <- hold$precip_cv
env_20m_qcs$precip_max <- hold$precip_max
env_20m_qcs$precip_mean <- hold$precip_mean
env_20m_qcs$precip_min <- hold$precip_min
env_20m_qcs$rsds_cv <- hold$rsds_cv
env_20m_qcs$rsds_max <- hold$rsds_max
env_20m_qcs$rsds_mean <- hold$rsds_mean
env_20m_qcs$rsds_min <- hold$rsds_min
env_20m_qcs$temp_air_cv <- hold$temp_air_cv
env_20m_qcs$temp_air_max <- hold$temp_air_max
env_20m_qcs$temp_air_mean <- hold$temp_air_mean
env_20m_qcs$temp_air_min <- hold$temp_air_min

env_20m_qcs <- env_20m_qcs %>% filter(do_5m_min_bccm >= 100, temp_s_min_bccm >= 0, temp_s_min_nep36 >= 0, !is.na(NH4_5m_mean_bccm)) # areas where BCCM doesn't make good predictions

env_20m_qcs$region <- "Queen Charlotte Strait"
env_20m_qcs$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_qcs$substrate]

summary(env_20m_qcs)

save(env_20m_qcs, file = "code/output_data/regional_prediction_datasets/prediction_QCS.RData")
ggplot(env_20m_qcs, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()



#west coast vancouver island
env_20m_wcvi <- as.data.frame(bathy_wcvi, xy=TRUE)
names(env_20m_wcvi) <- c("X_m", "Y_m", "depth")
env_20m_wcvi <- env_20m_wcvi %>% filter(depth <= max_depth, depth >= min_depth)
env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
WCVI_layers <- c(rei_wcvi, slope_wcvi, substrate_wcvi, freshwater_wcvi)
hold <- terra::extract(x = WCVI_layers, y = env_20m_wcvi_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_wcvi$rei <- hold$rei
env_20m_wcvi$substrate <- hold$substrate
env_20m_wcvi$slope <- hold$slope
env_20m_wcvi$freshwater <- hold$freshwater
hold <- extract(x = tidal_index_all, y = env_20m_wcvi_sf)
env_20m_wcvi$tidal <- hold$tidal
hold <- extract(x = culmulative_effects, y = env_20m_wcvi_sf)
env_20m_wcvi$cul_eff <- hold$cul_eff
env_20m_wcvi <- filter(env_20m_wcvi, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_wcvi_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_wcvi$NH4_5m_mean_bccm <- hold$NH4_5m_mean_bccm
env_20m_wcvi$NH4_5m_mean_nep36 <- hold$NH4_5m_mean_nep36
env_20m_wcvi$NO3_5m_mean_bccm <- hold$NO3_5m_mean_bccm
env_20m_wcvi$NO3_5m_mean_nep36 <- hold$NO3_5m_mean_nep36
env_20m_wcvi$salt_5m_mean_bccm <- hold$salt_5m_mean_bccm
env_20m_wcvi$salt_5m_mean_nep36 <- hold$salt_5m_mean_nep36
env_20m_wcvi$salt_5m_min_bccm <- hold$salt_5m_min_bccm
env_20m_wcvi$salt_5m_min_nep36 <- hold$salt_5m_min_nep36
env_20m_wcvi$salt_5m_cv_bccm <- hold$salt_5m_cv_bccm
env_20m_wcvi$salt_5m_cv_nep36 <- hold$salt_5m_cv_nep36
env_20m_wcvi$PAR_5m_mean_bccm <- hold$PAR_5m_mean_bccm
env_20m_wcvi$PAR_5m_mean_nep36 <- hold$PAR_5m_mean_nep36
env_20m_wcvi$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_wcvi$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_wcvi$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_wcvi$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_wcvi$temp_s_mean_bccm <- hold$temp_s_mean_bccm
env_20m_wcvi$temp_s_mean_nep36 <- hold$temp_s_mean_nep36
env_20m_wcvi$temp_s_max_bccm <- hold$temp_s_max_bccm
env_20m_wcvi$temp_s_max_nep36 <- hold$temp_s_max_nep36
env_20m_wcvi$temp_s_min_bccm <- hold$temp_s_min_bccm
env_20m_wcvi$temp_s_min_nep36 <- hold$temp_s_min_nep36
env_20m_wcvi$temp_s_cv_bccm <- hold$temp_s_cv_bccm
env_20m_wcvi$temp_s_cv_nep36 <- hold$temp_s_cv_nep36
env_20m_wcvi$temp_s_diff_bccm <- hold$temp_s_diff_bccm
env_20m_wcvi$temp_s_diff_nep36 <- hold$temp_s_diff_nep36
env_20m_wcvi$temp_5m_mean_bccm <- hold$temp_5m_mean_bccm
env_20m_wcvi$temp_5m_mean_nep36 <- hold$temp_5m_mean_nep36
env_20m_wcvi$temp_5m_max_bccm <- hold$temp_5m_max_bccm
env_20m_wcvi$temp_5m_max_nep36 <- hold$temp_5m_max_nep36
env_20m_wcvi$temp_5m_min_bccm <- hold$temp_5m_min_bccm
env_20m_wcvi$temp_5m_min_nep36 <- hold$temp_5m_min_nep36
env_20m_wcvi$temp_5m_cv_bccm <- hold$temp_5m_cv_bccm
env_20m_wcvi$temp_5m_cv_nep36 <- hold$temp_5m_cv_nep36
env_20m_wcvi$temp_5m_diff_bccm <- hold$temp_5m_diff_bccm
env_20m_wcvi$temp_5m_diff_nep36 <- hold$temp_5m_diff_nep36
env_20m_wcvi$do_5m_mean_bccm <- hold$do_5m_mean_bccm
env_20m_wcvi$do_5m_mean_nep36 <- hold$do_5m_mean_nep36
env_20m_wcvi$do_5m_min_bccm <- hold$do_5m_min_bccm
env_20m_wcvi$do_5m_min_nep36 <- hold$do_5m_min_nep36
env_20m_wcvi$precip_cv <- hold$precip_cv
env_20m_wcvi$precip_max <- hold$precip_max
env_20m_wcvi$precip_mean <- hold$precip_mean
env_20m_wcvi$precip_min <- hold$precip_min
env_20m_wcvi$rsds_cv <- hold$rsds_cv
env_20m_wcvi$rsds_max <- hold$rsds_max
env_20m_wcvi$rsds_mean <- hold$rsds_mean
env_20m_wcvi$rsds_min <- hold$rsds_min
env_20m_wcvi$temp_air_cv <- hold$temp_air_cv
env_20m_wcvi$temp_air_max <- hold$temp_air_max
env_20m_wcvi$temp_air_mean <- hold$temp_air_mean
env_20m_wcvi$temp_air_min <- hold$temp_air_min

env_20m_wcvi <- env_20m_wcvi %>% filter(do_5m_min_bccm >= 100, temp_s_min_bccm >= 0, temp_s_min_nep36 >= 0, !is.na(NH4_5m_mean_bccm)) # areas where BCCM doesn't make good predictions

env_20m_wcvi$region <- "West Coast Vancouver Island"
env_20m_wcvi$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_wcvi$substrate]

summary(env_20m_wcvi)

save(env_20m_wcvi, file = "code/output_data/regional_prediction_datasets/prediction_WCVI.RData")
ggplot(env_20m_wcvi, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()


#salish sea
env_20m_ss <- as.data.frame(bathy_ss, xy=TRUE)
names(env_20m_ss) <- c("X_m", "Y_m", "depth")
env_20m_ss <- env_20m_ss %>% filter(depth <= max_depth, depth >= min_depth)
env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
SS_layers <- c(rei_ss, slope_ss, substrate_ss, freshwater_ss)
hold <- terra::extract(x = SS_layers, y = env_20m_ss_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_ss$rei <- hold$rei
env_20m_ss$substrate <- hold$substrate
env_20m_ss$slope <- hold$slope
env_20m_ss$freshwater <- hold$freshwater
hold <- extract(x = tidal_index_all, y = env_20m_ss_sf)
env_20m_ss$tidal <- hold$tidal
hold <- extract(x = culmulative_effects, y = env_20m_ss_sf)
env_20m_ss$cul_eff <- hold$cul_eff
env_20m_ss <- filter(env_20m_ss, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_ss_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_ss$NH4_5m_mean_bccm <- hold$NH4_5m_mean_bccm
env_20m_ss$NH4_5m_mean_nep36 <- hold$NH4_5m_mean_nep36
env_20m_ss$NO3_5m_mean_bccm <- hold$NO3_5m_mean_bccm
env_20m_ss$NO3_5m_mean_nep36 <- hold$NO3_5m_mean_nep36
env_20m_ss$salt_5m_mean_bccm <- hold$salt_5m_mean_bccm
env_20m_ss$salt_5m_mean_nep36 <- hold$salt_5m_mean_nep36
env_20m_ss$salt_5m_min_bccm <- hold$salt_5m_min_bccm
env_20m_ss$salt_5m_min_nep36 <- hold$salt_5m_min_nep36
env_20m_ss$salt_5m_cv_bccm <- hold$salt_5m_cv_bccm
env_20m_ss$salt_5m_cv_nep36 <- hold$salt_5m_cv_nep36
env_20m_ss$PAR_5m_mean_bccm <- hold$PAR_5m_mean_bccm
env_20m_ss$PAR_5m_mean_nep36 <- hold$PAR_5m_mean_nep36
env_20m_ss$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_ss$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_ss$PAR_5m_min_bccm <- hold$PAR_5m_min_bccm
env_20m_ss$PAR_5m_min_nep36 <- hold$PAR_5m_min_nep36
env_20m_ss$temp_s_mean_bccm <- hold$temp_s_mean_bccm
env_20m_ss$temp_s_mean_nep36 <- hold$temp_s_mean_nep36
env_20m_ss$temp_s_max_bccm <- hold$temp_s_max_bccm
env_20m_ss$temp_s_max_nep36 <- hold$temp_s_max_nep36
env_20m_ss$temp_s_min_bccm <- hold$temp_s_min_bccm
env_20m_ss$temp_s_min_nep36 <- hold$temp_s_min_nep36
env_20m_ss$temp_s_cv_bccm <- hold$temp_s_cv_bccm
env_20m_ss$temp_s_cv_nep36 <- hold$temp_s_cv_nep36
env_20m_ss$temp_s_diff_bccm <- hold$temp_s_diff_bccm
env_20m_ss$temp_s_diff_nep36 <- hold$temp_s_diff_nep36
env_20m_ss$temp_5m_mean_bccm <- hold$temp_5m_mean_bccm
env_20m_ss$temp_5m_mean_nep36 <- hold$temp_5m_mean_nep36
env_20m_ss$temp_5m_max_bccm <- hold$temp_5m_max_bccm
env_20m_ss$temp_5m_max_nep36 <- hold$temp_5m_max_nep36
env_20m_ss$temp_5m_min_bccm <- hold$temp_5m_min_bccm
env_20m_ss$temp_5m_min_nep36 <- hold$temp_5m_min_nep36
env_20m_ss$temp_5m_cv_bccm <- hold$temp_5m_cv_bccm
env_20m_ss$temp_5m_cv_nep36 <- hold$temp_5m_cv_nep36
env_20m_ss$temp_5m_diff_bccm <- hold$temp_5m_diff_bccm
env_20m_ss$temp_5m_diff_nep36 <- hold$temp_5m_diff_nep36
env_20m_ss$do_5m_mean_bccm <- hold$do_5m_mean_bccm
env_20m_ss$do_5m_mean_nep36 <- hold$do_5m_mean_nep36
env_20m_ss$do_5m_min_bccm <- hold$do_5m_min_bccm
env_20m_ss$do_5m_min_nep36 <- hold$do_5m_min_nep36
env_20m_ss$precip_cv <- hold$precip_cv
env_20m_ss$precip_max <- hold$precip_max
env_20m_ss$precip_mean <- hold$precip_mean
env_20m_ss$precip_min <- hold$precip_min
env_20m_ss$rsds_cv <- hold$rsds_cv
env_20m_ss$rsds_max <- hold$rsds_max
env_20m_ss$rsds_mean <- hold$rsds_mean
env_20m_ss$rsds_min <- hold$rsds_min
env_20m_ss$temp_air_cv <- hold$temp_air_cv
env_20m_ss$temp_air_max <- hold$temp_air_max
env_20m_ss$temp_air_mean <- hold$temp_air_mean
env_20m_ss$temp_air_min <- hold$temp_air_min

env_20m_ss <- env_20m_ss %>% filter(do_5m_min_bccm >= 100, temp_s_min_bccm >= 0, temp_s_min_nep36 >= 0, !is.na(NH4_5m_mean_bccm)) # areas where BCCM doesn't make good predictions

env_20m_ss$region <- "Salish Sea"
env_20m_ss$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_ss$substrate]
summary(env_20m_ss)

rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "SpatRaster"])

save(env_20m_ss, file = "code/output_data/regional_prediction_datasets/prediction_SS.RData")

ggplot(env_20m_ss, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()

load("code/output_data/regional_prediction_datasets/prediction_SS.RData")
load("code/output_data/regional_prediction_datasets/prediction_NCC.RData")
load("code/output_data/regional_prediction_datasets/prediction_WCVI.RData")
load("code/output_data/regional_prediction_datasets/prediction_QCS.RData")
load("code/output_data/regional_prediction_datasets/prediction_HG.RData")


#combine all and filter
env_20m_all <- bind_rows(env_20m_hg, env_20m_ncc, env_20m_qcs, env_20m_wcvi,env_20m_ss)
summary(env_20m_all)

# ggplot(env_20m_all, aes(x = X_m, y = Y_m, color = substrate))+
#   geom_point(size = 0.01, pch = 15)+
#   theme_bw()+
#   coord_equal()

#need to think about this in reference to how much we want to extrapolate, will depend on future conditions
env_20m_all <- env_20m_all %>%
  filter(#freshwater < quantile(seagrass_data$freshwater, probs = 0.99),  # we have no surveys in freshwater areas, so better to exclude
         salt_5m_mean_bccm > quantile(seagrass_data$saltmean_bccm, probs = 0.001), # same as above
         !is.na(cul_eff),
         !is.na(NH4_5m_mean_nep36), #temperature > quantile(seagrass_data$temperature, probs = 0.001), # want to allow extrapolation into higher temperatures
         #temperature < quantile(seagrass_data$temperature, probs = 0.999),
         tidal < quantile(seagrass_data$tidal, probs = 0.9999), # remove high current areas where we can't sample, neither species observed in high values in quadrats 
         rei < quantile(seagrass_data$rei, probs = 0.99999))   # remove high exposure areas where we can't sample? PH found in moderately high but not

summary(env_20m_all)

#scale environmental variables function####
scale_fun_ref <- function(x, reference){
  (x  - mean(reference)) / sd(reference)
}

env_20m_all <- env_20m_all %>%
  mutate(depth_stnd = scale_fun_ref(depth, reference = seagrass_data$depth),
         rei_stnd = scale_fun_ref(rei, reference = seagrass_data$rei),
         rei_sqrt_stnd = scale_fun_ref(sqrt(rei), reference = sqrt(seagrass_data$rei)),
         tidal_stnd = scale_fun_ref(tidal, reference = seagrass_data$tidal),
         tidal_sqrt_stnd = scale_fun_ref(sqrt(tidal), reference = sqrt(seagrass_data$tidal)),
         freshwater_stnd = scale_fun_ref(freshwater, reference = seagrass_data$freshwater),
         freshwater_sqrt_stnd = scale_fun_ref(sqrt(freshwater), reference = sqrt(seagrass_data$freshwater)),
         slope_stnd = scale_fun_ref(slope, reference = seagrass_data$slope),
         slope_sqrt_stnd = scale_fun_ref(sqrt(slope), reference = sqrt(seagrass_data$slope)),
         NH4_bccm_stnd = scale_fun_ref(NH4_5m_mean_bccm, reference = seagrass_data$NH4_bccm),
         NH4_nep_stnd = scale_fun_ref(NH4_5m_mean_nep36, reference = seagrass_data$NH4_nep),
         NO3_bccm_stnd = scale_fun_ref(NO3_5m_mean_bccm, reference = seagrass_data$NO3_bccm),
         NO3_nep_stnd = scale_fun_ref(NO3_5m_mean_nep36, reference = seagrass_data$NO3_nep),
         saltmean_bccm_stnd = scale_fun_ref(salt_5m_mean_bccm, reference = seagrass_data$saltmean_bccm),
         saltmean_nep_stnd = scale_fun_ref(salt_5m_mean_nep36, reference = seagrass_data$saltmean_nep),
         saltmean_bccm_sq_stnd = scale_fun_ref((salt_5m_mean_bccm)^2, reference = (seagrass_data$saltmean_bccm)^2),
         saltmean_nep_sq_stnd = scale_fun_ref((salt_5m_mean_nep36)^2, reference = (seagrass_data$saltmean_nep)^2),
         saltmin_bccm_stnd = scale_fun_ref(salt_5m_min_bccm, reference = seagrass_data$saltmin_bccm),
         saltmin_nep_stnd = scale_fun_ref(salt_5m_min_nep36, reference = seagrass_data$saltmin_nep),
         saltmin_bccm_sq_stnd = scale_fun_ref((salt_5m_min_bccm)^2, reference = (seagrass_data$saltmin_bccm)^2),
         saltmin_nep_sq_stnd = scale_fun_ref((salt_5m_min_nep36)^2, reference = (seagrass_data$saltmin_nep)^2),
         saltcv_bccm_stnd = scale_fun_ref(salt_5m_cv_bccm, reference = seagrass_data$saltcv_bccm),
         saltcv_nep_stnd = scale_fun_ref(salt_5m_cv_nep36, reference = seagrass_data$saltcv_nep),
         PARmean_bccm_stnd = scale_fun_ref(PAR_5m_mean_bccm, reference = seagrass_data$PARmean_bccm),
         PARmean_nep_stnd = scale_fun_ref(PAR_5m_mean_nep36, reference = seagrass_data$PARmean_nep),
         PARmin_bccm_stnd = scale_fun_ref(PAR_5m_min_bccm, reference = seagrass_data$PARmin_bccm),
         PARmin_nep_stnd = scale_fun_ref(PAR_5m_min_nep36, reference = seagrass_data$PARmin_nep),
         surftempmean_bccm_stnd = scale_fun_ref(temp_s_mean_bccm, reference = seagrass_data$surftempmean_bccm),
         surftempmean_nep_stnd = scale_fun_ref(temp_s_mean_nep36, reference = seagrass_data$surftempmean_nep),
         surftempmin_bccm_stnd = scale_fun_ref(temp_s_min_bccm, reference = seagrass_data$surftempmin_bccm),
         surftempmin_nep_stnd = scale_fun_ref(temp_s_min_nep36, reference = seagrass_data$surftempmin_nep),
         surftempmax_bccm_stnd = scale_fun_ref(temp_s_max_bccm, reference = seagrass_data$surftempmax_bccm),
         surftempmax_nep_stnd = scale_fun_ref(temp_s_max_nep36, reference = seagrass_data$surftempmax_nep),
         surftempcv_bccm_stnd = scale_fun_ref(temp_s_cv_bccm, reference = seagrass_data$surftempcv_bccm),
         surftempcv_nep_stnd = scale_fun_ref(temp_s_cv_nep36, reference = seagrass_data$surftempcv_nep),
         surftempdiff_bccm_stnd = scale_fun_ref(temp_s_diff_bccm, reference = seagrass_data$surftempdiff_bccm),
         surftempdiff_nep_stnd = scale_fun_ref(temp_s_diff_nep36, reference = seagrass_data$surftempdiff_nep),
         tempmean_bccm_stnd = scale_fun_ref(temp_5m_mean_bccm, reference = seagrass_data$tempmean_bccm),
         tempmean_nep_stnd = scale_fun_ref(temp_5m_mean_nep36, reference = seagrass_data$tempmean_nep),
         tempmin_bccm_stnd = scale_fun_ref(temp_5m_min_bccm, reference = seagrass_data$tempmin_bccm),
         tempmin_nep_stnd = scale_fun_ref(temp_5m_min_nep36, reference = seagrass_data$tempmin_nep),
         tempmax_bccm_stnd = scale_fun_ref(temp_5m_max_bccm, reference = seagrass_data$tempmax_bccm),
         tempmax_nep_stnd = scale_fun_ref(temp_5m_max_nep36, reference = seagrass_data$tempmax_nep),
         tempcv_bccm_stnd = scale_fun_ref(temp_5m_cv_bccm, reference = seagrass_data$tempcv_bccm),
         tempcv_nep_stnd = scale_fun_ref(temp_5m_cv_nep36, reference = seagrass_data$tempcv_nep),
         tempdiff_bccm_stnd = scale_fun_ref(temp_5m_diff_bccm, reference = seagrass_data$tempdiff_bccm),
         tempdiff_nep_stnd = scale_fun_ref(temp_5m_diff_nep36, reference = seagrass_data$tempdiff_nep),
         DOmean_bccm_stnd = scale_fun_ref(do_5m_mean_bccm, reference = seagrass_data$DOmean_bccm),
         DOmean_nep_stnd = scale_fun_ref(do_5m_mean_nep36, reference = seagrass_data$DOmean_nep),
         DOmin_bccm_stnd = scale_fun_ref(do_5m_min_bccm, reference = seagrass_data$DOmin_bccm),
         DOmin_nep_stnd = scale_fun_ref(do_5m_min_nep36, reference = seagrass_data$DOmin_nep),
         airtempcv_stnd = scale_fun_ref(temp_air_cv, reference = seagrass_data$airtempcv), 
         airtempmax_stnd = scale_fun_ref(temp_air_max, reference = seagrass_data$airtempmax), 
         airtempmean_stnd = scale_fun_ref(temp_air_mean, reference = seagrass_data$airtempmean), 
         airtempmin_stnd = scale_fun_ref(temp_air_min, reference = seagrass_data$airtempmin),
         prcv_stnd = scale_fun_ref(precip_cv, reference = seagrass_data$prcv), 
         prmax_stnd = scale_fun_ref(precip_max, reference = seagrass_data$prmax), 
         prmin_stnd = scale_fun_ref(precip_min, reference = seagrass_data$prmin), 
         prmean_stnd = scale_fun_ref(precip_mean, reference = seagrass_data$prmean), 
         rsdscv_stnd = scale_fun_ref(rsds_cv, reference = seagrass_data$rsdscv), 
         rsdsmax_stnd = scale_fun_ref(rsds_max, reference = seagrass_data$rsdsmax), 
         rsdsmean_stnd = scale_fun_ref(rsds_mean, reference = seagrass_data$rsdsmean), 
         rsdsmin_stnd = scale_fun_ref(rsds_min, reference = seagrass_data$rsdsmin),
         cul_eff_stnd = scale_fun_ref(cul_eff, reference = seagrass_data$cul_eff)) %>%
  mutate(ID = 1:nrow(env_20m_all),
         HKey = factor("new"),
         X = X_m/1000,
         Y = Y_m/1000,
         substrate = factor(substrate))

summary(env_20m_all)
 
env_20m_all <- env_20m_all %>%
  filter(!is.na(freshwater_sqrt_stnd)) 

#env_20m_all_sf <- st_as_sf(env_20m_all, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
  

#### Compare the sampled predictor space to that of the total hindcast climatology
#2013-2023 hindcast climatology without depth and substrate and turn into a dataframe

hindcast_predictor_data <- env_20m_all %>% select(rei_stnd, tidal_stnd, slope_stnd, NH4_bccm_stnd, NO3_bccm_stnd, saltmin_bccm_stnd, saltcv_bccm_stnd, PARmin_bccm_stnd, surftempmean_bccm_stnd, surftempmin_bccm_stnd, 
                                                  surftempmax_bccm_stnd, surftempcv_bccm_stnd, tempmean_bccm_stnd, tempmin_bccm_stnd, tempmax_bccm_stnd, tempcv_bccm_stnd, DOmin_bccm_stnd, DOmean_bccm_stnd, freshwater_sqrt_stnd,
                                                  airtempmean_stnd, prmean_stnd, rsdsmean_stnd)

# transect predictor data
transect_predictor_data <- seagrass_data %>% select(rei_stnd, tidal_stnd, slope_stnd, NH4_bccm_stnd, NO3_bccm_stnd, saltmin_bccm_stnd, saltcv_bccm_stnd, PARmin_bccm_stnd, surftempmean_bccm_stnd, surftempmin_bccm_stnd, 
                                                    surftempmax_bccm_stnd, surftempcv_bccm_stnd, tempmean_bccm_stnd, tempmin_bccm_stnd, tempmax_bccm_stnd, tempcv_bccm_stnd, DOmin_bccm_stnd, DOmean_bccm_stnd, freshwater_sqrt_stnd,
                                                    airtempmean_stnd, prmean_stnd, rsdsmean_stnd)
#transect_predictor_data <- seagrass_data %>% select(rei_stnd, tidal_stnd, slope_stnd, NH4_stnd, NO3_stnd, saltmean_stnd, saltmin_stnd, PARmean_stnd, surftempmean_stnd, surftempmin_stnd, surftempmax_stnd, surftempcv_stnd,surftempdiff_stnd, tempmean_stnd, tempmin_stnd, tempmax_stnd, tempcv_stnd, tempdiff_stnd, DOmean_stnd)

#Create a PCA for the total predictor space of the hindcast
Hindcast_PCA <- princomp(hindcast_predictor_data)
summary(Hindcast_PCA)

#Generate a scree plot of the importance of each PC
factoextra::fviz_eig(Hindcast_PCA, addlabels = TRUE)
#Generate a quick biplot of the first two PCs in PCA-space
fviz_pca_var(Hindcast_PCA, col.var = "black")
#Generate a plot of Cos2 to represent how well each variable is represented in the PCA's first two PCs
fviz_cos2(Hindcast_PCA, choice = "var", axes = 1:2)
#Combine the two previous plots 
fviz_pca_var(Hindcast_PCA, col.var = "cos2", gradient.cols = viridis(n=200), repel = TRUE)

#Separate out the PCA scores
Hindcast_Scores <- as.data.frame(Hindcast_PCA$scores)
Hindcast_Scores$Type <- "Hindcast Predictor Space"

#Project the environmental conditions of the transect surveys onto the PCA space
Transect_Scores <- as.data.frame(predict(Hindcast_PCA, transect_predictor_data))
Transect_Scores$Type <- "Survey Predictor Space"

#Rowbind all scores into a single dataframe for plotting
All_Scores <- rbind(Hindcast_Scores, Transect_Scores)

#Create a function that creates a frequency plot for the values of the study area and transect surveys for a specified principle component
Density_Plot <- function(All_Scores, comp, title, x, y, show.legend = TRUE){
  Dense <- ggplot(All_Scores) + 
    geom_density(aes(x = All_Scores[, comp], fill = Type), alpha = 0.5, show.legend = show.legend) + 
    scale_fill_manual(values = c("#f28e2b", "#76b7b2")) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(title) +
    ylab("Frequency") +
    theme_classic() +
    guides(x.sec = "axis", y.sec = "axis") +
    theme(
      axis.line.y.left = element_line(size = 0.75, color = "black"),
      axis.line.y.right = element_line(size = 0.75, color = "black"),
      axis.line.x.top = element_line(size = 0.75, color = "black"),
      axis.line.x.bottom = element_line(size = 0.75, color = "black"),
      axis.ticks.x.top = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.x.top = element_blank(),
      axis.text.y.right = element_blank(),
      legend.title = element_blank(),
      legend.position = if (show.legend) c(x, y) else "none",  # inside first plot
      legend.text = element_text(size = 10),    # moderate size
      legend.key.size = unit(0.8, "cm"),        # balanced key size
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.background = element_blank()       # no box
    )
  return(Dense)
}

# Create plots (legend only in the first one)
p1 <- Density_Plot(All_Scores, 1, "Principal Component 1 Value (39%)", 0.3, 0.8, show.legend = TRUE)  # moved legend left
p2 <- Density_Plot(All_Scores, 2, "Principal Component 2 Value (19%)", 0.2, 0.81, show.legend = FALSE)
p3 <- Density_Plot(All_Scores, 3, "Principal Component 3 Value (15%)", 0.75, 0.81, show.legend = FALSE)
p4 <- Density_Plot(All_Scores, 4, "Principal Component 4 Value (9%)", 0.25, 0.81, show.legend = FALSE)

# Arrange in 2x2 grid
PCA_dense_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)

PCA_dense_plot
ggsave(filename = "./figures/pre-analysis/Predictor_Hindcast_Transect_PCA_Density_Plot.png", plot = PCA_dense_plot, scale = 1, height = 10, width = 12)

#Create a PCA plot of the first two principal components of our hindcast and survey predictor space with no points but 95% elipses
PCA_Plot_Ellipses <- ggplot(All_Scores) + 
  scale_colour_manual( values = c("#f28e2b", "#76b7b2")) + 
  xlab("Principal Component 1 (39%)") + 
  ylab("Principal Component 2 (19%)") + 
  theme_classic() + 
  guides(x.sec = "axis", y.sec = "axis")+ 
  theme(axis.line.y.left = element_line(size = 0.75, color = "black"), 
        axis.line.y.right = element_line(size = 0.75, color = "black"), 
        axis.line.x.top = element_line(size = 0.75, color = "black"), 
        axis.line.x.bottom = element_line(size = 0.75, color = "black"),
        axis.ticks.x.top = element_blank(), 
        axis.ticks.y.right = element_blank(), 
        axis.text.x.top = element_blank(), 
        axis.text.y.right = element_blank())+ 
  theme(legend.position = "bottom", 
        legend.justification = c("left","top"), 
        legend.box.just = "right") + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(x = Comp.1, y = Comp.2, color = Type, group = Type, fill = Type)) + 
  scale_fill_manual( values = c("#f28e2b", "#76b7b2")) + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size= 18), 
        legend.background = element_blank())

PCA_Plot_Ellipses
ggsave("./figures/pre-analysis/Predictor_Hindcast_Transect_PCA_Plot_Ellipses.png", plot = PCA_Plot_Ellipses, width = 10, height = 12)

#hindcast_MESS<- modEvA::MESS(transect_predictor_data, hindcast_predictor_data)
# refer to elith et al for more details on MESS. measures the similarity of any given point to a reference set of points, with respect to the chosen predictor variables. It reports the closeness of the point to the distribution of reference points, gives negative values for dissimilar points and maps these values across the whole prediction region
hindcast_MESS<- predicts::mess(x = hindcast_predictor_data, v= transect_predictor_data, full = FALSE)

env_20m_all <- env_20m_all %>% cbind(hindcast_MESS)

mess_raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, mess)
mess_raster_hg <- rast(x = mess_raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_hg, file.path("./raster/mess_predictions_hg.tif"), overwrite=TRUE)

mess_raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, mess)
mess_raster_ss <- rast(x = mess_raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_ss, file.path("./raster/mess_predictions_ss.tif"), overwrite=TRUE)

mess_raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, mess)
mess_raster_wcvi <- rast(x = mess_raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_wcvi, file.path("./raster/mess_predictions_wcvi.tif"), overwrite=TRUE)

mess_raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, mess)
mess_raster_ncc <- rast(x = mess_raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_ncc, file.path("./raster/mess_predictions_ncc.tif"), overwrite=TRUE)

mess_raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, mess)
mess_raster_qcs <- rast(x = mess_raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_qcs, file.path("./raster/mess_predictions_qcs.tif"), overwrite=TRUE)



# see if there is a difference between if we only use 1993-2012 data compared to the full 2013-2023 dataset

transect_predictor_data_1993_2009 <- seagrass_data %>% 
  filter(Year< 2010) %>%
  select(rei_stnd, tidal_stnd, slope_stnd, NH4_bccm_stnd, NO3_bccm_stnd, saltmin_bccm_stnd, saltcv_bccm_stnd, PARmin_bccm_stnd, surftempmean_bccm_stnd, surftempmin_bccm_stnd, 
         surftempmax_bccm_stnd, surftempcv_bccm_stnd, tempmean_bccm_stnd, tempmin_bccm_stnd, tempmax_bccm_stnd, tempcv_bccm_stnd, DOmin_bccm_stnd, DOmean_bccm_stnd, freshwater_sqrt_stnd,
         airtempmean_stnd, prmean_stnd, rsdsmean_stnd )

temporal_validation_MESS<- predicts::mess(x = hindcast_predictor_data, v= transect_predictor_data_1993_2009, full = FALSE)
temporal_validation_MESS <- temporal_validation_MESS %>%
  rename(mess_tv = mess)
env_20m_all <- env_20m_all %>% cbind(temporal_validation_MESS)

env_20m_all_sub <- env_20m_all %>% mutate(mess_col = ifelse(mess < -0.01, "Below", "Above"),
                                          mess_tv_col = ifelse(mess_tv < -0.01, "Below", "Above"))
cols<- c("Above" = "grey", "Below" = "red")
mess_plot<-ggplot(env_20m_all_sub)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=mess_col, fill=mess_col, width=20,height=20))+
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_sf(expand = FALSE)+
  ylab("")+
  xlab("") 
mess_plot
ggsave("./figures/pre-analysis/hindcast_mess.png", height = 6, width = 6)


temporal_val_mess_plot<-ggplot(env_20m_all_sub)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=mess_tv_col, fill = mess_tv_col, width=20,height=20))+
  #scale_colour_gradient2() +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_sf(expand = FALSE)+
  ylab("")+
  xlab("") 
temporal_val_mess_plot
ggsave("./figures/pre-analysis/hindcast_mess_tv.png", height = 6, width = 6)

mess_raster_hg_tv <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, mess_tv)
mess_raster_hg_tv <- rast(x = mess_raster_hg_tv %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_hg_tv, file.path("./raster/mess_predictions_hg_tv.tif"), overwrite=TRUE)

mess_raster_ss_tv <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, mess_tv)
mess_raster_ss_tv <- rast(x = mess_raster_ss_tv %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_ss_tv, file.path("./raster/mess_predictions_ss_tv.tif"), overwrite=TRUE)

mess_raster_wcvi_tv <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, mess_tv)
mess_raster_wcvi_tv <- rast(x = mess_raster_wcvi_tv %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_wcvi_tv, file.path("./raster/mess_predictions_wcvi_tv.tif"), overwrite=TRUE)

mess_raster_ncc_tv <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, mess_tv)
mess_raster_ncc_tv <- rast(x = mess_raster_ncc_tv %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_ncc_tv, file.path("./raster/mess_predictions_ncc_tv.tif"), overwrite=TRUE)

mess_raster_qcs_tv <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, mess_tv)
mess_raster_qcs_tv <- rast(x = mess_raster_qcs_tv %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(mess_raster_qcs_tv, file.path("./raster/mess_predictions_qcs_tv.tif"), overwrite=TRUE)


#save outputs####
save(env_20m_all, file = "code/output_data/prediction_model_inputs.RData")
