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

#load seagrass data
load("code/output_data/seagrass_model_inputs.RData")

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

substrate_hg <- rast("raw_data/substrate_20m/hg_20m.tif")
substrate_ncc <- rast("raw_data/substrate_20m/ncc_20m.tif")
substrate_qcs <- rast("raw_data/substrate_20m/qcs_20m.tif")
substrate_wcvi <- rast("raw_data/substrate_20m/wcvi_20m.tif")
substrate_ss <- rast("raw_data/substrate_20m/sog_20m.tif")
names(substrate_hg) <- names(substrate_ncc) <- names(substrate_qcs) <- names(substrate_wcvi) <- names(substrate_ss) <- "substrate"
crs(substrate_hg) <- crs(substrate_ncc) <- crs(substrate_qcs) <- crs(substrate_wcvi) <- crs(substrate_ss) <- "EPSG:3005"

hindcast2013_2023 <- terra::rast("code/output_data/processed_ocean_variables/Predictor_Hindcast_Climatologies_2013-2023.tif")
#hindcast1993_2023 <- terra::rast("code/output_data/processed_ocean_variables/Predictor_Hindcast_Climatologies_1993-2023.tif")


####make prediction data####
max_depth <- quantile(seagrass_data$depth, probs = c(0.99))
min_depth <- quantile(seagrass_data$depth, probs = c(0.01))

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
env_20m_hg <- filter(env_20m_hg, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_hg_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_hg$NH4_5m_mean <- hold$NH4_5m_mean
env_20m_hg$NO3_5m_mean <- hold$NO3_5m_mean
env_20m_hg$salt_5m_mean <- hold$salt_5m_mean
env_20m_hg$salt_5m_min <- hold$salt_5m_min
env_20m_hg$PAR_5m_mean <- hold$PAR_5m_mean
env_20m_hg$PAR_5m_min <- hold$PAR_5m_min
env_20m_hg$PAR_5m_max <- hold$PAR_5m_max
env_20m_hg$temp_s_mean <- hold$temp_s_mean
env_20m_hg$temp_s_max <- hold$temp_s_max
env_20m_hg$temp_s_min <- hold$temp_s_min
env_20m_hg$temp_5m_mean <- hold$temp_5m_mean
env_20m_hg$temp_5m_max <- hold$temp_5m_max
env_20m_hg$temp_5m_min <- hold$temp_5m_min
env_20m_hg$do_5m_mean <- hold$do_5m_mean
env_20m_hg$do_5m_min <- hold$do_5m_min
env_20m_hg <- env_20m_hg %>% filter(do_5m_min >= 100, temp_s_min >= 0, !is.na(NH4_5m_mean)) # areas where BCCM doesn't make good predictions

env_20m_hg$region <- "Haida Gwaii"
env_20m_hg$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_hg$substrate]
summary(env_20m_hg)

#save(env_20m_hg, file = "./data/prediction_HG.RData")

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
env_20m_ncc <- filter(env_20m_ncc, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_ncc_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_ncc$NH4_5m_mean <- hold$NH4_5m_mean
env_20m_ncc$NO3_5m_mean <- hold$NO3_5m_mean
env_20m_ncc$salt_5m_mean <- hold$salt_5m_mean
env_20m_ncc$salt_5m_min <- hold$salt_5m_min
env_20m_ncc$PAR_5m_mean <- hold$PAR_5m_mean
env_20m_ncc$PAR_5m_min <- hold$PAR_5m_min
env_20m_ncc$PAR_5m_max <- hold$PAR_5m_max
env_20m_ncc$temp_s_mean <- hold$temp_s_mean
env_20m_ncc$temp_s_max <- hold$temp_s_max
env_20m_ncc$temp_s_min <- hold$temp_s_min
env_20m_ncc$temp_5m_mean <- hold$temp_5m_mean
env_20m_ncc$temp_5m_max <- hold$temp_5m_max
env_20m_ncc$temp_5m_min <- hold$temp_5m_min
env_20m_ncc$do_5m_mean <- hold$do_5m_mean
env_20m_ncc$do_5m_min <- hold$do_5m_min
env_20m_ncc <- env_20m_ncc %>% filter(do_5m_min >= 100, temp_s_min >= 0, !is.na(NH4_5m_mean)) # areas where BCCM doesn't make good predictions

env_20m_ncc$region <- "North Central Coast"
env_20m_ncc$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_ncc$substrate]

summary(env_20m_ncc)

#save(env_20m_ncc, file = "./data/prediction_NCC.RData")
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
env_20m_qcs <- filter(env_20m_qcs, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_qcs_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_qcs$NH4_5m_mean <- hold$NH4_5m_mean
env_20m_qcs$NO3_5m_mean <- hold$NO3_5m_mean
env_20m_qcs$salt_5m_mean <- hold$salt_5m_mean
env_20m_qcs$salt_5m_min <- hold$salt_5m_min
env_20m_qcs$PAR_5m_mean <- hold$PAR_5m_mean
env_20m_qcs$PAR_5m_min <- hold$PAR_5m_min
env_20m_qcs$PAR_5m_max <- hold$PAR_5m_max
env_20m_qcs$temp_s_mean <- hold$temp_s_mean
env_20m_qcs$temp_s_max <- hold$temp_s_max
env_20m_qcs$temp_s_min <- hold$temp_s_min
env_20m_qcs$temp_5m_mean <- hold$temp_5m_mean
env_20m_qcs$temp_5m_max <- hold$temp_5m_max
env_20m_qcs$temp_5m_min <- hold$temp_5m_min
env_20m_qcs$do_5m_mean <- hold$do_5m_mean
env_20m_qcs$do_5m_min <- hold$do_5m_min
env_20m_qcs <- env_20m_qcs %>% filter(do_5m_min >= 100, temp_s_min >= 0, !is.na(NH4_5m_mean)) # areas where BCCM doesn't make good predictions
env_20m_qcs$region <- "Queen Charlotte Strait"
env_20m_qcs$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_qcs$substrate]

summary(env_20m_qcs)

#save(env_20m_qcs, file = "./data/prediction_QCS.RData")
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
env_20m_wcvi <- filter(env_20m_wcvi, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_wcvi_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_wcvi$NH4_5m_mean <- hold$NH4_5m_mean
env_20m_wcvi$NO3_5m_mean <- hold$NO3_5m_mean
env_20m_wcvi$salt_5m_mean <- hold$salt_5m_mean
env_20m_wcvi$salt_5m_min <- hold$salt_5m_min
env_20m_wcvi$PAR_5m_mean <- hold$PAR_5m_mean
env_20m_wcvi$PAR_5m_min <- hold$PAR_5m_min
env_20m_wcvi$PAR_5m_max <- hold$PAR_5m_max
env_20m_wcvi$temp_s_mean <- hold$temp_s_mean
env_20m_wcvi$temp_s_max <- hold$temp_s_max
env_20m_wcvi$temp_s_min <- hold$temp_s_min
env_20m_wcvi$temp_5m_mean <- hold$temp_5m_mean
env_20m_wcvi$temp_5m_max <- hold$temp_5m_max
env_20m_wcvi$temp_5m_min <- hold$temp_5m_min
env_20m_wcvi$do_5m_mean <- hold$do_5m_mean
env_20m_wcvi$do_5m_min <- hold$do_5m_min
env_20m_wcvi <- env_20m_wcvi %>% filter(do_5m_min >= 100, temp_s_min >= 0, !is.na(NH4_5m_mean)) # areas where BCCM doesn't make good predictions
env_20m_wcvi$region <- "West Coast Vancouver Island"
env_20m_wcvi$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_wcvi$substrate]

summary(env_20m_wcvi)

#save(env_20m_wcvi, file = "./data/prediction_WCVI.RData")
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
env_20m_ss <- filter(env_20m_ss, !is.na(rei), !is.na(substrate), !is.na(tidal))
env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = hindcast2013_2023, y = env_20m_ss_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
env_20m_ss$NH4_5m_mean <- hold$NH4_5m_mean
env_20m_ss$NO3_5m_mean <- hold$NO3_5m_mean
env_20m_ss$salt_5m_mean <- hold$salt_5m_mean
env_20m_ss$salt_5m_min <- hold$salt_5m_min
env_20m_ss$PAR_5m_mean <- hold$PAR_5m_mean
env_20m_ss$PAR_5m_min <- hold$PAR_5m_min
env_20m_ss$PAR_5m_max <- hold$PAR_5m_max
env_20m_ss$temp_s_mean <- hold$temp_s_mean
env_20m_ss$temp_s_max <- hold$temp_s_max
env_20m_ss$temp_s_min <- hold$temp_s_min
env_20m_ss$temp_5m_mean <- hold$temp_5m_mean
env_20m_ss$temp_5m_max <- hold$temp_5m_max
env_20m_ss$temp_5m_min <- hold$temp_5m_min
env_20m_ss$do_5m_mean <- hold$do_5m_mean
env_20m_ss$do_5m_min <- hold$do_5m_min
env_20m_ss <- env_20m_ss %>% filter(do_5m_min >= 100, temp_s_min >= 0, !is.na(NH4_5m_mean)) # areas where BCCM doesn't make good predictions
env_20m_ss$region <- "Salish Sea"
env_20m_ss$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_ss$substrate]
summary(env_20m_ss)

rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "SpatRaster"])

#save(env_20m_ss, file = "./data/prediction_SS.RData")

ggplot(env_20m_ss, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()


#combine all and filter
env_20m_all <- bind_rows(env_20m_hg, env_20m_ncc, env_20m_qcs, env_20m_wcvi,env_20m_ss)
summary(env_20m_all)

# ggplot(env_20m_all, aes(x = X_m, y = Y_m, color = substrate))+
#   geom_point(size = 0.01, pch = 15)+
#   theme_bw()+
#   coord_equal()

#need to think about this in reference to how much we want to extrapolate, will depend on future conditions
env_20m_all <- env_20m_all %>%
  filter(freshwater < quantile(seagrass_data$freshwater, probs = 0.99),  # we have no surveys in freshwater areas, so better to exclude
         salt_5m_mean > quantile(seagrass_data$saltmean, probs = 0.001), # same as above
         #temperature > quantile(seagrass_data$temperature, probs = 0.001), # want to allow extrapolation into higher temperatures
         #temperature < quantile(seagrass_data$temperature, probs = 0.999),
         tidal < quantile(seagrass_data$tidal, probs = 0.9999), # remove high current areas where we can't sample, neither species observed in high values in quadrats 
         rei < quantile(seagrass_data$rei, probs = 0.99999))   # remove high exposure areas where we can't sample? PH found in moderately high but not
# goes from 15.1 million prediction cells to 14.3 million

summary(env_20m_all)

#scale environmental variables function####
scale_fun_ref <- function(x, reference){
  (x  - mean(reference)) / sd(reference)
}

# 
# 
# 
# predict_grid <- env.grid %>% 
#   mutate(depth_ln_stnd = scale_fun_ref(log(depth), reference = log(all_sets_wide$depth_m)),
#          salinity_stnd = scale_fun_ref(salinity, reference = all_sets_wide$salinity),
#          salinity_range_stnd = scale_fun_ref(salinity_range, reference = all_sets_wide$salinity_range),
#          tidal_ln_stnd = scale_fun_ref(log(tidal), reference = log(all_sets_wide$tidal)),
#          circ_stnd = scale_fun_ref(circulation, reference = all_sets_wide$circulation),
#          BPI_stnd = scale_fun_ref(BBPI, reference = all_sets_wide$BBPI),
#          rocky_stnd = scale_fun_ref(rocky, reference = all_sets_wide$rocky), 
#          muddy_stnd = scale_fun_ref(muddy, reference = all_sets_wide$muddy))

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
         NH4_stnd = scale_fun_ref(NH4_5m_mean, reference = seagrass_data$NH4),
         NO3_stnd = scale_fun_ref(NO3_5m_mean, reference = seagrass_data$NO3),
         saltmean_stnd = scale_fun_ref(salt_5m_mean, reference = seagrass_data$saltmean),
         saltmean_sq_stnd = scale_fun_ref((salt_5m_mean)^2, reference = (seagrass_data$saltmean)^2),
         saltmin_stnd = scale_fun_ref(salt_5m_min, reference = seagrass_data$saltmin),
         saltmin_sq_stnd = scale_fun_ref((salt_5m_min)^2, reference = (seagrass_data$saltmin)^2),
         PARmean_stnd = scale_fun_ref(PAR_5m_mean, reference = seagrass_data$PARmean),
         PARmin_stnd = scale_fun_ref(PAR_5m_min, reference = seagrass_data$PARmin),
         PARmax_stnd = scale_fun_ref(PAR_5m_max, reference = seagrass_data$PARmax),
         surftempmean_stnd = scale_fun_ref(temp_s_mean, reference = seagrass_data$surftempmean),
         surftempmin_stnd = scale_fun_ref(temp_s_min, reference = seagrass_data$surftempmin),
         surftempmax_stnd = scale_fun_ref(temp_s_max, reference = seagrass_data$surftempmax),
         tempmean_stnd = scale_fun_ref(temp_5m_mean, reference = seagrass_data$tempmean),
         tempmin_stnd = scale_fun_ref(temp_5m_min, reference = seagrass_data$tempmin),
         tempmax_stnd = scale_fun_ref(temp_5m_max, reference = seagrass_data$tempmax),
         DOmean_stnd = scale_fun_ref(do_5m_mean, reference = seagrass_data$DOmean),
         DOmin_stnd = scale_fun_ref(do_5m_min, reference = seagrass_data$DOmin)) %>%
  mutate(ID = 1:nrow(env_20m_all),
         HKey = factor("new"),
         X = X_m/1000,
         Y = Y_m/1000,
         substrate = factor(substrate))

summary(env_20m_all)
  
#save outputs####
save(env_20m_all, file = "code/output_data/prediction_model_inputs.RData")

  

#### Compare the sampled predictor space to that of the total hindcast climatology
#2013-2023 hindcast climatology without depth and substrate and turn into a dataframe
hindcast_predictor_data <- env_20m_all %>% select(rei_stnd, tidal_stnd, freshwater_stnd, slope_stnd, NH4_stnd, NO3_stnd, saltmean_stnd, saltmin_stnd,
                                                  PARmean_stnd, PARmin_stnd, tempmin_stnd, tempmax_stnd, DOmean_stnd, DOmin_stnd)

# transect predictor data
transect_predictor_data <- seagrass_data %>% select(rei_stnd, tidal_stnd, freshwater_stnd, slope_stnd, NH4_stnd, NO3_stnd, saltmean_stnd, saltmin_stnd,
                                                    PARmean_stnd, PARmin_stnd, tempmin_stnd, tempmax_stnd, DOmean_stnd, DOmin_stnd)

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
Density_Plot <- function(All_Scores, comp, title, x, y){
  Dense <- ggplot(All_Scores) + 
    geom_density(aes(x = All_Scores[,comp], fill = Type), alpha = 0.5) + 
    scale_fill_manual(values = c("#404788FF", "#55C667FF"), name = "Feature Space" )  + 
    scale_y_continuous(expand = c(0, 0)) +
    xlab(title) +
    ylab("Frequency") + 
    theme_classic()   +
    guides(x.sec = "axis", y.sec = "axis")+
    theme(axis.line.y.left = element_line(size = 0.75, color = "black"),
          axis.line.y.right = element_line(size = 0.75, color = "black"),
          axis.line.x.top = element_line(size = 0.75, color = "black"),
          axis.line.x.bottom = element_line(size = 0.75, color = "black"),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank())  +
    theme(legend.position = c(x,y)) +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 8),
          legend.background = element_blank())
  return(Dense)
}

#PCA_dense_plot <-  gridExtra::grid.arrange(Dense1, Dense2, nrow = 2)
#Plot grid a density plot for the first four principal components
PCA_dense_plot <- gridExtra::grid.arrange(Density_Plot(All_Scores, 1, "Principal Component 1 Value (39%)", 0.75, 0.81), 
                                     Density_Plot(All_Scores, 2, "Principal Component 2 Value (19%)", 0.25, 0.81), 
                                     Density_Plot(All_Scores, 3, "Principal Component 3 Value (15%)", 0.75, 0.81), 
                                     Density_Plot(All_Scores, 4, "Principal Component 4 Value (9%)", 0.78, 0.81), nrow = 2)
PCA_dense_plot
ggsave(filename = "./figures/pre-analysis/Predictor_Hindcast_Transect_PCA_Density_Plot.png", plot = PCA_dense_plot, scale = 1, height = 10, width = 12)

#Create a PCA plot of the first two principal components of our hindcast and survey predictor space with no points but 95% elipses
PCA_Plot_Ellipses <- ggplot(All_Scores) + 
  scale_colour_manual( values = c("#404788FF", "#55C667FF")) + 
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
  scale_fill_manual( values = c("#404788FF", "#55C667FF")) + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size= 18), 
        legend.background = element_blank())

PCA_Plot_Ellipses
ggsave("./figures/pre-analysis/Predictor_Hindcast_Transect_PCA_Plot_Ellipses.png", plot = PCA_Plot_Ellipses, width = 10, height = 12)
