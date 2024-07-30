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

temperature_hg <- rast("raw_data/envlayers-20m-hg/temp_mean_summer.tif")
temperature_ncc <- rast("raw_data/envlayers-20m-ncc/temp_mean_summer.tif")
temperature_qcs <- rast("raw_data/envlayers-20m-qcs/temp_mean_summer.tif")
temperature_wcvi <- rast("raw_data/envlayers-20m-wcvi/temp_mean_summer.tif")
temperature_ss <- rast("raw_data/envlayers-20m-shelfsalishsea/temp_mean_summer.tif")

salinity_hg <- rast("raw_data/envlayers-20m-hg/salinity_mean_summer.tif")
salinity_ncc <- rast("raw_data/envlayers-20m-ncc/salinity_mean_summer.tif")
salinity_qcs <- rast("raw_data/envlayers-20m-qcs/salinity_mean_summer.tif")
salinity_wcvi <- rast("raw_data/envlayers-20m-wcvi/salinity_mean_summer.tif")
salinity_ss <- rast("raw_data/envlayers-20m-shelfsalishsea/salinity_mean_summer.tif")

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
names(tidal_all)<-"tidal_mean_summer"
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

####make prediction data####
scale <- 100 #scale of environmental predictions
max_depth <- quantile(seagrass_data$depth, probs = c(0.99))
min_depth <- quantile(seagrass_data$depth, probs = c(0.01))

#haida gwaii
env_20m_hg_all <- as.data.frame(bathy_hg, xy=TRUE)
names(env_20m_hg_all) <- c("X_m", "Y_m", "depth")
#env_20m_hg_all <- env_20m_hg_all %>% filter(depth < max_depth, depth > min_depth)

all_y <- unique(env_20m_hg_all$Y_m)[order(unique(env_20m_hg_all$Y_m))]
select_y <- all_y[seq(1, length(all_y), by = scale/20)]
all_x <- unique(env_20m_hg_all$X_m)[order(unique(env_20m_hg_all$X_m))]
select_x <- all_x[seq(1, length(all_x), by = scale/20)]

env_20m_hg <- env_20m_hg_all %>%
  filter(Y_m %in% select_y, X_m %in% select_x)

env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m","Y_m"), crs = "EPSG:3005")
hold <- extract(x = rei_hg, y = env_20m_hg_sf)
env_20m_hg$rei <- hold$rei_hg
hold <- extract(x = substrate_hg, y = env_20m_hg_sf)
env_20m_hg$substrate <- hold$substrate
env_20m_hg <- filter(env_20m_hg, !is.na(rei), !is.na(substrate))
env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = slope_hg, y = env_20m_hg_sf)
env_20m_hg$slope <- hold$slope
hold <- extract(x = salinity_hg, y = env_20m_hg_sf)
env_20m_hg$salinity <- hold$salinity_mean_summer
hold <- extract(x = temperature_hg, y = env_20m_hg_sf)
env_20m_hg$temperature <- hold$temp_mean_summer
hold <- extract(x = tidal_index_all, y = env_20m_hg_sf)
env_20m_hg$tidal <- hold$tidal_mean_summer
hold <- extract(x = freshwater_hg, y = env_20m_hg_sf)
env_20m_hg$freshwater <- hold$freshwater
env_20m_hg$region <- "Haida Gwaii"
env_20m_hg$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_hg$substrate]
summary(env_20m_hg)

#save(env_20m_hg, file = "./data/prediction_HG.RData")

ggplot(env_20m_hg, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()

#north central coast
env_20m_ncc_all <- as.data.frame(bathy_ncc, xy=TRUE)
names(env_20m_ncc_all) <- c("X_m", "Y_m", "depth")
#env_20m_ncc_all <- env_20m_ncc_all %>% filter(depth < max_depth, depth > min_depth)

all_y <- unique(env_20m_ncc_all$Y_m)[order(unique(env_20m_ncc_all$Y_m))]
select_y <- all_y[seq(1, length(all_y), by = scale/20)]
all_x <- unique(env_20m_ncc_all$X_m)[order(unique(env_20m_ncc_all$X_m))]
select_x <- all_x[seq(1, length(all_x), by = scale/20)]

env_20m_ncc <- env_20m_ncc_all %>%
  filter(Y_m %in% select_y, X_m %in% select_x)

env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = rei_ncc, y = env_20m_ncc_sf)
env_20m_ncc$rei <- hold$rei_ncc
hold <- extract(x = substrate_ncc, y = env_20m_ncc_sf)
env_20m_ncc$substrate <- hold$substrate
env_20m_ncc <- filter(env_20m_ncc, !is.na(rei), !is.na(substrate))
env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = slope_ncc, y = env_20m_ncc_sf)
env_20m_ncc$slope <- hold$slope
hold <- extract(x = salinity_ncc, y = env_20m_ncc_sf)
env_20m_ncc$salinity <- hold$salinity_mean_summer
hold <- extract(x = temperature_ncc, y = env_20m_ncc_sf)
env_20m_ncc$temperature <- hold$temp_mean_summer
hold <- extract(x = tidal_index_all, y = env_20m_ncc_sf)
env_20m_ncc$tidal <- hold$tidal_mean_summer
hold <- extract(x = freshwater_ncc, y = env_20m_ncc_sf)
env_20m_ncc$freshwater <- hold$freshwater
env_20m_ncc$region <- "North Central Coast"
env_20m_ncc$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_ncc$substrate]

summary(env_20m_ncc)

#save(env_20m_ncc, file = "./data/prediction_NCC.RData")
ggplot(env_20m_ncc, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()

#queen charlotte strait
env_20m_qcs_all <- as.data.frame(bathy_qcs, xy=TRUE)
names(env_20m_qcs_all) <- c("X_m", "Y_m", "depth")
#env_20m_qcs_all <- env_20m_qcs_all %>% filter(depth < max_depth, depth > min_depth)

all_y <- unique(env_20m_qcs_all$Y_m)[order(unique(env_20m_qcs_all$Y_m))]
select_y <- all_y[seq(1, length(all_y), by = scale/20)]
all_x <- unique(env_20m_qcs_all$X_m)[order(unique(env_20m_qcs_all$X_m))]
select_x <- all_x[seq(1, length(all_x), by = scale/20)]

env_20m_qcs <- env_20m_qcs_all %>%
  filter(Y_m %in% select_y, X_m %in% select_x)

env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = rei_qcs, y = env_20m_qcs_sf)
env_20m_qcs$rei <- hold$rei_qcs
hold <- extract(x = substrate_qcs, y = env_20m_qcs_sf)
env_20m_qcs$substrate <- hold$substrate
env_20m_qcs <- filter(env_20m_qcs, !is.na(rei), !is.na(substrate))
env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = slope_qcs, y = env_20m_qcs_sf)
env_20m_qcs$slope <- hold$slope
hold <- extract(x = salinity_qcs, y = env_20m_qcs_sf)
env_20m_qcs$salinity <- hold$salinity_mean_summer
hold <- extract(x = temperature_qcs, y = env_20m_qcs_sf)
env_20m_qcs$temperature <- hold$temp_mean_summer
hold <- extract(x = tidal_index_all, y = env_20m_qcs_sf)
env_20m_qcs$tidal <- hold$tidal_mean_summer
hold <- extract(x = freshwater_qcs, y = env_20m_qcs_sf)
env_20m_qcs$freshwater <- hold$freshwater
env_20m_qcs$region <- "Queen Charlotte Strait"
env_20m_qcs$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_qcs$substrate]

summary(env_20m_qcs)

#save(env_20m_qcs, file = "./data/prediction_QCS.RData")
ggplot(env_20m_qcs, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()



#west coast vancouver island
env_20m_wcvi_all <- as.data.frame(bathy_wcvi, xy=TRUE)
names(env_20m_wcvi_all) <- c("X_m", "Y_m", "depth")
#env_20m_wcvi_all <- env_20m_wcvi_all %>% filter(depth < max_depth, depth > min_depth)

all_y <- unique(env_20m_wcvi_all$Y_m)[order(unique(env_20m_wcvi_all$Y_m))]
select_y <- all_y[seq(1, length(all_y), by = scale/20)]
all_x <- unique(env_20m_wcvi_all$X_m)[order(unique(env_20m_wcvi_all$X_m))]
select_x <- all_x[seq(1, length(all_x), by = scale/20)]

env_20m_wcvi <- env_20m_wcvi_all %>%
  filter(Y_m %in% select_y, X_m %in% select_x)

env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = rei_wcvi, y = env_20m_wcvi_sf)
env_20m_wcvi$rei <- hold$rei_wcvi
hold <- extract(x = substrate_wcvi, y = env_20m_wcvi_sf)
env_20m_wcvi$substrate <- hold$substrate
env_20m_wcvi <- filter(env_20m_wcvi, !is.na(rei), !is.na(substrate))
env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = slope_wcvi, y = env_20m_wcvi_sf)
env_20m_wcvi$slope <- hold$slope
hold <- extract(x = salinity_wcvi, y = env_20m_wcvi_sf)
env_20m_wcvi$salinity <- hold$salinity_mean_summer
hold <- extract(x = temperature_wcvi, y = env_20m_wcvi_sf)
env_20m_wcvi$temperature <- hold$temp_mean_summer
hold <- extract(x = tidal_index_all, y = env_20m_wcvi_sf)
env_20m_wcvi$tidal <- hold$tidal_mean_summer
hold <- extract(x = freshwater_wcvi, y = env_20m_wcvi_sf)
env_20m_wcvi$freshwater <- hold$freshwater
env_20m_wcvi$region <- "West Coast Vancouver Island"
env_20m_wcvi$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_wcvi$substrate]

summary(env_20m_wcvi)

#save(env_20m_wcvi, file = "./data/prediction_WCVI.RData")
ggplot(env_20m_wcvi, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()


#salish sea
env_20m_ss_all <- as.data.frame(bathy_ss, xy=TRUE)
names(env_20m_ss_all) <- c("X_m", "Y_m", "depth")
#env_20m_ss_all <- env_20m_ss_all %>% filter(depth < max_depth, depth > min_depth)

all_y <- unique(env_20m_ss_all$Y_m)[order(unique(env_20m_ss_all$Y_m))]
select_y <- all_y[seq(1, length(all_y), by = scale/20)]
all_x <- unique(env_20m_ss_all$X_m)[order(unique(env_20m_ss_all$X_m))]
select_x <- all_x[seq(1, length(all_x), by = scale/20)]

env_20m_ss <- env_20m_ss_all %>%
  filter(Y_m %in% select_y, X_m %in% select_x)

env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = rei_ss, y = env_20m_ss_sf)
env_20m_ss$rei <- hold$rei_sog
hold <- extract(x = substrate_ss, y = env_20m_ss_sf)
env_20m_ss$substrate <- hold$substrate
env_20m_ss <- filter(env_20m_ss, !is.na(rei), !is.na(substrate))
env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m", "Y_m"), crs = "EPSG:3005")
hold <- extract(x = slope_ss, y = env_20m_ss_sf)
env_20m_ss$slope <- hold$slope
hold <- extract(x = salinity_ss, y = env_20m_ss_sf)
env_20m_ss$salinity <- hold$salinity_mean_summer
hold <- extract(x = temperature_ss, y = env_20m_ss_sf)
env_20m_ss$temperature <- hold$temp_mean_summer
hold <- extract(x = tidal_index_all, y = env_20m_ss_sf)
env_20m_ss$tidal <- hold$tidal_mean_summer
hold <- extract(x = freshwater_ss, y = env_20m_ss_sf)
env_20m_ss$freshwater <- hold$freshwater
env_20m_ss$region <- "Salish Sea"
env_20m_ss$substrate <- c("Rock", "Mixed", "Sand", "Mud")[env_20m_ss$substrate]
summary(env_20m_ss)

#save(env_20m_ss, file = "./data/prediction_SS.RData")

ggplot(env_20m_ss, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()


#combine all and filter
env_20m_all <- bind_rows(env_20m_hg, env_20m_ncc, env_20m_qcs, env_20m_wcvi,env_20m_ss)
summary(env_20m_all)

ggplot(env_20m_all, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()

#need to think about this in reference to how much we want to extrapolate, will depend on future conditions
# env_20m_all <- env_20m_all %>%
#   filter(freshwater < quantile(seagrass_data$freshwater, probs = 0.99),
#          salinity > quantile(seagrass_data$salinity, probs = 0.001),
#          #temperature > quantile(seagrass_data$temperature, probs = 0.001),
#          #temperature < quantile(seagrass_data$temperature, probs = 0.999),
#          #tidal < quantile(seagrass_data$tidal, probs = 0.999),
#          tidal >= 0,
#          rei < quantile(seagrass_data$rei, probs = 0.999),
#          !is.na(freshwater))


summary(env_20m_all)

#standardize prediction data####
scale_ref <- function(x, ref){
  mean_ref <- mean(ref)
  sd_ref <- sd(ref)
  scaled_x <- (x - mean_ref) / sd_ref
  return(scaled_x)
}

env_20m_all <- env_20m_all %>%
  mutate(tidal_sqrt = sqrt(tidal), rei_sqrt = sqrt(rei)) %>%
  mutate(depth_stnd = scale_ref(depth, seagrass_data$depth),
         slope_stnd = scale_ref(slope, seagrass_data$slope),
         rei_stnd = scale_ref(rei, seagrass_data$rei),
         rei_sqrt_stnd = scale_ref(rei_sqrt, seagrass_data$rei_sqrt),
         temperature_stnd = scale_ref(temperature, seagrass_data$temperature),
         salinity_stnd = scale_ref(salinity, seagrass_data$salinity),
         tidal_sqrt_stnd = scale_ref(tidal_sqrt, seagrass_data$tidal_sqrt)
  ) %>%
  mutate(ID = 1:nrow(env_20m_all),
         HKey = factor("new"),
         X = X_m/1000,
         Y = Y_m/1000,
         substrate = factor(substrate))
summary(env_20m_all)

ggplot(env_20m_all, aes(x = X_m, y = Y_m, color = substrate))+
  geom_point(size = 0.01, pch = 15)+
  theme_bw()+
  coord_equal()


#save outputs####
save(env_20m_all, file = "code/output_data/prediction_model_inputs.RData")
