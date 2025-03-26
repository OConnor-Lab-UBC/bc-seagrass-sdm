###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
#
# Objective:
# ---------
# Select only observations in SDM prediction area ( <20m depth) add in predictor covariates. 
# Where we don't have substrate and slope observations from divers, add in modelled data
# Remove observations where no depth observation from divers. 
#
###############################################################################
# Load packages
library(sf)
library(tidyverse)
library(terra)
library(GGally)
library(blockCV)
library(reproducible)
library(caret)

# coastline
coastline_full <- st_read("raw_data/CHS_HWL2015_Coastline.gdb", layer = "Line_CHS_Pacific_HWL_2015_5028437")
coastline <- coastline_full %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_crop(st_bbox(coastline_full) + c(0, 100000, -10000, -80000)) %>%
  st_transform(crs = "EPSG:3005")

#load_data####
#load("code/output_data/seagrass_data_spatialized.RData")
# 400,597 quadrats

load("code/output_data/seagrass_data_spatialized_aggregated.RData")
spat <- spat %>%
  rename(Slope = mean_slope,
         CorDepthM = mean_CorDepthM)
# 95,428 obs at 20 m resolution

# check covariate observations
summary(spat$Slope) 
summary(spat$CorDepthM) 
sum(is.na(spat$Substrate)) 

#remove obs with no depth observations (if using non aggregated)
#spat <- spat %>%  filter(!is.na(CorDepthM))

# make sf object
spatialised_sf<- spat %>% st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

#read in 20m rasters####
rei_all <- terra::vrt(c("raw_data/REI/rei_hg.tif", "raw_data/REI/rei_ncc.tif", "raw_data/REI/rei_qcs.tif", "raw_data/REI/rei_sog.tif", "raw_data/REI/rei_wcvi.tif"), "rei.vrt", overwrite=T)
#plot(rei_all)

slope_all <- terra::vrt(c("raw_data/envlayers-20m-hg/slope.tif", "raw_data/envlayers-20m-ncc/slope.tif", "raw_data/envlayers-20m-qcs/slope.tif", "raw_data/envlayers-20m-wcvi/slope.tif", "raw_data/envlayers-20m-shelfsalishsea/slope.tif"), "slope.vrt", overwrite=T)
#plot(slope_all)

substrate_all <- terra::vrt(c("raw_data/substrate_20m/hg_20m.tif", "raw_data/substrate_20m/ncc_20m.tif", "raw_data/substrate_20m/qcs_20m.tif", "raw_data/substrate_20m/sog_20m.tif", "raw_data/substrate_20m/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"
#plot(substrate_all)

freshwater_all <- vrt(c("raw_data/freshwater-index/hg_freshwater_index.tif", "raw_data/freshwater-index/ncc_freshwater_index.tif", "raw_data/freshwater-index/qcs_freshwater_index.tif", "raw_data/freshwater-index/wcvi_freshwater_index.tif", "raw_data/freshwater-index/salish_sea_freshwater_index.tif"), "freshwater.vrt", overwrite=T)
#plot(freshwater_all)

tidal_all <- vrt(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))
names(tidal_all)<-"tidal_current_index"
#change to index 0-1 scale
tidal_index_all <- tidal_all/(maxFn(tidal_all))
crs(tidal_index_all) <- "EPSG:3005"

#read in oceanographic data
hindcast1993_2002 <- terra::rast("code/output_data/processed_ocean_variables/Predictor_Hindcast_Climatologies_1993-2002.tif")
hindcast2003_2012 <- terra::rast("code/output_data/processed_ocean_variables/Predictor_Hindcast_Climatologies_2003-2012.tif")
hindcast2013_2023 <- terra::rast("code/output_data/processed_ocean_variables/Predictor_Hindcast_Climatologies_2013-2023.tif")

# read in cumulative effects
culmulative_effects <- terra::rast("code/output_data/culmulative_effects_all_20m.tif")


#extract from predictor layers
culeff_extract <- terra::extract(culmulative_effects, spatialised_sf)
summary(culeff_extract$focal_mean)
spatialised_sf$cul_eff <- culeff_extract$focal_mean

rei_extract <- terra::extract(rei_all, spatialised_sf)
summary(rei_extract$rei)
spatialised_sf$rei <- rei_extract$rei

tidal_extract <- terra::extract(tidal_index_all, spatialised_sf)
summary(tidal_extract$tidal_current_index)
spatialised_sf$tidal <- tidal_extract$tidal_current_index

freshwater_extract <- terra::extract(freshwater_all, spatialised_sf)
summary(freshwater_extract$freshwater)
spatialised_sf$freshwater <- freshwater_extract$freshwater

slope_extract <- terra::extract(slope_all, spatialised_sf)
summary(slope_extract$slope)
spatialised_sf$slope_mod <- slope_extract$slope

substrate_extract <- terra::extract(substrate_all, spatialised_sf)
summary(substrate_extract$substrate)
spatialised_sf$substrate_mod <- substrate_extract$substrate
spatialised_sf$substrate_mod <- c("Rock", "Mixed", "Sand", "Mud")[spatialised_sf$substrate_mod]

spatialised_sf <- spatialised_sf %>%
  mutate(slope = ifelse(is.na(Slope)|Slope > 90, slope_mod, Slope),
         substrate = ifelse(is.na(Substrate), substrate_mod, Substrate)) %>%
  rename(depth = CorDepthM) %>%
  select(-c(Slope, Substrate, slope_mod, substrate_mod))

spatialised_sf$slope<- round(spatialised_sf$slope, digits = 0)

spatialised_sf <- spatialised_sf %>%
  filter(!is.na(rei), !is.na(substrate), !is.na(tidal), !is.na(cul_eff))

#ocean layers have to be by decade slice
spatialised_sf_1993_2002 <- spatialised_sf %>% filter(Year < 2003)
spatialised_sf_2003_2012 <- spatialised_sf %>% filter(Year > 2002 & Year < 2013)
spatialised_sf_2013_2023 <- spatialised_sf %>% filter(Year > 2012)

#Extract ocean predictor values 
oceanvars_1993_2002 <- terra::extract(x = hindcast1993_2002, y = spatialised_sf_1993_2002, fun = "mean", touches = TRUE, bind = TRUE) %>% 
  terra::as.data.frame() 
oceanvars_2003_2012 <- terra::extract(x = hindcast2003_2012, y = spatialised_sf_2003_2012, fun = "mean", touches = TRUE, bind = TRUE) %>% 
  terra::as.data.frame() 
oceanvars_2013_2023 <- terra::extract(x = hindcast2013_2023, y = spatialised_sf_2013_2023, fun = "mean", touches = TRUE, bind = TRUE) %>% 
  terra::as.data.frame() 

oceanvars_allyears <- rbind(oceanvars_1993_2002, oceanvars_2003_2012, oceanvars_2013_2023) %>%
  rename(NH4 = NH4_5m_mean, NO3 = NO3_5m_mean, saltmean = salt_5m_mean, saltmin = salt_5m_min, saltcv = salt_5m_cv, PARmean = PAR_5m_mean, PARmin = PAR_5m_min,
         PARmax = PAR_5m_max, surftempmean = temp_s_mean, surftempmax = temp_s_max, surftempmin = temp_s_min, surftempcv = temp_s_cv, surftempdiff = temp_s_diff, tempmean = temp_5m_mean, 
         tempmax = temp_5m_max, tempmin = temp_5m_min, tempcv= temp_5m_cv, tempdiff = temp_5m_diff, DOmean = do_5m_mean, DOmin = do_5m_min) %>%
  select(ID, NH4, NO3, saltmean, saltmin, saltcv, PARmean, PARmin, PARmax, surftempmean, surftempmax, surftempmin, surftempcv, surftempdiff, tempmean, tempmax, tempmin, 
         tempcv, tempdiff, DOmean, DOmin)

spatialised_sf <- dplyr::full_join(spatialised_sf, oceanvars_allyears, by=c("ID")) %>% filter(!is.na(NH4))
spatialised_sf$mean_PerCovZO<- round(spatialised_sf$mean_PerCovZO, digits = 0)
#end with 91,420 obs

scale_fun <- function(x){
  (x  - mean(x)) / sd(x)
}

spatialised_sf <- spatialised_sf %>% 
  mutate(depth_stnd = scale_fun(depth),
         rei_stnd = scale_fun(rei),
         rei_sqrt_stnd = scale_fun(sqrt(rei)),
         tidal_stnd = scale_fun(tidal),
         tidal_sqrt_stnd = scale_fun(sqrt(tidal)),
         freshwater_stnd = scale_fun(freshwater),
         freshwater_sqrt_stnd = scale_fun(sqrt(freshwater)),
         slope_stnd = scale_fun(slope),
         slope_sqrt_stnd = scale_fun(sqrt(slope)),
         NH4_stnd = scale_fun(NH4),
         NO3_stnd = scale_fun(NO3),
         saltmean_stnd = scale_fun(saltmean),
         saltmean_sq_stnd = scale_fun((saltmean)^2),
         saltmin_stnd = scale_fun(saltmin),
         saltmin_sq_stnd = scale_fun((saltmin)^2),
         saltcv_stnd = scale_fun(saltcv),
         PARmean_stnd = scale_fun(PARmean),
         PARmin_stnd = scale_fun(PARmin),
         PARmax_stnd = scale_fun(PARmax),
         surftempmean_stnd = scale_fun(surftempmean),
         surftempmin_stnd = scale_fun(surftempmin),
         surftempmax_stnd = scale_fun(surftempmax),
         surftempcv_stnd = scale_fun(surftempcv),
         surftempdiff_stnd = scale_fun(surftempdiff),
         tempmean_stnd = scale_fun(tempmean),
         tempmin_stnd = scale_fun(tempmin),
         tempmax_stnd = scale_fun(tempmax),
         tempcv_stnd = scale_fun(tempcv),
         tempdiff_stnd = scale_fun(tempdiff),
         DOmean_stnd = scale_fun(DOmean),
         DOmin_stnd = scale_fun(DOmin),
         cul_eff_stnd = scale_fun(cul_eff)) 

#ggpairs(spatialised_sf %>% dplyr::select(depth_stnd:DOmin_stnd) %>% st_set_geometry(NULL))
# keep sqrt rei, tidal, freshwater. Maybe sqrt slope. Don't need sqrt of NO3 and NH4
# log salinity make it worse, squaring it makes it marginally better, will include it to next stage but might not be worth keeping it
# high correlation between salt mean and min, pick one
# high correlation between DO mean and min, pick one
# the PARs are all moderately correlated so pick one. 
# temperatures are all correlated

#### Create spatial blocks for CV 

# spatial clustering
#make seperate fold for eelgrass to test versus one that could be used for both species
sp_blocks_eelgrass <- cv_spatial(x = spatialised_sf,
                        column = NULL,
                        k = 10, # 10 or more is recommended best practice Yates et al 2023
                        hexagon = FALSE,
                        size = 40000,  #in meters, matern range is 32km for eelgrass and 71km for surfgrass
                        selection = "random",
                        biomod2 = FALSE,
                        seed = 42, # to ensure reproducibility
                        plot = FALSE)

sp_blocks_seagrass <- cv_spatial(x = spatialised_sf,
                                 column = NULL,
                                 k = 10, # 10 or more is recommended best practice Yates et al 2023
                                 hexagon = FALSE,
                                 size = 75000,  #in meters, matern range is 32km for eelgrass and 71km for surfgrass
                                 selection = "random",
                                 biomod2 = FALSE,
                                 seed = 42, # to ensure reproducibility
                                 plot = FALSE)
spatialised_sf$fold_eelgrass <- sp_blocks_eelgrass$folds_ids
spatialised_sf$fold_seagrass <- sp_blocks_seagrass$folds_ids
#spatialised_sf$fold[is.na(spatialised_sf$fold)] <- 10
table(spatialised_sf$fold_eelgrass)
table(spatialised_sf$fold_seagrass)

#check that there are presences in each fold
ph <- spatialised_sf %>% filter(PH ==1) 
unique(ph$fold_seagrass)
zo <- spatialised_sf %>% filter(ZO ==1) 
unique(zo$fold_seagrass)
unique(zo$fold_eelgrass)

cv_plot_eelgrass <- cv_plot(cv = sp_blocks_eelgrass, #  blockCV object
        x = spatialised_sf, # sample points
        r = rei_all) #  raster background
cv_plot_eelgrass
ggsave("./figures/pre-analysis/spatial_blocks_eelgrass_type1.png", height = 6, width = 6)

sp_blocks_e <- ggplot(spatialised_sf)+
  geom_sf(aes(color = factor(fold_eelgrass)))+
  geom_sf(data = coastline)
sp_blocks_e
ggsave("./figures/pre-analysis/spatial_blocks_eelgrass_type2.png", height = 6, width = 6)

cv_plot_seagrass <- cv_plot(cv = sp_blocks_seagrass, #  blockCV object
                            x = spatialised_sf, # sample points
                            r = rei_all) #  raster background
cv_plot_seagrass
ggsave("./figures/pre-analysis/spatial_blocks_seagrass_type1.png", height = 6, width = 6)

sp_blocks_ss <- ggplot(spatialised_sf)+
  geom_sf(aes(color = factor(fold_seagrass)))+
  geom_sf(data = coastline)
sp_blocks_ss
ggsave("./figures/pre-analysis/spatial_blocks_seagrass_type2.png", height = 6, width = 6)

cv_eelgrass <- list()
for ( f in 1:10 ){
  foldname <- paste0("fold",f)
  cv_eelgrass[[foldname]][['train']] <- which(sp_blocks_eelgrass$folds_ids!=f)
  cv_eelgrass[[foldname]][['test']] <- which(sp_blocks_eelgrass$folds_ids==f)
}

#return blockpoly, foldID and CV list
cv_list_eelgrass<-list( foldID=sp_blocks_eelgrass$foldID, blockpolys=sp_blocks_eelgrass$blocks, cv=cv_eelgrass) 

cv_seagrass <- list()
for ( f in 1:10 ){
  foldname <- paste0("fold",f)
  cv_seagrass[[foldname]][['train']] <- which(sp_blocks_seagrass$folds_ids!=f)
  cv_seagrass[[foldname]][['test']] <- which(sp_blocks_seagrass$folds_ids==f)
}

#return blockpoly, foldID and CV list
cv_list_seagrass<-list( foldID=sp_blocks_seagrass$foldID, blockpolys=sp_blocks_seagrass$blocks, cv=cv_seagrass) 

#prepare data for model####
seagrass_data <- spatialised_sf %>% st_drop_geometry()
XY <- as.data.frame(st_coordinates(spatialised_sf))/1000
seagrass_data$X <- XY$X
seagrass_data$Y <- XY$Y
seagrass_data$X_m <- XY$X*1000
seagrass_data$Y_m <- XY$Y*1000

seagrass_data_long <- seagrass_data %>%
  select(HKey, ID, fold_eelgrass, fold_seagrass, Year, substrate, depth_stnd:cul_eff_stnd, ZO, PH, mean_PerCovZO, X:Y_m) %>%
  gather(key = species, value = presence, ZO:PH) %>%
  mutate(presence = ifelse(presence > 1, 1, presence)) %>%
  mutate(HKey = factor(HKey), substrate = factor(substrate))

#save outputs####
save(seagrass_data_long, seagrass_data, coastline, cv_list_eelgrass, cv_list_seagrass, file = "code/output_data/seagrass_model_inputs.RData")

