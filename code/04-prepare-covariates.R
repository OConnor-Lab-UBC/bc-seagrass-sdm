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
# 99,397 obs at 20 m resolution

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

substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"
#plot(substrate_all)

freshwater_all <- vrt(c("raw_data/freshwater-index/hg_freshwater_index.tif", "raw_data/freshwater-index/ncc_freshwater_index.tif", "raw_data/freshwater-index/qcs_freshwater_index.tif", "raw_data/freshwater-index/wcvi_freshwater_index.tif", "raw_data/freshwater-index/salish_sea_freshwater_index.tif"), "freshwater.vrt", overwrite=T)
#plot(freshwater_all)

bathy_all <- vrt(c("raw_data/envlayers-20m-hg//bathymetry.tif", "raw_data/envlayers-20m-ncc/bathymetry.tif", "raw_data/envlayers-20m-qcs/bathymetry.tif", "raw_data/envlayers-20m-wcvi/bathymetry.tif", "raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif"), "bathy.vrt", overwrite=T)


tidal_all <- vrt(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))
names(tidal_all)<-"tidal_current_index"
#change to index 0-1 scale
tidal_index_all <- tidal_all/(maxFn(tidal_all))
crs(tidal_index_all) <- "EPSG:3005"

#read in oceanographic data
# Set folder path
folder_1993_2002 <- "code/output_data/processed_ocean_variables/years_1993-2002"
folder_2003_2012 <- "code/output_data/processed_ocean_variables/years_2003-2012"
folder_2013_2023 <- "code/output_data/processed_ocean_variables/years_2013-2023"

files_1993_2002 <- list.files(folder_1993_2002, pattern = "\\.tif$", full.names = TRUE)
files_2003_2012 <- list.files(folder_2003_2012, pattern = "\\.tif$", full.names = TRUE)
files_2013_2023 <- list.files(folder_2013_2023, pattern = "\\.tif$", full.names = TRUE)

hindcast1993_2002 <- terra::rast(files_1993_2002)
hindcast2003_2012 <- terra::rast(files_2003_2012)
hindcast2013_2023 <- terra::rast(files_2013_2023)

# read in cumulative effects
culmulative_effects <- terra::rast("code/output_data/processed_ocean_variables/culmulative_effects_all_20m.tif")


#extract from predictor layers
culeff_extract <- terra::extract(culmulative_effects, spatialised_sf)
summary(culeff_extract$focal_mean)
spatialised_sf$cul_eff <- culeff_extract$focal_mean

bathy_extract <- terra::extract(bathy_all, spatialised_sf)
summary(bathy_extract$bathy)
spatialised_sf$bathy <- bathy_extract$bathy

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
         substrate = ifelse(is.na(Substrate), substrate_mod, Substrate),
         depth = ifelse(is.na(CorDepthM), bathy, CorDepthM)) %>%
  select(-c(Slope, Substrate, slope_mod, substrate_mod, CorDepthM, bathy))

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
  rename(DOmean_bccm = do_5m_mean_bccm, DOmean_nep = do_5m_mean_nep36, DOmin_bccm = do_5m_min_bccm, DOmin_nep = do_5m_min_nep36,
         NH4_bccm = NH4_5m_mean_bccm, NH4_nep = NH4_5m_mean_nep36, NO3_bccm = NO3_5m_mean_bccm, NO3_nep = NO3_5m_mean_nep36,
         PARmean_bccm = PAR_5m_mean_bccm, PARmean_nep = PAR_5m_mean_nep36, PARmin_bccm = PAR_5m_min_bccm, PARmin_nep = PAR_5m_min_nep36, PARmax_bccm = PAR_5m_max_bccm, PARmax_nep = PAR_5m_max_nep36,
         prcv = precip_cv, prmax = precip_max, prmin = precip_min, prmean = precip_mean,
         rsdscv = rsds_cv, rsdsmax = rsds_max, rsdsmean = rsds_mean, rsdsmin = rsds_min,
         saltmean_bccm = salt_5m_mean_bccm, saltmean_nep = salt_5m_mean_nep36, saltmin_bccm = salt_5m_min_bccm, saltmin_nep = salt_5m_min_nep36 , saltcv_bccm = salt_5m_cv_bccm,  saltcv_nep = salt_5m_cv_nep36,
         surftempmean_bccm = temp_s_mean_bccm, surftempmean_nep = temp_s_mean_nep36, surftempmax_bccm = temp_s_max_bccm, surftempmax_nep = temp_s_max_nep36, surftempmin_bccm = temp_s_min_bccm, surftempmin_nep = temp_s_min_nep36, surftempcv_bccm = temp_s_cv_bccm, surftempcv_nep = temp_s_cv_nep36, surftempdiff_bccm = temp_s_diff_bccm, surftempdiff_nep = temp_s_diff_nep36,
         tempmean_bccm = temp_5m_mean_bccm, tempmean_nep = temp_5m_mean_nep36, tempmax_bccm = temp_5m_max_bccm, tempmax_nep = temp_5m_max_nep36, tempmin_bccm = temp_5m_min_bccm, tempmin_nep = temp_5m_min_nep36, tempcv_bccm = temp_5m_cv_bccm, tempcv_nep = temp_5m_cv_nep36, tempdiff_bccm = temp_5m_diff_bccm, tempdiff_nep = temp_5m_diff_nep36,
         airtempcv = temp_air_cv, airtempmax = temp_air_max, airtempmean = temp_air_mean, airtempmin = temp_air_min) %>%
  select(ID, NH4_bccm, NH4_nep, NO3_bccm, NO3_nep, saltmean_bccm, saltmean_nep, saltmin_bccm, saltmin_nep, saltcv_bccm, saltcv_nep,
         PARmean_bccm, PARmean_nep, PARmin_bccm, PARmin_nep, PARmax_bccm, PARmax_nep, surftempmean_bccm, surftempmean_nep, 
         surftempmax_bccm, surftempmax_nep, surftempmin_bccm, surftempmin_nep, surftempcv_bccm, surftempcv_nep, surftempdiff_bccm, surftempdiff_nep,
         tempmean_bccm, tempmean_nep, tempmax_bccm, tempmax_nep, tempmin_bccm, tempmin_nep, tempcv_bccm, tempcv_nep, tempdiff_bccm, tempdiff_nep,
         airtempcv, airtempmax, airtempmean, airtempmin,
         DOmean_bccm, DOmean_nep, DOmin_bccm, DOmin_nep, prcv, prmax, prmin, prmean, rsdscv, rsdsmax, rsdsmean, rsdsmin)

spatialised_sf <- dplyr::full_join(spatialised_sf, oceanvars_allyears, by=c("ID")) %>% filter(!is.na(NH4_bccm), !is.na(NH4_nep))
spatialised_sf$mean_PerCovZO<- round(spatialised_sf$mean_PerCovZO, digits = 0)
#end with 92,226 obs

summary(spatialised_sf)

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
         NH4_bccm_stnd = scale_fun(NH4_bccm),
         NH4_nep_stnd = scale_fun(NH4_nep),
         NO3_bccm_stnd = scale_fun(NO3_bccm),
         NO3_nep_stnd = scale_fun(NO3_nep),
         saltmean_bccm_stnd = scale_fun(saltmean_bccm),
         saltmean_nep_stnd = scale_fun(saltmean_nep),
         saltmean_bccm_sq_stnd = scale_fun((saltmean_bccm)^2),
         saltmean_nep_sq_stnd = scale_fun((saltmean_nep)^2),
         saltmin_bccm_stnd = scale_fun(saltmin_bccm),
         saltmin_nep_stnd = scale_fun(saltmin_nep),
         saltmin_bccm_sq_stnd = scale_fun((saltmin_bccm)^2),
         saltmin_nep_sq_stnd = scale_fun((saltmin_nep)^2),
         saltcv_bccm_stnd = scale_fun(saltcv_bccm),
         saltcv_nep_stnd = scale_fun(saltcv_nep),
         PARmean_bccm_stnd = scale_fun(PARmean_bccm),
         PARmean_nep_stnd = scale_fun(PARmean_nep),
         PARmin_bccm_stnd = scale_fun(PARmin_bccm),
         PARmin_nep_stnd = scale_fun(PARmin_nep),
         PARmax_bccm_stnd = scale_fun(PARmax_bccm),
         PARmax_nep_stnd = scale_fun(PARmax_nep),
         surftempmean_bccm_stnd = scale_fun(surftempmean_bccm),
         surftempmean_nep_stnd = scale_fun(surftempmean_nep),
         surftempmin_bccm_stnd = scale_fun(surftempmin_bccm),
         surftempmin_nep_stnd = scale_fun(surftempmin_nep),
         surftempmax_bccm_stnd = scale_fun(surftempmax_bccm),
         surftempmax_nep_stnd = scale_fun(surftempmax_nep),
         surftempcv_bccm_stnd = scale_fun(surftempcv_bccm),
         surftempcv_nep_stnd = scale_fun(surftempcv_nep),
         surftempdiff_bccm_stnd = scale_fun(surftempdiff_bccm),
         surftempdiff_nep_stnd = scale_fun(surftempdiff_nep),
         tempmean_bccm_stnd = scale_fun(tempmean_bccm),
         tempmean_nep_stnd = scale_fun(tempmean_nep),
         tempmin_bccm_stnd = scale_fun(tempmin_bccm),
         tempmin_nep_stnd = scale_fun(tempmin_nep),
         tempmax_bccm_stnd = scale_fun(tempmax_bccm),
         tempmax_nep_stnd = scale_fun(tempmax_nep),
         tempcv_bccm_stnd = scale_fun(tempcv_bccm),
         tempcv_nep_stnd = scale_fun(tempcv_nep),
         tempdiff_bccm_stnd = scale_fun(tempdiff_bccm),
         tempdiff_nep_stnd = scale_fun(tempdiff_nep),
         DOmean_bccm_stnd = scale_fun(DOmean_bccm),
         DOmean_nep_stnd = scale_fun(DOmean_nep),
         DOmin_bccm_stnd = scale_fun(DOmin_bccm),
         DOmin_nep_stnd = scale_fun(DOmin_nep),
         airtempcv_stnd = scale_fun(airtempcv), 
         airtempmax_stnd = scale_fun(airtempmax), 
         airtempmean_stnd = scale_fun(airtempmean), 
         airtempmin_stnd = scale_fun(airtempmin),
         prcv_stnd = scale_fun(prcv), 
         prmax_stnd = scale_fun(prmax), 
         prmin_stnd = scale_fun(prmin), 
         prmean_stnd = scale_fun(prmean),
         rsdscv_stnd = scale_fun(rsdscv), 
         rsdsmax_stnd = scale_fun(rsdsmax), 
         rsdsmean_stnd = scale_fun(rsdsmean), 
         rsdsmin_stnd = scale_fun(rsdsmin),
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

