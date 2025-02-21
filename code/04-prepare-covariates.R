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


#extract from environmental layers
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
  filter(!is.na(rei), !is.na(substrate), !is.na(tidal))

#ocean layers have to be by decadal slice
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
  rename(NH4 = NH4_5m_mean, NO3 = NO3_5m_mean, saltmean = salt_5m_mean, saltmin = salt_5m_min, PARmean = PAR_5m_mean, PARmin = PAR_5m_min,
         PARmax = PAR_5m_max, surftempmean = temp_s_mean, surftempmax = temp_s_max, surftempmin = temp_s_min, tempmean = temp_5m_mean, 
         tempmax = temp_5m_max, tempmin = temp_5m_min, DOmean = do_5m_mean, DOmin = do_5m_min) %>%
  select(ID, NH4, NO3, saltmean, saltmin, PARmean, PARmin, PARmax, surftempmean, surftempmax, surftempmin, tempmean, tempmax, tempmin, 
         DOmean, DOmin)

spatialised_sf <- dplyr::full_join(spatialised_sf, oceanvars_allyears, by=c("ID")) %>% filter(!is.na(NH4))
#end with 91,420 obs

#START HERE WITH EDITS
#scale environmental variables####
scale_fun_ref <- function(x, reference){
  (x  - mean(reference)) / sd(reference)
}

scale_fun <- function(x){
  (x  - mean(x)) / sd(x)
}

all_sets_wide <- all_sets_wide %>% 
  mutate(depth_ln_stnd = scale_fun(log(depth_m)),
         salinity_stnd = scale_fun(salinity),
         salinity_range_stnd = scale_fun(salinity_range),
         tidal_ln_stnd = scale_fun(log(tidal)),
         circ_stnd = scale_fun(circulation),
         BPI_stnd = scale_fun(BBPI),
         rocky_stnd = scale_fun(rocky), 
         muddy_stnd = scale_fun(muddy))

ggpairs(all_sets_wide %>% dplyr::select(depth_ln_stnd:muddy_stnd) %>% st_set_geometry(NULL))

#standardize environmental prediction grid####
predict_grid <- env.grid %>% 
  mutate(depth_ln_stnd = scale_fun_ref(log(depth), reference = log(all_sets_wide$depth_m)),
         salinity_stnd = scale_fun_ref(salinity, reference = all_sets_wide$salinity),
         salinity_range_stnd = scale_fun_ref(salinity_range, reference = all_sets_wide$salinity_range),
         tidal_ln_stnd = scale_fun_ref(log(tidal), reference = log(all_sets_wide$tidal)),
         circ_stnd = scale_fun_ref(circulation, reference = all_sets_wide$circulation),
         BPI_stnd = scale_fun_ref(BBPI, reference = all_sets_wide$BBPI),
         rocky_stnd = scale_fun_ref(rocky, reference = all_sets_wide$rocky), 
         muddy_stnd = scale_fun_ref(muddy, reference = all_sets_wide$muddy))

scale_fun <- function(x){
  (x  - mean(x)) / sd(x)
}

spatialised_sf %>%
  st_drop_geometry() %>%
  mutate(rei_ln = log(rei),
         rei_sqrt = sqrt(rei),
         tidal_sqrt = sqrt(tidal)) %>%
  select(depth, rei, rei_ln, rei_sqrt, slope, tidal, tidal_sqrt, salinity, temperature, freshwater) %>%
  ggpairs()


spatialised_sf <- spatialised_sf %>%
  mutate(rei_ln = log(rei),
         rei_sqrt = sqrt(rei), 
         tidal_sqrt = sqrt(tidal)) %>%
  mutate(depth_stnd = scale(depth), 
         rei_stnd = scale(rei),
         rei_ln_stnd = scale(rei_ln),
         rei_sqrt_stnd = scale(rei_sqrt), 
         slope_stnd = scale(slope),
         temperature_stnd = scale(temperature), 
         tidal_sqrt_stnd = scale(tidal_sqrt), 
         salinity_stnd = scale(salinity))


#### Create spatial blocks for CV  (blockCV package functions)

# spatial clustering
sp_blocks <- cv_spatial(x = spatialised_sf,
                        column = NULL,
                        k = 5,
                        hexagon = FALSE,
                        size = 25000,  
                        selection = "random",
                        biomod2 = FALSE,
                        seed = 42, # to ensure reproducibility
                        plot = FALSE)

spatialised_sf$fold <- sp_blocks$folds_ids
spatialised_sf$fold[is.na(spatialised_sf$fold)] <- 5
table(spatialised_sf$fold)

cv_plot(cv = sp_blocks, #  blockCV object
        x = spatialised_sf, # sample points
        r = rei_all) #  raster background

sp_blocks <- ggplot(spatialised_sf)+
  geom_sf(aes(color = factor(fold)))+
  geom_sf(data = coastline)
sp_blocks
ggsave("./figures/pre-analysis/spatial_blocks.pdf", height = 6, width = 6)


cv <- list()
for ( f in 1:5 ){
  foldname <- paste0("fold",f)
  cv[[foldname]][['train']] <- which(sp_blocks$folds_ids!=f)
  cv[[foldname]][['test']] <- which(sp_blocks$folds_ids==f)
}


#return blockpoly, foldID and CV list
cv_list<-list( foldID=sp_blocks$foldID, blockpolys=sp_blocks$blocks, cv=cv) 

#prepare data for model####
seagrass_data <- spatialised_sf %>% st_drop_geometry()
XY <- as.data.frame(st_coordinates(spatialised_sf))/1000
seagrass_data$X <- XY$X
seagrass_data$Y <- XY$Y
seagrass_data$X_m <- XY$X*1000
seagrass_data$Y_m <- XY$Y*1000

seagrass_data_long <- seagrass_data %>%
  select(HKey, ID, fold, Year, depth_stnd, slope_stnd, substrate, rei_stnd, rei_ln_stnd, rei_sqrt_stnd,
         temperature_stnd, tidal_sqrt_stnd, salinity_stnd, ZO, PH, X:Y_m) %>%
  gather(key = species, value = presence, ZO:PH) %>%
  mutate(presence = ifelse(presence > 1, 1, presence)) %>%
  mutate(HKey = factor(HKey), substrate = factor(substrate))

#save outputs####
save(seagrass_data_long, seagrass_data, coastline, cv_list, file = "code/output_data/seagrass_model_inputs.RData")

# # Convert to spdf
# spatialised.seagrass <- seagrass_data %>%
#   st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 
# 
# # export as shapefile
# # likely to have issues with attribute field names shortening
# spatialised.sf <- spatialised.sf %>% mutate_all(~replace(., is.na(.), -9999))
# 
# st_write(spatialised.sf, "code/output_data/seagrass_covariates.shp", append=FALSE)
