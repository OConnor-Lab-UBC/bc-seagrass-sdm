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

#load_data####
load("code/output_data/seagrass_data_spatialized.RData")
# 400,597 quadrats

# check covariate observations
summary(spatialised$Slope) #18,129 NAs
summary(spatialised$CorDepthM) #584 NAs
sum(is.na(spatialised$Substrate)) #1,941 NAs

spatialised <- spatialised %>%
  filter(!is.na(CorDepthM))




#read in 20m rasters####
rei_all <- terra::vrt(c("raw_data/REI/rei_hg.tif", "raw_data/REI/rei_ncc.tif", "raw_data/REI/rei_qcs.tif", "raw_data/REI/rei_sog.tif", "raw_data/REI/rei_wcvi.tif"), "rei.vrt", overwrite=T)
#plot(rei_all)

slope_all <- terra::vrt(c("raw_data/envlayers-20m-hg/slope.tif", "raw_data/envlayers-20m-ncc/slope.tif", "raw_data/envlayers-20m-qcs/slope.tif", "raw_data/envlayers-20m-wcvi/slope.tif", "raw_data/envlayers-20m-shelfsalishsea/slope.tif"), "slope.vrt", overwrite=T)
#plot(slope_all)

substrate_all <- terra::vrt(c("raw_data/substrate_20m/hg_20m.tif", "raw_data/substrate_20m/ncc_20m.tif", "raw_data/substrate_20m/qcs_20m.tif", "raw_data/substrate_20m/sog_20m.tif", "raw_data/substrate_20m/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"
#plot(substrate_all)

temperature_all <- vrt(c("raw_data/envlayers-20m-hg/temp_mean_summer.tif", "raw_data/envlayers-20m-ncc/temp_mean_summer.tif", "raw_data/envlayers-20m-qcs/temp_mean_summer.tif", "raw_data/envlayers-20m-wcvi/temp_mean_summer.tif", "raw_data/envlayers-20m-shelfsalishsea/temp_mean_summer.tif"), "temp_mean_summer.vrt", overwrite=T)
#plot(temperature_all)

salinity_all <- vrt(c("raw_data/envlayers-20m-hg/salinity_mean_summer.tif", "raw_data/envlayers-20m-ncc/salinity_mean_summer.tif", "raw_data/envlayers-20m-qcs/salinity_mean_summer.tif", "raw_data/envlayers-20m-wcvi/salinity_mean_summer.tif", "raw_data/envlayers-20m-shelfsalishsea/salinity_mean_summer.tif"), "salinity_mean_summer.vrt", overwrite=T)
#plot(salinity_all)

freshwater_all <- vrt(c("raw_data/freshwater-index/hg_freshwater_index.tif", "raw_data/freshwater-index/ncc_freshwater_index.tif", "raw_data/freshwater-index/qcs_freshwater_index.tif", "raw_data/freshwater-index/wcvi_freshwater_index.tif", "raw_data/freshwater-index/salish_sea_freshwater_index.tif"), "freshwater.vrt", overwrite=T)
#plot(freshwater_all)

#tidal_all <- vrt(c("raw_data/envlayers-20m-hg/tidal_mean_summer.tif", "raw_data/envlayers-20m-ncc/tidal_mean_summer.tif", "raw_data/envlayers-20m-qcs/tidal_mean_summer.tif", "raw_data/envlayers-20m-wcvi/tidal_mean_summer.tif", "raw_data/envlayers-20m-shelfsalishsea/tidal_mean_summer.tif"), "tidal_mean_summer.vrt", overwrite=T)
tidal_all <- vrt(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))
names(tidal_all)<-"tidal_mean_summer"
#change to index 0-1 scale
tidal_index_all <- tidal_all/(maxFn(tidal_all))
crs(tidal_index_all) <- "EPSG:3005"
#plot(tidal_index_all)

coastline_full <- st_read("raw_data/CHS_HWL2015_Coastline.gdb", layer = "Line_CHS_Pacific_HWL_2015_5028437")
coastline <- coastline_full %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_crop(st_bbox(coastline_full) + c(0, 100000, -10000, -80000)) %>%
  st_transform(crs = "EPSG:3005")

spatialised.sf<- spatialised %>%
  st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

ggplot(spatialised.sf)+
  geom_sf(data = coastline, size = 0.001)+
  geom_sf(color = "red")+
  theme_bw()+coord_sf(expand = FALSE)

#add extracted data to seagrass data####

rei_extract <- terra::extract(rei_all, spatialised.sf)
summary(rei_extract$rei)
spatialised.sf$rei <- rei_extract$rei

temperature_extract <- terra::extract(temperature_all, spatialised.sf)
summary(temperature_extract$temp_mean_summer)
spatialised.sf$temperature <- temperature_extract$temp_mean_summer

salinity_extract <- terra::extract(salinity_all, spatialised.sf)
summary(salinity_extract$salinity_mean_summer)
spatialised.sf$salinity <- salinity_extract$salinity_mean_summer

tidal_extract <- terra::extract(tidal_index_all, spatialised.sf)
summary(tidal_extract$tidal_mean_summer)
spatialised.sf$tidal <- tidal_extract$tidal_mean_summer

freshwater_extract <- terra::extract(freshwater_all, spatialised.sf)
summary(freshwater_extract$freshwater)
spatialised.sf$freshwater <- freshwater_extract$freshwater

slope_extract <- terra::extract(slope_all, spatialised.sf)
summary(slope_extract$slope)
spatialised.sf$slope_mod <- slope_extract$slope

substrate_extract <- terra::extract(substrate_all, spatialised.sf)
summary(substrate_extract$substrate)
spatialised.sf$substrate_mod <- substrate_extract$substrate
spatialised.sf$substrate_mod <- c("Rock", "Mixed", "Sand", "Mud")[spatialised.sf$substrate_mod]

spatialised.sf <- spatialised.sf %>%
  mutate(slope = ifelse(is.na(Slope)|Slope > 90, slope_mod, Slope),
         substrate = ifelse(is.na(Substrate), substrate_mod, Substrate)) %>%
  rename(depth = CorDepthM) %>%
  select(-c(Slope, Substrate, slope_mod, substrate_mod))


spatialised.sf$slope<- round(spatialised.sf$slope, digits = 0)

spatialised.sf <- spatialised.sf %>%
  filter(!is.na(rei), !is.na(substrate), !is.na(temperature), !is.na(tidal))

#scale environmental data####

scale_fun <- function(x){
  (x  - mean(x)) / sd(x)
}

# spatialised.sf %>%
#   st_drop_geometry() %>%
#   mutate(rei_ln = log(rei),
#          rei_sqrt = sqrt(rei),
#          tidal_sqrt = sqrt(tidal)) %>%
#   select(depth, rei, rei_ln, rei_sqrt, slope, tidal, tidal_sqrt, salinity, temperature, freshwater) %>%
#   ggpairs()


spatialised.sf <- spatialised.sf %>%
  mutate(rei_sqrt = sqrt(rei), 
         tidal_sqrt = sqrt(tidal)) %>%
  mutate(depth_stnd = scale(depth), 
         rei_stnd = scale(rei),
         rei_sqrt_stnd = scale(rei_sqrt), 
         slope_stnd = scale(slope),
         temperature_stnd = scale(temperature), 
         tidal_sqrt_stnd = scale(tidal_sqrt), 
         salinity_stnd = scale(salinity))


#### Create spatial blocks for CV  (blockCV package functions)

# spatial clustering
sp.blocks <- cv_spatial(x = spatialised.sf,
                        column = NULL,
                        k = 5,
                        hexagon = FALSE,
                        size = 25000,  
                        selection = "random",
                        biomod2 = FALSE,
                        seed = 42, # to ensure reproducibility
                        plot = FALSE)

spatialised.sf$fold <- sp.blocks$folds_ids
spatialised.sf$fold[is.na(spatialised.sf$fold)] <- 5
table(spatialised.sf$fold)

cv_plot(cv = sp.blocks, #  blockCV object
        x = spatialised.sf, # sample points
        r = rei_all) #  raster background

sp_blocks <- ggplot(spatialised.sf)+
  geom_sf(aes(color = factor(fold)))+
  geom_sf(data = coastline)
sp_blocks
ggsave("./figures/pre-analysis/spatial_blocks.pdf", height = 6, width = 6)


cv <- list()
for ( f in 1:5 ){
  foldname <- paste0("fold",f)
  cv[[foldname]][['train']] <- which(sp.blocks$folds_ids!=f)
  cv[[foldname]][['test']] <- which(sp.blocks$folds_ids==f)
}


#return blockpoly, foldID and CV list
cv.list<-list( foldID=sp.blocks$foldID, blockpolys=sp.blocks$blocks, cv=cv) 

#prepare data for model####
seagrass_data <- spatialised.sf %>% st_drop_geometry()
XY <- as.data.frame(st_coordinates(spatialised.sf))/1000
seagrass_data$X <- XY$X
seagrass_data$Y <- XY$Y
seagrass_data$X_m <- XY$X*1000
seagrass_data$Y_m <- XY$Y*1000

seagrass_data_long <- seagrass_data %>%
  select(HKey, ID, fold, Year, depth_stnd, slope_stnd, substrate, rei_stnd, rei_sqrt_stnd,
         temperature_stnd, tidal_sqrt_stnd, salinity_stnd, ZO, PH, X:Y_m) %>%
  gather(key = species, value = presence, ZO:PH) %>%
  mutate(presence = ifelse(presence > 1, 1, presence)) %>%
  mutate(HKey = factor(HKey), substrate = factor(substrate))

#save outputs####
save(seagrass_data_long, seagrass_data, coastline, cv.list, file = "code/output_data/seagrass_model_inputs.RData")

