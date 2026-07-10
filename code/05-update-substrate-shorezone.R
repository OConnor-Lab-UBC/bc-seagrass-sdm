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
# it is known that the substrate model over predicts rock and does not capture exposed sand or pocket beaches well
# Shorezone was used to replace the intertidal and 0 to 5 m depths in the susbtrate model in areas visually observed with shorezone to have sand
# In ArcPro shorezone coastal classes 27 (wide sand beach), 28 (sand flat), 30 (narrow sand beach),
# 16 (ramp with wide sand beach), and 17 (platform with wide sand beach) unit lines were used to select bottom patches depth bands intertidal and 0-5 ribbons. 
# Bottom patches polygons were only included that were within 100 m from the shorezone unit lines. Historical shorezone uses an outdated high water line that did not match the bottom patches so manual removal of not needed polygons was necessary
# classes 16 and 17 required a bit more consideration and some rock will be included in the substrate. Only polygons were retained that visually (from the satellite imagery) had over 50% sand versus rock
# we also included the class 29 (mud flat) to replace the susbtrate layer for mud, in case mud is not predicted well in some areas
###############################################################################

# load packages
library(sf)
library(tidyverse)
library(terra)

#load substrate layer
substrate_hg <- rast("raw_data/substrate_20m/hg_20m.tif")
substrate_ncc <- rast("raw_data/substrate_20m/ncc_20m.tif")
substrate_qcs <- rast("raw_data/substrate_20m/qcs_20m.tif")
substrate_wcvi <- rast("raw_data/substrate_20m/wcvi_20m.tif")
substrate_ss <- rast("raw_data/substrate_20m/sog_20m.tif")
names(substrate_hg) <- names(substrate_ncc) <- names(substrate_qcs) <- names(substrate_wcvi) <- names(substrate_ss) <- "substrate"
crs(substrate_hg) <- crs(substrate_ncc) <- crs(substrate_qcs) <- crs(substrate_wcvi) <- crs(substrate_ss) <- "EPSG:3005"

substrate<-mosaic(substrate_hg, substrate_ncc, substrate_qcs, substrate_wcvi, substrate_ss, fun = "max")
names(substrate) <- "substrate"
crs(substrate) <- "EPSG:3005"

sand_vector_27_28_30 <- vect("raw_data/substrate_20m/shorezone_sand_class27_28_30_replace_intertidal_to_5m.shp")
sand_vector_16_17<- vect("raw_data/substrate_20m/shorezone_sand_class16_17_replace_intertidal_to_5m.shp")
mud_vector_29<- vect("raw_data/substrate_20m/shorezone_mud_class29_replace_intertidal_to_5m.shp")
sand_raster_27_28_30 <- rasterize(sand_vector_27_28_30, substrate, fun='count')
sand_raster_27_28_30 <- ifel(!is.na(sand_raster_27_28_30), 3, NA)
sand_raster_16_17 <- rasterize(sand_vector_16_17, substrate, fun='count')
sand_raster_16_17 <- ifel(!is.na(sand_raster_16_17), 3, NA)
mud_raster <- rasterize(mud_vector_29, substrate, fun='count')
mud_raster <- ifel(!is.na(mud_raster), 4, NA)

substrate_sandupdate <- ifel(!is.na(sand_raster_27_28_30), sand_raster_27_28_30, substrate)
substrate_sandupdate2 <- ifel(!is.na(sand_raster_16_17), sand_raster_16_17, substrate_sandupdate)
substrate_update <- ifel(!is.na(mud_raster), mud_raster, substrate_sandupdate2)

substrate_update_hg<- terra::crop(substrate_update, substrate_hg)
substrate_update_hg <- resample(substrate_update_hg, substrate_hg, method = "near")

substrate_update_ncc<- terra::crop(substrate_update, substrate_ncc)
substrate_update_ncc <- resample(substrate_update_ncc, substrate_ncc, method = "near")

substrate_update_qcs<- terra::crop(substrate_update, substrate_qcs)
substrate_update_qcs <- resample(substrate_update_qcs, substrate_qcs, method = "near")

substrate_update_wcvi<- terra::crop(substrate_update, substrate_wcvi)
substrate_update_wcvi <- resample(substrate_update_wcvi, substrate_wcvi, method = "near")

substrate_update_ss<- terra::crop(substrate_update, substrate_ss)
substrate_update_ss <- resample(substrate_update_ss, substrate_ss, method = "near")

writeRaster(substrate_update_hg, file.path("raw_data/substrate_20m/updated/hg_20m.tif"), overwrite=TRUE)
writeRaster(substrate_update_ncc, file.path("raw_data/substrate_20m/updated/ncc_20m.tif"), overwrite=TRUE)
writeRaster(substrate_update_qcs, file.path("raw_data/substrate_20m/updated/qcs_20m.tif"), overwrite=TRUE)
writeRaster(substrate_update_wcvi, file.path("raw_data/substrate_20m/updated/wcvi_20m.tif"), overwrite=TRUE)
writeRaster(substrate_update_ss, file.path("raw_data/substrate_20m/updated/sog_20m.tif"), overwrite=TRUE)





