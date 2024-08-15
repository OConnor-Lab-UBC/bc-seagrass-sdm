
library(dplyr)
library(terra)
library(sf)
#library(sp)
library(rnaturalearth)
library(ggplot2)
library(gstat)
library(MultiscaleDTM)
#install.packages("devtools")
#devtools::install_github("bio-oracle/biooracler")
library(biooracler)
library(stars)

#Create empty raster
BaseRast <- rast(extent=ext(474232.8235, 1304032.8235, 314708.7366, 1250588.7366), crs = "EPSG:3005", resolution = 200)


dataset_id_baseline <- BioOracle_vars[grep("baseline", BioOracle_vars$dataset_id), ]

#baseline constraints
time = c('2001-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
latitude = c(48, 56)
longitude = c(-134,-122)
# latitude = c(10, 20)
# longitude = c(120, 130)


constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

dir ="predictor_layers/baseline_2000_2020/"


#thetao (Ocean temperature - Bottom)
variables = c("thetao_max", "thetao_min")

data_id <- dataset_id_baseline[grep("thetao", dataset_id_baseline$dataset_id),]


baseline <- download_layers("thetao_baseline_2000_2019_depthmean", variables, constraints, fmt = "raster")

baseline <- subst(baseline, NA, -9999)
baseline <- terra::project(baseline,BaseRast,method="near")
baseline <- subst(baseline, -9999, NA)

writeRaster(baseline$thetao_max_1, file="code/output_data/baseline_2000_2020/bt_max_2000_baseline.tif",overwrite=TRUE)

writeRaster(baseline$thetao_max_2, file="predictor_layers/baseline_2000_2020/bt_max_2010_baseline.tif",overwrite=TRUE)

writeRaster(baseline$thetao_min_1, file="predictor_layers/baseline_2000_2020/bt_min_2000_baseline.tif",overwrite=TRUE)

writeRaster(baseline$thetao_min_2, file="predictor_layers/baseline_2000_2020/bt_min_2010_baseline.tif",overwrite=TRUE)

