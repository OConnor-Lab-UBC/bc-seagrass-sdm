###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# Generate independent datasets to validate SDMs. Primaraly NETforce for eelgrass and ShoreZone for surfgrass
#
###############################################################################
# Eelgrass
# raster of netforce eelgrass polygon, point and line data: want full netforce data, years prior to sdm obs (1974-1992), years from obs to fit sdm (1993-2023)

# Surfgrass
# raster of shorezone surfgrass, make assumption that surfgrass is found near high water line that is used by shorezone. 
# shorezone surfgrass data: selected lines from 2014-2024 (could exclude cammano if jsut want up to 2023) and only kept records with along line occurence >25% (SURF_L). = 310 records. Buffered line by 40 m to capture where surfgrass likely
# still underrepresented where surfgrass is likely to be on shoreline so selected bottom patches depth ribbons within 100 m from shoreline unit. Manually went through and selected appropriate polygons. Combined buffered unit length and BoPs for final surfgrass validation layer. 
# went through gbif records and most were i naturalist so went back to origional photos to verify they were not drift
# also used Matt's ROV points. Used the lat long that we closest to shroe as this is likely where the surfgrass was observed

#load packages####
library(sf)
library(tidyverse)
library(terra)

#### Load modelling functions ####
source("code/modelling-functions.R")

#load netforce eelgrass  data
poly<- vect("raw_data/netforce/netforce_eelgrass_BC_polygons.shp")
#only use data that have the spatial accuracy we need
poly <- poly[poly$QC_Score >= 3, ]

lines1<- vect("raw_data/netforce/netforce_eelgrass_BC_line.shp")
names(lines1)[names(lines1) == "YEAR"] <- "Year"
lines1<- lines1[lines1$QC_Score >= 3, ]

lines2<- vect("raw_data/netforce/netforce_eelgrass_BC_line_Transects.shp")
lines2 <- lines2[lines2$QC_Score >= 3, ]
lines<-rbind(lines1, lines2)

points<- vect("raw_data/netforce/netforce_eelgrass_BC_points.shp")
names(points)[names(points) == "QC_score"] <- "QC_Score"
points <- points[points$QC_Score >= 3, ]

# template raster
template_rast <- rast(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))

# Identify all unique years
all_years <- sort(unique(c(poly$Year, lines1$Year, lines2$Year, points$Year)))
all_years <- all_years[all_years != "-999"]
# r was crashing so had to write each year directly to disk and start over a few times
# years_to_remove <- c("1974", "1975", "1976", "1977", "1978", "1979", "1980", 
#                      "1981", "1995", "1996", "1997", "1999", "2000", "2002", 
#                      "2003", "2004", "2005", "2006", "2007", "2008", "2009",
#                      "2010", "2011", "2012", "2013", "2014", "2015", "2016", 
#                      "2017", "2018", "2019", "2020", "2021", "2022", "2023",
#                      "2024")
# 
# # Filter all_years to exclude the unwanted years
# all_years <- all_years[!all_years %in% years_to_remove]

# Initialize a list to store rasters
eelgrass_rasters <- list()

# Loop through each year and rasterize
for (yr in all_years) {
  cat("Processing year:", yr, "\n")
  r <- rasterize_eelgrass_year(
    year = yr,
    poly_data = poly,
    line_data = lines,
    point_data = points,
    template_rast = template_rast
    )
  # Store in list
  eelgrass_rasters[[as.character(yr)]] <- r
  
  # Save raster to file 
  writeRaster(r, filename = paste0("code/output_data/independent_validation/eelgrass_netforce_raster_", yr, ".tif"), 
              overwrite = TRUE)
}



# Subset the list to only the years of interest
prediction_years <- 2013:2023
pre_prediction_years <-1974:2012

# Create full file paths for prediction years
prediction_year_rasters <- paste0("code/output_data/independent_validation/eelgrass_netforce_raster_", prediction_years, ".tif")

# Load first raster
r <- rast(prediction_year_rasters[1])
presence_count <- ifel(!is.na(r) & r > 0, 1, 0)  # Count 1 if presence, 0 if absence or NA

# Loop through remaining rasters
for (i in 2:length(prediction_year_rasters)) {
  r <- rast(prediction_year_rasters[i])
  presence_count <- presence_count + ifel(!is.na(r) & r > 0, 1, 0)
  gc()
}

# Set cells to NA if they were NA in *all* years
presence_count[presence_count == 0] <- NA

writeRaster(presence_count, "code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif", overwrite = TRUE)

# Create full file paths for pre-prediction years
preprediction_year_rasters <- paste0("code/output_data/independent_validation/eelgrass_netforce_raster_", pre_prediction_years, ".tif")

# Load first raster
r <- rast(preprediction_year_rasters[1])
presence_count <- ifel(!is.na(r) & r > 0, 1, 0)  # Count 1 if presence, 0 if absence 

# Loop through remaining rasters
for (i in 2:length(preprediction_year_rasters)) {
  r <- rast(preprediction_year_rasters[i])
  presence_count <- presence_count + ifel(!is.na(r) & r > 0, 1, 0)
  gc()
}

# Set cells to NA if they were NA in *all* years
presence_count[presence_count == 0] <- NA

writeRaster(presence_count, "code/output_data/independent_validation/BCeelgrass_netforce_1974_2012.tif", overwrite = TRUE)



#### phyllospadix, there is so little data with overlap will not be separating by years and doing a count, just a yes or no

#load shorezone surfgrass  data
poly_sg<- vect("raw_data/shorezone/phyllospadix2014_2024_shorezone_BOP_depthribbons.shp")

# load Matt's rov data 2020-2023 (though removed most 2020 data as there was only one deep lat and long and not the shore one)
points_sg1<- vect("raw_data/Baum/ROV_Phyllospadix_obs.shp")
#load gbif data
points_sg2<- vect("raw_data/gbif/phyllospadix_gbif.shp")

# template raster
template_rast <- rast(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))

# Rasterize all with value = 1
r_poly  <- rasterize(poly_sg, template_rast, field = 1, background = NA, touches = TRUE)
r_pts1  <- rasterize(points_sg1, template_rast, field = 1, background = NA, touches = TRUE)
r_pts2  <- rasterize(points_sg2, template_rast, field = 1, background = NA, touches = TRUE)

# Combine them using cover()
surfgrass_combined <- cover(r_pts1, r_poly)
surfgrass_combined <- cover(r_pts2, surfgrass_combined)

writeRaster(surfgrass_combined, "code/output_data/independent_validation/surfgrass_validation_raster_2013_2024.tif", overwrite = TRUE)
