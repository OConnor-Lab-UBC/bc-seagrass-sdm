###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# Generate independent datasets to validate SDMs. NETforce for eelgrass and ShoreZone for surfgrass
#
###############################################################################

# raster of netforce eelgrass polygon, point and line data: want full netforce data, years prior to sdm obs (1974-1992), years from obs to fit sdm (1993-2023)
# raster of shorezone surfgrass, make assumption that surfgrass is found near high water line that is used by shorezone. Could make polygons that extend down to 0 (ask Kayleigh she did with bottom patches)

# what about GBIF, Inat? reeflife
# Matt's ROV surfgrass points

#load packages####
library(sf)
library(tidyverse)
library(terra)

#load netforce data
poly<- vect("raw_data/netforce/netforce_eelgrass_BC_polygons.shp")
#only use data that have the spatial accuracy we need
poly <- poly[poly$QC_Score >= 3, ]
#poly_subset_post2012 <- poly[poly$Year >= 2013 & poly$Year <= 2023, ]

lines1<- vect("raw_data/netforce/netforce_eelgrass_BC_line.shp")
names(lines1)[names(lines1) == "YEAR"] <- "Year"
lines1<- lines1[lines1$QC_Score >= 3, ]
#lines1_subset_post2012 <- lines1[lines1$Year >= 2013 & lines1$Year <= 2023, ]

lines2<- vect("raw_data/netforce/netforce_eelgrass_BC_line_Transects.shp")
lines2 <- lines2[lines2$QC_Score >= 3, ]
#lines2_subset_post2012 <- lines2[lines2$Year >= 2013 & lines2$Year <= 2023, ]

lines<-rbind(lines1, lines2)

points<- vect("raw_data/netforce/netforce_eelgrass_BC_points.shp")
names(points)[names(points) == "QC_score"] <- "QC_Score"
points <- points[points$QC_Score >= 3, ]
#points_subset_post2012 <- points[points$Year >= 2013 & points$Year <= 2023, ]

# template raster
template_rast <- rast(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))

# Identify all unique years
all_years <- sort(unique(c(poly$Year, lines1$Year, lines2$Year, points$Year)))
all_years <- all_years[all_years != "-999"] 

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
    template_rast = template_rast)
  # Store in list
  eelgrass_rasters[[as.character(yr)]] <- r
}

save(eelgrass_rasters, file = "code/output_data/eelgrass_netforce_rasters.RData")



# Subset the list to only the years of interest, Combine the yearly rasters into a SpatRaster stack
netforce_start_year <- 1974
presurvey_start_year <- 1992
survey_start_year <- 1993
survey_end_year <- 2023
prediction_start_year <- 2013
preprediction_start_year <-2012
presurvey_years_1974_1992 <- as.character(seq(netforce_start_year, presurvey_start_year))
survey_years_1993_2012 <- as.character(seq(survey_start_year, preprediction_start_year))
predicted_years_2013_2023 <- as.character(seq(prediction_start_year, survey_end_year))

# Sum across years (ignoring NA, only counting 1s)
# Set cells that were NA in all years back to NA
# Save rasters

# selected_stack_1974_2023 <- rast(eelgrass_rasters)
# eelgrass_years_observed_1974_2023 <- sum(selected_stack_1974_2023, na.rm = TRUE)
# eelgrass_years_observed_1974_2023[sum(is.na(selected_stack_1974_2023)) == nlyr(selected_stack_1974_2023)] <- NA
# writeRaster(eelgrass_years_observed_1974_2023, "netforce_1974_2023.tif", overwrite = TRUE)
# remove(selected_stack_1974_2023, eelgrass_years_observed_1974_2023)

selected_rasters_1974_1992 <- eelgrass_rasters[names(eelgrass_rasters) %in% as.character(presurvey_years_1974_1992)]
selected_stack_1974_1992 <- rast(selected_rasters_1974_1992)
eelgrass_years_observed_1974_1992 <- sum(selected_stack_1974_1992, na.rm = TRUE)
eelgrass_years_observed_1974_1992[sum(is.na(selected_stack_1974_1992)) == nlyr(selected_stack_1974_1992)] <- NA
writeRaster(eelgrass_years_observed_1974_1992, "netforce_1974_1992.tif", overwrite = TRUE)
remove(selected_rasters_1974_1992, selected_stack_1974_1992, eelgrass_years_observed_1974_1992)

selected_rasters_1993_2012 <- eelgrass_rasters[names(eelgrass_rasters) %in% as.character(survey_years_1993_2012)]
selected_stack_1993_2012 <- rast(selected_rasters_1993_2012)
eelgrass_years_observed_1993_2012 <- sum(selected_stack_1993_2012, na.rm = TRUE)
eelgrass_years_observed_1993_2012[sum(is.na(selected_stack_1993_2012)) == nlyr(selected_stack_1993_2012)] <- NA
writeRaster(eelgrass_years_observed_1993_2012, "netforce_1993_2012.tif", overwrite = TRUE)
remove(selected_rasters_1993_2012, selected_stack_1993_2012, eelgrass_years_observed_1993_2012)

selected_rasters_2013_2023 <- eelgrass_rasters[names(eelgrass_rasters) %in% as.character(predicted_years_2013_2023)]
selected_stack_2013_2023 <- rast(selected_rasters_2013_2023)
eelgrass_years_observed_2013_2023 <- sum(selected_stack_2013_2023, na.rm = TRUE)
eelgrass_years_observed_2013_2023[sum(is.na(selected_stack_2013_2023)) == nlyr(selected_stack_2013_2023)] <- NA
writeRaster(eelgrass_years_observed_2013_2023, "netforce_2013_2023.tif", overwrite = TRUE)
remove(selected_rasters_2013_2023, selected_stack_2013_2023, eelgrass_years_observed_2013_2023)

