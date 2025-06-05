###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# validate SDMs with independent datasets. NETforce for eelgrass and ShoreZone for surfgrass
#
###############################################################################

# raster of netforce
# compare overlap between netforce and different sdm thresholds. 

#load packages####
library(sf)
library(tidyverse)
library(terra)

#load netforce data
poly<- vect("raw_data/netforce/netforce_eelgrass_BC_polygons.shp")
poly_subset <- poly[poly$Year >= 2013 & poly$Year <= 2023, ]
poly_subset <- poly_subset[poly_subset$QC_Score >= 3, ]
lines1<- vect("raw_data/netforce/netforce_eelgrass_BC_line.shp")
names(lines1)[names(lines1) == "YEAR"] <- "Year"
lines1_subset <- lines1[lines1$Year >= 2013 & lines1$Year <= 2023, ]
lines1_subset <- lines1_subset[lines1_subset$QC_Score >= 3, ]
lines2<- vect("raw_data/netforce/netforce_eelgrass_BC_line_Transects.shp")
lines2_subset <- lines2[lines2$Year >= 2013 & lines2$Year <= 2023, ]
lines2_subset <- lines2_subset[lines2_subset$QC_Score >= 3, ]
points<- vect("raw_data/netforce/netforce_eelgrass_BC_points.shp")
names(points)[names(points) == "QC_score"] <- "QC_Score"
points_subset <- points[points$Year >= 2013 & points$Year <= 2023, ]
points_subset <- points_subset[points_subset$QC_Score >= 3, ]

#excluding lines 2 as all data is less than 2010
geoms <- list(points = points_subset)
lines = lines1_subset, poly = poly_subset) #, 
# Unique years across all types
all_years <- unique(unlist(lapply(geoms, \(g) unique(g$Year))))


# template raster
template_rast <- rast(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))

# Rasterize each geometry type for each year
yearly_rasters_points <- list()

for (yr in all_years) {
  yearly_stack <- list()
  
  for (type in names(geoms)) {
    g <- geoms[[type]]
    g_year <- g[g$Year == yr, ]
    if (nrow(g_year) > 0) {
      r <- rasterize(g_year, template_rast, field = 1, background = NA)
      r[!is.na(r)] <- 1
      yearly_stack[[type]] <- r
    }
  }
  
  # Merge the geometry-type rasters for this year (cell = 1 if touched by any)
  if (length(yearly_stack) == 1) {
    # Only one geometry type had data for this year
    r_combined <- yearly_stack[[1]]
  } else if (length(yearly_stack) > 1) {
    # Combine multiple rasters using cover()
    r_combined <- do.call(cover, yearly_stack)
  } else {
    # No geometry types had data for this year — skip to next year
    next
  }
  # Set touched cells to 1
  r_combined[!is.na(r_combined)] <- 1
  yearly_rasters_points[[as.character(yr)]] <- r_combined
}

# Sum across all years
r_stack <- rast(yearly_rasters)
r_years <- app(r_stack, sum, na.rm = TRUE)



yearly_rasters <- list()

for (yr in all_years) {
  yearly_stack <- list()
  
  for (type in names(geoms)) {
    g <- geoms[[type]]
    g_year <- g[g$year == yr, ]
    
    if (!is.null(g_year) && nrow(g_year) > 0) {
      r <- try(rasterize(g_year, template, field = 1, background = NA), silent = TRUE)
      if (!inherits(r, "try-error")) {
        r[!is.na(r)] <- 1
        yearly_stack[[type]] <- r
      }
    }
  }
  
  # Defensive logic
  if (length(yearly_stack) == 1) {
    r_combined <- yearly_stack[[1]]
  } else if (length(yearly_stack) > 1) {
    # Remove any NULLs (shouldn’t happen, but just in case)
    non_null_stack <- yearly_stack[!sapply(yearly_stack, is.null)]
    if (length(non_null_stack) >= 2) {
      r_combined <- do.call(cover, non_null_stack)
    } else {
      r_combined <- non_null_stack[[1]]
    }
  } else {
    next  # Skip this year; no data
  }
  
  # Set cell values to 1
  r_combined[!is.na(r_combined)] <- 1
  yearly_rasters[[as.character(yr)]] <- r_combined
}
