###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# extract field validation data points 
#
###############################################################################


#### load packages####
library(DBI)
library(odbc)
library(reshape2)
library(sf)
library(tidyverse)

#### Load util.R gfdata functions ####
source("code/r-sql-link-functions.R")

#### Get field validation dive survey data ####

# Connect to mdb
sdm_mdb <- mdb_connection("../../../Field validation/Database/SDM_Field_Validation2025UpdatedLocations.accdb")

# Load queries
sdm_sql <- readLines("code/sql/get-sdm-records.sql")
sdm_sql <- paste(sdm_sql, collapse = "\n ")


# Run queries
sdm <- DBI::dbGetQuery( sdm_mdb, sdm_sql )

# Convert to sf object (specify coordinate columns and CRS)
sdmpoints_sf <- st_as_sf(sdm, coords = c("LonDeep", "LatDeep"), crs = 4326)  # WGS84

# Check
print(sdmpoints_sf)

# Write to shapefile
st_write(sdmpoints_sf, "code/output_data/field_validation/2025divesites.shp", delete_dsn = TRUE)
