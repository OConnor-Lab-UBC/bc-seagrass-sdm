###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# Get independent data from GBIF
#
###############################################################################

# download occurrence data from GBIF using the 'rgbif' package. The data is then cleaned up. The script has been modified but is based on the tutorial at this link: https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/

# Libraries -----------------------------------------
library(rgbif)              #for downloading data from GBIF
library(CoordinateCleaner)  #helps clean coordinate data
library(maps)               #for creating maps
library(dplyr)              #used for data frame management, e.g. filter(), select(), group_by()
library(readr)              #to save df as a csv and read csv files in

# Downloading data for zostera and phyllospadix in BC
gbif_data <- occ_download(
  pred_in("taxonKey", c(2863967, 2863952, 2863958, 2863954, 2863957)), 
  pred_within("POLYGON((-133.0 48.0, -133.0 56.0, -122.0 56.0, -122.0 48.0,  -133.0 48.0))"),
  pred("hasCoordinate", TRUE))
    
d <- occ_download_get('0025972-250525065834625') %>%
  occ_download_import()

#filter to get more accurate data in canada of presence only data
d_filtered <- d %>%
  filter(occurrenceStatus == "PRESENT",
         countryCode == "CA",
         hasGeospatialIssues == "FALSE",
         coordinateUncertaintyInMeters < 21)

# cleaning up coordinates using the package CoordinateCleaner
# Identify Invalid Lat/Long Coordinates using cc_val
myspecies_clean <- cc_val(d_filtered,
                          lon = "decimalLongitude",
                          lat = "decimalLatitude",
                          value = "clean",
                          verbose = TRUE)

# map the cleaned occurrence data:
map("world", xlim = range(myspecies_clean$decimalLongitude), ylim = range(myspecies_clean$decimalLatitude))  # if the map doesn't appear right at first, run this command again
points(myspecies_clean[ , c("decimalLongitude", "decimalLatitude")], col = "red", pch = ".")


# map the cleaned occurrence records with a different colour on top of the raw ones:
points(myspecies_clean[ , c("decimalLongitude", "decimalLatitude")], pch = 20, cex = 0.5, col = "turquoise")

myspecies_clean$lon <- as.double(myspecies_clean$decimalLongitude) ## cast lon from char to double
myspecies_clean$lat <- as.double(myspecies_clean$decimalLatitude) 

data_sf <- sf::st_as_sf(myspecies_clean, coords = c("decimalLongitude", "decimalLatitude"), remove = FALSE, crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3005))

sf::st_write(data_sf, "raw_data/gbif/gbifextract.shp", overwrite = TRUE)
