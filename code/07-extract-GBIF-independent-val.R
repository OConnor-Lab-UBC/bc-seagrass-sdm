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

# Downloading data 
myspecies <- c("Zostera marina L.", "Phyllospadix Hook.", "Phyllospadix scouleri Hook.", "Phyllospadix torreyi S.Watson", "Phyllospadix serrulatus Rupr. ex Asch.")

gbif_data <- occ_download(scientificName = myspecies, hasCoordinate = TRUE, limit = 100000, decimalLongitude = "-133, -122", decimalLatitude = "48, 56")  
gbif_data

# get the DOI for citing these data properly:
gbif_citation(gbif_data)  # unfortunately it is more complicated to obtain with R a proper citation for a dataset with multiple species. To get a DOI for these data, download the dataset directly from www.gbif.org and then import the .csv to R. It is very important to properly cite the data sources! GBIF is not a source, just a repository for many people who put in very hard work to collect these data and make them available

# check how the data are organized:
names(gbif_data)
names(gbif_data[[myspecies[1]]])
names(gbif_data[[myspecies[1]]]$meta)
names(gbif_data[[myspecies[1]]]$data)

# create and fill a list with only the 'data' section for each species:
# removed: "depth", "depthAccuracy", "lifestage",  "sampleSizeUnit","samplingEffort", "coordinatePrecision", "fieldNotes", "samplingProtocol" because they were causing errors in R which was saying that they didn't exist.
myspecies_coords_list <- vector("list", length(myspecies))
names(myspecies_coords_list) <- myspecies
for (s in myspecies) {
  coords <- gbif_data[[s]]$data[ , c("country", "stateProvince", "decimalLatitude", "decimalLongitude", "issues", "protocol", "occurrenceStatus", "genus", "species", "acceptedScientificName", "coordinateUncertaintyInMeters","dateIdentified", "year", "month", "day", "datasetName", "recordedBy", "geodeticDatum", "habitat", "individualCount",  "elevation", "elevationAccuracy", "higherGeography", "georeferenceProtocol")]
  myspecies_coords_list[[s]] <- data.frame(species = s, coords)
}
lapply(myspecies_coords_list, head)

# collapse the list into a data frame:
myspecies_coords <- as.data.frame(do.call(rbind, myspecies_coords_list), row.names = 1:sum(sapply(myspecies_coords_list, nrow)))
head(myspecies_coords)
tail(myspecies_coords)

# filter country column to Canada only
myspecies_coords <- myspecies_coords %>% filter(country == "Canada")

# if you want to remove records of absence or zero-abundance (if any):
# I did not run this part of the script because I wanted to keep records of absence
names(myspecies_coords)
sort(unique(myspecies_coords$individualCount))  # notice if some points correspond to zero abundance
sort(unique(myspecies_coords$occurrenceStatus))  # check for different indications of "absent", which could be in different languages! and remember that R is case-sensitive
absence_rows <- which(myspecies_coords$individualCount == 0 | myspecies_coords$occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"))
length(absence_rows)
if (length(absence_rows) > 0) {
  myspecies_coords <- myspecies_coords[-absence_rows, ]
}

# cleaning up coordinates using the package CoordinateCleaner

# Identify Invalid Lat/Long Coordinates using cc_val
myspecies_clean <- cc_val(myspecies_coords,
                          lon = "decimalLongitude",
                          lat = "decimalLatitude",
                          value = "clean",
                          verbose = TRUE)

# Geographic Cleaning of Coordinates from Biological Collections using clean_coordinates
flagged <- clean_coordinates(myspecies_clean, seas_buffer = 5000, species = "scientificName")
summary(flagged)

# export data frame as a csv
write_csv(flagged, "C:\\Users\\laush\\Documents\\UVic\\Baum Lab\\Blue Carbon Canada\\salt_marsh_atlantic_arctic\\GBIF_clean_withTestResults.csv")

# map the cleaned occurrence data:
map("world", xlim = range(myspecies_clean$decimalLongitude), ylim = range(myspecies_clean$decimalLatitude))  # if the map doesn't appear right at first, run this command again
points(myspecies_clean[ , c("decimalLongitude", "decimalLatitude")], col = "red", pch = ".")
# possible erroneous points e.g. on the Equator (lat and lon = 0) should have disappeared now

#this function doesn't work
# also eliminate presences with reported coordinate uncertainty (location error, spatial resolution) larger than 5 km (5000 m):
myspecies_clean <- coord_uncertain(myspecies_clean, coorduncertainityLimit = 5000)
nrow(myspecies_clean)
# but note that this will only get rid of records where coordinate uncertainty is adequately reported, which may not always be the case! Careful mapping and visual inspection is necessary

# map the cleaned occurrence records with a different colour on top of the raw ones:
points(myspecies_clean[ , c("decimalLongitude", "decimalLatitude")], pch = 20, cex = 0.5, col = "turquoise")

myspecies_clean$lon <- as.double(myspecies_clean$decimalLongitude) ## cast lon from char to double
myspecies_clean$lat <- as.double(myspecies_clean$decimalLatitude) 

data_sf <- sf::st_as_sf(myspecies_clean, coords = c("lon", "lat"), remove = FALSE, crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3005))

sf::st_write(data_sf, "raw_data/gbif/gbifextract.shp")
