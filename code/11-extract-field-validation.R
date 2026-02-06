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
library(units)
library(terra)

#### Load functions ####
source("code/modelling-functions.R")

sub.cat <- read.csv( "./lookup-tbls/SubstrateCategories.csv", header=T, sep="," )

#---------------------------------------------------------------------#
#### Get SDM field validation dive survey data ####

# Connect to mdb
sdm_mdb <- mdb_connection("//ent.dfo-mpo.ca/dfo-mpo/GROUP/PAC/Reg_Shares/EOS/ES/MSEA/Programs/IMRP/Field_validation/Database/SDM_Field_Validation2025Final.accdb")
# Load queries -- HEADER
sdm_sql <- readLines("code/sql/get-sdm-records.sql")
sdm_sql <- paste(sdm_sql, collapse = "\n ")
# Run queries
dat <- DBI::dbGetQuery( sdm_mdb, sdm_sql )
# Source
dat$Source <- "SDM validation"

length(unique(dat$HKey))

# Load queries -- BELT
sdm_sql <- readLines("code/sql/get-sdm-belt-records.sql")
sdm_sql <- paste(sdm_sql, collapse = "\n ")
# Run queries
belt <- DBI::dbGetQuery( sdm_mdb, sdm_sql )

length(unique(belt$HKey))

diff <- anti_join(dat,belt, by = "HKey")
diff

#### Get centroid of transects to match to predictions ####

#Find sites with matching deep & shallow lat/lon values == cliffs
# Find rows with matching deep and shallow lat/lon values
singlePts_df1 <- dat %>%
  filter((LatDeep == LatShallow) & (LonDeep == LonShallow)) %>% 
  select(HKey, LatDeep, LonDeep) %>% 
  rename(Lat = LatDeep, Lon = LonDeep)

# Remove filtered rows dataset
remove <- singlePts_df1$HKey
df <- dat[!(dat$HKey %in% remove),]


#Create shallow points and deep points 
# Select the deep and shallow columns 
deep <- df %>% 
  select(HKey, LatDeep, LonDeep) %>% 
  rename(Lat = LatDeep, Lon = LonDeep)
head(deep) 
shallow <- df %>%
  select(HKey, LatShallow, LonShallow)%>% 
  rename(Lat = LatShallow, Lon = LonShallow)
head(shallow) 
deep2shallow <- rbind(deep, shallow)
head(deep2shallow)


# Find HKey with shallow coords == 0 or NA --- drop camera or aborted sites
HKey_deepOnly <- shallow %>% 
  filter(is.na(Lat) | Lat == 0 | is.na(Lon) | Lon == 0) %>% 
  pull(HKey)

# Remove these rows from deep2shallow, only complete pairs
deep2shallow <- deep2shallow %>% 
  filter(!HKey %in% HKey_deepOnly)
deep2shallow <- deep2shallow[order(deep2shallow$HKey),]
head(deep2shallow)

# Create new df with these single point records using deep lat/lon values
singlePts_df2 <- deep %>% 
  filter(HKey %in% HKey_deepOnly)
head(singlePts_df2)

# Combine both single points df (drop camera + cliff locations) 
singlePts_df <- rbind(singlePts_df1, singlePts_df2)
rm(singlePts_df1,singlePts_df2, deep, shallow)

# Convert deep2shallow df to sf object 
points_sf <- st_as_sf(deep2shallow, coords = c("Lon", "Lat"), crs = 4326) # CRS 4326 for WGS84
plot(points_sf)

# 3. Create transect lines 
# Create lines from points
transectlines_sf <- points_sf %>%
  group_by(HKey) %>%
  summarise(do_union = FALSE) %>% # do_union = FALSE is important for retaining individual lines
  st_cast("MULTILINESTRING")
plot(transectlines_sf)

# Search for data entry errors - transect lines that are too long and remove --- this code should not be needed once QAQC complete
transectlines_sf$length <- st_length(transectlines_sf)
head(transectlines_sf[order(transectlines_sf$length, decreasing = TRUE), ], 20)
# Remove the unit value
transectlines_sf$length <- drop_units(transectlines_sf$length)

# Create recordset of really long transects for QA/QC 
long_transectlines_sf <- transectlines_sf %>%
  filter(length > 45)
plot(long_transectlines_sf) # These are all sand only MSBIDS sites

transectlines_sf <- transectlines_sf %>% select(-length)

# 4. Create centroid points 
# Calculate the centroid points for complete transects
centroid_point <- st_centroid(transectlines_sf)
print(centroid_point)

# Create points object from single points df
singlePts_sf <- st_as_sf(singlePts_df, coords = c("Lon", "Lat"), crs = 4326) # CRS 4326 for WGS84
plot(singlePts_sf)

# Combine centroids and single points
allpts_sf <- rbind(centroid_point, singlePts_sf)
plot(allpts_sf)

#add back in relative abundance and other info
allpts_sf <- allpts_sf %>%
  left_join(dat, by = "HKey")

# Convert to the appropriate projection 
validation_sf <- st_transform(allpts_sf, crs = 3005)


####extract predictions from sdms####
#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"

bathy_all <- terra::vrt(c("raw_data/envlayers-20m-hg/bathymetry.tif", "raw_data/envlayers-20m-ncc/bathymetry.tif", "raw_data/envlayers-20m-qcs/bathymetry.tif", "raw_data/envlayers-20m-wcvi/bathymetry.tif", "raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif"), "bathy.vrt", overwrite=T)
slope_all <- terra::vrt(c("raw_data/envlayers-20m-hg/slope.tif", "raw_data/envlayers-20m-ncc/slope.tif", "raw_data/envlayers-20m-qcs/slope.tif", "raw_data/envlayers-20m-wcvi/slope.tif", "raw_data/envlayers-20m-shelfsalishsea/slope.tif"), "slope.vrt", overwrite=T)

# get predictions (tried loading from mosaics earlier and didn't work)
parent_dir <- "raster/eelgrass"
folders <- list.dirs(parent_dir, recursive = FALSE)

# load and mosaic eelgrass SDM rasters from 4 different SDM predictions
mosaics <- list()

for (f in folders) {
  tif_files <- list.files(
    f,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  # split files
  se_files  <- tif_files[grepl("se", basename(tif_files), ignore.case = TRUE)]
  pred_files <- tif_files[!grepl("se", basename(tif_files), ignore.case = TRUE)]
  
  folder_name <- basename(f)
  mosaics[[folder_name]] <- list()
  
  # mosaic prediction files (no se)
  if (length(pred_files) > 0) {
    r_pred <- lapply(pred_files, rast)
    mosaics[[folder_name]]$pred <- do.call(mosaic, r_pred)
  }
  
  #mosaic SE files
  if (length(se_files) > 0) {
    r_se <- lapply(se_files, rast)
    mosaics[[folder_name]]$se <- do.call(mosaic, r_se)
  }
}


mosaics_flat <- unlist(mosaics, recursive = FALSE)
sdm_resampled <- lapply(mosaics_flat, function(p) {
  
  lapply(p, function(r) {
    
    r <- terra::writeRaster(r, tempfile(fileext = ".tif"), overwrite = TRUE)
    
    if (!terra::compareGeom(r, bathy_all, crs = TRUE, stopOnError = FALSE)) {
      r <- terra::project(r, bathy_all)
    }
    
    terra::crop(
      terra::resample(r, bathy_all, method = "bilinear"),
      bathy_all
    )
  })
})


rstack <- c(sdm_resampled[[1]], sdm_resampled[[2]], sdm_resampled[[3]], sdm_resampled[[4]], sdm_resampled[[5]], sdm_resampled[[6]], sdm_resampled[[7]], sdm_resampled[[8]])
names(rstack) <- c("eelgrass_bccm_nospatial_pred", "eelgrass_bccm_nospatial_se", "eelgrass_bccm_spatial_pred", "eelgrass_bccm_spatial_se", "eelgrass_nep_nospatial_pred", "eelgrass_nep_nospatial_pred", "eelgrass_nep_spatial_pred", "eelgrass_nep_spatial_se")

rstack <- terra::rast(rstack)


bathy_extract <- terra::extract(bathy_all, validation_sf)
summary(bathy_extract$bathy)
validation_sf$bathy_mod <- bathy_extract$bathy

substrate_extract <- terra::extract(substrate_all, validation_sf)
summary(substrate_extract$substrate)
#validation_sf$substrate_num <- substrate_extract$substrate
validation_sf$substrate_mod <- substrate_extract$substrate
validation_sf$substrate_mod <- c("Rock", "Mixed", "Sand", "Mud")[validation_sf$substrate_mod]

slope_extract <- terra::extract(slope_all, validation_sf)
summary(slope_extract$slope)
validation_sf$slope_mod <- slope_extract$slope

sdms_extract <- terra::extract(rstack, validation_sf)
summary(sdms_extract)
validation_sf <- cbind(validation_sf, sdms_extract)

# remove NA in bathy
validation_sf <- validation_sf %>% filter(!is.na(bathy_mod))

# remove rows only when none of the sdms have predictions
validation_sf <- validation_sf %>%
  filter(!if_all(c(eelgrass_bccm_nospatial_pred,), is.na)) %>%
  mutate(across(
    where(is.character),
    ~ if_else(.x == "None", "Absent", .x)
  ))
# resulted in 38 sites removed. 


#### extract belt information and join with site and predictions ####
# Convert all logical columns to 1's & 0's
belt <- belt %>% mutate(across(where(is.logical), as.numeric))

# add in substrate classification for each belt
belt <- belt %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) %>%
  select(-RMSM.cat,-BType1, -BType2, -Btype.description)

# Sum of p/a observations for species at each site (sum = 0 to 12 observations)
df_sum_species <- belt %>% 
  group_by(HKey) %>% 
  summarise(across(c(RSU,GSU,PSU,RSC,PT,PH,ZM), sum)  ) 
df_sum_species

# Presence/Absence for each target species at each site
num_cols <- which(sapply(df_sum_species, is.numeric))[-1] 
df_pa <- df_sum_species %>%
  mutate(
    across(
      all_of(names(df_sum_species)[num_cols]),
      ~ pmin(., 1)
    )
  )
df_pa



# Frequency of p/a observations for species at each site
# AND average depth at sites
df_depth_spp_freq <- belt %>% 
  group_by(HKey) %>% 
  summarise(
    avgCorDepth_obs = mean(CorDepthM, na.rm = TRUE),
    slope_obs = {
      depth_2R  <- CorDepthM[BeltTransect == "2R"][1]
      depth_22R <- CorDepthM[BeltTransect == "22R"][1]
      depth_18R <- CorDepthM[BeltTransect == "18R"][1]
      
      slope_rad <- if (!is.na(depth_22R) & !is.na(depth_2R)) {
        atan2((depth_2R - depth_22R), 20)
      } else if (!is.na(depth_18R) & !is.na(depth_2R)) {
        atan2((depth_2R - depth_18R), 20)
      } else {
        NA_real_
      }
      
      if (!is.na(slope_rad)) {
        round(abs(slope_rad) * 180 / pi, digits = 0)
      } else {
        NA_real_
      }
    },
    substrate_mode_obs = Mode(Substrate),
    num_belts = n(),
    RSU_freq = sum(RSU == TRUE)/n(),
    GSU_freq = sum(GSU == TRUE)/n(),
    PSU_freq = sum(PSU == TRUE)/n(),
    RSC_freq = sum(RSC == TRUE)/n(),
    PT_freq = sum(PT == TRUE)/n(),
    PH_freq = sum(PH == TRUE)/n(),
    ZM_freq = sum(ZM == TRUE)/n(),
  ) 
df_depth_spp_freq

obs_belts <- df_depth_spp_freq %>% 
  left_join(df_pa, by = "HKey")


# Combine datasets
validation_sf <- left_join(validation_sf, obs_belts, by="HKey")
summary(validation_sf)
validation_sf <- validation_sf %>% select(-RA_RSU, -RA_GSU, -RA_PSU, -RA_RSC, -PC_PT, -RSU_freq, -RSU, -GSU_freq, -GSU, -PSU_freq, -PSU, -RSC_freq, -RSC, -PT_freq, -PT)
save(validation_sf, file = "code/output_data/field_validation/validation_dataset.RData")


