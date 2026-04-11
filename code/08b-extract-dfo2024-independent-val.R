###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
#
# Overview:
# Processing scripts that query and standardize zostera, phyllospadix, percent cover, depth, and substrate observations from Pacific shellfish databases through SQL Server and BHM database through MS Access. 
# Data range is 2024 only
# Calculate slope based off quadrat size, quadrat skipping, and depth observations
# Clean data to remove where divers may have misidentified eelgrass and surfgrass for each other (based on substrate observations)
# For more information of survey protocols see Data-dictionary.txt 


# Data use:
# As per the Policy for Scientific Data, DFO scientific data are a public resource and subject to full and open access within two years of being acquired or generated. 
# As per the Application for Access to Shellfish Data Holdings, data collected two years prior to today should not be included in data pulls to allow the Principle Investigator a two year period to analyze and report findings. 
# Exceptions to the use of data collected in the last two years may be made if the user has met the conditions outlined in the Application for Access or has discussed the use of the data with the Principle Investigator for the survey

#Script built in R 4.4.0

# Requirements:
# r-sql-link-functions.R 
# Access to shellfish SQL server (VPN connection and access permitted from DFO Shellfish Data Unit)
# Access to benthic habitat mapping access database on DFO Spatial Datasets drive (VPN connection and access permitted from MSEA section)
###############################################################################



#### load packages####
library(DBI)
library(odbc)
library(reshape2)
library(sf)
library(tidyverse)
library(terra)
library(matrixStats)
library(geosphere)
library(Hmisc)

#### Load util.R gfdata functions ####
source("code/r-sql-link-functions.R")

#### Read in substrate category table ####
sub.cat <- read.csv( "./lookup-tbls/SubstrateCategories.csv", header=T, sep="," )

coastline_full <- st_read("raw_data/CHS_HWL2015_Coastline.gdb", layer = "Line_CHS_Pacific_HWL_2015_5028437")
boundary <- coastline_full %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_crop(st_bbox(coastline_full) + c(0, 100000, -10000, -80000)) %>%
  st_transform(crs = "EPSG:3005") 

# Load bathy raster
bathy <- vrt("raw_data/Bathymetry/Coastwide_20m_mosaic.tif")

#### Outputs ####
# Subfolders for outputs
outdir <- file.path("code", "output_data")
# Create the main directory, and subfolders
dir.create( outdir, recursive=TRUE )

# final field names for dive surveys
fnames <-  c( "Type", "Source", "Survey", "HKey", "Method",
              "Year","Month", "Day", "LatDeep","LonDeep",
              "LatShallow","LonShallow", "Transect_length", "Quadrat", "CorDepthM","Substrate", "Slope", "PerCovZO",
              "PH", "ZO", "Identification")  

#---------------------------------------------------------------------#
#### Get GSU_bio (green sea urchin) dive survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2005 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)

## Extract data from SQL Server
# Read SQL query
gsu_sql <- readLines("code/sql/2024/get-gsu-records.sql")
# Load query
gsu_sql <- paste(gsu_sql, collapse = "\n ")
# Run query
gsu_dat <- DBI::dbGetQuery( sf_db_connection(), gsu_sql )

## Check for duplicates of species records from a single finishing event
dup_ind <- paste(gsu_dat$HKey, gsu_dat$Transect, gsu_dat$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( gsu_dat[inds,] )
  # Remove duplicate records
  gsu_dat <- gsu_dat[!duplicated(gsu_dat),]
}

## Calculate slope using arc tangent method, make assumption that slope is same on Quadrat 1 as Quadrat 2 (since we don't have starting depth at bottom)
gsu_dat <- gsu_dat %>%
  mutate(StartDepth = ifelse(Quadrat==1, NA, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff, Quadrat_distance)) %>%
  mutate(Slope = ifelse(Quadrat==1, lead(Slope, n=1), Slope)) %>%
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
gsu_dat$Slope <- round(abs(gsu_dat$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
gsu_dat <- gsu_dat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#get percent cover of ZO (not possible to get percent cover of PH as outher algae would be found alongside PH, for ZO algae not as likely to also be present, i.e. ZO usually in monoculture)
gsu_dat <- gsu_dat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, UnderstoryPct)) 

gsu_dat$PerCovZO[gsu_dat$PerCovZO > 100] <- 100

## Melt all invert and algae species into species column
#gsu_dat <- melt(gsu_dat, measure.vars=sp_pa, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
gsu_dat$Type <- "Research"
# Source
gsu_dat$Source <- "GSU_bio"
# GSU doesn't have a consistently filled transect field in database
gsu_dat$Transect <- NA
# Method
gsu_dat$Method <- "Dive"

## remove where phyllo and zostera are likely mis identified based on substrate type
gsu_dat <- gsu_dat %>%
  mutate(Identification = case_when((ZO == 1 & Btype.description == "Bedrock dominant" &(Substrate3 != 7 & Substrate3 != 9) %>% replace_na(TRUE)) ~ "Remove",
                                    Substrate != "Rock" & Substrate != "Mixed" & PH == 1 & Btype.description != "Sand/shell, some rock" ~ "Remove",
                                    .default = "Keep"))

# Rename
gsu_dat <- gsu_dat[fnames]

# Order rows
gsu_dat <- gsu_dat[order(gsu_dat$HKey, gsu_dat$Year, gsu_dat$Transect, gsu_dat$Quadrat),]



#---------------------------------------------------------------------#
#### Get Cuke_bio (California sea cucumber) dive survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2000 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)

## Extract data from SQL Server
# Read SQL query
rsc_sql <- readLines("code/sql/2024/get-cuke-records.sql")
# Load query
rsc_sql <- paste(rsc_sql, collapse = "\n ")
# Run query
rsc_dat <- DBI::dbGetQuery( sf_db_connection(), rsc_sql )

## Check for duplicates of species records from a single finishing event
dup_ind <- paste(rsc_dat$HKey, rsc_dat$Transect, rsc_dat$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( rsc_dat[inds,] )
  # Remove duplicate records
  rsc_dat <- rsc_dat[!duplicated(rsc_dat),]
}

#Years >2016 have Quad0 so know bottom depths, surveys prior we often have no start depth of first quadrat of each transect but making assumption it has similar slope as quadrat after it. 
rsc_dat <- rsc_dat %>%
  mutate(StartDepth = case_when(Year>2016 & Quadrat==0 ~ CorDepthM,
                                Year>2016 & Quadrat!=0 ~ lag(CorDepthM, n=1),
                                Year<2017 & Quadrat==1 ~ NA,
                                Year<2017 & Quadrat!=1 ~ lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff,5)) %>%
  mutate(Slope = ifelse(Year<2017 & Quadrat==1, lead(Slope, n=1), Slope)) %>%
  filter(Quadrat!=0) %>% # Need to remove quad 0 from each transect not all transects have Quad zero and it only has depth recorded
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
rsc_dat$Slope <- round(abs(rsc_dat$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
rsc_dat <- rsc_dat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#get percent cover of ZO (not possible to get percent cover of PH as outher algae would be found alongside PH, for ZO algae not as likely to also be present, i.e. ZO usually in monoculture)
rsc_dat <- rsc_dat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, PctCover)) 

rsc_dat$PerCovZO[rsc_dat$PerCovZO > 100] <- 100

## Melt all invert and algae species into species column
#rsc_dat <- melt(rsc_dat, measure.vars=sp_pa, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
rsc_dat$Type <- "Research"
# Source
rsc_dat$Source <- "Cuke_bio"
# Method
rsc_dat$Method <- "Dive"

## remove where phyllo and zostera are likely mis identified based on substrate type
rsc_dat <- rsc_dat %>%
  mutate(Identification = case_when((ZO == 1 & Btype.description == "Bedrock dominant" &(Substrate3 != 7 & Substrate3 != 9 & Substrate3 != 10) %>% replace_na(TRUE)) ~ "Remove",
                                    Substrate != "Rock" & Substrate != "Mixed" & PH == 1 & Btype.description != "Sand/shell, some rock" & Substrate3 != 3 ~ "Remove",
                                    .default = "Keep"))

# Rename
rsc_dat <- rsc_dat[fnames]

# Order rows
rsc_dat <- rsc_dat[order(rsc_dat$HKey, rsc_dat$Year, rsc_dat$Transect, rsc_dat$Quadrat),]


#---------------------------------------------------------------------#
#### Get Multispecies_bio (multi-species) dive survey data ####
# P/A observations: PH, ZO 
# These DFO surveys began in 2016
# SQL code has filtered out Start Lat and Long are not null and CorDepth is not null
# Transect length is rounded up to the nearest 25m (there is no record of transect length in the database, so quadrat skipping info was used)
# All species observations are converted to presence (1)/ absence (0)

## Extract data from SQL Server
# Read SQL query
multi_sql <- readLines("code/sql/2024/get-multispecies-records.sql")
# Load query
multi_sql <- paste(multi_sql, collapse = "\n ")
# Run query
multi_dat <- DBI::dbGetQuery( sf_db_connection(), multi_sql )

## Check for duplicates of species records from a single sampling event
dup_ind <- paste(multi_dat$HKey, multi_dat$Transect, multi_dat$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( multi_dat[inds,] )
  # Remove duplicate records
  multi_dat <- multi_dat[!duplicated(multi_dat),]
}


## Calculate slope using arc tangent method, make assumption that slope is same on Quadrat 1 as Quadrat 2 (since we don't have starting depth at bottom)
multi_dat <- multi_dat %>%
  mutate(StartDepth = ifelse(Quadrat==1, NA, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff,Quadrat_distance)) %>%
  mutate(Slope = ifelse(Quadrat==1, lead(Slope, n=1), Slope)) %>%
  select(-StartDepth,-EndDepth, -Elev.Diff)

#change to degrees
multi_dat$Slope <- round(abs(multi_dat$Slope) * 180/pi, digits = 0)


## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
multi_dat <- multi_dat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#get percent cover of ZO (not possible to get percent cover of PH as outher algae would be found alongside PH, for ZO algae not as likely to also be present, i.e. ZO usually in monoculture)
multi_dat <- multi_dat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, UnderstoryPct)) 

multi_dat$PerCovZO[multi_dat$PerCovZO > 100] <- 100

## Melt all invert and algae species into species column
#multi_dat <- melt(multi_dat, measure.vars=sp_pa, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
multi_dat$Type <- "Research"
# Source
multi_dat$Source <- "Multispecies_bio"
# No time type in multispecies db
multi_dat$Time_type <- NA
# Method
multi_dat$Method <- "Dive"

## remove where phyllo and zostera are likely mis identified based on substrate type
multi_dat <- multi_dat %>%
  mutate(Identification = case_when(ZO == 1 & Btype.description == "Bedrock dominant"  ~ "Remove",
                                    Substrate != "Rock" & Substrate != "Mixed" & PH == 1 & Btype.description != "Sand/shell, some rock" ~ "Remove",
                                    .default = "Keep"))

# Rename
multi_dat <- multi_dat[fnames]

# Order rows
multi_dat <- multi_dat[order(multi_dat$HKey, multi_dat$Year, multi_dat$Transect, multi_dat$Quadrat),]



#---------------------------------------------------------------------#
#### Get GDK_bio (geoduck) dive survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2005 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)
# This dataset has different spatialization requirements so is kept separate

## Extract data from SQL Server
# Read SQL query
gdk_sql <- readLines("code/sql/2024/get-geoduck-records.sql")
# Load query
gdk_sql <- paste(gdk_sql, collapse = "\n ")
# Run query
gdk_dat <- DBI::dbGetQuery( sf_db_connection(), gdk_sql )

## Check for duplicates of species records from a single finishing event
dup_ind <- paste(gdk_dat$HKey, gdk_dat$Transect, gdk_dat$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( gdk_dat[inds,] )
  # Remove duplicate records
  gdk_dat <- gdk_dat[!duplicated(gdk_dat),]
}

#Most transects have a Quad0 so know bottom depths, but some don't so need to find those
gdk_nozeroquad <- gdk_dat %>%
  group_by(HKey) %>%
  filter(!all(0 %in% Quadrat)) %>%
  ungroup %>%
  distinct(HKey) %>%
  pull(HKey)

#errors in database of transect_dist_from_start
gdk_dat[gdk_dat$HKey=="10491" & gdk_dat$Quadrat== 15, "Transect_dist_from_start"] <- 300
gdk_dat[gdk_dat$HKey=="10493" & gdk_dat$Quadrat== 15, "Transect_dist_from_start"] <- 215
gdk_dat[gdk_dat$HKey=="10541" & gdk_dat$Quadrat== 12, "Transect_dist_from_start"] <- 170
gdk_dat[gdk_dat$HKey=="10545" & gdk_dat$Quadrat== 11, "Transect_dist_from_start"] <- 205
gdk_dat[gdk_dat$HKey=="10545" & gdk_dat$Quadrat== 12, "Transect_dist_from_start"] <- 225
gdk_dat[gdk_dat$HKey=="10545" & gdk_dat$Quadrat== 13, "Transect_dist_from_start"] <- 245
gdk_dat[gdk_dat$HKey=="10545" & gdk_dat$Quadrat== 14, "Transect_dist_from_start"] <- 265
gdk_dat[gdk_dat$HKey=="10545" & gdk_dat$Quadrat== 15, "Transect_dist_from_start"] <- 285
gdk_dat[gdk_dat$HKey=="10545" & gdk_dat$Quadrat== 16, "Transect_dist_from_start"] <- 305

# Calculate slope using arc tangent method
gdk_dat <- gdk_dat %>%
  filter(!HKey %in% 13094:13123) %>% # these only have two quadrats in each transect
  mutate(StartDepth = ifelse(Quadrat==0 | (Quadrat==1 & HKey %in% gdk_nozeroquad), CorDepthM, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Quadrat_distance = ifelse(Quadrat==0 | (Quadrat==1 & HKey %in% gdk_nozeroquad), Transect_dist_from_start, Transect_dist_from_start - lag(Transect_dist_from_start, n = 1)),
         Slope = atan2(Elev.Diff, Quadrat_distance)) %>%
  filter(Quadrat!=0) %>% # Need to remove quad 0 from each transect not all transects have Quad zero and it only has depth recorded
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
gdk_dat$Slope <- round(abs(gdk_dat$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
gdk_dat <- gdk_dat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#add percent cover of ZO, this survey does not record percent cover so all absences are assigned 0 and presences assigned NA
gdk_dat <- gdk_dat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, NA)) 


## Melt all invert and algae species into species column
#gdk_dat <- melt(gdk_dat, measure.vars=sp_pa, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
gdk_dat$Type <- "Research"
# Source
gdk_dat$Source <- "GDK_bio"
# Method
gdk_dat$Method <- "Dive"
O <- NA

## remove where phyllo and zostera are likely mis identified based on substrate type
gdk_dat <- gdk_dat %>%
  mutate(Identification = case_when(ZO == 1 & Btype.description == "Bedrock dominant" &(Substrate3 != 7 & Substrate3 != 10) %>% replace_na(TRUE)  ~ "Remove",
                                    Substrate != "Rock" & Substrate != "Mixed" & PH == 1 ~ "Remove",
                                    .default = "Keep"))

# Rename
gdk_dat <- gdk_dat[fnames]

# Order rows
gdk_dat <- gdk_dat[order(gdk_dat$HKey, gdk_dat$Year, gdk_dat$Transect, gdk_dat$Quadrat),]




#### Get BHM (Benthic Habitat Mapping) dive survey data ####
# P/A observations: PH (phyllospadix), ZO (eelgrass) 
# All species observations are converted to presence (1)/ absence (0)

# Connect to mdb
bhm_mdb <- mdb_connection("//ent.dfo-mpo.ca/dfo-mpo/GROUP/PAC/Reg_Shares/EOS/ES/Spatial_Datasets/Dive_Surveys/Database/BHM_DiveSurveys_Jan2025.mdb")

# Load queries
bhm_sql <- readLines("code/sql/2024/get-bhm-records-quads.sql")
bhm_sql <- paste(bhm_sql, collapse = "\n ")
bhm_spp_sql <- readLines("code/sql/2024/get-bhm-records-spp.sql")
bhm_spp_sql <- paste(bhm_spp_sql, collapse = "\n ")

# Run queries
bhm_quads <- DBI::dbGetQuery( bhm_mdb, bhm_sql )
bhm_spp <- DBI::dbGetQuery( bhm_mdb, bhm_spp_sql )


#edit quads table to calculate quadat distance
bhm_quads <- bhm_quads %>%
  mutate(Quadrat_distance = case_when(Quadrat==1 ~ 5,#the slope of first quadrat on transect is always 5
                                      QuadratSkipping==1 ~ 5, #sampled every  quadrat
                                      QuadratSkipping==2 ~ 10, #sampled every other quadrat
                                      QuadratSkipping==3 ~ 15)) ##sampled every 3rd quadrat

## Edit species table ###
bhm_spp$SpType <- toupper(bhm_spp$SpType) # capitalize

#create unique species codes
bhm_spp$Species <- paste0(bhm_spp$SpType,"_",bhm_spp$SpeciesCode)

bhm_spp$Species[bhm_spp$Species=="A_ZO"] <- "ZO" 
bhm_spp$Species[bhm_spp$Species=="A_PH"] <- "PH" 

sppQuad <- reshape2::dcast( bhm_spp, HKey+Quadrat~Species, fun=length, value.var = "Species" )
sppQuad <- sppQuad %>% select(HKey, Quadrat, ZO, PH)

#join species data to quads
bhm_dat <- bhm_quads %>% 
  left_join(sppQuad, by = join_by(HKey, Quadrat)) %>% 
  mutate(ZO = ifelse(is.na(ZO), 0, ZO), 
         PH = ifelse(is.na(PH), 0, PH))

#change to zero/one
bhm_dat$ZO[bhm_dat$ZO > 0] <- 1 
bhm_dat$PH[bhm_dat$PH > 0] <- 1 

# Calculate slope using arc tangent method 
bhm_dat <- bhm_dat %>%
  mutate(StartDepth = ifelse(Quadrat==0, CorDepthM, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff, Quadrat_distance)) %>%
  filter(Quadrat!=0) %>% # Need to remove quad 0 from each transect 
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
bhm_dat$Slope <- round(abs(bhm_dat$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
bhm_dat <- bhm_dat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#add percent cover of ZO
bhm_dat <- bhm_dat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, UnderstoryPct)) 

bhm_dat$PerCovZO[bhm_dat$PerCovZO > 100] <- 100

## Melt all invert and algae species into species column
#bhm_dat <- melt(bhm_dat, measure.vars=sp_pa, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
bhm_dat$Type <- "Research"
# Source
bhm_dat$Source <- "BHM"
# Method
bhm_dat$Method <- "Dive"

## remove where phyllo and zostera are likely mis identified based on substrate type
bhm_dat <- bhm_dat %>%
  mutate(Identification = case_when(ZO == 1 & Btype.description == "Bedrock dominant" &(Substrate3 != 7 & Substrate3 != 9 & Substrate3 != 10) %>% replace_na(TRUE)  ~ "Remove",
                                    Substrate != "Rock" & Substrate != "Mixed" & PH == 1 & Btype.description != "Sand/shell, some rock" &(Substrate3 != 2 & Substrate3 != 4) %>% replace_na(TRUE) ~ "Remove",
                                    .default = "Keep"))
# Rename
bhm_dat <- bhm_dat[fnames]

# Order rows
bhm_dat <- bhm_dat[order(bhm_dat$HKey, bhm_dat$Year, bhm_dat$Transect, bhm_dat$Quadrat),]



####Combine datasets ####
dat <- bind_rows(gsu_dat, rsc_dat, multi_dat, gdk_dat, bhm_dat)

#need to make new HKey incase there is overlap in HKey between source types
dat <- dat %>% 
  mutate (HKey = paste0(Source,"_",HKey))

# when ZO has value of 1 and percent cover has value of 0, need to change to NA
dat$PerCovZO[dat$PerCovZO == 0 & dat$ZO == 1] <- NA

# remove ZO presence observations in high intertidal, likely Z. japonica
dat$Identification[dat$ZO == 1 & dat$CorDepthM < -2] <- "Remove"

# make shapefile of dat to further explore any issues 

# Convert to spdf and export
# create spatial points
dat_sf <- dat %>%
  mutate(Latitude = case_when(!is.na(LatDeep) ~ LatDeep,
                              is.na(LatDeep) ~ LatShallow),
         Longitude = case_when(!is.na(LonDeep) ~ LonDeep,
                               is.na(LonDeep) ~ LonShallow)) %>%
  filter(!is.na(Latitude))  %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

# export as shapefile
# likely to have issues with attribute field names shortening
st_write(dat_sf, "code/output_data/processed_observations/dat_transect_2024.shp", append=FALSE)

#----------------------------------------------------------------------------#
# Clean input data
# Find missing start and end latitude and longitude values

#Find number of quadrats per HKey
qCnt <- dat %>% 
  group_by(HKey) %>%
  dplyr::summarize( qCnt = length(Quadrat))

# Add number of quadrats to df
dat <- left_join(dat, qCnt, by="HKey")

# Find transects which don't have x,y values for either the start or end and remove
missing_both <- which( apply(is.na(dat[c("LonShallow","LatShallow","LonDeep","LatDeep")]), 1, all) )
no_xy <- dat[missing_both,]

#find transects that only have one xy coodinate per transect
missing_startorend <- which( !complete.cases(dat[, c("LonShallow","LatShallow","LonDeep","LatDeep")]) )
single_xy <- dat[c(missing_startorend),]
unique(single_xy$Source) 

#find transects that have same lat and lon 
same_startend <- which( paste0(dat$LonDeep, dat$LatDeep) == paste0(dat$LonShallow, dat$LatShallow))
same_xy <- dat[c(same_startend),]

#make cliffs dataset
cliffs <- dat[c(missing_startorend, same_startend),]
cliffs <- cliffs%>% filter(Transect_length < 41) #this is the cliffs data set to be spatialized separately and retain all quadrats
unique(cliffs$HKey) #3 transects are cliffs

#clean up dat
dat_for_trans <- dat[-c(missing_startorend,same_startend, missing_both),] #remove those from dat, took out  becasue we spatialize some of these. Need to confirm this right to do

single_xy <- single_xy%>% filter(Transect_length > 40) #remove cliffs from single_xy
unique(single_xy$HKey) #0 transects have single xy
same_xy <- same_xy%>% filter(Transect_length > 40) #remove cliffs from same_xy 
unique(same_xy$HKey) #2 transects have same xy

# #likely can only spatialize if have deep end coordinate and transect length is <200m (otherwise other points of land may confuse spatialization)
# single_xy_tospatialize <- single_xy %>% 
#   filter(Transect_length < 201,
#          !is.na(LatDeep),
#          !is.na(LonDeep))
# unique(single_xy_tospatialize$HKey) #1666 transects have single xy that can be spatialized

# #get seperate dataframe of where we can extract some quadrats at the one point
# single_xy_topoint <- single_xy %>% 
#   filter(Transect_length > 200|Transect_length < 201 & is.na(LatDeep)) %>%
#   rbind(same_xy) # add in same_xy transects
# unique(single_xy_topoint$HKey) #396 transects have single xy we can get quadrats from one point

# Create transect dataset with filtered dat dataset
trans <- dat_for_trans[!duplicated(dat_for_trans$HKey), 
                       c("Survey","Source", "Year","Month","Day","HKey",
                         "LonShallow","LatShallow","LonDeep","LatDeep")]

# # Create transect dataset with filtered single_xy to spatialize dataset
# single_xy_trans <- single_xy_tospatialize[!duplicated(single_xy_tospatialize$HKey),
#                                           c("Survey","Source","Year","Month","Day","HKey",
#                                             "LonShallow","LatShallow","LonDeep","LatDeep")]


#----------------------------------------------------------------------------#
# Create spatial lines from shallow to deep transect points
# extended by 'dist_extend' from the deep point
# extended by distance to coastline or 'threshold' from shallow point, threshold decided to be 50m
# (shallow point tends not to reach shoreline)

# For transects with a unique shallow and deep position,
# Calculate distance from shallow point to deep point
shallow.sf <- trans %>% st_as_sf(coords = c("LonShallow", "LatShallow"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")
deep.sf <- trans %>% st_as_sf(coords = c("LonDeep", "LatDeep"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

s2e <- st_distance(shallow.sf, deep.sf, by_element = TRUE)
trans$length <- as.integer(s2e)

# Remove transects with start and end points greater than this distance apart.
# dist_transect apart because either the deep or shallow position is likely wrong
trans <- trans %>%
  mutate(Keep = case_when(length<800 & Source=="GSU_bio" ~ 'Yes',
                          length<400 & Source=="Cuke_bio" ~ 'Yes',
                          length<2000 & Source=="GDK_bio" ~ 'Yes',
                          length<350 & Source=="BHM" ~ 'Yes',
                          length<150 & Source=="Multispecies_bio" ~ 'Yes',
                          TRUE ~ 'No')) %>%
  filter(Keep == "Yes") %>%
  select(-Keep)


#Convert shallow points to spatial class with dataset now filtered to remove transects with likely incorrect xy
shallow.sf<- trans %>% st_as_sf(coords = c("LonShallow", "LatShallow"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

nearest <- st_distance(shallow.sf, boundary)
nearest_min <- matrixStats::rowMins(nearest)
trans$distance <- as.integer(nearest_min)

# If distance to coastline is greater than threshold, use threshold, which is 50m
trans$distance[ trans$distance > 50 ] <- 50

# Create featureID and re-assign row names
trans$featureID <- 1:nrow(trans)
row.names(trans) <- trans$featureID

# Distance to extend deep transect point toward deeper depths
dist_extend <- 50

# convert to lines extended past points along transect bearing
lines.list <- list()
for (i in 1:nrow(trans)){
  # get deep, shallow and distance from shallow point to coastline
  deep <- trans[i, c("LonDeep","LatDeep")]
  shallow <- trans[i, c("LonShallow","LatShallow")]
  d <- trans[i, "distance"]
  # if deep and shallow points differ create lines
  if ( ! all(deep == shallow) ){
    # calculate bearing from deep to shallow
    b <- finalBearing(deep, shallow)
    # extend point from shallow (nearshore)
    e.shallow <- destPoint(p=shallow, b=b, d=d) # extend by the nearest distance to shore
    # extend point from deep
    e.deep <- destPoint(p=deep, b=b-180, d=dist_extend) # extend 50m past deep point
    # create line from extended points
    l <- rbind( as.matrix( e.deep ),
                as.matrix( e.shallow ))
    l <- l %>% as.data.frame %>% 
      st_as_sf(coords = c("lon","lat"), crs = "EPSG:4326") %>%
      st_transform(crs = "EPSG:3005") %>%
      group_by() %>% 
      dplyr::summarize() %>%
      st_cast("LINESTRING")
    # hkey
    hkey <- as.character( trans[i, "HKey"] )
    # add to list
    lines.list[[hkey]] <- st_as_sfc(l)
  }
}

# create spatial lines and spatial lines dataframe
sl<- lines.list %>%
  bind_rows() %>%
  pivot_longer(cols = everything(), names_to = as.character("HKey"), values_to = "Lines") %>%
  st_as_sf(crs = "EPSG:3005")


#----------------------------------------------------------------------------#
# Generate points along the transect at a certain distance set by 'dist_points'

# Generate points
spdf.list <- list()
for( i in 1:nrow(sl)){
  # Get Hkey
  hkey <- trans$HKey[trans$featureID == i]
  # Get length of line
  linelength <- st_length(sl[i,])
  # Number of points for breaks at 'dist_points' distance
  npts <- ceiling( linelength / 20 ) # 20 becasue 20 m bathy layer
  # Set seed for consistent results
  set.seed(42)
  # Regularly sample points along line
  pts <- st_line_sample(sl[i,], n = npts, type = "regular") 
  ckey <- as.character(hkey)
  spdf.list[[ckey]] <-
    data.frame(st_coordinates(st_as_sf(pts, crs =	"EPSG:3005")),
               featureID = rep(i, npts),
               HKey = rep(hkey, npts)) 
}

# bind all points together
spdf <- do.call("rbind", spdf.list)
spdf<- spdf %>%
  st_as_sf(coords = c("X", "Y"), crs = 	"EPSG:3005")


##############----------------------------------------------------------------------------#
# Add points with only a single x,y at end to spatialize
# Create featureID and re-assign row names

# singlexy.sf <- single_xy_trans %>% st_as_sf(coords = c("LonDeep", "LatDeep"), crs = "EPSG:4326") %>%
#   st_transform(crs = "EPSG:3005")
# 
# # create spatial lines and remove transects likely not correct due to incorrect xy by distance
# sl2<- singlexy.sf %>% 
#   mutate(
#     Lines = st_nearest_points(geometry, boundary[st_nearest_feature(geometry, boundary, pairwise = TRUE),], pairwise = TRUE),
#     closest_point = st_cast(Lines, 'POINT')[seq(2, nrow(.)*2, 2)],
#     length = as.integer(st_length(Lines)))
# 
# #filter transect database by distance
# single_xy_trans$length <- sl2$length
# 
# single_xy_trans <- single_xy_trans %>%
#   mutate(Keep = case_when(length<800 & Source=="GSU_bio" ~ 'Yes',
#                           length<400 & Source=="Cuke_bio" ~ 'Yes',
#                           length<2000 & Source=="GDK_bio" ~ 'Yes',
#                           length<350 & Source=="BHM" ~ 'Yes',
#                           length<150 & Source=="Multispecies_bio" ~ 'Yes',
#                           TRUE ~ 'No')) %>%
#   filter(Keep == "Yes") %>%
#   select(-Keep)
# 
# #filter spatial lines dataset
# sl2<- sl2 %>%  
#   mutate(Keep = case_when(length<800 & Source=="GSU_bio" ~ 'Yes',
#                           length<400 & Source=="Cuke_bio" ~ 'Yes',
#                           length<2000 & Source=="GDK_bio" ~ 'Yes',
#                           length<350 & Source=="BHM" ~ 'Yes',
#                           length<150 & Source=="Multispecies_bio" ~ 'Yes',
#                           TRUE ~ 'No')) %>%
#   filter(Keep == "Yes") %>%
#   select(HKey, Lines) %>%
#   st_drop_geometry() %>%
#   st_as_sf(crs = "EPSG:3005")
# 
# single_xy_trans$featureID <- 1:nrow(single_xy_trans)
# row.names(single_xy_trans) <- single_xy_trans$featureID

# Generate points along the transect at a certain distance set by 'dist_points'
# Generate points
# spdf.list2 <- list()
# for( i in 1:nrow(sl2)){
#   # Get Hkey
#   hkey <- single_xy_trans$HKey[single_xy_trans$featureID == i]
#   # Get length of line
#   linelength <- st_length(sl2[i,])
#   # Number of points for breaks at 'dist_points' distance
#   npts <- ceiling( linelength / 20 ) # 20 becasue 20 m bathy layer
#   # Set seed for consistent results
#   set.seed(42)
#   # Regularly sample points along line
#   pts <- st_line_sample(sl2[i,], n = npts, type = "regular") 
#   ckey <- as.character(hkey)
#   spdf.list2[[ckey]] <-
#     data.frame(st_coordinates(st_as_sf(pts, crs =	"EPSG:3005")),
#                featureID = rep(i, npts),
#                HKey = rep(hkey, npts)) 
# }

# bind all points together
# spdf2 <- do.call("rbind", spdf.list2)
# spdf2<- spdf2 %>%
#   st_as_sf(coords = c("X", "Y"), crs = 	"EPSG:3005")
# 
# # extract depth from bathy raster
# extractdepth3 <- terra::extract(x = bathy, y = spdf2)
# names(extractdepth3) <-c("ID", "bathy")
# 
# # combine with spdf
# ptsdat2<- cbind(st_drop_geometry(spdf2), st_coordinates(st_as_sf(spdf2, coords = c("x", "y"), crs = "EPSG:3005")), bathy=extractdepth3$bathy)



#----------------------------------------------------------------------------#
# Extract depth from bathymetry raster at spdf points

# extract depth from bathy raster
extractdepth <- terra::extract(x = bathy, y = spdf)
names(extractdepth) <-c("ID", "bathy")

# combine with spdf
ptsdat<- cbind(st_drop_geometry(spdf), st_coordinates(st_as_sf(spdf, coords = c("x", "y"), crs = "EPSG:3005")), bathy=extractdepth$bathy)

# ptsdat <- ptsdat %>% rbind(ptsdat2)

#----------------------------------------------------------------------------#
# Merge points with quadrats using closest depth match
# If loop throws error increase tol in find.matches (doesn't work when set to Inf)

# Remove HKeys not in ptsdat
qdat <- dat[dat$HKey %in% unique(ptsdat$HKey),]


# Match pts to quadrats by depth
matchSpatial <- function( x ){  
  # require
  require(Hmisc)
  # subset by hkey
  quad <- qdat[qdat$HKey == x,]
  pts <- ptsdat[ptsdat$HKey == x,]
  # match
  matchdepth <- find.matches(quad$CorDepthM, pts$bathy, tol=1000, maxmatch=1)
  # merge quad and pts based off matchdepth
  mdat <- data.frame( quad, pts[matchdepth$matches, c("X", "Y", "bathy")],
                      ID = paste(quad$HKey, matchdepth$matches, sep="_"))
  # Calculate difference between quadrat depth and bathy
  mdat$depthdiff <- abs( mdat$CorDepthM - mdat$bathy )
  # return
  return( mdat )
}

# Run getDist in parallel
cl <- parallel::makeCluster( parallel::detectCores() - 1 )
## make variables available to cluster
parallel::clusterExport( cl, varlist=c("qdat","ptsdat") )
# run on cluster
spatialised.list <- parallel::parLapply( cl, unique(ptsdat$HKey), matchSpatial )
## stop cluster
parallel::stopCluster( cl )

# bind data.frames together
spatialised <- do.call("rbind", spatialised.list)

# Remove quadrats where bathy is NA becasue outside extent of bathymetry
spatialised <- spatialised[ which(complete.cases(spatialised$bathy)), ]



#### cliffs dataset
# get lat lon of one point
cliffs<- cliffs %>%
  mutate(LatDeep = case_when(!is.na(LatDeep) ~ LatDeep,
                             is.na(LatDeep) ~ LatShallow),
         LonDeep = case_when(!is.na(LonDeep) ~ LonDeep,
                             is.na(LonDeep) ~ LonShallow))%>%
  filter(!is.na(LatDeep))

cliffs.sf<- cliffs %>% 
  st_as_sf(coords = c("LonDeep", "LatDeep"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

cliffs.spdf<- cbind(cliffs, st_coordinates(st_as_sf(cliffs.sf, coords = c("x", "y"), crs = "EPSG:3005")))
cliffs.spdf$bathy <- NA
cliffs.spdf$ID <- paste(cliffs.spdf$HKey, "1", sep="_")
cliffs.spdf$depthdiff <- NA

spatialised<- spatialised%>%
  rbind(cliffs.spdf) # add in cliffs to spatialized

#### sites single xy that are just one point
# get lat lon of one point
# single_xy_topoint<- single_xy_topoint %>%
#   mutate(LatDeep = case_when(!is.na(LatDeep) ~ LatDeep,
#                              is.na(LatDeep) ~ LatShallow),
#          LonDeep = case_when(!is.na(LonDeep) ~ LonDeep,
#                              is.na(LonDeep) ~ LonShallow))%>%
#   filter(!is.na(LatDeep))
# 
# singleptxy.sf <- single_xy_topoint %>% st_as_sf(coords = c("LonDeep", "LatDeep"), crs = 	"EPSG:4326") %>%
#   st_transform(crs = "EPSG:3005")
# 
# # extract depth from bathy raster
# extractdepth2 <- terra::extract(x = bathy, y = singleptxy.sf)
# names(extractdepth2) <-c("ID", "bathy")
# 
# # combine with singleptxy.sf
# singleptxy.spdf<- cbind(single_xy_topoint, st_coordinates(st_as_sf(singleptxy.sf, coords = c("x", "y"), crs = "EPSG:3005")), bathy=extractdepth2$bathy)
# 
# singleptxy.spdf$ID <- paste(singleptxy.spdf$HKey, "1", sep="_")
# singleptxy.spdf$depthdiff <- abs( singleptxy.spdf$CorDepthM - singleptxy.spdf$bathy)
# 
# #filter out records with depth diff >5m
# singleptxy.spdf <- singleptxy.spdf %>% 
#   filter(depthdiff<= 5)
# 
# spatialised<- spatialised%>%
#   rbind(singleptxy.spdf) # add in single xy point to spatialized

#----------------------------------------------------------------------------#
# Clean up - Remove points that are likely incorrect

spatialised<- spatialised %>%
  filter(Identification == "Keep")


# check depth diff
summary(spatialised$depthdiff)

# Combine transect and quadrat to get unique quadrat identifier
spatialised$QID <- paste(spatialised$HKey, spatialised$Quadrat, sep="_")

#not including this one for now as diver depths are more accurate than bathy raster
# Remove quadrats with depth difference > 5
# between corrected quadrat depth and depth from bathymetry raster
#spatialised <- spatialised[which(spatialised$depthdiff <= 5),]

# Remove NAs in x and y's (usually due to NAs present in CorDepthM)
spatialised <- spatialised[complete.cases(spatialised[, c("X","Y")]),]

#Mean number of quadrats aggregated into a single spatial point
nquads <- aggregate( Quadrat ~ ID, data= spatialised, function(x) length(unique(x)))
cat( round( mean(nquads$Quadrat), 1 ), "quadrats per spatial point on average",
     "(range", min(nquads$Quadrat), "to", max(nquads$Quadrat), "quadrats)","\n\n")
names(nquads)[2] <- "NumQuadrats"

# Quadrats retained
#select 
spatialised<-spatialised %>% select("Survey","Year","Month","Day","HKey","ID" , "X", "Y",
                                    "LonDeep","LatDeep","LonShallow","LatShallow", "CorDepthM", "Slope", "Substrate", "PerCovZO", "PH", "ZO")

# Convert to spdf
spatialised.spdf <- spatialised %>%
  st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 


#----------------------------------------------------------------------------#

# Aggregate presence / absence by spatial points
#
# # Transect attributes
att <- spatialised[!duplicated(spatialised$ID),]
att <- att[c("Survey","Year","Month","Day","HKey","ID" ,"X","Y",
             "LonDeep","LatDeep","LonShallow","LatShallow")]

## Quadrat attributes - Mean depth from quadrats aggregated to spatialised points
mean_att <- aggregate( . ~ ID, mean, data = spatialised[c("ID", "CorDepthM", "Slope", "PerCovZO", "PH", "ZO")], na.action = na.pass)
names(mean_att) <- c("ID", "mean_CorDepthM", "mean_slope", "mean_PerCovZO", "PH", "ZO")
#
##Quadrat attributes - most common substrate
Mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mode_att <- aggregate( . ~ ID, Mode, data = spatialised[c("ID", "Substrate")], na.action = na.pass)

## Add back attributes
att <- merge(att, mean_att, by="ID")
spat <- merge(att, mode_att, by="ID")
#spat <- merge(spat, nquads, by="ID")
spat <- spat[order(spat$HKey),]

## Ensure presence/absence
spat$PH[spat$PH > 0] <- 1
spat$ZO[spat$ZO > 0] <- 1




#JUST KEEP WHERE VALUES ARE PRESENCES????????


# Convert to spdf and export
spat.spdf <- spat %>%
  st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005")

# export as shapefile
# likely to have issues with attribute field names shortening
st_write(spat.spdf, "code/output_data/processed_observations/SpatializedQuadrats_aggregated_2024only.shp", append=FALSE)

save(spat, file="code/output_data/processed_observations/seagrass_data_spatialized_aggregated_2024only.RData")
