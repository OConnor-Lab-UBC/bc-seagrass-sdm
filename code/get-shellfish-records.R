###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC seagrass SDM
#
# Overview:
# Processing scripts that query and standardize zostera, phyllospadix, depth, and substrate observations from Pacific shellfish databases through SQL Server

# Data use:
# As per the Policy for Scientific Data, DFO scientific data are a public resource and subject to full and open access within two years of being acquired or generated. 
# As per the Application for Access to Shellfish Data Holdings, data collected two years prior to today should not be included in data pulls to allow the Principle Investigator a two year period to analyze and report findings. 
# Exceptions to the use of data collected in the last two years may be made if the user has met the conditions outlined in the Application for Access or has discussed the use of the data with the Principle Investigator for the survey

# Requirements:
# r-sql-link-functions.R 
# Access to shellfish SQL server (VPN connection and access permitted from DFO Shellfish Data Unit)

###############################################################################

#### load packages####
library(DBI)
library(odbc)
library(reshape2)
library(sf)
library(tidyverse)

#### Load util.R gfdata functions ####
source("code/r-sql-link-functions.R")

#### Read in substrate category table ####
sub.cat <- read.csv( "./lookup-tbls/SubstrateCategories.csv", header=T, sep="," )

#### Outputs ####
# Subfolders for outputs
outdir <- file.path("DataSources", "Data")
# Create the main directory, and subfolders
dir.create( outdir, recursive=TRUE )
# final field names for dive surveys
fnames <-  c( "Type", "Source", "Survey", "HKey", "Method",
              "Year","Month", "Major_stat_area_code", "Stat_subarea_code", "LatDeep","LonDeep",
              "LatShallow","LonShallow", "Transect_length", "Quadrat", "CorDepthM","Substrate", "Slope",
              "SpNum", "Species")  


#---------------------------------------------------------------------#
#### Get RSU_bio (red sea urchin) presence/absence survey data #### 

# P/A observations: PH (phyloospadix), ZO (eelgrass) 
# SQL code has filtered out years 2005 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)

##Extract data from SQL Server
# read SQL query
rsu_sql <- readLines("code/sql/get-rsu-records.sql")
# Load query
rsu_sql <- paste(rsu_sql, collapse = "\n ")
# Run query
rsu_queried <- DBI::dbGetQuery( sf_db_connection(), rsu_sql )

## Check for duplicates of species records 
dup_ind <- paste(rsu_queried$HKey, rsu_queried$Transect, rsu_queried$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( rsu_queried[inds,] )
  # Remove duplicate records
  rsu_queried <- rsu_queried[!duplicated(rsu_queried),]
}


slope <- rsu_queried %>%
  select(HKey, Quadrat, CorDepthM) 

# Set start depth as the recorded depth from the previous quadrat, unless the quadrat is number 0
slope$StartDepth <- ifelse(slope$Quadrat==1, slope$CorDepthM, lag(slope$CorDepthM, n=1))
slope$EndDepth <- slope$CorDepthM

# Calculate slope for each quadrat using the arc-tangent
slope$Elev.Diff <- slope$StartDepth-slope$EndDepth
slope$Slope <- atan2(slope$Elev.Diff,5) 
slope$Slope <- round(slope$Slope, digits=4)

#For records with slope greater than 1.0 or less than -1.0 set slope to either 1.0 or -1.0 to correct for typo or instances where the swell likely exacerbated the difference between two height recordings 
slope$Slope[slope$Slope > 1.0] <- 1.0
slope$Slope[slope$Slope < -1.0] <- -1.0 

BHM_data$Slope <- abs(BHM_data$Slope) * 180/piBHM_data$Slope <- abs(BHM_data$Slope) * 180/pi




## Melt all species into species column
# list p/a species. 
sp_pa_rsu <- c("PH", "ZO")
# melt table for p/a species
rsu_dat <- melt(rsu_queried, measure.vars=sp_pa_rsu, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
rsu_dat$Type <- "Research"
# Source
rsu_dat$Source <- "RSU_bio"
# Method
rsu_dat$Method <- "Dive"
# Rename
rsu_dat <- rsu_dat[fnames]

# Order rows
rsu_dat <- rsu_dat[order(rsu_dat$HKey, rsu_dat$Year, rsu_dat$Transect, rsu_dat$Quadrat),]

## Save files
save(rsu_dat, file="code/data/rsu_db_data.RData")
write.csv(rsu_dat, "code/data/rsu_db_quadrat.csv", row.names = F)

#---------------------------------------------------------------------#
#### Get GSU_bio (green sea urchin) presence/absence survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2005 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)

## Extract data from SQL Server
# Read SQL query
gsu_sql <- readLines("code/sql/get-gsu-records.sql")
# Load query
gsu_sql <- paste(gsu_sql, collapse = "\n ")
# Run query
gsu_queried <- DBI::dbGetQuery( sf_db_connection(), gsu_sql )

## Check for duplicates of species records from a single finishing event
dup_ind <- paste(gsu_queried$HKey, gsu_queried$Transect, gsu_queried$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( gsu_queried[inds,] )
  # Remove duplicate records
  gsu_queried <- gsu_queried[!duplicated(gsu_queried),]
}

## Melt all invert and algae species into species column
# List p/a species. 
sp_pa_gsu <- c("PH", "ZO")
# melt table for p/a species
gsu_dat <- melt(gsu_queried, measure.vars=sp_pa_gsu, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
gsu_dat$Type <- "Research"
# Source
gsu_dat$Source <- "GSU_bio"
# GSU doesn't have a consistently filled transect field in database
gsu_dat$Transect <- ""
# Method
gsu_dat$Method <- "Dive"
# Rename
gsu_dat <- gsu_dat[fnames]

# Order rows
gsu_dat <- gsu_dat[order(gsu_dat$HKey, gsu_dat$Year, gsu_dat$Transect, gsu_dat$Quadrat),]

## Save files
save(gsu_dat, file="code/data/gsu_db_data.RData")
write.csv(gsu_dat, "code/data/gsu_db_quadrat.csv", row.names = F)


#---------------------------------------------------------------------#
#### Get Cuke_bio (California sea cucumber) survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2005 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)

## Extract data from SQL Server
# Read SQL query
rsc_sql <- readLines("code/sql/get-cuke-records.sql")
# Load query
rsc_sql <- paste(rsc_sql, collapse = "\n ")
# Run query
rsc_queried <- DBI::dbGetQuery( sf_db_connection(), rsc_sql )

## Check for duplicates of species records from a single finishing event
dup_ind <- paste(rsc_queried$HKey, rsc_queried$Transect, rsc_queried$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( rsc_queried[inds,] )
  # Remove duplicate records
  rsc_queried <- rsc_queried[!duplicated(rsc_queried),]
}

## Melt all invert and algae species into species column
# list p/a species
sp_pa_rsc <- c("PH", "ZO")
# melt table for p/a species
rsc_dat <- melt(rsc_queried, measure.vars=sp_pa_rsc, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
rsc_dat$Type <- "Research"
# Source
rsc_dat$Source <- "Cuke_bio"
# Method
rsc_dat$Method <- "Dive"
# Rename
rsc_dat <- rsc_dat[fnames]

# Order rows
rsc_dat <- rsc_dat[order(rsc_dat$HKey, rsc_dat$Year, rsc_dat$Transect, rsc_dat$Quadrat),]

## Save files
save(rsc_dat, file="code/data/rsc_db_data.RData")
write.csv(rsc_dat, "code/data/rsc_db_quadrat.csv", row.names = F)

#---------------------------------------------------------------------#
#### Get Multispecies_bio (multi-species) dive survey data ####
# P/A observations: PH, ZO 
# These DFO surveys began in 2016
# SQL code has filtered out Start Lat and Long are not null and CorDepth is not null
# Transect length is rounded up to the nearest 25m (there is no record of transect length in the database, so quadrat skipping info was used)
# All species observations are converted to presence (1)/ absence (0)

## Extract data from SQL Server
# Read SQL query
multi_sql <- readLines("code/sql/get-multispecies-records.sql")
# Load query
multi_sql <- paste(multi_sql, collapse = "\n ")
# Run query
multi_queried <- DBI::dbGetQuery( sf_db_connection(), multi_sql )

## Check for duplicates of species records from a single sampling event
dup_ind <- paste(multi_queried$HKey, multi_queried$Transect, multi_queried$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( multi_queried[inds,] )
  # Remove duplicate records
  multi_queried <- multi_queried[!duplicated(multi_queried),]
}


## Calculate slope using arc tangent method
multi_queried <- multi_queried %>%
  mutate(StartDepth = ifelse(Quadrat==1, CorDepthM, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff,Quadrat_distance)) %>%
  select(-StartDepth,-EndDepth, -Elev.Diff)
  
#For records with slope greater than 1.0 or less than -1.0 set slope to either 1.0 or -1.0 to correct for typo or instances where the swell likely exacerbated the difference between two height recordings 
multi_queried$Slope[multi_queried$Slope > 1.0] <- 1.0
multi_queried$Slope[multi_queried$Slope < -1.0] <- -1.0 

#change to degrees
multi_queried$Slope <- round(abs(multi_queried$Slope) * 180/pi, digits = 0)


## Calculate the substrate that represents > 50% for each quad

# Match substrateID to substrate category
multi_queried <- multi_queried %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 
  

## Melt all invert and algae species into species column
# list p/a species
sp_pa_multi <- c("PH", "ZO")
# melt table for p/a species
multi_dat <- melt(multi_queried, measure.vars=sp_pa_multi, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
multi_dat$Type <- "Research"
# Source
multi_dat$Source <- "Multispecies_bio"
# No time type in multispecies db
multi_dat$Time_type <- ""
# Method
multi_dat$Method <- "Dive"

# Rename
multi_dat <- multi_dat[fnames]

# Order rows
multi_dat <- multi_dat[order(multi_dat$HKey, multi_dat$Year, multi_dat$Transect, multi_dat$Quadrat),]

## Save files
save(multi_dat, file="code/data/multispecies_db_data.RData")
write.csv(multi_dat, "code/data/multispecies_db_quadrat.csv", row.names = F)



#---------------------------------------------------------------------#
#### Get ABL_bio dive survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out Start Lat and Long are not null 
# All species observations are converted to presence (1)/ absence (0)
# Zostera would be recorded if observed on Ab surveys but Ab surveys are typically in the wrong habitat for Zostera.
# Phyllospadix when observed on Ab surveys are recorded. It is in shallow exposed areas so divers may not be able to work as shallow as the Phyllospadix band depending on conditions. 
# This dataset has different spatialization requirements so is kept separate


# Load data
dat <- read.csv("code/ab-data-restricted/Shellfish_Bio_Abalone_HeadersDensity_2005-2021.csv",header=T, fileEncoding="UTF-8-BOM", stringsAsFactors = F)
algaedat <- read.csv("code/ab-data-restricted/Shellfish_Bio_Abalone_Habitat_2005-2021.csv",header=T, fileEncoding="UTF-8-BOM", stringsAsFactors = F)

# Only one lat/long per site, remove sites with no la/long
dat <- dat %>%
  mutate (Latitude = coalesce(LatDeep, LatShallow), Longitude = coalesce(LonDeep, LonShallow))%>%
  tidyr::drop_na (Latitude)

dat<- dat %>% 
  group_by(HKey) %>% 
  mutate(CorDepthM = mean(CorDepthM))%>%
  ungroup() %>%
  distinct(HKey, .keep_all = TRUE) %>%
  select(HKey, SurveyName, Year, Month, Day, Time, Time_type, Major_stat_area_code, Stat_subarea_code, Latitude, Longitude, TransectNum, Transect_length, CorDepthM)%>%
  rename(Survey = SurveyName, Transect = TransectNum, Bathymetry = CorDepthM)

# add PH column to algae data
algaedat$PH <- 0
algaedat$PH[algaedat$CanopySpecies1 == "PH"] <- 1
algaedat$PH[algaedat$CanopySpecies2 == "PH"] <- 1
algaedat$PH[algaedat$UndStySpecies1 == "PH"] <- 1
algaedat$PH[algaedat$UndStySpecies2 == "PH"] <- 1
algaedat$PH[algaedat$UndStySpecies3 == "PH"] <- 1
algaedat$PH[algaedat$UndStySpecies4 == "PH"] <- 1
algaedat$PH[algaedat$UndStySpecies5 == "PH"] <- 1
algaedat$PH[algaedat$TurfSpecies1 == "PH"] <- 1
algaedat$PH[algaedat$TurfSpecies2 == "PH"] <- 1

algaedat$ZO <- 0
algaedat$ZO[algaedat$CanopySpecies1 == "ZO"] <- 1
algaedat$ZO[algaedat$CanopySpecies2 == "ZO"] <- 1
algaedat$ZO[algaedat$UndStySpecies1 == "ZO"] <- 1
algaedat$ZO[algaedat$UndStySpecies2 == "ZO"] <- 1
algaedat$ZO[algaedat$UndStySpecies3 == "ZO"] <- 1
algaedat$ZO[algaedat$UndStySpecies4 == "ZO"] <- 1
algaedat$ZO[algaedat$UndStySpecies5 == "ZO"] <- 1
algaedat$ZO[algaedat$TurfSpecies1 == "ZO"] <- 1
algaedat$ZO[algaedat$TurfSpecies2 == "ZO"] <- 1

algaedat<- algaedat %>% 
  group_by(HKey) %>% 
  mutate(PH = sum(PH), ZO = sum(ZO))%>%
  ungroup() %>%
  distinct(HKey, .keep_all = TRUE) %>%
  select(HKey, PH, ZO)

#change to 0 and 1
algaedat$PH[algaedat$PH>0]<-1
algaedat$ZO[algaedat$ZO>0]<-1

dat <- merge(dat, algaedat, by = "HKey")

## Create new fields and rename
# Type
dat$Type <- "Research"
# Source
dat$Source <- "ABL_bio"
# Method
dat$Method <- "Dive"

# Convert to spdf and export
# create spatial points
dat_sf <- dat %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

# export as shapefile
# likely to have issues with attribute field names shortening
st_write(dat_sf, "code/data/ABLBio_quadrat.shp") 


#---------------------------------------------------------------------#
#### Get GDK_bio (geoduck) survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2005 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)
# This dataset has different spatialization requirements so is kept separate

## Extract data from SQL Server
# Read SQL query
gdk_sql <- readLines("code/sql/get-geoduck-records.sql")
# Load query
gdk_sql <- paste(gdk_sql, collapse = "\n ")
# Run query
gdk_queried <- DBI::dbGetQuery( sf_db_connection(), gdk_sql )

## Check for duplicates of species records from a single finishing event
dup_ind <- paste(gdk_queried$HKey, gdk_queried$Transect, gdk_queried$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( gdk_queried[inds,] )
  # Remove duplicate records
  gdk_queried <- gdk_queried[!duplicated(gdk_queried),]
}

## Melt all invert and algae species into species column
# list p/a species
sp_pa_gdk <- c("PH", "ZO")
# melt table for p/a species
gdk_dat <- melt(gdk_queried, measure.vars=sp_pa_gdk, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
gdk_dat$Type <- "Research"
# Source
gdk_dat$Source <- "GDK_bio"
# Method
gdk_dat$Method <- "Dive"

# Rename
gdk_dat <- gdk_dat[fnames]

# Order rows
gdk_dat <- gdk_dat[order(gdk_dat$HKey, gdk_dat$Year, gdk_dat$Transect, gdk_dat$Quadrat),]

## Save files
save(gdk_dat, file="code/data/gdk_db_data.RData")
write.csv(gdk_dat, "code/data/gdk_db_quadrat.csv", row.names = F)
