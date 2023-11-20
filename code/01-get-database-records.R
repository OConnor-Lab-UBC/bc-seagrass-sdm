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
outdir <- file.path("code", "output_data")
# Create the main directory, and subfolders
dir.create( outdir, recursive=TRUE )
# final field names for dive surveys
fnames <-  c( "Type", "Source", "Survey", "HKey", "Method",
              "Year","Month", "Major_stat_area_code", "Stat_subarea_code", "LatDeep","LonDeep",
              "LatShallow","LonShallow", "Transect_length", "Quadrat", "CorDepthM","Substrate", "Slope",
              "SpNum", "Species")  


#---------------------------------------------------------------------#
#### Get RSU_bio (red sea urchin) presence/absence survey data #### 

# P/A observations: PH (phyllospadix), ZO (eelgrass) 
# SQL code has filtered out years 2000 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)
# Most data pre 2005 only has a starting lat/long and not an ending. May not working to use this earlier data

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

rsu_queried %>% 
  group_by(Survey) %>%
  summarise(count= n_distinct(Location))

# Complicated regarding skip patterns, see data dictionary for details
rsu_queried <- rsu_queried %>%
  mutate(Quadrat_distance = case_when(Survey=="RES-H" | Survey == "RES-T" | Survey == "RBB" | Survey == "RBS" ~ 2, #sampled every other quadrat
                                      Survey=="RES-P" ~ NA)) #no metadata available

## Calculate slope using arc tangent method, make assumption that slope is same on Quadrat 1 as Quadrat 2 (since we don't have starting depth at bottom)
rsu_queried <- rsu_queried %>%
  mutate(StartDepth = ifelse(Quadrat==1 | Survey=="RES-P", NA, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff, Quadrat_distance)) %>%
  mutate(Slope = ifelse(Quadrat==1, lead(Slope, n=1), Slope)) #%>%
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
rsu_queried$Slope2 <- round(abs(rsu_queried$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
rsu_queried <- rsu_queried %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

## Melt all species into species column
# list p/a species. 
sp_pa_rsu <- c("PH", "ZO")
# melt table for p/a species
rsu_dat <- melt(rsu_queried, measure.vars=sp_pa_rsu, value.name = "SpNum", variable.name = "Species")

#remove Survey=="RES-P" as don't know their collection methods for algae or much other metadata
rsu_dat <-rsu_dat %>%
  filter(Survey!="RES-P")

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
save(rsu_dat, file="code/output_data/rsu_db_data.RData")
write.csv(rsu_dat, "code/output_data/rsu_db_quadrat.csv", row.names = F)

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

## Calculate slope using arc tangent method, make assumption that slope is same on Quadrat 1 as Quadrat 2 (since we don't have starting depth at bottom)
gsu_queried <- gsu_queried %>%
  mutate(StartDepth = ifelse(Quadrat==1, NA, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff, Quadrat_distance)) %>%
  mutate(Slope = ifelse(Quadrat==1, lead(Slope, n=1), Slope)) %>%
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
gsu_queried$Slope <- round(abs(gsu_queried$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
gsu_queried <- gsu_queried %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

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
save(gsu_dat, file="code/output_data/gsu_db_data.RData")
write.csv(gsu_dat, "code/output_data/gsu_db_quadrat.csv", row.names = F)


#---------------------------------------------------------------------#
#### Get Cuke_bio (California sea cucumber) survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2000 to present and Start Lat and Long are not null and CorDepth is not null
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

#Years >2016 have Quad0 so know bottom depths, surveys prior we have no start depth of first quadrat of each transect but making assumption it has similar slope as quadrat after it. 
rsc_queried <- rsc_queried %>%
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
rsc_queried$Slope <- round(abs(rsc_queried$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
rsc_queried <- rsc_queried %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

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
save(rsc_dat, file="code/output_data/rsc_db_data.RData")
write.csv(rsc_dat, "code/output_data/rsc_db_quadrat.csv", row.names = F)

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


## Calculate slope using arc tangent method, make assumption that slope is same on Quadrat 1 as Quadrat 2 (since we don't have starting depth at bottom)
multi_queried <- multi_queried %>%
  mutate(StartDepth = ifelse(Quadrat==1, NA, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff,Quadrat_distance)) %>%
  mutate(Slope = ifelse(Quadrat==1, lead(Slope, n=1), Slope)) %>%
  select(-StartDepth,-EndDepth, -Elev.Diff)
  
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
save(multi_dat, file="code/output_data/multispecies_db_data.RData")
write.csv(multi_dat, "code/output_data/multispecies_db_quadrat.csv", row.names = F)



#---------------------------------------------------------------------#
#### Get ABL_bio (abalone) dive survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out Start Lat and Long are not null 
# All species observations are converted to presence (1)/ absence (0)
# Zostera would be recorded if observed on Ab surveys but Ab surveys are typically in the wrong habitat for Zostera.
# Phyllospadix when observed on Ab surveys are recorded. It is in shallow exposed areas so divers may not be able to work as shallow as the Phyllospadix band depending on conditions. 
# This dataset has different spatialization requirements so is kept separate, Transects follows zig zag pattern and transect is assumed to be about a 20x20m area
# Looking at 95% CI of the depths at quadrats within a transect most quadrats are within 2m (2 SE) of the mean, and a few are within 3m. This makes it more reasonable to amalgamate all quadrats from one transect
# cannot get slope from this survey as transect is zig zag

# Load data
ABLdat <- read.csv("raw_data/ab-data-restricted/Shellfish_Bio_Abalone_HeadersDensity_Substrate_2005-2021.csv",header=T, fileEncoding="UTF-8-BOM", stringsAsFactors = F)
algaedat <- read.csv("raw_data/ab-data-restricted/Shellfish_Bio_Abalone_Habitat_2005-2021.csv",header=T, fileEncoding="UTF-8-BOM", stringsAsFactors = F)

# Only one lat/long per site, remove sites with no lat/long
ABLdat <- ABLdat %>%
  mutate (Latitude = coalesce(LatDeep, LatShallow), Longitude = coalesce(LonDeep, LonShallow))%>%
  tidyr::drop_na (Latitude)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
ABLdat <- ABLdat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#Assign substrate to each quadrat based on most prevalent substrate at site
ABLsubdat<- ABLdat %>% 
  group_by(HKey) %>%
  summarize (Substrate = names(which.max(table(Substrate))))

ABLdat<- ABLdat %>% 
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

ABLdat <- merge(ABLdat, algaedat, by = "HKey")
ABLdat <- merge(ABLdat, ABLsubdat, by = "HKey")

## Create new fields and rename
# Type
ABLdat$Type <- "Research"
# Source
ABLdat$Source <- "ABL_bio"
# Method
ABLdat$Method <- "Dive"
#Slope is not possible to calculate from this survey, so will need to use modelled slope
ABLdat$Slope <- "NA"

# Convert to spdf and export
# create spatial points
ABLdat_sf <- ABLdat %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

# export as shapefile
# likely to have issues with attribute field names shortening
st_write(ABLdat_sf, "code/output_data/ABLBio_quadrat.shp") 


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

#Most transects have a Quad0 so know bottom depths, but some don't so need to find those
gdk_nozeroquad <- gdk_queried %>%
  group_by(HKey) %>%
  filter(!all(0 %in% Quadrat)) %>%
  ungroup %>%
  distinct(HKey) %>%
  pull(HKey)

#errors in database of transect_dist_from_start
gdk_queried[gdk_queried$HKey=="10491" & gdk_queried$Quadrat== 15, "Transect_dist_from_start"] <- 300
gdk_queried[gdk_queried$HKey=="10493" & gdk_queried$Quadrat== 15, "Transect_dist_from_start"] <- 215
gdk_queried[gdk_queried$HKey=="10541" & gdk_queried$Quadrat== 12, "Transect_dist_from_start"] <- 170
gdk_queried[gdk_queried$HKey=="10545" & gdk_queried$Quadrat== 11, "Transect_dist_from_start"] <- 205
gdk_queried[gdk_queried$HKey=="10545" & gdk_queried$Quadrat== 12, "Transect_dist_from_start"] <- 225
gdk_queried[gdk_queried$HKey=="10545" & gdk_queried$Quadrat== 13, "Transect_dist_from_start"] <- 245
gdk_queried[gdk_queried$HKey=="10545" & gdk_queried$Quadrat== 14, "Transect_dist_from_start"] <- 265
gdk_queried[gdk_queried$HKey=="10545" & gdk_queried$Quadrat== 15, "Transect_dist_from_start"] <- 285
gdk_queried[gdk_queried$HKey=="10545" & gdk_queried$Quadrat== 16, "Transect_dist_from_start"] <- 305

gdk_queried <- gdk_queried %>%
  filter(!HKey %in% 13094:13123) %>% # these only have two quadrats in each transect
  mutate(StartDepth = ifelse(Quadrat==0 | (Quadrat==1 & HKey %in% gdk_nozeroquad), CorDepthM, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Quadrat_distance = ifelse(Quadrat==0 | (Quadrat==1 & HKey %in% gdk_nozeroquad), Transect_dist_from_start, Transect_dist_from_start - lag(Transect_dist_from_start, n = 1)),
         Slope = atan2(Elev.Diff, Quadrat_distance)) %>%
  filter(Quadrat!=0) %>% # Need to remove quad 0 from each transect not all transects have Quad zero and it only has depth recorded
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
gdk_queried$Slope <- round(abs(gdk_queried$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
gdk_queried <- gdk_queried %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

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
save(gdk_dat, file="code/output_data/gdk_db_data.RData")
write.csv(gdk_dat, "code/output_data/gdk_db_quadrat.csv", row.names = F)
