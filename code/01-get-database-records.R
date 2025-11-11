###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
#
# Overview:
# Processing scripts that query and standardize zostera, phyllospadix, percent cover, depth, and substrate observations from Pacific shellfish databases through SQL Server and BHM database through MS Access. 
# Data range is 1993-2023 to correspond with ROMs hindcast (1993-2021), NEP36, and Salish Sea Cast
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
              "Year","Month", "Day", "LatDeep","LonDeep",
              "LatShallow","LonShallow", "Transect_length", "Quadrat", "CorDepthM","Substrate", "Slope", "PerCovZO",
              "PH", "ZO", "Identification")  


#---------------------------------------------------------------------#
#### Get RSU_bio (red sea urchin) dive survey data #### 

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
rsu_dat <- DBI::dbGetQuery( sf_db_connection(), rsu_sql )

## Check for duplicates of species records 
dup_ind <- paste(rsu_dat$HKey, rsu_dat$Transect, rsu_dat$Quadrat, sep="-")
if ( any (duplicated(dup_ind)) ){
  inds <- c(which(duplicated(dup_ind, fromLast = FALSE)), 
            which(duplicated(dup_ind, fromLast = TRUE))) %>% sort()
  # show duplicated records
  print( rsu_dat[inds,] )
  # Remove duplicate records
  rsu_dat <- rsu_dat[!duplicated(rsu_dat),]
}

# Complicated regarding skip patterns, see data dictionary for details
rsu_dat <- rsu_dat %>%
  mutate(Quadrat_distance = case_when(Survey=="RES-H" | Survey == "RES-T" | Survey == "RBB" | Survey == "RBS" ~ 2, #sampled every other quadrat
                                      Survey=="RES-P" ~ NA)) #no metadata available

## Calculate slope using arc tangent method, make assumption that slope is same on Quadrat 1 as Quadrat 2 (since we don't have starting depth at bottom)
rsu_dat <- rsu_dat %>%
  mutate(StartDepth = ifelse(Quadrat==1 | Survey=="RES-P", NA, lag(CorDepthM, n=1)),
         EndDepth = CorDepthM,
         Elev.Diff = StartDepth - EndDepth,
         Slope = atan2(Elev.Diff, Quadrat_distance)) %>%
  mutate(Slope = ifelse(Quadrat==1, lead(Slope, n=1), Slope)) %>%
  select(-StartDepth,-EndDepth, -Elev.Diff)  

#change to degrees
rsu_dat$Slope <- round(abs(rsu_dat$Slope) * 180/pi, digits = 0)

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
rsu_dat <- rsu_dat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

## Melt all species into species column
# rsu_dat <- melt(rsu_dat, measure.vars=sp_pa, value.name = "SpNum", variable.name = "Species")

#get percent cover of ZO (not possible to get percent cover of PH as outher algae would be found alongside PH, for ZO algae not as likely to also be present, i.e. ZO usually in monoculture)
rsu_dat <- rsu_dat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, UnderstoryPct)) 

rsu_dat$PerCovZO[rsu_dat$PerCovZO > 100] <- 100

#remove Survey=="RES-P" as don't know their collection methods for algae or much other metadata and have not gotten first nation approval to share data
rsu_dat <-rsu_dat %>% filter(Survey!="RES-P")

## Create new fields and rename
# Type
rsu_dat$Type <- "Research"
# Source
rsu_dat$Source <- "RSU_bio"
# Method
rsu_dat$Method <- "Dive"

## remove where phyllo and zostera are likely mis identified based on substrate type
rsu_dat <- rsu_dat %>%
  mutate(Identification = case_when((ZO == 1 & Btype.description == "Bedrock dominant" &(Substrate3 != 7 & Substrate3 != 9) %>% replace_na(TRUE)) ~ "Remove",
                                    Substrate != "Rock" & Substrate != "Mixed" & PH == 1 & Btype.description != "Sand/shell, some rock" ~ "Remove",
                                    .default = "Keep"))

# Rename
rsu_dat <- rsu_dat[fnames]

# Order rows
rsu_dat <- rsu_dat[order(rsu_dat$HKey, rsu_dat$Year, rsu_dat$Transect, rsu_dat$Quadrat),]


#---------------------------------------------------------------------#
#### Get GSU_bio (green sea urchin) dive survey data ####
# P/A observations: PH, ZO 
# SQL code has filtered out years 2005 to present and Start Lat and Long are not null and CorDepth is not null
# All species observations are converted to presence (1)/ absence (0)

## Extract data from SQL Server
# Read SQL query
gsu_sql <- readLines("code/sql/get-gsu-records.sql")
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
rsc_sql <- readLines("code/sql/get-cuke-records.sql")
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
multi_sql <- readLines("code/sql/get-multispecies-records.sql")
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
ABLdat <- read.csv("raw_data/ab-data-restricted/New/AbaloneDensitySubstrate1978-2021.csv",header=T, fileEncoding="UTF-8-BOM", stringsAsFactors = F)
algaedat <- read.csv("raw_data/ab-data-restricted/New/AbaloneHabitat2002-2021.csv",header=T, fileEncoding="UTF-8-BOM", stringsAsFactors = F)
algaedat_earlier <- read.csv("raw_data/ab-data-restricted/New/AbaloneDensitySubstrateAlgae1995-2005.csv",header=T, fileEncoding="UTF-8-BOM", stringsAsFactors = F)

# Only one lat/long per site, remove sites with no lat/long
ABLdat <- ABLdat %>%
  mutate (Latitude = coalesce(LatDeep, LatShallow), Longitude = coalesce(LonDeep, LonShallow))%>%
  tidyr::drop_na (Latitude)

ABLdat<- ABLdat %>%
  filter(Year > 2002 & SurveyDesign < 4) #algae data not collected on quadrat level before this. Other survey designs are not clear

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
ABLdat <- ABLdat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#Assign substrate to each HKey based on most prevalent substrate at site, going to test not agregating quadrats. if doesnt work use this code
# ABLsubdat<- ABLdat %>%
#   filter(!is.na(Substrate)) %>%
#   group_by(HKey) %>%
#   summarize (Substrate = names(which.max(table(Substrate))))

ABLdat<- ABLdat %>% 
  # group_by(HKey) %>%
  # mutate(CorDepthM = mean(CorDepthM))%>%
  # ungroup() %>%
  # distinct(HKey, .keep_all = TRUE) %>%
  select(SurveyName, Year, Month, Day, HKey, Quadrat, LonDeep, LatDeep, LonShallow, LatShallow, CorDepthM, Latitude, Longitude, Substrate)%>%
  rename(Survey = SurveyName)


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
algaedat$PH[algaedat$CanopySpecies1 == "PK"] <- 1
algaedat$PH[algaedat$CanopySpecies2 == "PK"] <- 1
algaedat$PH[algaedat$UndStySpecies1 == "PK"] <- 1
algaedat$PH[algaedat$UndStySpecies2 == "PK"] <- 1
algaedat$PH[algaedat$UndStySpecies3 == "PK"] <- 1
algaedat$PH[algaedat$UndStySpecies4 == "PK"] <- 1
algaedat$PH[algaedat$UndStySpecies5 == "PK"] <- 1
algaedat$PH[algaedat$TurfSpecies1 == "PK"] <- 1
algaedat$PH[algaedat$TurfSpecies2 == "PK"] <- 1

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
algaedat$ZO[algaedat$CanopySpecies1 == "ZM"] <- 1
algaedat$ZO[algaedat$CanopySpecies2 == "ZM"] <- 1
algaedat$ZO[algaedat$UndStySpecies1 == "ZM"] <- 1
algaedat$ZO[algaedat$UndStySpecies2 == "ZM"] <- 1
algaedat$ZO[algaedat$UndStySpecies3 == "ZM"] <- 1
algaedat$ZO[algaedat$UndStySpecies4 == "ZM"] <- 1
algaedat$ZO[algaedat$UndStySpecies5 == "ZM"] <- 1
algaedat$ZO[algaedat$TurfSpecies1 == "ZM"] <- 1
algaedat$ZO[algaedat$TurfSpecies2 == "ZM"] <- 1

algaedat<- algaedat %>% 
  # group_by(HKey) %>% 
  # mutate(PH = sum(PH), ZO = sum(ZO))%>%
  # ungroup() %>%
  # distinct(HKey, .keep_all = TRUE) %>%
  select(HKey, QuadratNum, PH, ZO, UnderstoryPct)%>%
  rename(Quadrat = QuadratNum)

# earlier data
algaedat_earlier <- algaedat_earlier %>%
  select(-c(PT, ZO, NT, MA)) 

algaedat_earlier$PH <- 0
algaedat_earlier$PH[algaedat_earlier$Algae1 == "PH"] <- 1
algaedat_earlier$PH[algaedat_earlier$Algae2 == "PH"] <- 1
algaedat_earlier$PH[algaedat_earlier$Algae3 == "PH"] <- 1
algaedat_earlier$PH[algaedat_earlier$Algae1 == "PK"] <- 1
algaedat_earlier$PH[algaedat_earlier$Algae2 == "PK"] <- 1
algaedat_earlier$PH[algaedat_earlier$Algae3 == "PK"] <- 1

algaedat_earlier$ZO <- 0
algaedat_earlier$ZO[algaedat_earlier$Algae1 == "ZO"] <- 1
algaedat_earlier$ZO[algaedat_earlier$Algae2 == "ZO"] <- 1
algaedat_earlier$ZO[algaedat_earlier$Algae3 == "ZO"] <- 1
algaedat_earlier$ZO[algaedat_earlier$Algae1 == "ZM"] <- 1
algaedat_earlier$ZO[algaedat_earlier$Algae2 == "ZM"] <- 1
algaedat_earlier$ZO[algaedat_earlier$Algae3 == "ZM"] <- 1

algaedat_earlier<- algaedat_earlier %>% 
  # group_by(HKey) %>% 
  # mutate(PH = sum(PH), ZO = sum(ZO))%>%
  mutate(UnderstoryPct = NA) %>%
  # ungroup() %>%
  # distinct(HKey, .keep_all = TRUE) %>%
  select(HKey, Quadrat, PH, ZO, UnderstoryPct)

algaedat<- algaedat %>% 
  rbind(algaedat_earlier)

#change to 0 and 1
# algaedat$PH[algaedat$PH>0]<-1
# algaedat$ZO[algaedat$ZO>0]<-1


ABLdat <- merge(ABLdat, algaedat, by = c("HKey", "Quadrat"))
# ABLdat <- merge(ABLdat, ABLsubdat, by = "HKey")

#add percent cover of ZO
ABLdat <- ABLdat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, UnderstoryPct)) 

ABLdat$PerCovZO[ABLdat$PerCovZO > 100] <- 100

## Create new fields and rename
# Source
ABLdat$Source <- "ABL_bio"
#Slope is not possible to calculate from this survey, so will need to use modelled slope
ABLdat$Slope <- NA

ABLdat <- ABLdat %>% 
  mutate (HKey = paste0(Source,"_",HKey))

ABLdat$ID <- paste(ABLdat$HKey, "1", sep="_")

## remove where zostera are likely mis identified based on substrate type
ABLdat <- ABLdat %>%
  filter(!(ZO == 1 & Substrate == "Rock"))

ABLdat_sf <- ABLdat %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

ABL.spdf<- cbind(ABLdat, st_coordinates(st_as_sf(ABLdat_sf, coords = c("x", "y"), crs = "EPSG:3005")))

ABL.spdf<-ABL.spdf %>% select("Survey","Year","Month","Day","HKey","ID" , "X", "Y",
                              "LonDeep","LatDeep","LonShallow","LatShallow", "CorDepthM", "Slope", "Substrate", "PerCovZO", "PH", "ZO")


## Save files
save(ABL.spdf, file="code/output_data/abl_seagrass_data.RData")

# export as shapefile
# likely to have issues with attribute field names shortening
ABLsp_sf <- ABLdat %>%
  mutate(Latitude = case_when(!is.na(LatDeep) ~ LatDeep,
                              is.na(LatDeep) ~ LatShallow),
         Longitude = case_when(!is.na(LonDeep) ~ LonDeep,
                               is.na(LonDeep) ~ LonShallow)) %>%
  filter(!is.na(Latitude))  %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")


st_write(ABLsp_sf, "code/output_data/processed_observations/abl_transect.shp", append=FALSE)


#---------------------------------------------------------------------#
#### Get GDK_bio (geoduck) dive survey data ####
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
bhm_mdb <- mdb_connection("//ent.dfo-mpo.ca/dfo-mpo/GROUP/PAC/Reg_Shares/EOS/ES/Spatial_Datasets/Dive_Surveys/Database/BHM_DiveSurveys_CLEAN_July2024.mdb")

# Load queries
bhm_sql <- readLines("code/sql/get-bhm-records-quads.sql")
bhm_sql <- paste(bhm_sql, collapse = "\n ")
bhm_spp_sql <- readLines("code/sql/get-bhm-records-spp.sql")
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




#### Get Hilo (Gulf Islands High and low current) dive survey data ####

#load data
hilo_dat <- read.csv("raw_data/hilo/HiLo_quadrats.csv",header=T, stringsAsFactors = F)

# change lat long into decimal degrees
hilo_dat <- hilo_dat %>%
  mutate(LonShallow = (-1)*abs(Long_deg_s+(Long_min_s/60)),
         LatShallow = (Lat_deg_s+(Lat_min_s/60)),
         LonDeep = (-1)*abs(Long_deg_d+(Long_min_d/60)),
         LatDeep = (Lat_deg_d+(Lat_min_d/60)))

## Calculate the substrate that represents > 50% for each quad
# Match substrateID to substrate category
hilo_dat <- hilo_dat %>% 
  dplyr::left_join(sub.cat, by=c("Substrate1", "Substrate2")) %>%
  rename (Substrate = RMSM.Nme) 

#add percent cover of ZO, this survey does not record percent cover so all absences are assigned 0 and presences assigned NA
hilo_dat <- hilo_dat %>%
  mutate(PerCovZO = ifelse(ZO==0, 0, NA)) 

## Melt all invert and algae species into species column
#hilo_dat <- melt(hilo_dat, measure.vars=sp_pa, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
hilo_dat$Type <- "Research"
# Source
hilo_dat$Source <- "MSEA"
# Method
hilo_dat$Method <- "Dive"
# Slope is not possible to determine from the method of this survey because don't know distance between quadrats
hilo_dat$Slope <- NA
# transect length is not possible to determine from the method of this survey 
hilo_dat$Transect_length <- NA


## remove where phyllo and zostera are likely mis identified based on substrate type
hilo_dat <- hilo_dat %>%
  mutate(Identification = case_when(ZO == 1 & Btype.description == "Bedrock dominant"  ~ "Remove",
                                    Substrate != "Rock" & Substrate != "Mixed" & PH == 1 ~ "Remove",
                                    .default = "Keep"))

# Rename
hilo_dat <- hilo_dat[fnames]
# Order rows
hilo_dat <- hilo_dat[order(hilo_dat$HKey, hilo_dat$Transect, hilo_dat$Quadrat),]

####Combine datasets ####
dat <- bind_rows(rsu_dat, gsu_dat, rsc_dat, multi_dat, gdk_dat, bhm_dat, hilo_dat)

#need to make new HKey incase there is overlap in HKey between source types
dat <- dat %>% 
  mutate (HKey = paste0(Source,"_",HKey))

# when ZO has value of 1 and percent cover has value of 0, need to change to NA
dat$PerCovZO[dat$PerCovZO == 0 & dat$ZO == 1] <- NA

# remove ZO presence observations in high intertidal, likely Z. japonica
dat$Identification[dat$ZO == 1 & dat$CorDepthM < -2] <- "Remove"


## Save files
save(dat, file="code/output_data/processed_observations/seagrass_data.RData")





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
st_write(dat_sf, "code/output_data/processed_observations/dat_transect.shp", append=FALSE)



