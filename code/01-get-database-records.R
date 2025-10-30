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
source("r-sql-link-functions.R")

#### Read in substrate category table ####
sub.cat <- read.csv( "./lookup-tbls/SubstrateCategories.csv", header=T, sep="," )

#### Outputs ####
# Subfolders for outputs
outdir <- file.path("output_data")
# Create the main directory, and subfolders
dir.create( outdir, recursive=TRUE )

# final field names for dive surveys
fnames <-  c( "Type", "Source", "Survey", "HKey", "Method",
              "Year","Month", "Day", "LatDeep","LonDeep",
              "LatShallow","LonShallow")  


#### Get BHM (Benthic Habitat Mapping) dive survey data ####
# P/A observations: PH (phyllospadix), ZO (eelgrass) 
# All species observations are converted to presence (1)/ absence (0)

# Connect to mdb
bhm_mdb <- mdb_connection("//ent.dfo-mpo.ca/dfo-mpo/GROUP/PAC/Reg_Shares/EOS/ES/Spatial_Datasets/Dive_Surveys/Database/BHM_DiveSurveys_Jan2025.mdb")

# Load queries
bhm_sql <- readLines("sql/get-bhm-records-quads.sql")
bhm_sql <- paste(bhm_sql, collapse = "\n ")
bhm_spp_sql <- readLines("sql/get-bhm-records-spp.sql")
bhm_spp_sql <- paste(bhm_spp_sql, collapse = "\n ")

# Run queries
bhm_quads <- DBI::dbGetQuery( bhm_mdb, bhm_sql )
bhm_spp <- DBI::dbGetQuery( bhm_mdb, bhm_spp_sql )

bhm_quads <- bhm_quads %>% distinct(HKey, .keep_all = TRUE)

bhm_quads <- bhm_quads %>%
  mutate(
    lat = if_else(!is.na(LatShallow), LatShallow, LatDeep),
    long = if_else(!is.na(LonShallow), LonShallow, LonDeep)
  )

bhm_quads<- bhm_quads %>% select (c(Survey, HKey, Year, Month, Day, lat, long))

library(sf)

# Convert to sf object
bhm_quads_sf <- st_as_sf(
  bhm_quads,
  coords = c("long", "lat"),   # column names for coordinates
  crs = 4326                   # WGS84 coordinate reference system
)

# Save to shapefile
st_write(bhm_quads_sf, "bhm_sites.shp", delete_dsn = TRUE)



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
save(dat, file="code/output_data/seagrass_data.RData")





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
st_write(dat_sf, "code/output_data/dat_transect.shp", append=FALSE)



