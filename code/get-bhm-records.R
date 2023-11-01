###############################################################################
#
# Authors:      Ashley Park & Jessica Nephin
# Affiliation:  Fisheries and Oceans Canada (DFO)
# Group:        Marine Spatial Ecology and Analysis
# Location:     Institute of Ocean Sciences
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      Query DFO Databases
#
# Overview:
# Processing scripts that query invertebrate and algae/marine plant observations from MSEA Benthic Habitat Mapping database and standardizes it to be processed through the spatialize scripts (https://gitlab.com/dfo-msea/spatialise-dive-quadrats)

# Data use:
# As per the Policy for Scientific Data, DFO scientific data are a public resource and subject to full and open access within two years of being acquired or generated. 
# As per the Application for Access to Shellfish Data Holdings, data collected two years prior to today should not be included in data pulls to allow the Principle Investigator a two year period to analyze and report findings. 
# Exceptions to the use of data collected in the last two years may be made if the user has met the conditions outlined in the Application for Access or has discussed the use of the data with the Principle Investigator for the survey

# Requirements:
# r-sql-link-functions.R 
# Access to benthic habitat mapping access database on Spatial Datasets drive (VPN connection and access permitted from MSEA section)

###############################################################################

# Update: Access query Now works with 64 bit   Warning MS Access (.mdb) database queries only work when running 32bit R
# Check which version of R is being used and reset if necessary
#Sys.getenv("R_ARCH")   # "/i386" 32 bit R
# To reset: Select Tools menu | Global Options... | R Version:

#### load packages####
library(DBI)
library(dplyr)
library(tidyverse)
library(reshape2)

#### Load util.R gfdata functions ####
source("DataSources/r-sql-link-functions.R")

#### Outputs ####
# Subfolders for outputs
outdir <- file.path("DataSources", "Data")
# Create the main directory, and subfolders
dir.create( outdir, recursive=TRUE )

# Connect to mdb
bhm_mdb <- mdb_connection("//dcbcpbsna01a/Spatial_Datasets/Dive_Surveys/Database/BHM_DiveSurveys_CURRENT_March2022.mdb")

# Load query
bhm_sql <- readLines("DataSources/sql/get-bhm-records.sql")
bhm_sql <- paste(bhm_sql, collapse = "\n ")
# Run query
bhm_dat <- DBI::dbGetQuery( bhm_mdb, bhm_sql )

## Edit Observations table ###
bhm_dat$SpeciesCode <- as.character(bhm_dat$Species)
bhm_dat$SpType <- toupper(bhm_dat$SpType) # capitalize


# These species were either only recorded in 2013 or recorded as a species group in later surveys
# Comments explain what observations are being recoded to and counts of dive sites with presence
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="AB" & bhm_dat$SpType=="I"] <- "OB"  # Acorn barnacle = Barnacle - other --- chk yrs
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="BB" & bhm_dat$SpType=="I"] <- "BRF" # Blue encrusting bryozoan(2) = Bryozoan flat
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="BS" & bhm_dat$SpType=="I"] <- "SPE" # Boot sponge = Sponge erect
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="TS" & bhm_dat$SpType=="I"] <- "SPE" # Trumpet sponge = Sponge erect
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="GS" & bhm_dat$SpType=="I"] <- "SPE" # Glove sponge = Sponge erect
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="HS" & bhm_dat$SpType=="I"] <- "SH"  # Humpback shrimp(4) = Shrimp
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="NS" & bhm_dat$SpType=="I"] <- "SH"  # Coonstripe shrimp(13) = Shrimp
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="SN" & bhm_dat$SpType=="I"] <- "SH"  # Spot prawn(2) = Shrimp
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="SW" & bhm_dat$SpType=="I"] <- "SP"  # Sea whip(1) = sea pen
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="FT" & bhm_dat$SpType=="I"] <- "SL"  # Chelyosoma productum(1) = tunicate - solitary - other

# Non-target species
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="LP" & bhm_dat$SpType=="I"] <- "nontarget"  # Six-armed star(6) = " "
#bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="BC" & bhm_dat$SpType=="I"] <- "nontarget"  # Butter clam (Saxidomus gigantea)(5) = ?
#bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="CO" & bhm_dat$SpType=="I"] <- "nontarget"  # Cockle (Clinocardium nuttallii)(7) = ?
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="XH" & bhm_dat$SpType=="I"] <- "nontarget"  # Stevens' hermit crab(5) = ?
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="KL" & bhm_dat$SpType=="I"] <- "nontarget"  # Limpets(4) - removed
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="TW" & bhm_dat$SpType=="I"] <- "nontarget"  # Calcareous tubeworm(36) - removed b/c ubiquitous
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="DA" & bhm_dat$SpType=="I"] <- "nontarget"  # Nudibranch (Dirona albolineata) - removed
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="OD" & bhm_dat$SpType=="I"] <- "nontarget"  # Nudibranch (Dirona pellucita) - removed
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="CN" & bhm_dat$SpType=="I"] <- "nontarget"  # Nudibranch (Hermissenda crassicornis) - removed
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="RT" & bhm_dat$SpType=="I"] <- "nontarget"  # Brown turban --- chk yrs
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="PC" & bhm_dat$SpType=="I"] <- "nontarget"  # Psolus --- chk yrs
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="SY" & bhm_dat$SpType=="I"] <- "nontarget"  # Spiny red star --- chk yrs
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="WU" & bhm_dat$SpType=="I"] <- "nontarget"  # White urchin --- chk yrs
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="AG" & bhm_dat$SpType=="A"] <- "nontarget"  # Agarum # OR is this a typo --- Red turban snail?
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="AL" & bhm_dat$SpType=="A"] <- "nontarget"  # Alaria 
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="BO" & bhm_dat$SpType=="A"] <- "nontarget"  # Bossiella 
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="DE" & bhm_dat$SpType=="A"] <- "nontarget"  # Desmarestia 
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="EN" & bhm_dat$SpType=="A"] <- "nontarget"  # Encrusting 
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="LA" & bhm_dat$SpType=="A"] <- "nontarget"  # Laminaria 
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="OP" & bhm_dat$SpType=="A"] <- "nontarget"  # Opuntella 
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="BP" & bhm_dat$SpType=="A"] <- "nontarget"  # Botryocladia pseudodichotoma
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="DG" & bhm_dat$SpType=="A"] <- "nontarget"  # Derbesia marina

# Species codes that are likely meant to be algae not inverts and vise versa
bhm_dat$SpType[bhm_dat$SpeciesCode=="BF" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="BH" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="BO" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="CL" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="DB" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="GI" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="GR" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="LA" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="LB" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="PH" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="DG" & bhm_dat$SpType=="I"] <- "A"
bhm_dat$SpType[bhm_dat$SpeciesCode=="BA" & bhm_dat$SpType=="A"] <- "I" # Brooding anemone
bhm_dat$SpType[bhm_dat$SpeciesCode=="DC" & bhm_dat$SpType=="A"] <- "I" # Dungeness crab
bhm_dat$SpType[bhm_dat$SpeciesCode=="HF" & bhm_dat$SpType=="A"] <- "I" # Hydrocoral
bhm_dat$SpType[bhm_dat$SpeciesCode=="OB" & bhm_dat$SpType=="A"] <- "I" # Other barnacle
bhm_dat$SpType[bhm_dat$SpeciesCode=="PU" & bhm_dat$SpType=="A"] <- "I" # Purple urchin 
bhm_dat$SpType[bhm_dat$SpeciesCode=="RS" & bhm_dat$SpType=="A"] <- "I" # Rock scallop
bhm_dat$SpType[bhm_dat$SpeciesCode=="SI" & bhm_dat$SpType=="A"] <- "I" # Swimming anemone

# Typos that were removed, do not correspond with any code
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="GSC" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="JL" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="PY" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="TH" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="UA" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="YZ" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="MO" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="PG" & bhm_dat$SpType=="I"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="AN" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="AT" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="BR" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="DR" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="FY" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="LD" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="LM" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="NE" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="RG" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="RM" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="SY" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="*" & bhm_dat$SpType=="A"] <- "nontarget"
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="SU" & bhm_dat$SpType=="A"] <- "nontarget"

# Typo that was not feasibly observed at depth 
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="OR"& bhm_dat$SpType=="I"] <- "nontarget" # Pacific Oyster

# SpeciesCode code changes occurred during the duration of the BHM surveys
bhm_dat$SpeciesCode[bhm_dat$SpeciesCode=="LB" & bhm_dat$SpType=="A"] <- "LG"  # former name Laminaria bongardiana 

# Remove nontarget species codes
bhm_dat <- filter( bhm_dat, SpeciesCode!="nontarget" )

# Trim white space around Species code
bhm_dat$SpeciesCode <- trimws(bhm_dat$SpeciesCode)

# Remove rows with blank entry for Species code
bhm_dat = bhm_dat[!(is.na(bhm_dat$SpeciesCode) | bhm_dat$SpeciesCode==""), ]

#create unique species codes
bhm_dat$Species <- paste0(bhm_dat$SpType,"_",bhm_dat$SpeciesCode)

bhm_dat$Species[bhm_dat$Species=="I_ABL"] <- "ABL" 
bhm_dat$Species[bhm_dat$Species=="I_RSU"] <- "RSU" 
bhm_dat$Species[bhm_dat$Species=="I_PU"] <- "PSU" 
bhm_dat$Species[bhm_dat$Species=="I_GSU"] <- "GSU" 
bhm_dat$Species[bhm_dat$Species=="I_GDK"] <- "GDK" 
bhm_dat$Species[bhm_dat$Species=="I_HSC"] <- "HSC" 
bhm_dat$Species[bhm_dat$Species=="I_PO"] <- "PO" 
bhm_dat$Species[bhm_dat$Species=="I_PYC"] <- "PYC" 
bhm_dat$Species[bhm_dat$Species=="I_RSC"] <- "RSC" 
bhm_dat$Species[bhm_dat$Species=="A_ZO"] <- "ZO" 
bhm_dat$Species[bhm_dat$Species=="A_PT"] <- "PT" 
bhm_dat$Species[bhm_dat$Species=="A_MA"] <- "MA" 
bhm_dat$Species[bhm_dat$Species=="A_NT"] <- "NT" 

bhm_dat <- dplyr::select( bhm_dat, Survey, HKey, Year, Month, Day, LatDeep, LonDeep, LatShallow, LonShallow, Quadrat, CorDepthM, Species)

## Save files
save(bhm_dat, file="DataSources/Data/bhm_db_data.RData")
write.csv(bhm_dat, "DataSources/Data/bhm_db_quadrat.csv", row.names = F)


