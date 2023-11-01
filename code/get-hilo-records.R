#### load packages####
library(dplyr)
library(reshape2)

#### Outputs ####
# Subfolders for outputs
outdir <- file.path("DataSources", "Data")
# Create the main directory, and subfolders
dir.create( outdir, recursive=TRUE )
# final field names for dive surveys
fnames <-  c( "Type", "Source", "Survey", "HKey", "Method",
              "Year","Month", "Day","LatDeep","LonDeep",
              "LatShallow","LonShallow", "Transect", "Quadrat", 
              "CorDepthM","SpNum", "Species")  


dat <- read.csv("DataSources/sql/HiLo_quadrats.csv",header=T, stringsAsFactors = F)


# change lat long into decimal degrees
dat$LonShallow <- (-1)*abs(dat$Long_deg_s+(dat$Long_min_s/60))
dat$LatShallow <- (dat$Lat_deg_s+(dat$Lat_min_s/60))
dat$LonDeep <- (-1)*abs(dat$Long_deg_d+(dat$Long_min_d/60))
dat$LatDeep <- (dat$Lat_deg_d+(dat$Lat_min_d/60))

## Melt all invert and algae species into species column
# list p/a species. 
sp_pa_hilo <- c("RSU","GSU", "PSU", "GDK","PT", "ZO", "MA","NT","ABL","MU","PO","PYN")
# melt table for p/a species
hilo_dat <- melt(dat, measure.vars=sp_pa_hilo, value.name = "SpNum", variable.name = "Species")

## Create new fields and rename
# Type
hilo_dat$Type <- "Research"
# Source
hilo_dat$Source <- "MSEA"
# Method
hilo_dat$Method <- "Dive"
# Rename
hilo_dat <- hilo_dat[fnames]
# Order rows
hilo_dat <- hilo_dat[order(hilo_dat$HKey, hilo_dat$Transect, hilo_dat$Quadrat),]


## Save files
save(hilo_dat, file="DataSources/Data/hilo_db_quadrat.RData")
write.csv(hilo_dat, "DataSources/Data/hilo_db_quadrat.csv", row.names = F)
