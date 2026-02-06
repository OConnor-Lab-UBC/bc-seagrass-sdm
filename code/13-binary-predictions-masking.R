###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# take best model, apply binary threshold, mask with remote sensed marsh and tidal flats to deal with overpredictions
# add mapped eelgrass, run through field validation again to assess how well masked layer is doing and if deals with overprediciton in intertidal
#compare to human anthropogenic: shoreline modification and culmulative impacts 


# TO DO
# decide on best threshold, from model and whether to modify it based on field validation
###############################################################################

library(terra)

# --- Load rasters ---
eelgrass_20m <- rast("eelgrass_20m.tif")   # values 0–1
marsh_10m    <- rast("marsh_10m.tif")      # classes 1–5

# --- 1. Threshold eelgrass raster to binary (presence / absence) ---
# define threshold (adjust as needed)
thr <- 0.5

eel_bin <- eelgrass_20m >= thr
eel_bin <- classify(eel_bin, cbind(NA, 0))  # ensure FALSE = 0, TRUE = 1
names(eel_bin) <- "eelgrass_bin"

# --- 2. Convert marsh raster to binary (classes 1–3 = marsh) ---
marsh_bin_10m <- marsh_10m %in% c(1, 2, 3)
marsh_bin_10m <- classify(marsh_bin_10m, cbind(NA, 0))
names(marsh_bin_10m) <- "marsh_bin"

# --- 3. Resample marsh raster to 20 m (to match eelgrass raster) ---
# Use max so any marsh presence within a 20 m cell is retained
marsh_bin_20m <- aggregate(
  marsh_bin_10m,
  fact = 2,        # 10 m → 20 m
  fun  = max,
  na.rm = TRUE
)

# Align exactly to eelgrass grid (important)
marsh_bin_20m <- resample(marsh_bin_20m, eel_bin, method = "near")

# --- 4. Remove eelgrass presence where marsh is present ---
eel_final <- eel_bin
eel_final[marsh_bin_20m == 1] <- 0

names(eel_final) <- "eelgrass_final"