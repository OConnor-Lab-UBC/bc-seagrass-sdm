###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# Compare independent data sets to sdms
# when adding multiple different sdms will want to look at areas correctly classified (when creating sf object just have different columns be different models)
#
###############################################################################

#load packages####
library(sf)
library(tidyverse)
library(terra)
library(ggplot2)

#### Eelgrass

#load netforce eelgrass  data 2013-2023
eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))

# load and mosaic eelgrass SDM rasters
hg <- rast("raster/eelgrass_predictions_hg.tif")
ncc <- rast("raster/eelgrass_predictions_ncc.tif")
qcs <-rast("raster/eelgrass_predictions_qcs.tif")
wcvi <- rast("raster/eelgrass_predictions_wcvi.tif")
ss <-rast("raster/eelgrass_predictions_ss.tif")

eelgrass_sdm<- mosaic(hg, ncc, qcs, wcvi, ss, fun = "mean")

#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"


# Ensure rasters align
eelgrass_sdm_resampled <- resample(eelgrass_sdm, eelgrass_indep, method = "bilinear")

eelgrass_stack <- c(eelgrass_indep, eelgrass_sdm_resampled)
names(eelgrass_stack) <- c("obs", "mod")

r_points <- as.points(eelgrass_stack, na.rm = TRUE)
eelgrass_sf<- r_points %>% st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

#extract substrate data to see if that is causing mismatch
substrate_extract <- terra::extract(substrate_all, eelgrass_sf)
eelgrass_sf$substrate <- substrate_extract$substrate
eelgrass_sf$substrate <- c("Rock", "Mixed", "Sand", "Mud")[eelgrass_sf$substrate]

eelgrass_sf$rock_group <- ifelse(eelgrass_sf$substrate == "Rock", "Rock", "Not Rock")

# plot modelled probability by rock group
boxplot(mod ~ rock_group, data = eelgrass_sf,
        xlab = "Modelled Substrate Type",
        ylab = "Modelled Probability")

ggplot(eelgrass_sf, aes(x = as.factor(obs), y = mod, fill = rock_group)) +
  geom_boxplot() +
  labs(title = "Modelled Probability by Observation Frequency and Substrate Type",
       x = "Years Observed",
       y = "Modelled Probability") +
  scale_fill_manual(values = c("Rock" = "grey", "Not Rock" = "green")) +
  theme_minimal()

# should remove areas with modelled rock substrate from comparison as the substrate model is likely wrong in these areas
eelgrass_sf_norock <- eelgrass_sf %>% filter(substrate != "Rock")

boxplot(mod ~ obs, data = eelgrass_sf_norock,
        #main = "Modelled Probability by Observation Frequency",
        xlab = "Years with Observed Eelgrass",
        ylab = "Modelled Probability")

# to help estimate a threshold should use high confidence locations - ie where eelgrass has been observed more than once
high_confidence <- eelgrass_sf_norock[eelgrass_sf_norock$obs >= 2, ]
threshold <- quantile(high_confidence$mod, 0.10)
print(threshold) # threshold is 0.124, if using all obs and not just ones with high confidence than threshold is 0.14


model_binary <- eelgrass_sdm_resampled >= threshold

false_negatives <- (eelgrass_indep >= 1) & (model_binary == 0)
fn_raster <- classify(false_negatives, cbind(0, NA))

plot(fn_raster, main = "False Negatives (Missed Observed Eelgrass)", col = "red") 

# Evaluate model scoring with correlation
spearman_cor<- cor(eelgrass_sf_norock$obs, eelgrass_sf_norock$mod, method = "spearman")
# -0.03034606  A Spearman correlation of –0.0303 suggests that higher modelled probabilities do not correspond in any meaningful way to more years of observed presence.

# modify when have different models to compare; this includes the areas with rock, which should end up equally out across models. i just didnt want to include rock in developing the threshold
eelgrass_sf$obs_bin <- ifelse(eelgrass_sf$obs > 0, 1, 0)
eelgrass_sf$mod1_bin <- ifelse(eelgrass_sf$mod >= threshold, 1, 0)

cell_area <- 20 * 20 #(area in each cell is 400 m2)

#identify true positives
eelgrass_sf$TP1 <- with(eelgrass_sf, obs_bin == 1 & mod1_bin == 1)
area_tp1 <- sum(eelgrass_sf$TP1, na.rm = TRUE) * cell_area  # in m²
area_tp1_km2 <- area_tp1 / 1e6 # 36km 2

# % of observed area correctly predicted
n_obs_cells <- sum(eelgrass_sf$obs_bin == 1, na.rm = TRUE)
obs_area <- n_obs_cells * cell_area

percent_tp1 <- area_tp1 / obs_area * 100  #58%

save(eelgrass_sf, threshold, eelgrass_sf_norock, spearman_cor, file = "code/output_data/independent_validation/eelgrass_independent_validation.RData")





#### Surfgrass

#load surfgrass independent data
surfgrass_indep <- rast(c("code/output_data/independent_validation/surfgrass_validation_raster_2013_2024.tif"))

# load and mosaic surfgrass SDM rasters
hg <- rast("raster/surfgrass_predictions_hg.tif")
ncc <- rast("raster/surfgrass_predictions_ncc.tif")
qcs <-rast("raster/surfgrass_predictions_qcs.tif")
wcvi <- rast("raster/surfgrass_predictions_wcvi.tif")
ss <-rast("raster/surfgrass_predictions_ss.tif")

surfgrass_sdm<- mosaic(hg, ncc, qcs, wcvi, ss, fun = "mean")

# Ensure rasters align
surfgrass_sdm_resampled <- resample(surfgrass_sdm, surfgrass_indep, method = "bilinear")

surfgrass_stack <- c(surfgrass_indep, surfgrass_sdm_resampled)
names(surfgrass_stack) <- c("obs", "mod")

r_points <- as.points(surfgrass_stack, na.rm = TRUE)
surfgrass_sf<- r_points %>% st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

#extract substrate data to see if that is causing mismatch
substrate_extract <- terra::extract(substrate_all, surfgrass_sf)
surfgrass_sf$substrate <- substrate_extract$substrate
surfgrass_sf$substrate <- c("Rock", "Mixed", "Sand", "Mud")[surfgrass_sf$substrate]

surfgrass_sf$soft_group <- ifelse(surfgrass_sf$substrate == "Mud", "Soft", "Not soft")
surfgrass_sf$soft_group <- ifelse(surfgrass_sf$substrate == "Sand", "Soft", "Not soft")
# note can have surfgrass on soft substrate, but that is not the norm

# plot modelled probability by substrate group
boxplot(mod ~ soft_group, data = surfgrass_sf,
        xlab = "Modelled Substrate Type",
        ylab = "Modelled Probability")

# should remove areas with modelled soft substrate from comparison 
surfgrass_sf_nosoft <- surfgrass_sf %>% filter(soft_group == "Not soft")

threshold <- quantile(surfgrass_sf_nosoft$mod, 0.10)
print(threshold) # threshold is 0.017 (which is low!). Hesitent to use this as a threshold as the independent validation data is not the best from shorezone. Can use it to help model selection though

model_binary <- surfgrass_sdm_resampled >= threshold

false_negatives <- (surfgrass_indep >= 1) & (model_binary == 0)
fn_raster <- classify(false_negatives, cbind(0, NA))

plot(fn_raster, main = "False Negatives (Missed Observed Surfgrass)", col = "red") 

# modify when have different models to compare; this includes the areas with rock, which should end up equally out across models. i just didnt want to include rock in developing the threshold
surfgrass_sf$obs_bin <- ifelse(surfgrass_sf$obs > 0, 1, 0)
surfgrass_sf$mod1_bin <- ifelse(surfgrass_sf$mod >= threshold, 1, 0)

cell_area <- 20 * 20 #(area in each cell is 400 m2)

#identify true positives
surfgrass_sf$TP1 <- with(surfgrass_sf, obs_bin == 1 & mod1_bin == 1)
area_tp1 <- sum(surfgrass_sf$TP1, na.rm = TRUE) * cell_area  # in m²
area_tp1_km2 <- area_tp1 / 1e6 # 22km 2

# % of observed area correctly predicted
n_obs_cells <- sum(surfgrass_sf$obs_bin == 1, na.rm = TRUE)
obs_area <- n_obs_cells * cell_area

percent_tp1 <- area_tp1 / obs_area * 100  #84%

save(surfgrass_sf, threshold, surfgrass_sf_nosoft,  file = "code/output_data/independent_validation/surfgrass_independent_validation.RData")

