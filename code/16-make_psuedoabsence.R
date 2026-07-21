###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# Generate psuedo absences based on seagrass presences in netforce data all years and dfo dive data
#
###############################################################################

library(sf)
library(terra)
library(dplyr)
library(purrr)
library(tidyr)
library(Metrics)
library(spThin)


source("code/modelling-functions.R")
template_rast <- rast("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif")



#eelgrass
eelgrass_domain <- terra::vrt(c("raster/eelgrass/xgb_nep/eelgrass_predictions_hg_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ncc_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_qcs_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ss_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_wcvi_xgb_nep.tif"), "eelgrass_xgb_nep.vrt", overwrite=T)   # values 0–1

terra::values(eelgrass_domain) <- ifelse(!is.na(terra::values(eelgrass_domain)), 1, NA)

# load independent eelgrass rasters
eelgrass_indep <- rast("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif")
values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1
eelgrass_indep2 <- rast("code/output_data/independent_validation/BCeelgrass_netforce_1974_2012.tif")
values(eelgrass_indep2)[values(eelgrass_indep2) >= 1] <- 1
dfo2024<-vect("code/output_data/processed_observations/SpatializedQuadrats_aggregated_2024only_ZO.shp")

dfo2024_rast <- terra::rasterize(
  x = dfo2024,
  y = template_rast,
  field = "ZO",
  touches = TRUE
)

#Peng et al 2026 remote sense global seagrass dataset https://www.nature.com/articles/s41586-026-10704-3
peng<-vect("raw_data/Pengetal2026/bc_seagrass_2023_2024.shp")
peng_proj <- project(peng, template_rast)

peng_rast <- terra::rasterize(
  x = peng_proj,
  y = template_rast,
  field = "b1",
  touches = TRUE
)


# load model inputs
load("code/output_data/seagrass_model_inputs.RData")
data <- seagrass_data_long %>%
  filter(species == "ZO", presence == 1) %>%
  select(X_m, Y_m) %>%
  distinct()
data_sf <- st_as_sf(data, coords = c("X_m", "Y_m"), crs = 3005)
data_sf$pres <- 1
data_vect <- terra::vect(data_sf)
training_pres_rast <- terra::rasterize(
  x = data_vect,
  y = template_rast,
  field = "pres",
  touches = TRUE
)
training_pres_rast <- ifel(!is.na(training_pres_rast), 1, NA)
freq(training_pres_rast)

#create dataset to inform psuedoabsence creation to avoid areas where eelgrass has been/or likely has been observed at any point in time
combined_exclusion <- terra::ifel(
  !is.na(training_pres_rast) | !is.na(eelgrass_indep) | !is.na(eelgrass_indep2) | !is.na(dfo2024_rast) | !is.na(peng_rast),
  1,
  NA
)

freq(combined_exclusion)

# make sure predictions and independent datasets are aligned and crop to eelgrass predictions
combined_exclusion_resample<- terra::crop(
  terra::resample(combined_exclusion, eelgrass_domain, method = "bilinear"),
  eelgrass_domain)
writeRaster(combined_exclusion_resample, "code/output_data/independent_validation/combined_exclusion.tif", overwrite = TRUE)




#### need to make presence dataset from netforce 2013-2023 data that is thinned 
# convert eelgrass presence cells to points
# from looking previously there is no difference in the predictions if eelgrass has been observed more times so just change all to 1. A Spearman correlation of –0.0303 suggests that higher modelled probabilities do not correspond in any meaningful way to more years of observed presence.
values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1
names(eelgrass_indep)<-"obs"

eelgrass_pts <- terra::as.points(
  eelgrass_indep,
  na.rm = TRUE
)

coords <- terra::crds(eelgrass_pts)

eelgrass_sf <- st_as_sf(
  data.frame(x = coords[,1],
             y = coords[,2]),
  coords = c("x","y"),
  crs = terra::crs(eelgrass_indep)
)

#randomize points
set.seed(101)
eelgrass_sf <- eelgrass_sf[sample(nrow(eelgrass_sf)), ]

min_dist <- 100  # metres

selected <- vector("logical", nrow(eelgrass_sf))
selected[1] <- TRUE

selected_xy <- st_coordinates(eelgrass_sf[1, ])

for (i in 2:nrow(eelgrass_sf)) {
  
  pt_xy <- st_coordinates(eelgrass_sf[i, ])
  
  d <- sqrt(
    (selected_xy[,1] - pt_xy[1,1])^2 +
      (selected_xy[,2] - pt_xy[1,2])^2
  )
  
  if (all(d >= min_dist)) {
    selected[i] <- TRUE
    selected_xy <- rbind(selected_xy, pt_xy)
  }
}

eelgrass_thinned <- eelgrass_sf[selected, ]
save(eelgrass_thinned, file = "code/output_data/independent_validation/thinned_netforce_2013_2023_eelgrassobs.RData")
st_write(eelgrass_thinned, "code/output_data/independent_validation/thinned_netforce_2013_2023_eelgrassobs.shp", append=FALSE)
# there were a few weird points on land and in the middle of the pacific so deleted them in the shapefile, though they will be eliminated when extracting rpedictions


candidate_rast <- make_pseudoabsence_candidate_rast(
  domain_rast = eelgrass_domain,
  exclusion_rast = combined_exclusion_resample,
  buffer_m = 100,
  filename = "code/output_data/independent_validation/pseudoabsence_candidate_cells.tif",
  overwrite = TRUE
)

terra::global(!is.na(candidate_rast), "sum", na.rm = TRUE)


#sample 10 pseudoabsence sets

n_pres <- nrow(eelgrass_thinned)

pa_list <- vector("list", 10)

for (i in seq_along(pa_list)) {
  pa_list[[i]] <- sample_pseudoabsences_min_dist(
    candidate_rast = candidate_rast,
    n_pa = n_pres,
    min_dist_m = 100,
    seed = 100 + i
  )
  pa_list[[i]]$pa_df$pa_set <- i
}
pseudoabsences_all <- do.call(
  rbind,
  lapply(pa_list, function(x) x$pa_df)
)
save(
  pseudoabsences_all,
  file = "code/output_data/independent_validation/pseudoabsences_thinned_netforce_2013_2023.RData"
)

pseudoabsences_all_sf <- st_as_sf(
  pseudoabsences_all,
  coords = c("x", "y"),
  crs = terra::crs(eelgrass_domain)
)

st_write(
  pseudoabsences_all_sf,
  "code/output_data/independent_validation/pseudoabsences_all.shp",
  delete_dsn = TRUE
)



#surfgrass
surfgrass_domain <- terra::vrt(c("raster/surfgrass/xgb_nep/surfgrass_predictions_hg_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_ncc_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_qcs_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_ss_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_wcvi_xgb_nep.tif"), "surfgrass_xgb_nep.vrt", overwrite=T)   # values 0–1

terra::values(surfgrass_domain) <- ifelse(!is.na(terra::values(surfgrass_domain)), 1, NA)

# load independent surfgrass rasters
surfgrass_indep <- rast("code/output_data/independent_validation/surfgrass_validation_raster_2013_2024.tif")
values(surfgrass_indep)[values(surfgrass_indep) >= 1] <- 1
dfo2024<-vect("code/output_data/processed_observations/SpatializedQuadrats_aggregated_2024only_PH.shp")

dfo2024_rast <- terra::rasterize(
  x = dfo2024,
  y = template_rast,
  field = "PH",
  touches = TRUE
)

# load model inputs
load("code/output_data/seagrass_model_inputs.RData")
data <- seagrass_data_long %>%
  filter(species == "PH", presence == 1) %>%
  select(X_m, Y_m) %>%
  distinct()
data_sf <- st_as_sf(data, coords = c("X_m", "Y_m"), crs = 3005)
data_sf$pres <- 1
data_vect <- terra::vect(data_sf)
training_pres_rast <- terra::rasterize(
  x = data_vect,
  y = template_rast,
  field = "pres",
  touches = TRUE
)
training_pres_rast <- ifel(!is.na(training_pres_rast), 1, NA)
freq(training_pres_rast)

#create dataset to inform psuedoabsence creation to avoid areas where surfgrass has been/or likely has been observed at any point in time
combined_exclusion <- terra::ifel(
  !is.na(training_pres_rast) | !is.na(surfgrass_indep) | !is.na(dfo2024_rast) ,
  1,
  NA
)

freq(combined_exclusion)

# make sure predictions and independent datasets are aligned and crop to surfgrass predictions
combined_exclusion_resample<- terra::crop(
  terra::resample(combined_exclusion, surfgrass_domain, method = "bilinear"),
  surfgrass_domain)
writeRaster(combined_exclusion_resample, "code/output_data/independent_validation/combined_exclusion_surfgrass.tif", overwrite = TRUE)




#### need to make presence dataset from netforce 2013-2023 data that is thinned 
# convert surfgrass presence cells to points
names(surfgrass_indep)<-"obs"

surfgrass_pts <- terra::as.points(
  surfgrass_indep,
  na.rm = TRUE
)

coords <- terra::crds(surfgrass_pts)

surfgrass_sf <- st_as_sf(
  data.frame(x = coords[,1],
             y = coords[,2]),
  coords = c("x","y"),
  crs = terra::crs(surfgrass_indep)
)

#randomize points
set.seed(101)
surfgrass_sf <- surfgrass_sf[sample(nrow(surfgrass_sf)), ]

min_dist <- 100  # metres

selected <- vector("logical", nrow(surfgrass_sf))
selected[1] <- TRUE

selected_xy <- st_coordinates(surfgrass_sf[1, ])

for (i in 2:nrow(surfgrass_sf)) {
  
  pt_xy <- st_coordinates(surfgrass_sf[i, ])
  
  d <- sqrt(
    (selected_xy[,1] - pt_xy[1,1])^2 +
      (selected_xy[,2] - pt_xy[1,2])^2
  )
  
  if (all(d >= min_dist)) {
    selected[i] <- TRUE
    selected_xy <- rbind(selected_xy, pt_xy)
  }
}

surfgrass_thinned <- surfgrass_sf[selected, ]
save(surfgrass_thinned, file = "code/output_data/independent_validation/thinned_2013_2024_surfgrassobs.RData")
st_write(surfgrass_thinned, "code/output_data/independent_validation/thinned_2013_2024_surfgrassobs.shp", append=FALSE)


candidate_rast <- make_pseudoabsence_candidate_rast(
  domain_rast = surfgrass_domain,
  exclusion_rast = combined_exclusion_resample,
  buffer_m = 100,
  filename = "code/output_data/independent_validation/pseudoabsence_candidate_cells_surfgrass.tif",
  overwrite = TRUE
)

terra::global(!is.na(candidate_rast), "sum", na.rm = TRUE)


#sample 10 pseudoabsence sets

n_pres <- nrow(surfgrass_thinned)

pa_list <- vector("list", 10)

for (i in seq_along(pa_list)) {
  pa_list[[i]] <- sample_pseudoabsences_min_dist(
    candidate_rast = candidate_rast,
    n_pa = n_pres,
    min_dist_m = 100,
    seed = 100 + i
  )
  pa_list[[i]]$pa_df$pa_set <- i
}
pseudoabsences_all <- do.call(
  rbind,
  lapply(pa_list, function(x) x$pa_df)
)
save(
  pseudoabsences_all,
  file = "code/output_data/independent_validation/pseudoabsences_thinned_2013_2023_surfgrass.RData"
)

pseudoabsences_all_sf <- st_as_sf(
  pseudoabsences_all,
  coords = c("x", "y"),
  crs = terra::crs(surfgrass_domain)
)

st_write(
  pseudoabsences_all_sf,
  "code/output_data/independent_validation/pseudoabsences_all_surfgrass.shp",
  delete_dsn = TRUE
)
