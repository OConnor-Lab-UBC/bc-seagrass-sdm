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


###############################################################################

library(terra)
library(dplyr)
library(purrr)
library(tibble)
library(caret)

#marsh_10m    <- rast("raw_data/XiuchengMaps/Pacific_cover_2024_Sentinel2.tif")      # classes 1–5
# 1-5: tidal marsh, tidal flats, open water, water-land interface, low-elevation uplands

# --- Convert marsh raster to binary (classes 1–2 = marsh&tidalflats&upland) ---
# marsh_bin_10m <- ifel(marsh_10m %in% c(1,2,4,5), 1, 0)
# 
# # --- Resample marsh raster to 20 m (to match eelgrass raster) ---
# # Use min so along edges we make sure we are not losing areas that could potentially be eelgrass, precautionary principle
# marsh_bin_20m <- aggregate(marsh_bin_10m,fact = 2,      # 10 m → 20 m
#   fun  = max,
#   na.rm = TRUE
# )
# 
# # Align exactly to eelgrass grid (important)
# marsh_bin_20m <- project(marsh_bin_20m, eel_bin, method="near")
# #save(marsh_bin_20m, file = "raw_data/XiuchengMaps/marsh_mask_20m.RData")
# 
# writeRaster(marsh_bin_20m, file.path("raw_data/XiuchengMaps/marsh_mask_20m.tif"), overwrite = TRUE)

marsh_bin_20m <- rast("raw_data/XiuchengMaps/marsh_mask_20m.tif")   

rei_all <- terra::vrt(c("raw_data/REI/rei_hg.tif", "raw_data/REI/rei_ncc.tif", "raw_data/REI/rei_qcs.tif", "raw_data/REI/rei_sog.tif", "raw_data/REI/rei_wcvi.tif"), "rei.vrt", overwrite=T)
zero_vec <- vect("raw_data/substrate_20m/exposed_shorezone_bottompatches_sand_class16_17_27_28_30_plusextra.shp")
toodeep_vec <- vect("raw_data/substrate_20m/bathy_updated_toodeep_mask.shp")
zero_vec_surfgrass <- vect("raw_data/substrate_20m/exposed_shorezone_bottompatches_sand_class16_17_27_28_30_plusextra_surfgrass.shp")


#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"


# eelgrass
load("./code/output_data/model_results/combined_metrics_eelgrass_4_validations.RData")
# for eelgras XGBOOST nep is best model. tss threshold from cv is 0.031

# --- Load rasters ---
eelgrass_20m <- terra::vrt(c("raster/eelgrass/xgb_nep/eelgrass_predictions_hg_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ncc_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_qcs_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ss_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_wcvi_xgb_nep.tif"), "eelgrass_xgb_nep.vrt", overwrite=T)   # values 0–1
eelgrass_20m_glmm_spatial_nep <- terra::vrt(c("raster/eelgrass/nep_spatial/eelgrass_predictions_hg_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_ncc_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_qcs_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_ss_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_wcvi_nep_spatial.tif"), "eelgrass_nep_spatial.vrt", overwrite=T)   # values 0–1



# --- 1. Threshold eelgrass raster to binary (presence / absence) ---
# define threshold (adjust as needed)
thr <- 0.031
thr1 <- 0.037

eel_bin <- ifel(eelgrass_20m >= thr, 1, 0)
eel_bin_nepspatial <- ifel(eelgrass_20m_glmm_spatial_nep >= thr1, 1, 0)
writeRaster(eel_bin, file.path("raster/eelgrass/eelgrass_predictions_xgb_nep_binary_notmasked.tif"), overwrite = TRUE)
writeRaster(eel_bin_nepspatial, file.path("raster/eelgrass/eelgrass_predictions_nepspatial_binary_notmasked.tif"), overwrite = TRUE)


# --- 4. Remove eelgrass presence where marsh, tidal flats, upland and beach is present ---
eel_final <- ifel(marsh_bin_20m == 1, 0, eel_bin)

# remove eelgrass predictions in areas with rei > 0.20
rei_aligned <- project(rei_all, eel_final, method="near")

eel_final <- ifel(rei_aligned > 0.2, 0, eel_final)


# remove areas with high exposure sand
zero_vec <- project(zero_vec, eel_final)
zero_ras <- rasterize(zero_vec, eel_final, field = 1, background = NA)

eel_final <- ifel(!is.na(zero_ras), 0, eel_final)

#remove areas with incorrect bathy based on updated CHS products
toodeep_vec <- project(toodeep_vec, eel_final)
toodeep_ras <- rasterize(toodeep_vec, eel_final, field = 1, background = NA)

eel_final <- ifel(!is.na(toodeep_ras), 0, eel_final)

# remove eelgrass predictions in areas that are rock
substrate_aligned <- project(substrate_all, eel_final, method="near")

#surf_final <- ifelse(!is.na(substrate_aligned) & substrate_aligned %in% c(3, 4), 0, surf_final)
eel_final[!is.na(substrate_aligned) & (substrate_aligned %in% c(1))] <- 0


# Add in areas where eelgrass has been mapped. 
#load netforce eelgrass  data 2013-2023
eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))

values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1

eelgrass_indep <- resample(eelgrass_indep, eel_final, method = "near")

eel_final_plus_mapped <- ifel(eelgrass_indep == 1, 1, eel_final)

n1 <- global(eel_final_plus_mapped == 1, "sum", na.rm = TRUE)
n1

area_m2 <- n1 * 400
area_m2
# 583987600 m2

area_km2 <- area_m2 / 1e6
area_km2
# 583.9876 km2

area_ha <- area_m2 / 10000
area_ha
# 58398.76 hectares

writeRaster(eel_final_plus_mapped, file.path("raster/eelgrass/eelgrass_predictions.tif"), overwrite = TRUE)

#eel_final_plus_mapped <- rast("raster/eelgrass/eelgrass_predictions.tif")   


# surfgrass
load("./code/output_data/model_results/combined_metrics_surfgrass_4_validations.RData")
# for surfgrass bccm spatial is best model. tss threshold from cv is 0.014

# --- Load rasters ---
surfgrass_20m <- terra::vrt(c("raster/surfgrass/bccm_spatial/surfgrass_predictions_hg_bccm_spatial.tif", "raster/surfgrass/bccm_spatial/surfgrass_predictions_ncc_bccm_spatial.tif", "raster/surfgrass/bccm_spatial/surfgrass_predictions_qcs_bccm_spatial.tif", "raster/surfgrass/bccm_spatial/surfgrass_predictions_ss_bccm_spatial.tif", "raster/surfgrass/bccm_spatial/surfgrass_predictions_wcvi_bccm_spatial.tif"), "surfgrass_bccm_spatial.vrt", overwrite=T)   # values 0–1

# --- 1. Threshold surfgrass raster to binary (presence / absence) ---
# define threshold (adjust as needed)
thr <- 0.014

surf_bin <- ifel(surfgrass_20m >= thr, 1, 0)

# --- 4. Remove surfgrass presence where marsh, tidal flats, upland and beach is present ---
surf_final <- ifel(marsh_bin_20m == 1, 0, surf_bin)

# remove surfgrass predictions in areas sand and mud (sometimes you get this but it is not a very common occurance and sand is getting overpredicted)
substrate_aligned <- project(substrate_all, surf_final, method="near")

#surf_final <- ifelse(!is.na(substrate_aligned) & substrate_aligned %in% c(3, 4), 0, surf_final)
surf_final[!is.na(substrate_aligned) & (substrate_aligned %in% c(3, 4))] <- 0

# remove areas with high exposure sand
zero_vec_surfgrass <- project(zero_vec_surfgrass, surf_final)
zero_ras_surfgrass <- rasterize(zero_vec_surfgrass, surf_final, field = 1, background = NA)

surf_final <- ifel(!is.na(zero_ras_surfgrass), 0, surf_final)

#remove areas with incorrect bathy based on updated CHS products
toodeep_vec <- project(toodeep_vec, surf_final)
toodeep_ras <- rasterize(toodeep_vec, surf_final, field = 1, background = NA)

surf_final <- ifel(!is.na(toodeep_ras), 0, surf_final)

n1 <- global(surf_final == 1, "sum", na.rm = TRUE)
n1

area_m2 <- n1 * 400
area_m2
# 1047827600 m2

area_km2 <- area_m2 / 1e6
area_km2
# 1047.828 km2

area_ha <- area_m2 / 10000
area_ha
# 104782.8

writeRaster(surf_final, file.path("raster/surfgrass/surfgrass_predictions.tif"), overwrite = TRUE)















# compare predictions to field validation dataset to see if overprediction into intertidal areas and exposed sandy areas has been adequately dealt with
load("code/output_data/field_validation/validation_dataset.RData")

preds_extract <- terra::extract(eel_final_plus_mapped, validation_sf)
summary(preds_extract$eelgrass_xgb_nep)
validation_sf$final_eelgrass_preds <- preds_extract$eelgrass_xgb_nep

validation_sf <- validation_sf %>%
  relocate(PH_freq, PH, ZM_freq, ZM, .after = PC_PH)

# need to also drop sites that after reviewing the database should not have been included as they were aborted or vis was too bad to see species
validation_sf <- validation_sf %>%
  filter(!HKey %in% c(74, 84, 85, 161))

# noticed error in database that have alerted to Sandie so remove this once it is changed
validation_sf$PC_ZM[validation_sf$HKey == "125"] <- "51-75"
validation_sf$PC_PH[validation_sf$HKey == "125"] <- "0"

# here are some more errors that need to be addressed in database still
validation_sf$PC_ZM[validation_sf$HKey == "160"] <- "0"
validation_sf$PC_ZM[validation_sf$HKey == "536"] <- "0"
validation_sf$PC_ZM[validation_sf$HKey == "585"] <- "0"
validation_sf$PC_PH[validation_sf$HKey == "160"] <- "0"
validation_sf$PC_PH[validation_sf$HKey == "536"] <- "0"
validation_sf$PC_PH[validation_sf$HKey == "588"] <- "26-50"
validation_sf$ZM[validation_sf$HKey == "35"] <- "0"
validation_sf$ZM[validation_sf$HKey == "521"] <- "0"
validation_sf$PH[validation_sf$HKey == "521"] <- "1"
validation_sf$ZM[validation_sf$HKey == "545"] <- "0"
validation_sf$PH[validation_sf$HKey == "545"] <- "1"

# determine difference between observed and modelled environmental layers
# if either observed or modelled layers were mixed these were always retained as we are just looking at predominant substrate and most sites are actually mixed in reality
validation_sf <- validation_sf %>%
  mutate(Depth_diff = abs(bathy_mod - avgCorDepth_obs),
         Slope_diff = abs(slope_mod - slope_obs),
         Substrate_diff =
           if_else(
             substrate_mod == "Mixed" | substrate_mode_obs == "Mixed", TRUE,
             substrate_mod == substrate_mode_obs))

validation_sf <- validation_sf %>% select (HKey, Region, Survey, Visibility, High.vegetation, avgCorDepth_obs, bathy_mod,
                                           Depth_diff, slope_obs, slope_mod, Slope_diff, substrate_mode_obs, 
                                           substrate_mod, Substrate_diff, 
                                           ZM, ZM_freq, PC_ZM, final_eelgrass_preds, 
                                           PH, PH_freq, PC_PH, geometry)
summary(validation_sf)

pred_col <- "final_eelgrass_preds"
obs_pa   <- "ZM"     # 0/1 observed
obs_pc   <- "PC_ZM"  # percent cover bins

# ---- Prep df ----
df <- validation_sf %>%
  mutate(
    PC_mid = case_when(
      PC_ZM == "0" ~ 0,
      PC_ZM == "1-25" ~ 12.5,
      PC_ZM == "26-50" ~ 38,
      PC_ZM == "51-75" ~ 63,
      PC_ZM == "76-100" ~ 88
    ),
    Presence = if_else(PC_ZM == "0", 0, 1)
  )

# ---- Correlations (binary preds vs ordinal cover midpoint) ----
spearman <- cor.test(df[[pred_col]], df$PC_mid, method = "spearman")
kendall  <- cor.test(df[[pred_col]], df$PC_mid, method = "kendall")

cor_results <- tibble(
  model        = pred_col,
  spearman_rho = unname(spearman$estimate),
  spearman_p   = spearman$p.value,
  kendall_tau  = unname(kendall$estimate),
  kendall_p    = kendall$p.value
)

cor_results


# ---- Confusion matrix (NO thresholding needed) ----
cm <- confusionMatrix(
  factor(df[[pred_col]], levels = c(0, 1)),
  factor(df[[obs_pa]], levels = c(0, 1)),
  positive = "1"
)

cm


# ---- Summary metrics ----
sens <- cm$byClass["Sensitivity"]
spec <- cm$byClass["Specificity"]

metrics <- tibble(
  model = pred_col,
  accuracy = cm$overall["Accuracy"],
  kappa    = cm$overall["Kappa"],
  sensitivity = sens,
  specificity = spec,
  bal_accuracy = cm$byClass["Balanced Accuracy"],
  precision = cm$byClass["Pos Pred Value"],
  npv       = cm$byClass["Neg Pred Value"],
  F1        = cm$byClass["F1"],
  TSS       = sens + spec - 1
)

metrics


# observed presence rate
obs_prev <- mean(df$ZM == 1, na.rm = TRUE)

# predicted presence rate
pred_prev <- mean(df$final_eelgrass_preds == 1, na.rm = TRUE)

obs_prev
pred_prev

prev_table <- tibble(
  n_total = sum(!is.na(df$ZM) & !is.na(df$final_eelgrass_preds)),
  obs_presence_n = sum(df$ZM == 1, na.rm = TRUE),
  obs_presence_rate = mean(df$ZM == 1, na.rm = TRUE),
  pred_presence_n = sum(df$final_eelgrass_preds == 1, na.rm = TRUE),
  pred_presence_rate = mean(df$final_eelgrass_preds == 1, na.rm = TRUE)
)

prev_table
cm$table

#Strong at ruling out eelgrass (NPV 0.94, spec 0.89)
#Moderately good at detecting eelgrass (sens 0.64)
#tends to overpredict presence

TN <- 374
FP <- 48
FN <- 22
TP <- 39

metrics_extra <- tibble(
  n_total = TN + FP + FN + TP,
  
  obs_prev = (TP + FN) / (TN + FP + FN + TP),
  pred_prev = (TP + FP) / (TN + FP + FN + TP),
  
  sensitivity = TP / (TP + FN),
  specificity = TN / (TN + FP),
  
  FPR = FP / (FP + TN),
  FNR = FN / (FN + TP),
  
  precision = TP / (TP + FP),
  NPV = TN / (TN + FN),
  
  accuracy = (TP + TN) / (TN + FP + FN + TP),
  TSS = (TP / (TP + FN)) + (TN / (TN + FP)) - 1
)

metrics_extra

prev_plot_df <- tibble(
  type = c("Observed presence rate", "Predicted presence rate"),
  rate = c(obs_prev, pred_prev)
)

prev_plot_df

library(ggplot2)

ggplot(prev_plot_df, aes(x = type, y = rate)) +
  geom_col() +
  ylab("Proportion presence") +
  xlab("") +
  theme_minimal()




