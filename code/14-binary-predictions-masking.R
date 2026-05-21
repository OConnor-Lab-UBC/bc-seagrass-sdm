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
library(sf)
library(tidyr)
library(tibble)
library(ggplot2)

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
surfgrass_20m_glmm_nospatial_nep <- terra::vrt(c("raster/surfgrass/nep_nospatial/surfgrass_predictions_hg_nep_nospatial.tif", "raster/surfgrass/nep_nospatial/surfgrass_predictions_ncc_nep_nospatial.tif", "raster/surfgrass/nep_nospatial/surfgrass_predictions_qcs_nep_nospatial.tif", "raster/surfgrass/nep_nospatial/surfgrass_predictions_ss_nep_nospatial.tif", "raster/surfgrass/nep_nospatial/surfgrass_predictions_wcvi_nep_nospatial.tif"), "surfgrass_nep_nospatial.vrt", overwrite=T)   # values 0–1

# --- 1. Threshold surfgrass raster to binary (presence / absence) ---
# define threshold (adjust as needed)
thr <- 0.014
thr1 <- 0.014

surf_bin <- ifel(surfgrass_20m >= thr, 1, 0)
surf_bin_nepnospatial <- ifel(surfgrass_20m_glmm_nospatial_nep >= thr1, 1, 0)
writeRaster(surf_bin, file.path("raster/surfgrass/surfgrass_predictions_bccmspatial_binary_notmasked.tif"), overwrite = TRUE)
writeRaster(surf_bin_nepnospatial, file.path("raster/surfgrass/surfgrass_predictions_nepnospatial_binary_notmasked.tif"), overwrite = TRUE)


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












eel_masked<- rast("raster/eelgrass/eelgrass_predictions.tif")
eel_unmasked <-rast("raster/eelgrass/eelgrass_predictions_xgb_nep_binary_notmasked.tif")


# compare predictions to field validation dataset to see if overprediction into intertidal areas and exposed sandy areas has been adequately dealt with
load("code/output_data/field_validation/validation_dataset.RData")

masked_extract <- terra::extract(eel_masked, validation_sf)
unmasked_extract <- terra::extract(eel_unmasked, validation_sf)

validation_sf$pred_masked_binary <- masked_extract[[2]]
validation_sf$pred_unmasked_binary <- unmasked_extract[[2]]

validation_sf <- validation_sf %>%
  filter(!HKey %in% c(74, 84, 85, 161))

validation_sf <- validation_sf %>%
  relocate(PH_freq, PH, ZM_freq, ZM, .after = PC_PH)

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

# validation_sf <- validation_sf %>% select (HKey, Region, Survey, Visibility, High.vegetation, avgCorDepth_obs, bathy_mod,
#                                            Depth_diff, slope_obs, slope_mod, Slope_diff, substrate_mode_obs, 
#                                            substrate_mod, Substrate_diff, 
#                                            ZM, ZM_freq, PC_ZM, final_eelgrass_preds, 
#                                            PH, PH_freq, PC_PH, geometry)

# ---- Prep df ----
df <- validation_sf %>%
  mutate(
    Presence = as.integer(ZM),
    pred_masked_binary = as.integer(pred_masked_binary),
    pred_unmasked_binary = as.integer(pred_unmasked_binary)
  )




binary_validation_metrics <- function(data, obs_col, pred_col, model_name = pred_col) {
  
  obs <- data[[obs_col]]
  pred <- data[[pred_col]]
  
  keep <- !is.na(obs) & !is.na(pred)
  obs <- obs[keep]
  pred <- pred[keep]
  
  tp <- sum(pred == 1 & obs == 1)
  tn <- sum(pred == 0 & obs == 0)
  fp <- sum(pred == 1 & obs == 0)
  fn <- sum(pred == 0 & obs == 1)
  
  sensitivity <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  npv <- if ((tn + fn) == 0) NA_real_ else tn / (tn + fn)
  
  tibble(
    model = model_name,
    n = tp + tn + fp + fn,
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    observed_prevalence = mean(obs == 1),
    predicted_prevalence = mean(pred == 1),
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = sensitivity + specificity - 1,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    false_positive_rate = if ((fp + tn) == 0) NA_real_ else fp / (fp + tn),
    false_negative_rate = if ((fn + tp) == 0) NA_real_ else fn / (fn + tp),
    precision = precision,
    npv = npv,
    F1 = if (
      is.na(precision) | is.na(sensitivity) |
      (precision + sensitivity) == 0
    ) {
      NA_real_
    } else {
      2 * precision * sensitivity / (precision + sensitivity)
    }
  )
}


binary_pred_cols <- c(
  "pred_unmasked_binary",
  "pred_masked_binary"
)

binary_metrics <- map_dfr(binary_pred_cols, function(p) {
  binary_validation_metrics(
    data = df,
    obs_col = "Presence",
    pred_col = p,
    model_name = p
  )
})

binary_metrics


binary_metrics_clean <- binary_metrics %>%
  select(
    model,
    n,
    observed_prevalence,
    predicted_prevalence,
    sensitivity,
    specificity,
    TSS,
    balanced_accuracy,
    false_positive_rate,
    false_negative_rate,
    precision,
    npv,
    F1,
    TP, TN, FP, FN
  )

binary_metrics_clean

metric_change <- binary_metrics_clean %>%
  select(
    model,
    predicted_prevalence,
    sensitivity,
    specificity,
    TSS,
    balanced_accuracy,
    false_positive_rate,
    false_negative_rate,
    precision,
    npv,
    F1
  ) %>%
  pivot_longer(
    cols = -model,
    names_to = "metric",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = model,
    values_from = value
  ) %>%
  mutate(
    change_masked_minus_unmasked =
      pred_masked_binary - pred_unmasked_binary
  )

metric_change



df_outcomes <- df %>%
  mutate(site_id = row_number()) %>%
  pivot_longer(
    cols = all_of(binary_pred_cols),
    names_to = "model",
    values_to = "prediction"
  ) %>%
  mutate(
    outcome = case_when(
      Presence == 1 & prediction == 1 ~ "True presence",
      Presence == 1 & prediction == 0 ~ "False absence",
      Presence == 0 & prediction == 1 ~ "False presence",
      Presence == 0 & prediction == 0 ~ "True absence",
      TRUE ~ NA_character_
    ),
    outcome = factor(
      outcome,
      levels = c(
        "True presence",
        "False absence",
        "False presence",
        "True absence"
      )
    )
  )


df_outcomes %>%
  count(model, outcome)

# have to do with with bathy_mod because so many of the intertidal sites could not observe depth

df_outcomes <- df_outcomes %>%
  mutate(
    depth_zone = case_when(
      bathy_mod < 0 ~ "Intertidal / exposed",
      bathy_mod >= 0 & bathy_mod < 2 ~ "Very shallow subtidal",
      bathy_mod >= 2 & bathy_mod < 5 ~ "Shallow subtidal",
      bathy_mod >= 5 ~ "Deeper subtidal",
      TRUE ~ NA_character_
    ),
    depth_zone = factor(
      depth_zone,
      levels = c(
        "Intertidal / exposed",
        "Very shallow subtidal",
        "Shallow subtidal",
        "Deeper subtidal"
      )
    )
  )


df_outcomes2 <- df_outcomes %>%
  mutate(
    depth_zone = case_when(
      bathy_mod < 0 ~ "Intertidal / exposed",
      bathy_mod >= 0 & bathy_mod < 2 ~ "Very shallow subtidal",
      bathy_mod >= 2 & bathy_mod < 5 ~ "Shallow subtidal",
      bathy_mod >= 5 ~ "Deeper subtidal",
      TRUE ~ NA_character_
    ),
    depth_zone = factor(
      depth_zone,
      levels = c(
        "Intertidal / exposed",
        "Very shallow subtidal",
        "Shallow subtidal",
        "Deeper subtidal"
      )
    ),
    substrate_depth = interaction(depth_zone, substrate_mode_obs, sep = " – ")
  )
df_outcomes2 <- sf::st_drop_geometry(df_outcomes2)

counts_cross <- df_outcomes2 %>%
  filter(!is.na(substrate_mode_obs), !is.na(depth_zone)) %>%
  count(substrate_depth) %>%
  rename(n_sites = n)

outcomes_cross <- df_outcomes2 %>%
  filter(!is.na(substrate_mode_obs), !is.na(depth_zone)) %>%
  count(model, substrate_depth, outcome) %>%
  group_by(model, substrate_depth) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  left_join(counts_cross, by = "substrate_depth")

outcomes_cross <- outcomes_cross %>%
  mutate(substrate_depth_lab = paste0(substrate_depth, "\n(n = ", n_sites, ")"))

ggplot(outcomes_cross, aes(x = substrate_depth_lab, y = prop, fill = outcome)) +
  geom_col() +
  facet_wrap(~ model) +
  coord_flip() +
  labs(
    x = "Depth zone × substrate type",
    y = "Proportion of validation sites",
    fill = "Prediction outcome",
    title = "Prediction outcomes by depth and substrate"
  ) +
  theme_minimal()


stress_test <- df_outcomes %>%
  mutate(
    sandy_site = substrate_mode_obs %in% c("Sand"),
    intertidal_or_exposed = avgCorDepth_obs < 0,
    stress_zone = case_when(
      intertidal_or_exposed & sandy_site ~ "Intertidal/exposed sand",
      intertidal_or_exposed & !sandy_site ~ "Intertidal/exposed non-sand",
      !intertidal_or_exposed & sandy_site ~ "Subtidal sand",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(model, stress_zone) %>%
  summarise(
    n = n(),
    observed_absences = sum(Presence == 0, na.rm = TRUE),
    false_positives = sum(Presence == 0 & prediction == 1, na.rm = TRUE),
    false_positive_rate = if_else(
      observed_absences == 0,
      NA_real_,
      false_positives / observed_absences
    ),
    predicted_prevalence = mean(prediction == 1, na.rm = TRUE),
    observed_prevalence = mean(Presence == 1, na.rm = TRUE),
    .groups = "drop"
  )

stress_test



binary_metrics_clean %>%
  select(model, sensitivity, specificity, TSS, false_positive_rate, false_negative_rate) %>%
  pivot_longer(
    cols = -model,
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(aes(x = model, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    x = "",
    y = "Metric value",
    title = "Validation metrics for unmasked and masked binary eelgrass maps"
  ) +
  theme_minimal()


ggplot(binary_metrics_clean, aes(x = model, y = TSS, fill = model)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "",
    y = "TSS",
    title = "TSS comparison between unmasked and masked binary eelgrass maps"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


binary_metrics_clean %>%
  select(
    model,
    n,
    predicted_prevalence,
    sensitivity,
    specificity,
    TSS,
    false_positive_rate,
    false_negative_rate,
    precision,
    npv,
    TP, TN, FP, FN
  )

stress_test











