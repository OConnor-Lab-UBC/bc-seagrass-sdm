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
load("./code/output_data/model_results/final_metrics_eelgrass.RData")
# for eelgras XGBOOST nep is best model. 
xgb_eelgrass<- final_metrics_eelgrass %>% filter (model == "XGBoost_nep") %>% select (threshold_spatial, threshold_field)
# do not want to make threshold based off independent because these are psudo absences and tempral becasue these maps are for current distirbution



# --- Load rasters ---
eelgrass_20m <- terra::vrt(c("raster/eelgrass/xgb_nep/eelgrass_predictions_hg_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ncc_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_qcs_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ss_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_wcvi_xgb_nep.tif"), "eelgrass_xgb_nep.vrt", overwrite=T)   # values 0–1
#eelgrass_20m_glmm_spatial_nep <- terra::vrt(c("raster/eelgrass/nep_spatial/eelgrass_predictions_hg_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_ncc_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_qcs_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_ss_nep_spatial.tif", "raster/eelgrass/nep_spatial/eelgrass_predictions_wcvi_nep_spatial.tif"), "eelgrass_nep_spatial.vrt", overwrite=T)   # values 0–1



# --- 1. Threshold eelgrass raster to binary (presence / absence) ---
# define threshold (adjust as needed)
thr<- mean(as.numeric(xgb_eelgrass[1, ]))


eel_bin <- ifel(eelgrass_20m >= thr, 1, 0)
#eel_bin_nepspatial <- ifel(eelgrass_20m_glmm_spatial_nep >= thr1, 1, 0)


n0 <- global(eel_bin == 1, "sum", na.rm = TRUE)
n0

area_m2_notmasked <- n0 * 400
area_m2_notmasked
# 1625356400 m2

area_km2_notmasked <- area_m2_notmasked / 1e6
area_km2_notmasked
# 1625.356 km2

area_ha_notmasked <- area_m2_notmasked / 10000
area_ha_notmasked
# 162535.6 hectares



writeRaster(eel_bin, file.path("raster/eelgrass/eelgrass_predictions_xgb_nep_binary_notmasked.tif"), overwrite = TRUE)
#writeRaster(eel_bin_nepspatial, file.path("raster/eelgrass/eelgrass_predictions_nepspatial_binary_notmasked.tif"), overwrite = TRUE)


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

eel_final[!is.na(substrate_aligned) & (substrate_aligned %in% c(1))] <- 0
writeRaster(eel_final, file.path("raster/eelgrass/eelgrass_predictions_nomapped.tif"), overwrite = TRUE)


# Add in areas where eelgrass has been mapped. 
#load netforce eelgrass  data 2013-2023
eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))
wrong_vec <- vect("code/output_data/independent_validation/netforce_mask.shp")

wrong_vec <- project(wrong_vec, eel_final)
wrong_ras <- rasterize(wrong_vec, eel_final, field = 1, background = NA)

values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1

eelgrass_indep <- resample(eelgrass_indep, eel_final, method = "near")

eel_final_plus_mapped <- ifel(eelgrass_indep == 1, 1, eel_final)

eel_final_plus_mapped <- ifel(!is.na(wrong_ras), 0, eel_final_plus_mapped)

n1 <- global(eel_final_plus_mapped == 1, "sum", na.rm = TRUE)
n1

area_m2 <- n1 * 400
area_m2
# 699017600 m2

area_km2 <- area_m2 / 1e6
area_km2
# 699.0176 km2

area_ha <- area_m2 / 10000
area_ha
# 69901.76 hectares

writeRaster(eel_final_plus_mapped, file.path("raster/eelgrass/eelgrass_predictions.tif"), overwrite = TRUE)



# surfgrass
load("./code/output_data/model_results/final_metrics_surfgrass.RData")
xgb_surfgrass<- final_metrics_surfgrass %>% filter (model == "XGBoost_nep") %>% select (threshold_spatial, threshold_field)

# --- Load rasters ---
surfgrass_20m <- terra::vrt(c("raster/surfgrass/xgb_nep/surfgrass_predictions_hg_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_ncc_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_qcs_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_ss_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_wcvi_xgb_nep.tif"), "surfgrass_xgb_nep.vrt", overwrite=T)   # values 0–1

# --- 1. Threshold surfgrass raster to binary (presence / absence) ---
# define threshold (adjust as needed)
thr<- mean(as.numeric(xgb_surfgrass[1, ]))

surf_bin <- ifel(surfgrass_20m >= thr, 1, 0)

n0 <- global(surf_bin == 1, "sum", na.rm = TRUE)
n0

area_m2_notmasked <- n0 * 400
area_m2_notmasked
# 947055600 m2

area_km2_notmasked <- area_m2_notmasked / 1e6
area_km2_notmasked
# 947.0556 km2

area_ha_notmasked <- area_m2_notmasked / 10000
area_ha_notmasked
# 94705.56 hectares

#surf_bin_nepnospatial <- ifel(surfgrass_20m_glmm_nospatial_nep >= thr1, 1, 0)
writeRaster(surf_bin, file.path("raster/surfgrass/surfgrass_predictions_xgb_nep_binary_notmasked.tif"), overwrite = TRUE)
#writeRaster(surf_bin_nepnospatial, file.path("raster/surfgrass/surfgrass_predictions_nepnospatial_binary_notmasked.tif"), overwrite = TRUE)


# --- 4. Remove surfgrass presence where marsh, tidal flats, upland and beach is present ---
#may not do this becasue intertidal overprediction was not a problem with surfgrass
#surf_final <- ifel(marsh_bin_20m == 1, 0, surf_bin)
surf_final <- surf_bin
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
# 813831600 m2

area_km2 <- area_m2 / 1e6
area_km2
# 813.8316 km2

area_ha <- area_m2 / 10000
area_ha
# 81383.16

writeRaster(surf_final, file.path("raster/surfgrass/surfgrass_predictions.tif"), overwrite = TRUE)









load("code/output_data/model_results/final_metrics_eelgrass_fieldnointertidal.RData")
load("code/output_data/model_results/final_metrics_eelgrass.RData")

final_metrics_eelgrass <- final_metrics_eelgrass %>% 
  select (model, auc_field, tjur_field, brier_field, logloss_field, tss_field, sensitivity_field, specificity_field) %>%
  filter(model == "XGBoost_nep") %>%
  mutate(
    model = recode(model,
                   "XGBoost_nep" = "XGBoost_nep (continuous & unrefined)")
  )
final_metrics_eelgrass_nointertidal <- final_metrics_eelgrass_nointertidal %>% 
  select (model, auc_field, tjur_field, brier_field, logloss_field, tss_field, sensitivity_field, specificity_field) %>%
  filter(model == "XGBoost_nep") %>%
  mutate(
    model = recode(model,
                   "XGBoost_nep" = "XGBoost_nep subtidal only (continuous & unrefined)"))

field_metrics_eelgrass <- rbind(final_metrics_eelgrass, final_metrics_eelgrass_nointertidal)

field_metrics_eelgrass <- field_metrics_eelgrass %>%
  rename(AUC = auc_field,
         TjurR2 = tjur_field,
         Brier = brier_field,
         Logloss = logloss_field,
         TSS = tss_field,
         Sensitivity = sensitivity_field,
         Specificity = specificity_field)



eel_binary_refined<- rast("raster/eelgrass/eelgrass_predictions.tif")
eel_binary_notrefined<- rast("raster/eelgrass/eelgrass_predictions_xgb_nep_binary_notmasked.tif")

# compare predictions to field validation dataset to see if overprediction into intertidal areas and exposed sandy areas has been adequately dealt with
load("code/output_data/field_validation/validation_dataset.RData")

validation_sf <- validation_sf %>%
  filter(!HKey %in% c(74, 84, 85, 161))

refined_extract <- terra::extract(eel_binary_refined, validation_sf)
notrefined_extract <- terra::extract(eel_binary_notrefined, validation_sf)

validation_sf$pred_refined_binary <- refined_extract[[2]]
validation_sf$pred_notrefined_binary <- notrefined_extract[[2]]

# ---- Prep df ----
df <- validation_sf %>%
  mutate(
    Presence = as.integer(ZM),
    pred_refined_binary = as.integer(pred_refined_binary),
    pred_notrefined_binary = as.integer(pred_notrefined_binary)
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
  "pred_refined_binary", "pred_notrefined_binary"
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
    sensitivity,
    specificity,
    TSS
  )

binary_metrics_clean <- binary_metrics_clean %>%
  rename(Sensitivity = sensitivity,
         Specificity = specificity)

field_metrics_eelgrass <- bind_rows(field_metrics_eelgrass, binary_metrics_clean)

field_metrics_eelgrass <- field_metrics_eelgrass %>%
  mutate(
  model = recode(model,
                 "pred_refined_binary" = "XGBoost_nep (binary & refined)",
                 "pred_notrefined_binary" = "XGBoost_nep (binary & unrefined)"),
  Species = "Eelgrass")






load("code/output_data/model_results/final_metrics_surfgrass.RData")
load("code/output_data/model_results/final_metrics_surfgrass_fieldnointertidal.RData")


final_metrics_surfgrass <- final_metrics_surfgrass %>% 
  select (model, auc_field, tjur_field, brier_field, logloss_field, tss_field, sensitivity_field, specificity_field) %>%
  filter(model == "XGBoost_nep") %>%
  mutate(
    model = recode(model,
                   "XGBoost_nep" = "XGBoost_nep (continuous & unrefined)")
  )

final_metrics_surfgrass <- final_metrics_surfgrass %>%
  rename(AUC = auc_field,
         TjurR2 = tjur_field,
         Brier = brier_field,
         Logloss = logloss_field,
         TSS = tss_field,
         Sensitivity = sensitivity_field,
         Specificity = specificity_field)

final_metrics_surfgrass_nointertidal <- final_metrics_surfgrass_nointertidal %>% 
  select (model, auc_field, tjur_field, brier_field, logloss_field, tss_field, sensitivity_field, specificity_field) %>%
  filter(model == "XGBoost_nep") %>%
  mutate(
    model = recode(model,
                   "XGBoost_nep" = "XGBoost_nep subtidal only (continuous & unrefined)"))

final_metrics_surfgrass_nointertidal <- final_metrics_surfgrass_nointertidal %>%
  rename(AUC = auc_field,
         TjurR2 = tjur_field,
         Brier = brier_field,
         Logloss = logloss_field,
         TSS = tss_field,
         Sensitivity = sensitivity_field,
         Specificity = specificity_field)

field_metrics_surfgrass <- rbind(final_metrics_surfgrass, final_metrics_surfgrass_nointertidal)

surf_binary_refined<- rast("raster/surfgrass/surfgrass_predictions.tif")
surf_binary_notrefined<- rast("raster/surfgrass/surfgrass_predictions_xgb_nep_binary_notmasked.tif")

refined_extract_surf <- terra::extract(surf_binary_refined, validation_sf)
notrefined_extract_surf <- terra::extract(surf_binary_notrefined, validation_sf)

validation_sf$pred_refined_binary_surf <- refined_extract_surf[[2]]
validation_sf$pred_notrefined_binary_surf <- notrefined_extract_surf[[2]]

# ---- Prep df ----
df_surf <- validation_sf %>%
  mutate(
    Presence = as.integer(PH),
    pred_refined_binary_surf = as.integer(pred_refined_binary_surf),
    pred_notrefined_binary_surf = as.integer(pred_notrefined_binary_surf)
  )


binary_pred_cols_surf <- c(
  "pred_refined_binary_surf", "pred_notrefined_binary_surf"
)

binary_metrics_surf <- map_dfr(binary_pred_cols_surf, function(p) {
  binary_validation_metrics(
    data = df_surf,
    obs_col = "Presence",
    pred_col = p,
    model_name = p
  )
})

binary_metrics_surf


binary_metrics_clean_surf <- binary_metrics_surf %>%
  select(
    model,
    sensitivity,
    specificity,
    TSS
  )

binary_metrics_clean_surf <- binary_metrics_clean_surf %>%
  rename(Sensitivity = sensitivity,
         Specificity = specificity)

field_metrics_surfgrass <- bind_rows(field_metrics_surfgrass, binary_metrics_clean_surf)

field_metrics_surfgrass <- field_metrics_surfgrass %>%
  mutate(
    model = recode(model,
                   "pred_refined_binary_surf" = "XGBoost_nep (binary & refined)",
                   "pred_notrefined_binary_surf" = "XGBoost_nep (binary & unrefined)"),
    Species = "Surfgrass")


field_metrics_seagrass<- bind_rows(field_metrics_eelgrass, field_metrics_surfgrass)


library(tidyverse)

metrics_long <- field_metrics_seagrass %>%
  pivot_longer(
    cols = AUC:Specificity,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c(
        "AUC",
        "TjurR2",
        "Brier",
        "Logloss",
        "TSS",
        "Sensitivity",
        "Specificity"
      )
    ),
    model = recode(
      model,
      "XGBoost_nep (continuous & unrefined)" = "Continuous",
      "XGBoost_nep subtidal only (continuous & unrefined)" = "Continuous (subtidal)",
      "XGBoost_nep (binary & unrefined)" = "Binary",
      "XGBoost_nep (binary & refined)" = "Binary (ecologically refined)"
    ),
    model = factor(
      model,
      levels = c(
        "Continuous",
        "Continuous (subtidal)",
        "Binary",
        "Binary (ecologically refined)"
      )
    )
  )

metrics_long <- metrics_long %>%
  mutate(
    Value_plot = case_when(
      Metric == "Brier" ~ 1 - Value,
      Metric == "Logloss" ~ 1 - Value,
      TRUE ~ Value
    ),
    Metric_plot = case_when(
      Metric == "Brier" ~ "1 - Brier",
      Metric == "Logloss" ~ "1 - Log loss",
      Metric == "TjurR2" ~ "Tjur R²",
      TRUE ~ Metric
    )
  )

metrics_long <- metrics_long %>%
  mutate(
    Metric_plot = factor(
      Metric_plot,
      levels = c(
        "AUC",
        "Tjur R²",
        "1 - Brier",
        "1 - Log loss",
        "TSS",
        "Sensitivity",
        "Specificity"
      )
    )
  )

eelgrass_plot <- metrics_long %>%
  filter(Species == "Eelgrass") %>%
  ggplot(aes(x = Value_plot,
             y = Metric_plot,
             colour = model)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.3),
             na.rm = TRUE) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  )  +
  labs(
    title = "Eelgrass",
    x = "Metric value",
    y = NULL,
    colour = NULL
  ) +
  scale_colour_manual(
    values = c(
      "Continuous" = "#0072B2",
      "Continuous (subtidal)" = "#56B4E9",
      "Binary" = "#E69F00",
      "Binary (ecologically refined)" = "#009E73"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    legend.position = "none"
  )
eelgrass_plot

surfgrass_plot <- metrics_long %>%
  filter(Species == "Surfgrass") %>%
  ggplot(aes(x = Value_plot,
             y = Metric_plot,
             colour = model)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.3),
             na.rm = TRUE) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  )  +
  labs(
    title = "Surfgrass",
    x = "Metric value",
    y = NULL,
    colour = NULL
  ) +
  scale_colour_manual(
    values = c(
      "Continuous" = "#0072B2",
      "Continuous (subtidal)" = "#56B4E9",
      "Binary" = "#E69F00",
      "Binary (ecologically refined)" = "#009E73"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    legend.position = "none"
  )
surfgrass_plot

library(cowplot)

legend <- get_legend(
  eelgrass_plot +
    theme(
      legend.position = "bottom"
    )
)

# Remove legends from both plots
eelgrass_nolegend <- eelgrass_plot +
  theme(legend.position = "none")

surfgrass_nolegend <- surfgrass_plot +
  theme(legend.position = "none")

# Combine plots + legend
combined_plot <- (
  eelgrass_nolegend / surfgrass_nolegend
) /
  wrap_elements(legend) +
  plot_layout(
    heights = c(3, 3, 0.3)
  )

combined_plot

ggsave(
  "figures/field_validation_alternative_metrics.png",
  combined_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
