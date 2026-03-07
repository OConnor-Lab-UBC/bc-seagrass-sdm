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
# models that does best across 4 forms of validation is sdmtmb bccm spatial

# tss from sdm tmb models was from full model, should get tss also from cv folds and the 

###############################################################################

library(terra)
library(dplyr)
library(purrr)
library(tibble)
library(caret)

load("code/output_data/model_results/eelgrass_eval_final_models.RData")
threshold_methods <- c("MaxSens+Spec", "MaxKappa", "MaxPCC")
model <- names(eval_results)
threshold_table <- do.call(rbind, lapply(model, function(q) {
  th_list <- eval_results[[q]]$thresholds
  th_list <- th_list[th_list$Method %in% threshold_methods, ]
  data.frame(
    Model = q,
    Threshold_Method = th_list$Method,
    Threshold = th_list$Predicted,
    stringsAsFactors = FALSE
  )
}))
# tss threshold for bccm spatial is 0.04, kappa is 0.22, PCC is 0.4
# from field validation 

# --- Load rasters ---
eelgrass_20m <- terra::vrt(c("raster/eelgrass/bccm_spatial/eelgrass_predictions_hg_bccm_spatial.tif", "raster/eelgrass/bccm_spatial/eelgrass_predictions_ncc_bccm_spatial.tif", "raster/eelgrass/bccm_spatial/eelgrass_predictions_qcs_bccm_spatial.tif", "raster/eelgrass/bccm_spatial/eelgrass_predictions_ss_bccm_spatial.tif", "raster/eelgrass/bccm_spatial/eelgrass_predictions_wcvi_bccm_spatial.tif"), "eelgrass_bccm_spatial.vrt", overwrite=T)   # values 0–1

#marsh_10m    <- rast("raw_data/XiuchengMaps/Pacific_cover_2024_Sentinel2.tif")      # classes 1–5
# 1-5: tidal marsh, tidal flats, open water, water-land interface, low-elevation uplands

# --- 1. Threshold eelgrass raster to binary (presence / absence) ---
# define threshold (adjust as needed)
thr <- 0.04

eel_bin <- ifel(eelgrass_20m >= thr, 1, 0)

# --- 2. Convert marsh raster to binary (classes 1–2 = marsh&tidalflats&upland) ---
# marsh_bin_10m <- ifel(marsh_10m %in% c(1,2,4,5), 1, 0)
# 
# # --- 3. Resample marsh raster to 20 m (to match eelgrass raster) ---
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


# --- 4. Remove eelgrass presence where marsh, tidal flats, upland and beach is present ---
eel_final <- ifel(marsh_bin_20m == 1, 0, eel_bin)

# remove eelgrass predictions in areas with rei > 0.20
rei_all <- terra::vrt(c("raw_data/REI/rei_hg.tif", "raw_data/REI/rei_ncc.tif", "raw_data/REI/rei_qcs.tif", "raw_data/REI/rei_sog.tif", "raw_data/REI/rei_wcvi.tif"), "rei.vrt", overwrite=T)
rei_aligned <- project(rei_all, eel_final, method="near")

eel_final <- ifel(rei_aligned > 0.2, 0, eel_final)


# remove areas with high exposure sand
zero_vec <- vect("raw_data/substrate_20m/exposed_shorezone_bottompatches_sand_class16_17_27_28_30_plusextra.shp")
zero_vec <- project(zero_vec, eel_final)
zero_ras <- rasterize(zero_vec, eel_final, field = 1, background = NA)

eel_final <- ifel(!is.na(zero_ras), 0, eel_final)

#remove areas with incorrect bathy based on updated CHS products
toodeep_vec <- vect("raw_data/substrate_20m/bathy_updated_toodeep_mask.shp")
toodeep_vec <- project(toodeep_vec, eel_final)
toodeep_ras <- rasterize(toodeep_vec, eel_final, field = 1, background = NA)

eel_final <- ifel(!is.na(toodeep_ras), 0, eel_final)



# Add in areas where eelgrass has been mapped. 
#load netforce eelgrass  data 2013-2023
eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))

# from looking previously there is no difference in the predictions if eelgrass has been observed more times so just change all to 1. A Spearman correlation of –0.0303 suggests that higher modelled probabilities do not correspond in any meaningful way to more years of observed presence.
values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1

eelgrass_indep <- resample(eelgrass_indep, eel_final, method = "near")

eel_final_plus_mapped <- ifel(eelgrass_indep == 1, 1, eel_final)



n1 <- global(eel_final_plus_mapped == 1, "sum", na.rm = TRUE)
n1

area_m2 <- n1 * 400
area_m2
# 784243600 m2

area_km2 <- area_m2 / 1e6
area_km2
# 784.2436

area_ha <- area_m2 / 10000
area_ha
# 78424.36 hectares

writeRaster(eel_final_plus_mapped, file.path("raster/eelgrass/eelgrass_predictions.tif"), overwrite = TRUE)

#eel_final_plus_mapped <- rast("raster/eelgrass/eelgrass_predictions.tif")   

# compare predictions to field validation dataset to see if overprediction into intertidal areas and exposed sandy areas has been adequetly dealt with
load("code/output_data/field_validation/validation_dataset.RData")

preds_extract <- terra::extract(eel_final_plus_mapped, validation_sf)
summary(preds_extract$eelgrass_bccm_spatial)
validation_sf$final_eelgrass_preds <- preds_extract$eelgrass_bccm_spatial

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
# Also in HKey 35 need to remove NT and UL as they were drift observations, and also remove PT observation as must be wrong, PC is 0 and there is no rock at site

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
                                           ZM, ZM_freq, PC_ZM, final_eelgrass_preds, eelgrass_bccm_spatial_se,
                                           PH, PH_freq, PC_PH, geometry) # need to add in surfgrass predictions
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

#Strong at ruling out eelgrass (NPV 0.94, spec 0.84)
#Moderately good at detecting eelgrass (sens 0.64)
#tends to overpredict presence

TN <- 355
FP <- 67
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
