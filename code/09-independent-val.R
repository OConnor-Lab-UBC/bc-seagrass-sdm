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
library(ecospat)
library(Metrics)
library(reshape2)

#### Eelgrass ####

#load netforce eelgrass  data 2013-2023
eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))

# from looking previously there is no difference in the predictions if eelgrass has been observered more times so just change all to 1. A Spearman correlation of –0.0303 suggests that higher modelled probabilities do not correspond in any meaningful way to more years of observed presence.
values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1
names(eelgrass_indep)<-"obs"

# load model thresholds created
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

# get predictions
parent_dir <- "raster/eelgrass"
folders <- list.dirs(parent_dir, recursive = FALSE)

# load and mosaic eelgrass SDM rasters from 4 different SDM predictions
mosaics <- list()

for (f in folders) {
  tif_files <- list.files(
    f,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  tif_files <- tif_files[!grepl("se", basename(tif_files))]
  r_list <- lapply(tif_files, rast)
  m <- do.call(mosaic, r_list)
  mosaics[[basename(f)]] <- m
}

#load substrate layer
substrate_all <- terra::vrt(c("raw_data/substrate_20m/updated/hg_20m.tif", "raw_data/substrate_20m/updated/ncc_20m.tif", "raw_data/substrate_20m/updated/qcs_20m.tif", "raw_data/substrate_20m/updated/sog_20m.tif", "raw_data/substrate_20m/updated/wcvi_20m.tif"), "substrate.vrt", overwrite=T)
names(substrate_all)<-"substrate"
crs(substrate_all) <- "EPSG:3005"


# Ensure rasters align, just take model predictions where there is observed data
# we can only compare in areas where eelgrass has been observed as there is no absence data from independent dataset

eelgrass_sdm_resampled <- lapply(mosaics, function(r) {
  if (!compareGeom(r, eelgrass_indep, crs = TRUE, stopOnError = FALSE)) {
    r <- project(r, eelgrass_indep)  }
  crop(
    resample(r, eelgrass_indep, method = "bilinear"),
    eelgrass_indep)})
plot(eelgrass_sdm_resampled[[1]])
observed_cells <- which(values(eelgrass_indep) > 0)
coords <- xyFromCell(eelgrass_indep, observed_cells)

eelgrass_stack <- c(eelgrass_sdm_resampled[[1]], eelgrass_sdm_resampled[[2]], eelgrass_sdm_resampled[[3]], eelgrass_sdm_resampled[[4]])
names(eelgrass_stack) <- c("bccm_nospatial", "bccm_spatial", "nep_nospatial", "nep_spatial")

r_points <- as.points(eelgrass_indep, na.rm = TRUE)
eelgrass_sf<- r_points %>% st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

prediction_extract <- terra::extract(eelgrass_stack, eelgrass_sf)
eelgrass_sf <- eelgrass_sf %>% bind_cols(prediction_extract)

eelgrass_sf <- eelgrass_sf %>%
  filter(if_all(everything(), ~ !is.na(.)))

#extract substrate data to see if that is causing mismatch
substrate_extract <- terra::extract(substrate_all, eelgrass_sf)
eelgrass_sf$substrate <- substrate_extract$substrate
eelgrass_sf$substrate <- c("Rock", "Mixed", "Sand", "Mud")[eelgrass_sf$substrate]

eelgrass_sf$rock_group <- ifelse(eelgrass_sf$substrate == "Rock", "Rock", "Not Rock")
# should remove areas with modelled rock substrate from comparison as the substrate model is likely wrong in these areas



####Compare between models####
models <- c("bccm_nospatial", "bccm_spatial", "nep_nospatial", "nep_spatial")

# Initialize results dataframe
results <- data.frame(
  Model = models,
  RMSE = NA,
  MAE = NA,
  Mean_Prediction = NA,
  Prop_Above_0.1 = NA,
  Prop_Above_0.3 = NA,
  Prop_Above_0.5 = NA
)

# Thresholds for proportion calculations
thresholds <- c(0.1, 0.3, 0.5)

# Calculate metrics for each model
for (i in seq_along(models)) {
  pred <- eelgrass_sf[[models[i]]]
  obs <- eelgrass_sf$obs  # should be 1 for all rows

  results$RMSE[i] <- rmse(obs, pred)
  results$MAE[i]  <- mae(obs, pred)
  results$Mean_Prediction[i] <- mean(pred)
  results$Prop_Above_0.1[i] <- mean(pred >= thresholds[1])
  results$Prop_Above_0.3[i] <- mean(pred >= thresholds[2])
  results$Prop_Above_0.5[i] <- mean(pred >= thresholds[3])
}

# Rank models by RMSE (lower is better)
results <- results %>%
  arrange(RMSE)

# Print results
print(results)
#Lower MAE = better model
#Lower RMSE = better model
#Mean prediciton, higher is better
#rank proportions, higher in each proportion better
#nep no spatial model comes out as best, but pretty marginal


# make a plot of predictions for each model
eelgrass_df <- sf::st_drop_geometry(eelgrass_sf) %>% select(-ID, -substrate)

df_melt <- melt(eelgrass_df, id.vars = c("obs", "rock_group"), variable.name = "Model", value.name = "Prediction")

ggplot(df_melt, aes(x = Model, y = Prediction, fill = rock_group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  ylab("Predicted Probability") +
  ggtitle("Comparison of Eelgrass Model Predictions by Substrate") +
  labs(fill = "Substrate") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )
#nep no spatial comes out as best

# false negatives
false_negatives_by_substrate <- lapply(models, function(m) {
  thr <- threshold[m]

  eelgrass_df %>%
    mutate(
      pred = as.numeric(.data[[m]]),
      false_negative = pred < thr
    ) %>%
    filter(!is.na(pred)) %>%
    group_by(rock_group) %>%
    summarise(
      Model = m,
      FN_Rate = mean(false_negative, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}) %>%
  bind_rows()

false_negatives_by_substrate

# Plot: FN rates by model and substrate
ggplot(false_negatives_by_substrate,
       aes(x = Model, y = FN_Rate, fill = rock_group)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Rock" = "#D95F02", "Not Rock" = "#1B9E77")) +
  labs(
    y = "False Negative Rate (Missed Observed Eelgrass)",
    x = "Model",
    fill = "Substrate",
    title = "False Negative Rates by Model and Substrate"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )
# none of the models preform any better when looking at false negatives


results_thresholds <- do.call(rbind, lapply(1:nrow(threshold_table), function(i) {
  row <- threshold_table[i, ]
  m <- row$Model
  thr <- row$Threshold
  method <- row$Threshold_Method

  eelgrass_df %>%
    mutate(
      pred = as.numeric(.data[[m]]),
      false_negative = pred < thr
    ) %>%
    filter(!is.na(pred)) %>%
    group_by(rock_group) %>%
    summarise(
      Model = m,
      Threshold_Method = method,
      FN_Rate = mean(false_negative, na.rm = TRUE),
      Pred_Suitable = mean(pred >= thr, na.rm = TRUE), # proportion predicted suitable
      n = n(),
      .groups = "drop"
    )
})) %>% bind_rows()

ggplot(results_thresholds,
       aes(x = Threshold_Method, y = FN_Rate, fill = rock_group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Model) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Rock" = "#D95F02", "Not Rock" = "#1B9E77")) +
  labs(
    y = "False Negative Rate",
    x = "Threshold Method",
    fill = "Substrate",
    title = "False Negative Rates by Threshold Method and Model"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(results_thresholds %>% filter(rock_group == "Not Rock"),
       aes(x = Pred_Suitable, y = FN_Rate, color = Threshold_Method, label = Model)) +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.01) +
  labs(
    x = "Proportion Predicted Suitable",
    y = "False Negative Rate (Not Rock)",
    color = "Threshold Method",
    title = "Model Performance Across Threshold Methods"
  ) +
  theme_bw()

#area retained
cell_area <- 20 * 20  # 400 m²

area_retained <- do.call(rbind, lapply(seq_len(nrow(threshold_table)), function(i) {
  
  m      <- threshold_table$Model[i]
  thr    <- threshold_table$Threshold[i]
  method <- threshold_table$Threshold_Method[i]
  
  pred_bin <- ifelse(eelgrass_df[[m]] >= thr, 1, 0)
  
  retained_cells <- eelgrass_df$obs == 1 & pred_bin == 1
  obs_cells <- eelgrass_df$obs == 1
  
  retained_area_m2 <- sum(retained_cells, na.rm = TRUE) * cell_area
  obs_area_m2 <- sum(obs_cells, na.rm = TRUE) * cell_area
  
  data.frame(
    Model = m,
    Threshold_Method = method,
    Threshold = thr,
    Retained_Area_m2 = retained_area_m2,
    Obs_Area_m2 = obs_area_m2,
    Percent_Retained = retained_area_m2 / obs_area_m2 * 100
  )
}))

ggplot(area_retained,
       aes(x = Threshold_Method, y = Percent_Retained, fill = Model)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    y = "Observed Eelgrass Area Retained After Thresholding",
    x = "Threshold Method",
    fill = "Model",
    title = "Retention of Observed Eelgrass by Model and Threshold"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Within a given threshold method, models are very similar, but choice of threshold really affects results
#Model choice had minimal influence on the proportion of observed eelgrass retained after thresholding, with differences between models consistently <5% for a given thresholding method. In contrast, threshold selection strongly influenced retained observed area, with differences exceeding 50% within individual models.

ggplot(area_retained,
       aes(x = Threshold_Method,
           y = Percent_Retained,
           group = Model,
           color = Model)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 3) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    x = "Threshold Method",
    y = "Observed Eelgrass Area Retained (%)",
    title = "Threshold Choice Dominates Model Differences",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )






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





models <- names(eval_results)
threshold_methods <- c("MaxSens+Spec", "MaxKappa", "MaxPCCC")
cell_area <- 20 * 20  # 400 m²


# 3️⃣ Drop geometry for calculations
eelgrass_df <- sf::st_drop_geometry(eelgrass_sf)

# 4️⃣ Compute RMSE, MAE, mean prediction, and %TP for each model x threshold
summary_metrics <- do.call(rbind, lapply(1:nrow(threshold_table), function(i) {
  row <- threshold_table[i, ]
  m <- row$Model
  thr <- row$Threshold
  method <- row$Threshold_Method
  
  pred <- eelgrass_df[[m]]
  obs <- eelgrass_df$obs
  
  obs_bin <- ifelse(obs > 0, 1, 0)
  pred_bin <- ifelse(pred >= thr, 1, 0)
  TP <- obs_bin == 1 & pred_bin == 1
  FN_rate <- 1 - mean(pred_bin, na.rm = TRUE)  # since predictions only where obs exist
  
  data.frame(
    Model = m,
    Threshold_Method = method,
    Threshold = thr,
    RMSE = rmse(obs, pred),
    MAE = mae(obs, pred),
    Mean_Prediction = mean(pred, na.rm = TRUE),
    Percent_TP = mean(TP, na.rm = TRUE) * 100,
    FN_rate = FN_rate,
    stringsAsFactors = FALSE
  )
}))

summary_metrics


substrate_metrics <- do.call(rbind, lapply(1:nrow(threshold_table), function(i) {
  row <- threshold_table[i, ]
  m <- row$Model
  thr <- row$Threshold
  method <- row$Threshold_Method
  
  eelgrass_df %>%
    mutate(pred_bin = ifelse(.data[[m]] >= thr, 1, 0),
           obs_bin = ifelse(obs > 0, 1, 0)) %>%
    group_by(rock_group) %>%
    summarise(
      FN_rate = mean(obs_bin == 1 & pred_bin == 0, na.rm = TRUE),
      Percent_TP = mean(obs_bin == 1 & pred_bin == 1, na.rm = TRUE) * 100,
      n = n(),
      Model = m,
      Threshold_Method = method,
      .groups = "drop"
    )
})) 
substrate_metrics

ggplot(substrate_metrics, aes(x = Threshold_Method, y = Percent_TP, fill = rock_group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Model) +
  labs(
    y = "Percent Observed Eelgrass Captured",
    x = "Threshold Method",
    fill = "Substrate",
    title = "Observed Eelgrass Coverage by Model, Threshold, and Substrate"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

