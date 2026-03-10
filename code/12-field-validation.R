###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# compare sdm predictions to sdm validations to choose best model
#
###############################################################################

#when looking at subtidal and intertidal data together nep spatial, nep no spatial, and bccm spatial are all comparable
# when looking at subtidal nep spatial and bccm spatial are the best preforming and very similar

#functions
source("code/modelling-functions.R")

# Load packages
library(sf)
library(tidyverse)
library(terra)
library(patchwork)
library(effsize)
library(ggpubr)
library(pROC)
library(forcats)
library(caret)
library(purrr)

#themes for figures
boxed_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title.y = element_text(margin = margin(r = 8))
  )

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
# NEP spatial TSS threshold is 0.04

# load sdm validation dataset
load("code/output_data/field_validation/validation_dataset.RData")
validation_sf <- validation_sf %>%
  relocate(PH_freq, PH, ZM_freq, ZM, .after = PC_PH)

# need to also drop sites that after reviewing the database should not have been included as they were aborted or vis was too bad to see species
validation_sf <- validation_sf %>%
  filter(!HKey %in% c(74, 84, 85, 161))

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
                                           ZM, ZM_freq, PC_ZM, eelgrass_bccm_nospatial_pred, eelgrass_bccm_nospatial_se, eelgrass_bccm_spatial_pred, eelgrass_bccm_spatial_se, eelgrass_nep_nospatial_pred, eelgrass_nep_nospatial_pred, eelgrass_nep_spatial_pred, eelgrass_nep_spatial_se,
                                           PH, PH_freq, PC_PH, geometry) # need to add in surfgrass predictions
summary(validation_sf)
# substrate update diff false is 130 and substrate diff is 175 false (shows that update with shorezone substrate data makes substrate more accurate)
# at least 26% of sites had substrate mismatch between modelled and obs during field validation study

#not going to remove sites at this time as this issue would persist across the whole province. good to know that we are having a substrate mismatch >25% of the time
# remove sites that we know will have wrong predictions becasue the environmental layers do not match reality
# validation_sf <- validation_sf %>%
#   filter(
#     (Depth_diff < 5 | is.na(Depth_diff)),
#     (Slope_diff < 20 | is.na(Slope_diff)), # this was just a guess, need to add something more thought out
#     #(Substrate_diff == TRUE | is.na(Substrate_diff)), #leaving this out for now as that means more get eliminated for zm and ph then necessary
#     (Substrate_diff == TRUE | is.na(Substrate_diff))
#   )
# 334 passed threshold 
# 3 sites >20 degrees slope difference, 40 sites has >5m depth difference, 112 sites had substrate mismatch in updated substrate


##ZM
#wilcox
zm_pa_results <- list(
  e_bccm_nospatial = evaluate_occurrence(validation_sf, observed_col = ZM, predicted_col = eelgrass_bccm_nospatial_pred),
  e_bccm_spatial = evaluate_occurrence(validation_sf, observed_col = ZM, predicted_col = eelgrass_bccm_spatial_pred),
  e_nep_nospatial = evaluate_occurrence(validation_sf, observed_col = ZM, predicted_col = eelgrass_nep_nospatial_pred),
  e_nep_spatial = evaluate_occurrence(validation_sf, observed_col = ZM, predicted_col = eelgrass_nep_spatial_pred)
)

model_comparison <- imap_dfr(zm_pa_results, function(res, model) {
  tibble(
    model = model,
    wilcox_p     = res$wilcox_test$p.value,
    cliffs_delta = res$cliffs_delta$estimate#,
    #median_pres  = res$summary_stats
   )
})

model_comparison %>%
  arrange(
    desc(abs(cliffs_delta)),   # strongest separation
    wilcox_p                   # strongest significance
  )
#NEP model no spatial and then nep spatial come out as best but not much difference to bccm spatial. bccm no spatial is not significant

#spearmans
df <- validation_sf %>%
  mutate(
    # numeric midpoints for correlation
    PC_mid = case_when(
      PC_ZM == "0" ~ 0,
      PC_ZM == "1-25" ~ 12.5,
      PC_ZM == "26-50" ~ 38,
      PC_ZM == "51-75" ~ 63,
      PC_ZM == "76-100" ~ 88
    ),
    # collapse to presence/absence for binary metrics
    Presence = if_else(PC_ZM == "0", 0, 1)
  )

pred_cols <- c(
  "eelgrass_bccm_nospatial_pred",
  "eelgrass_bccm_spatial_pred",
  "eelgrass_nep_nospatial_pred",
  "eelgrass_nep_spatial_pred"
)

obs_pc   <- "PC_mid"
obs_pa   <- "ZM"
threshold_default <- 0.05 # bccm spatial preds TSS threshold

cor_results <- map_dfr(pred_cols, function(p) {
  
  spearman <- cor.test(df[[p]], df[[obs_pc]], method = "spearman")
  kendall  <- cor.test(df[[p]], df[[obs_pc]], method = "kendall")
  
  tibble(
    model        = p,
    spearman_rho = unname(spearman$estimate),
    spearman_p   = spearman$p.value,
    kendall_tau  = unname(kendall$estimate),
    kendall_p    = kendall$p.value
  )
})

cor_results
#model nep nospatial shows the strongest monotonic relationship but bccm spatial and nep spatial similar. below is rho and pval for bccm spatial
rho <- 0.103
pval <- 0.0241

roc_results <- map_dfr(pred_cols, function(p) {
  
  roc_obj <- roc(df[[obs_pa]], df[[p]], quiet = TRUE)
  
  coords_df <- coords(
    roc_obj,
    x = "all",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  ) %>%
    mutate(TSS = sensitivity + specificity - 1)
  
  best <- coords_df %>% slice_max(TSS, n = 1)
  
  tibble(
    model         = p,
    auc           = as.numeric(auc(roc_obj)),
    tss_threshold = best$threshold,
    max_tss       = best$TSS
  )
})

roc_results

#threshold for bccm spatial
tss_threshold <- 0.0195


# nep spatial still best but AUC and TSS is horrible 0.596 and 0.203. BCCM spatial is 0.589 and 0.214
# very little difference between all 4 models. nep spatial has and oddly low tss threshold

confusion_results <- map(pred_cols, function(p) {
  
  thr <- roc_results %>%
    filter(model == p) %>%
    pull(tss_threshold)
  
  df_bin <- df %>%
    mutate(Pred_binary = if_else(.data[[p]] >= thr, 1, 0))
  
  confusionMatrix(
    factor(df_bin$Pred_binary),
    factor(df_bin[[obs_pa]]),
    positive = "1"
  )
})

names(confusion_results) <- pred_cols
confusion_results$eelgrass_bccm_spatial_pred
# this model is very sensitive, so catches all true presences but most predicted 1s are wrong. Terrible specificity, fails at identifying absences. 
# what about if we exclude intertidal sites becasue model not built with intertidal data


df_nointertidal <- df %>% filter(Survey != "Intertidal")
validation_sf_nointertidal <- validation_sf %>% filter(Survey != "Intertidal")


#wilcox
zm_pa_results_no_intertidal <- list(
  e_bccm_nospatial = evaluate_occurrence(validation_sf_nointertidal, observed_col = ZM, predicted_col = eelgrass_bccm_nospatial_pred),
  e_bccm_spatial = evaluate_occurrence(validation_sf_nointertidal, observed_col = ZM, predicted_col = eelgrass_bccm_spatial_pred),
  e_nep_nospatial = evaluate_occurrence(validation_sf_nointertidal, observed_col = ZM, predicted_col = eelgrass_nep_nospatial_pred),
  e_nep_spatial = evaluate_occurrence(validation_sf_nointertidal, observed_col = ZM, predicted_col = eelgrass_nep_spatial_pred)
)

model_comparison_no_intertidal <- imap_dfr(zm_pa_results_no_intertidal, function(res, model) {
  tibble(
    model = model,
    wilcox_p     = res$wilcox_test$p.value,
    cliffs_delta = res$cliffs_delta$estimate#,
    #median_pres  = res$summary_stats
  )
})

model_comparison_no_intertidal %>%
  arrange(
    desc(abs(cliffs_delta)),   # strongest separation
    wilcox_p                   # strongest significance
  )
#NEP model spatial and bccm spatial come out the same and good 

cor_nointertidal <- map_dfr(pred_cols, function(p) {
  
  spearman <- cor.test(df_nointertidal[[p]], df_nointertidal[[obs_pc]], method = "spearman")
  kendall  <- cor.test(df_nointertidal[[p]], df_nointertidal[[obs_pc]], method = "kendall")
  
  tibble(
    model        = p,
    spearman_rho = unname(spearman$estimate),
    spearman_p   = spearman$p.value,
    kendall_tau  = unname(kendall$estimate),
    kendall_p    = kendall$p.value
  )
})

cor_nointertidal
# nep spatial and bccm spatial best and comparable

rho_nointertidal <- 0.421 # this is for bccm spatial
pval_nointertidal <- 0.0000001 


roc_nointertidal <- map_dfr(pred_cols, function(p) {
  
  roc_obj <- roc(df_nointertidal[[obs_pa]], df_nointertidal[[p]], quiet = TRUE)
  
  coords_df <- coords(
    roc_obj,
    x = "all",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  ) %>%
    mutate(TSS = sensitivity + specificity - 1)
  
  best <- coords_df %>% slice_max(TSS, n = 1)
  
  tibble(
    model         = p,
    auc           = as.numeric(auc(roc_obj)),
    tss_threshold = best$threshold,
    max_tss       = best$TSS
  )
})


roc_nointertidal
# bccm spatial and nep spatial best AUC and TSS is but only BCCM spatial tss threshold is reasonable, unlike nep spatial threshold is too low
# bccm spatial auc = 0.855, TSS= 0.645

tss_threshold_nointertidal <- 0.0203 # this is for bccm spatial though



validation_sf <- validation_sf %>%
  mutate(
    Tidal_zone = if_else(Survey == "Intertidal", "Intertidal", "Subtidal"),
    Tidal_zone = factor(Tidal_zone, levels = c("Subtidal", "Intertidal"))
  )


presence_df <- validation_sf %>%
  mutate(
    Presence = factor(ifelse(ZM == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  )

abundance_df <- validation_sf %>%
  filter(PC_ZM %in% c("1-25", "26-50", "51-75", "76-100")) %>% 
  mutate(
    PC_ZM = factor(PC_ZM, levels = c("1-25", "26-50", "51-75", "76-100"))
    )

presence_counts <- presence_df %>%
  group_by(Presence, Tidal_zone) %>%
  summarise(n = n(), .groups = "drop")

abundance_counts <- abundance_df %>%
  group_by(PC_ZM, Tidal_zone) %>%
  summarise(n = n(), .groups = "drop")

wilcox_presence <- presence_df %>%
  group_by(Presence) %>%
  summarise(
    p = wilcox.test(
      eelgrass_bccm_spatial_pred ~ Tidal_zone
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  )

threshold_df <- data.frame(
  threshold = c(threshold_default, tss_threshold),
  type = c("Model threshold", "Field validated threshold")
)

# Determine overall y-axis limits
y_min <- min(c(presence_df$eelgrass_bccm_spatial_pred, abundance_df$eelgrass_bccm_spatial_pred), na.rm = TRUE)
y_max <- max(c(presence_df$eelgrass_bccm_spatial_pred, abundance_df$eelgrass_bccm_spatial_pred), na.rm = TRUE)

# Offset for stars and sample sizes
y_star_offset <- 0.01 * y_max
y_n_label <- y_max * 1.08  # for both plots


p_presence <- ggplot(
  presence_df,
  aes(
    x = Presence,
    y = eelgrass_bccm_spatial_pred,
    fill = Tidal_zone
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.alpha = 0.5
  ) +
  geom_hline(
    data = threshold_df,
    aes(yintercept = threshold, linetype = type),
    colour = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = presence_counts,
    aes(
      x = Presence,
      y = y_max * 1.05,
      label = paste0("n=", n)
    ),
    position = position_dodge2(width = 0.75, padding = 0.1),
    angle = 90,
    size = 3,
    vjust = 0,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("Subtidal" = "grey60", "Intertidal" = "grey85")
  ) +
  labs(
    x = "Presence",
    y = "Relative probability of occurrence",
    fill = "Site type"
  ) +
  boxed_theme +
  coord_cartesian(ylim = c(y_min, y_max * 1.12))


kw_abundance <- kruskal.test(
  eelgrass_bccm_nospatial_pred ~ PC_ZM,
  data = abundance_df
)

kw_label <- case_when(
  kw_abundance$p.value < 0.001 ~ "***",
  kw_abundance$p.value < 0.01  ~ "**",
  kw_abundance$p.value < 0.05  ~ "*",
  TRUE                         ~ "ns"
)



abundance_counts <- abundance_df %>%
  group_by(PC_ZM, Tidal_zone) %>%
  summarise(n = n(), .groups = "drop")

# Build plot
p_abundance <- ggplot(
  abundance_df,
  aes(x = PC_ZM, y = eelgrass_bccm_spatial_pred, fill = Tidal_zone)
) +
  geom_boxplot(
    position = position_dodge2(width = 0.75), 
    outlier.alpha = 0.5
  ) +
  # Add n labels above each box
  geom_text(
    data = abundance_counts,
    aes(
      x = PC_ZM,
      y = y_max * 1.05,
      label = paste0("n=", n)
    ),
    position = position_dodge2(width = 0.75, padding = 0.1),
    angle = 90,
    size = 3,
    vjust = 0,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("Subtidal" = "grey60", "Intertidal" = "grey85")
  ) +
  labs(
    x = "Percent cover",
    y = "Relative probability of occurrence",
    fill = "Site type"
  ) +
  boxed_theme +
  coord_cartesian(ylim = c(y_min, y_max * 1.12))

p_abundance <- p_abundance +
  labs(y = NULL) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

p_abundance <- p_abundance +
  coord_cartesian(ylim = c(y_min, y_max * 1.08))

p_presence <- p_presence +
  labs(y = "Relative probability of occurrence") +
  scale_linetype_manual(
    name = "Threshold",   
    values = c(
      "Model threshold" = "dashed",
      "Field validated threshold" = "solid"
    )
  )

# Combine plots
zm_plot <- (p_presence + p_abundance) +
  plot_layout(
    widths = c(1, 2),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

zm_plot <- zm_plot +
  plot_annotation(
    title = "Eelgrass",
    theme = theme(
      plot.title = element_text(
        hjust = 0,
        face = "bold",
        size = 14,
        margin = margin(b = 5, l = 5)
      )
    )
  )

zm_plot




confusion_results_no_intertidal <- map(pred_cols, function(p) {
  
  thr <- roc_results %>%
    filter(model == p) %>%
    pull(tss_threshold)
  
  df_bin_nointertidal <- df_nointertidal %>%
    mutate(Pred_binary = if_else(.data[[p]] >= thr, 1, 0))
  
  confusionMatrix(
    factor(df_bin_nointertidal$Pred_binary),
    factor(df_bin_nointertidal[[obs_pa]]),
    positive = "1"
  )
})

names(confusion_results_no_intertidal) <- pred_cols
confusion_results_no_intertidal$eelgrass_bccm_spatial_pred
# this model is still sensitive (0.88), so catches most true presence. Specificity = 0.7593 so 76% ob absences
# Pos Pred Value (precision) = 0.3607 So when the model predicts eelgrass, it’s right only 36% of the time.not great
# Kappa = 0.397 this model still has moderate skill







