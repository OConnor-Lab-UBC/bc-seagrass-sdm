###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# rank all models to find best model
# create plot of all metrics in on 0-1
###############################################################################
# AUC, Tjur, Specificity, Sensitivity range 0-1 with higher values better
# Brier ranges 0-1 with low values better
# Logloss ranges 0- infinity with low values better and values >1 very poor or highly overconfident incorrect predictions (will cap 0-1)
# TSS ranges from -1 to 1 but values <0 are worse than random so can cap it at 0-1


library(dplyr)
library(fmsb)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(cowplot)
library(patchwork)
library(grid)
library(tidyr)
library(ggplot2)
library(tibble)
library(tidyverse)


# # ---------- Percentile rank functions ----------
# # Higher = better
pr_pos <- function(x) {
  (rank(x, ties.method = "average") - 1) / (sum(!is.na(x)) - 1)
}

pr_neg <- function(x) {
  (rank(-x, ties.method = "average") - 1) / (sum(!is.na(x)) - 1)
}

my_cols <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100)
breaks <- seq(0, 1, length.out = 101)

# ---- Colors ----
model_cols <- c(
  "GLM_bccm" = "#66C2A5",
  "GLMM_bccm" = "#FC8D62",
  "GBM_bccm" = "#8DA0CB",
  "XGBoost_bccm" = "#E78AC3",
  "GLM_nep" = "#A6D854",
  "GLMM_nep" = "#FFD92F",
  "GBM_nep" = "#E5C494",
  "XGBoost_nep" = "#B3B3B3"
)

model_order <- c(
  "GLM_bccm", "GLMM_bccm", "GBM_bccm", "XGBoost_bccm",
  "GLM_nep", "GLMM_nep", "GBM_nep", "XGBoost_nep"
)

###eelgrass
load("code/output_data/model_results/final_metrics_eelgrass.RData")
load("code/output_data/model_results/final_metrics_eelgrass_fieldnointertidal.RData")

rank_table_eelgrass <- final_metrics_eelgrass %>%
  select(-threshold_spatial, -threshold_field, -threshold_temporal, -threshold_independent) %>%  # Exclude all columns with 'threshold' in their names
  rename(Model = model, 
         AUC_Spatial = auc_spatial,
         Tjur_Spatial = tjur_spatial,
         Brier_Spatial = brier_spatial,
         Logloss_Spatial = logloss_spatial,
         Sensitivity_Spatial = sensitivity_spatial,
         Specificity_Spatial = specificity_spatial,
         TSS_Spatial = tss_spatial,
         AUC_Temporal = auc_temporal,
         Tjur_Temporal = tjur_temporal,
         Brier_Temporal = brier_temporal,
         Logloss_Temporal = logloss_temporal,
         Sensitivity_Temporal = sensitivity_temporal,
         Specificity_Temporal = specificity_temporal,
         TSS_Temporal = tss_temporal,
         AUC_Independent = auc_independent,
         Tjur_Independent = tjur_independent,
         Brier_Independent = brier_independent,
         Logloss_Independent = logloss_independent,
         Sensitivity_Independent = sensitivity_independent,
         Specificity_Independent = specificity_independent,
         TSS_Independent = tss_independent,
         AUC_Targeted = auc_field,
         Tjur_Targeted = tjur_field,
         Brier_Targeted = brier_field,
         Logloss_Targeted = logloss_field,
         Sensitivity_Targeted = sensitivity_field,
         Specificity_Targeted = specificity_field,
         TSS_Targeted = tss_field) 

#maybe better to have as percentile ranks to compare more broadly across whole grid
rank_table_pr_eelgrass <- rank_table_eelgrass %>%
  mutate(
    
    # Spatial
    Spatial_AUC = pr_pos(AUC_Spatial),
    Spatial_Tjur = pr_pos(Tjur_Spatial),
    Spatial_Brier = pr_neg(Brier_Spatial),
    Spatial_Logloss = pr_neg(Logloss_Spatial),
    Spatial_Sensitivity = pr_pos(Sensitivity_Spatial),
    Spatial_Specificity = pr_pos(Specificity_Spatial),
    Spatial_TSS = pr_pos(TSS_Spatial),
    
    # Temporal
    Temporal_AUC = pr_pos(AUC_Temporal),
    Temporal_Tjur = pr_pos(Tjur_Temporal),
    Temporal_Brier = pr_neg(Brier_Temporal),
    Temporal_Logloss = pr_neg(Logloss_Temporal),
    Temporal_Sensitivity = pr_pos(Sensitivity_Temporal),
    Temporal_Specificity = pr_pos(Specificity_Temporal),
    Temporal_TSS = pr_pos(TSS_Temporal),
    
    # Independent
    Independent_AUC = pr_pos(AUC_Independent),
    Independent_Tjur = pr_pos(Tjur_Independent),
    Independent_Brier = pr_neg(Brier_Independent),
    Independent_Logloss = pr_neg(Logloss_Independent),
    Independent_Sensitivity = pr_pos(Sensitivity_Independent),
    Independent_Specificity = pr_pos(Specificity_Independent),
    Independent_TSS = pr_pos(TSS_Independent),
    
    # Targeted
    Targeted_AUC = pr_pos(AUC_Targeted),
    Targeted_Tjur = pr_pos(Tjur_Targeted),
    Targeted_Brier = pr_neg(Brier_Targeted),
    Targeted_Logloss = pr_neg(Logloss_Targeted),
    Targeted_Sensitivity = pr_pos(Sensitivity_Targeted),
    Targeted_Specificity = pr_pos(Specificity_Targeted),
    Targeted_TSS = pr_pos(TSS_Targeted)
  )

rank_table_pr_eelgrass <- rank_table_pr_eelgrass %>%
  mutate(
    Spatial = rowMeans(
      select(., starts_with("Spatial_")),
      na.rm = TRUE
    ),
    
    Temporal = rowMeans(
      select(., starts_with("Temporal_")),
      na.rm = TRUE
    ),
    
    Independent = rowMeans(
      select(., starts_with("Independent_")),
      na.rm = TRUE
    ),
    
    Targeted = rowMeans(
      select(., starts_with("Targeted_")),
      na.rm = TRUE
    ),
    
    FinalScore = rowMeans(
      select(., 
             starts_with("Spatial_"),
             starts_with("Temporal_"),
             starts_with("Independent_"),
             starts_with("Targeted_")),
      na.rm = TRUE
    )
  )

# ---- Build grouped dataset with spacers ----
# rank_table_pr <- rank_table_pr %>%
#   mutate(Model = factor(Model, levels = model_order)) %>%
#   arrange(Model)




plot_df_eelgrass <- rank_table_eelgrass %>%
  select(-matches("_sd_")) %>%
  pivot_longer(
    cols = -Model,
    names_to = c("Metric", "Validation"),
    names_pattern = "(.+)_(Spatial|Temporal|Independent|Targeted)",
    values_to = "Value"
  )


plot_sd_eelgrass <- rank_table_eelgrass %>%
  select(Model, matches("_sd_")) %>%
  pivot_longer(
    cols = -Model,
    names_to = c("Metric", "Validation"),
    names_pattern = "^(.*)_sd_(.*)$",
    values_to = "SD"
  ) %>%
  mutate(
    Metric = str_to_title(Metric),
    Metric = recode(
      Metric,
      Auc = "AUC",
      Tss = "TSS",
      Tjur = "Tjur",
      Logloss = "Logloss"
    ),
    Validation = str_to_title(Validation)
  )


# combine means + SD
plot_df_eelgrass <- plot_df_eelgrass %>%
  left_join(
    plot_sd_eelgrass,
    by = c("Model", "Metric", "Validation")
  ) %>%
  mutate(
    # transform metrics so higher = better
    PlotValue = case_when(
      Metric == "Brier"   ~ pmax(1 - Value, 0),
      Metric == "Logloss" ~ pmax(1 - Value, 0),
      Metric == "TSS"     ~ pmax(Value, 0),
      TRUE                ~ Value
    ),
    
    # SD does not change for 1-x transformations
    SD_plot = SD,
    
    xmin = pmax(PlotValue - SD_plot, 0),
    xmax = pmin(PlotValue + SD_plot, 1),
    
    Metric = factor(
      Metric,
      levels = c(
        "Brier",
        "Logloss",
        "TSS",
        "Specificity",
        "Sensitivity",
        "Tjur",
        "AUC"
      ),
      labels = c(
        "1 - Brier",
        "1 - Log Loss",
        "TSS",
        "Specificity",
        "Sensitivity",
        "Tjur R²",
        "AUC"
      )
    ),
    
    Validation = factor(
      Validation,
      levels = c("Spatial", "Temporal", "Independent", "Targeted")
    ),
    
    Model = factor(Model, levels = model_order)
  )

metric_plot_eelgrass <- ggplot(
  plot_df_eelgrass,
  aes(x = PlotValue, y = Metric, colour = Model)
) +
  geom_errorbar(
    aes(xmin = xmin, xmax = xmax),
    linewidth = 0.5,
    alpha = 0.7,
    na.rm = TRUE,
    position = position_dodge(width = 0.7)
  ) +
  geom_point(
    shape = 16,      # solid circle
    size = 1.5,
    alpha = 0.85,
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~Validation, nrow = 2) +
  scale_colour_manual(
    values = model_cols,
    breaks = model_order
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.1f", x)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Metric value",
    y = NULL,
    colour = "Model"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.5, "cm")
  ) +
  guides(
    colour = guide_legend(
      nrow = 4,
      byrow = TRUE
    )
  )

# Extract model legend BEFORE removing it
model_leg_eelgrass <- cowplot::get_legend(
  metric_plot_eelgrass +
    theme(legend.position = "bottom")
)

# Remove legend from the plot used in the figure
metric_plot_eelgrass_noleg <- metric_plot_eelgrass +
  theme(legend.position = "none")

eelgrass_tags <- data.frame(
  Validation = factor(
    c("Spatial", "Temporal", "Independent", "Targeted"),
    levels = c("Spatial", "Temporal", "Independent", "Targeted")
  ),
  Metric = factor("AUC", levels = levels(plot_df_eelgrass$Metric)),
  x = 0.02,
  lab = c("A", "B", "C", "D")
)

metric_plot_eelgrass_noleg <- metric_plot_eelgrass_noleg +
  geom_text(
    data = eelgrass_tags,
    aes(x = x, y = Metric, label = lab),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 5
  )

# Heat map
  heat_df_cat_eelgrass <- rank_table_pr_eelgrass %>%
    select(Model, Model, Spatial, Temporal, Independent, Targeted, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  
  # Rebuild heatmap legend (continuous percentile scale)
  heat_leg_plot_eelgrass <- ggplot(
    data.frame(
      x = seq(0, 1, length.out = 100),
      y = 1,
      val = seq(0, 1, length.out = 100)
    ),
    aes(x = x, y = y, fill = val)
  ) +
    geom_tile() +
    scale_fill_gradientn(
      colours = my_cols,
      limits = c(0, 1),
      name = "Relative performance percentile",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(5.5, "cm"),
        barheight = unit(0.35, "cm"),
        ticks = TRUE
      )
    ) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1")
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  heat_leg_eelgrass <- cowplot::get_legend(heat_leg_plot_eelgrass)
  
  
  
  # Create heatmaps
  heat_plot_eelgrass <- as.ggplot(function() {
    pheatmap(
      heat_df_cat_eelgrass,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      treeheight_row = 0,
      treeheight_col = 0,
      legend = FALSE,
      color = my_cols,
      breaks = breaks,
      cellwidth = 20,
      cellheight = 20,
      fontsize = 10,
      fontsize_row = 9,
      fontsize_col = 10
    )
  })

  heat_plot_eelgrass <- cowplot::ggdraw(heat_plot_eelgrass) +
    cowplot::draw_label(
      "E",
      x = 0.75,
      y = 0.17,
      hjust = 1,
      vjust = 0,
      fontface = "bold",
      size = 14
    )

  # Combine plot + heatmap
  combined_plot_eelgrass <- metric_plot_eelgrass_noleg + heat_plot_eelgrass +
    plot_layout(
      ncol = 2,
      widths = c(2.3, 1.2)
    )

  # For eelgrass  pick XGBOOST NEP

  
  
  
  ###surfgrass
  load("code/output_data/model_results/final_metrics_surfgrass.RData")
  load("code/output_data/model_results/final_metrics_surfgrass_fieldnointertidal.RData")
  
  rank_table_surfgrass <- final_metrics_surfgrass %>%
    select(-threshold_spatial, -threshold_field, -threshold_temporal, -threshold_independent) %>%  # Exclude all columns with 'threshold' in their names
    rename(Model = model, 
           AUC_Spatial = auc_spatial,
           Tjur_Spatial = tjur_spatial,
           Brier_Spatial = brier_spatial,
           Logloss_Spatial = logloss_spatial,
           Sensitivity_Spatial = sensitivity_spatial,
           Specificity_Spatial = specificity_spatial,
           TSS_Spatial = tss_spatial,
           AUC_Temporal = auc_temporal,
           Tjur_Temporal = tjur_temporal,
           Brier_Temporal = brier_temporal,
           Logloss_Temporal = logloss_temporal,
           Sensitivity_Temporal = sensitivity_temporal,
           Specificity_Temporal = specificity_temporal,
           TSS_Temporal = tss_temporal,
           AUC_Independent = auc_independent,
           Tjur_Independent = tjur_independent,
           Brier_Independent = brier_independent,
           Logloss_Independent = logloss_independent,
           Sensitivity_Independent = sensitivity_independent,
           Specificity_Independent = specificity_independent,
           TSS_Independent = tss_independent,
           AUC_Targeted = auc_field,
           Tjur_Targeted = tjur_field,
           Brier_Targeted = brier_field,
           Logloss_Targeted = logloss_field,
           Sensitivity_Targeted = sensitivity_field,
           Specificity_Targeted = specificity_field,
           TSS_Targeted = tss_field) 
  
  #maybe better to have as percentile ranks to compare more broadly across whole grid
  rank_table_pr_surfgrass <- rank_table_surfgrass %>%
    mutate(
      
      # Spatial
      Spatial_AUC = pr_pos(AUC_Spatial),
      Spatial_Tjur = pr_pos(Tjur_Spatial),
      Spatial_Brier = pr_neg(Brier_Spatial),
      Spatial_Logloss = pr_neg(Logloss_Spatial),
      Spatial_Sensitivity = pr_pos(Sensitivity_Spatial),
      Spatial_Specificity = pr_pos(Specificity_Spatial),
      Spatial_TSS = pr_pos(TSS_Spatial),
      
      # Temporal
      Temporal_AUC = pr_pos(AUC_Temporal),
      Temporal_Tjur = pr_pos(Tjur_Temporal),
      Temporal_Brier = pr_neg(Brier_Temporal),
      Temporal_Logloss = pr_neg(Logloss_Temporal),
      Temporal_Sensitivity = pr_pos(Sensitivity_Temporal),
      Temporal_Specificity = pr_pos(Specificity_Temporal),
      Temporal_TSS = pr_pos(TSS_Temporal),
      
      # Independent
      Independent_AUC = pr_pos(AUC_Independent),
      Independent_Tjur = pr_pos(Tjur_Independent),
      Independent_Brier = pr_neg(Brier_Independent),
      Independent_Logloss = pr_neg(Logloss_Independent),
      Independent_Sensitivity = pr_pos(Sensitivity_Independent),
      Independent_Specificity = pr_pos(Specificity_Independent),
      Independent_TSS = pr_pos(TSS_Independent),
      
      # Targeted
      Targeted_AUC = pr_pos(AUC_Targeted),
      Targeted_Tjur = pr_pos(Tjur_Targeted),
      Targeted_Brier = pr_neg(Brier_Targeted),
      Targeted_Logloss = pr_neg(Logloss_Targeted),
      Targeted_Sensitivity = pr_pos(Sensitivity_Targeted),
      Targeted_Specificity = pr_pos(Specificity_Targeted),
      Targeted_TSS = pr_pos(TSS_Targeted)
    )
  
  rank_table_pr_surfgrass <- rank_table_pr_surfgrass %>%
    mutate(
      Spatial = rowMeans(
        select(., starts_with("Spatial_")),
        na.rm = TRUE
      ),
      
      Temporal = rowMeans(
        select(., starts_with("Temporal_")),
        na.rm = TRUE
      ),
      
      Independent = rowMeans(
        select(., starts_with("Independent_")),
        na.rm = TRUE
      ),
      
      Targeted = rowMeans(
        select(., starts_with("Targeted_")),
        na.rm = TRUE
      ),
      
      FinalScore = rowMeans(
        select(., 
               starts_with("Spatial_"),
               starts_with("Temporal_"),
               starts_with("Independent_"),
               starts_with("Targeted_")),
        na.rm = TRUE
      )
    )
  
  # ---- Build grouped dataset with spacers ----
  # rank_table_pr <- rank_table_pr %>%
  #   mutate(Model = factor(Model, levels = model_order)) %>%
  #   arrange(Model)
  
  
  
  
  plot_df_surfgrass <- rank_table_surfgrass %>%
    select(-matches("_sd_")) %>%
    pivot_longer(
      cols = -Model,
      names_to = c("Metric", "Validation"),
      names_pattern = "(.+)_(Spatial|Temporal|Independent|Targeted)",
      values_to = "Value"
    )
  
  
  plot_sd_surfgrass <- rank_table_surfgrass %>%
    select(Model, matches("_sd_")) %>%
    pivot_longer(
      cols = -Model,
      names_to = c("Metric", "Validation"),
      names_pattern = "^(.*)_sd_(.*)$",
      values_to = "SD"
    ) %>%
    mutate(
      Metric = str_to_title(Metric),
      Metric = recode(
        Metric,
        Auc = "AUC",
        Tss = "TSS",
        Tjur = "Tjur",
        Logloss = "Logloss"
      ),
      Validation = str_to_title(Validation)
    )
  
  
  # combine means + SD
  plot_df_surfgrass <- plot_df_surfgrass %>%
    left_join(
      plot_sd_surfgrass,
      by = c("Model", "Metric", "Validation")
    ) %>%
    mutate(
      # transform metrics so higher = better
      PlotValue = case_when(
        Metric == "Brier"   ~ pmax(1 - Value, 0),
        Metric == "Logloss" ~ pmax(1 - Value, 0),
        Metric == "TSS"     ~ pmax(Value, 0),
        TRUE                ~ Value
      ),
      
      # SD does not change for 1-x transformations
      SD_plot = SD,
      
      xmin = pmax(PlotValue - SD_plot, 0),
      xmax = pmin(PlotValue + SD_plot, 1),
      
      Metric = factor(
        Metric,
        levels = c(
          "Brier",
          "Logloss",
          "TSS",
          "Specificity",
          "Sensitivity",
          "Tjur",
          "AUC"
        ),
        labels = c(
          "1 - Brier",
          "1 - Log Loss",
          "TSS",
          "Specificity",
          "Sensitivity",
          "Tjur R²",
          "AUC"
        )
      ),
      
      Validation = factor(
        Validation,
        levels = c("Spatial", "Temporal", "Independent", "Targeted")
      ),
      
      Model = factor(Model, levels = model_order)
    )
  
  metric_plot_surfgrass <- ggplot(
    plot_df_surfgrass,
    aes(x = PlotValue, y = Metric, colour = Model)
  ) +
    geom_errorbar(
      aes(xmin = xmin, xmax = xmax),
      linewidth = 0.5,
      alpha = 0.7,
      na.rm = TRUE,
      position = position_dodge(width = 0.7)
    ) +
    geom_point(
      shape = 16,      # solid circle
      size = 1.5,
      alpha = 0.85,
      position = position_dodge(width = 0.7)
    ) +
    facet_wrap(~Validation, nrow = 2) +
    scale_colour_manual(
      values = model_cols,
      breaks = model_order
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.1f", x)),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "Metric value",
      y = NULL,
      colour = "Model"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.4),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      legend.spacing.x = unit(0.5, "cm")
    ) +
    guides(
      colour = guide_legend(
        nrow = 4,
        byrow = TRUE
      )
    )
  
  # Extract model legend BEFORE removing it
  model_leg_surfgrass <- cowplot::get_legend(
    metric_plot_surfgrass +
      theme(legend.position = "bottom")
  )
  
  # Remove legend from the plot used in the figure
  metric_plot_surfgrass_noleg <- metric_plot_surfgrass +
    theme(legend.position = "none")
  
  surfgrass_tags <- data.frame(
    Validation = factor(
      c("Spatial", "Temporal", "Independent", "Targeted"),
      levels = c("Spatial", "Temporal", "Independent", "Targeted")
    ),
    Metric = factor("AUC", levels = levels(plot_df_surfgrass$Metric)),
    x = 0.02,
    lab = c("F", "G", "H", "I")
  )
  
  metric_plot_surfgrass_noleg <- metric_plot_surfgrass_noleg +
    geom_text(
      data = surfgrass_tags,
      aes(x = x, y = Metric, label = lab),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 5
    )
  
  # Also test if we make a score based on each of the 4 validaitons. 
  # Heat map
  heat_df_cat_surfgrass <- rank_table_pr_surfgrass %>%
    select(Model, Model, Spatial, Temporal, Independent, Targeted, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  
  # Rebuild heatmap legend (continuous percentile scale)
  heat_leg_plot_surfgrass <- ggplot(
    data.frame(
      x = seq(0, 1, length.out = 100),
      y = 1,
      val = seq(0, 1, length.out = 100)
    ),
    aes(x = x, y = y, fill = val)
  ) +
    geom_tile() +
    scale_fill_gradientn(
      colours = my_cols,
      limits = c(0, 1),
      name = "Relative performance percentile",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(5.5, "cm"),
        barheight = unit(0.35, "cm"),
        ticks = TRUE
      )
    ) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1")
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  heat_leg_surfgrass <- cowplot::get_legend(heat_leg_plot_surfgrass)
  
  # Create heatmaps
  heat_plot_surfgrass <- as.ggplot(function() {
    pheatmap(
      heat_df_cat_surfgrass,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      treeheight_row = 0,
      treeheight_col = 0,
      legend = FALSE,
      color = my_cols,
      breaks = breaks,
      cellwidth = 20,
      cellheight = 20,
      fontsize = 10,
      fontsize_row = 9,
      fontsize_col = 10
    )
  })
  
  
  
  heat_plot_surfgrass <- cowplot::ggdraw(heat_plot_surfgrass) +
    cowplot::draw_label(
      "J",
      x = 0.75,
      y = 0.17,
      hjust = 1,
      vjust = 0,
      fontface = "bold",
      size = 14
    )
  
  # Combine plot + heatmap
  combined_plot_surfgrass <- metric_plot_surfgrass_noleg + heat_plot_surfgrass +
    plot_layout(
      ncol = 2,
      widths = c(2.3, 1.2)
    )
  
  # ---- Combine legends using patchwork/cowplot wrapping ----
  legend_row <- wrap_elements(model_leg_surfgrass) + wrap_elements(heat_leg_surfgrass) +
    plot_layout(
      ncol = 2,
      widths = c(2, 2)
    )
  
  eelgrass_title <- ggplot() +
    annotate(
      "text",
      x = -0.15, y = 0.5,
      label = "Eelgrass",
      hjust = 0,
      size = 4,        # reduce from 6
      fontface = "bold"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
  
  surfgrass_title <- ggplot() +
    annotate(
      "text",
      x = -0.15, y = 0.5,
      label = "Surfgrass",
      hjust = 0,
      size = 4,
      fontface = "bold"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
  
  
  # ---- Final figure ----
  final_plot_seagrass <-
    eelgrass_title /
    combined_plot_eelgrass /
    surfgrass_title /
    combined_plot_surfgrass /
    legend_row +
    plot_layout(heights = c(0.06, 1, 0.06, 1, 0.4))
  final_plot_seagrass
  
  
  ggsave(
    filename = "figures/model_performance_seagrass.png",
    plot = final_plot_seagrass,
    width = 190,
    height = 250,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
