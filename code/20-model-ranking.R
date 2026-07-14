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
load("code/output_data/model_results/final_metrics_eelgrass.RData")
load("code/output_data/model_results/final_metrics_eelgrass_fieldnointertidal.RData")
#load("./code/output_data/model_results/combined_metrics_surfgrass_4_validations.RData")

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
rank_table <- final_metrics_eelgrass %>%
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
rank_table_pr <- rank_table %>%
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

rank_table_pr <- rank_table_pr %>%
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




plot_df <- rank_table %>%
  select(-matches("_sd_")) %>%
  pivot_longer(
    cols = -Model,
    names_to = c("Metric", "Validation"),
    names_pattern = "(.+)_(Spatial|Temporal|Independent|Targeted)",
    values_to = "Value"
  )


plot_sd <- rank_table %>%
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
plot_df <- plot_df %>%
  left_join(
    plot_sd,
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
  plot_df,
  aes(x = PlotValue, y = Metric, colour = Model)
) +
  geom_errorbar(
    aes(xmin = xmin, xmax = xmax),
    linewidth = 1,
    alpha = 0.7,
    na.rm = TRUE,
    position = position_dodge(width = 0.65)
  ) +
  geom_point(
    shape = 16,      # solid circle
    size = 2.5,
    alpha = 0.85,
    position = position_dodge(width = 0.65)
  ) +
  facet_wrap(~Validation, nrow = 2) +
  scale_colour_manual(
    values = model_cols,
    breaks = model_order
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
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
      nrow = 2,
      byrow = TRUE
    )
  )
metric_plot_eelgrass
# Extract model legend BEFORE removing it
model_leg <- cowplot::get_legend(
  metric_plot_eelgrass +
    theme(legend.position = "bottom")
)

# Remove legend from the plot used in the figure
metric_plot_eelgrass_noleg <- metric_plot_eelgrass +
  theme(legend.position = "none")


# Also test if we make a score based on each of the 4 validaitons. 
  # Heat map
  heat_df_cat <- rank_table_pr %>%
    select(Model, Model, Spatial, Temporal, Independent, Targeted, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  
  # Rebuild heatmap legend (continuous percentile scale)
  heat_leg_plot <- ggplot(
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
  
  heat_leg <- cowplot::get_legend(heat_leg_plot)
  
  # Create heatmaps
  heat_plot_eelgrass <- as.ggplot(function() {
    pheatmap(
      heat_df_cat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      treeheight_row = 0,
      treeheight_col = 0,
      legend = FALSE,
      color = my_cols,
      breaks = breaks,
      cellwidth = 70,
      cellheight = 35,
      fontsize = 12,
      fontsize_row = 11,
      fontsize_col = 12
    )
  })

 
  # ---- Combine the two legends ----
  legend_row <- plot_grid(
    model_leg,
    heat_leg,
    ncol = 2,
    rel_widths = c(1.4, 0.8)
  )
  
  
  # Combine plot + heatmap
  combined_plot_eelgrass <- metric_plot_eelgrass_noleg + heat_plot_eelgrass +
    plot_layout(
      widths = c(3, 1.5)
    )
  
  
  # ---- Final figure with legends underneath ----
  final_plot_eelgrass <- plot_grid(
    combined_plot_eelgrass,
    legend_row,
    ncol = 1,
    rel_heights = c(1, 0.15)
  )
  
  final_plot_eelgrass
  
  # For eelgrass  pick XGBOOST NEP
 