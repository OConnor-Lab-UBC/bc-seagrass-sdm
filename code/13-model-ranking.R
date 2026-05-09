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
###############################################################################

load("./code/output_data/model_results/combined_metrics_eelgrass_4_validations.RData")
load("./code/output_data/model_results/combined_metrics_surfgrass_4_validations.RData")

library(dplyr)
library(fmsb)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(cowplot)
library(patchwork)
library(grid)

my_cols <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100)
breaks <- seq(0, 1, length.out = 101)

model_order <- c(
  "GLMM_bccm", "GLMM_spatial_bccm", "GBM_bccm", "XGB_bccm",
  "GLMM_nep", "GLMM_spatial_nep", "GBM_nep", "XGB_nep"
)

###eelgrass
rank_table <- combined_metrics_all_eelgrass %>%
  select(-threshold_spatial, -threshold_field, -sensitivity_field, - specificity_field) %>%  # Exclude all columns with 'threshold' in their names
  rename(Model = model, 
         AUC_Spatial = auc_spatial,
         Tjur_Spatial = tjur_spatial,
         RMSE_Spatial = rmse_spatial,
         TSS_Spatial = tss_spatial,
         AUC_Temporal = auc_temporal,
         Tjur_Temporal = tjur_temporal,
         RMSE_Temporal = rmse_temporal,
         TSS_Temporal = tss_temporal,
         MPS = mps_independent,
         FPPS = fpps_independent,
         FNR = fnr_independent,
         CBI = cbi_independent,
         AUC_Field = auc_field,
         Tjur_Field = tjur_r2_field,
         RMSE_Field = rmse_field,
         TSS_Field = TSS_field,
         Sensitivity_FieldCV = sensitivity_field_cv,
         Specificity_FieldCV = specificity_field_cv,
         TSS_FieldCV = TSS_field_cv,
         Cliffs_delta = cliffs_delta_field) 
# Normalize metrics
norm_pos <- function(x) (x - min(x)) / (max(x) - min(x))
norm_neg <- function(x) (max(x) - x) / (max(x) - min(x))
# Scale metrics and create a final composite score
rank_table_scaled <- rank_table %>%
  mutate(
    CV_AUC = norm_pos(AUC_Spatial),
    CV_Tjur = norm_pos(Tjur_Spatial),
    CV_RMSE = norm_neg(RMSE_Spatial),
    CV_TSS = norm_pos(TSS_Spatial),
    
    TEMP_AUC = norm_pos(AUC_Temporal),
    TEMP_Tjur = norm_pos(Tjur_Temporal),
    TEMP_RMSE = norm_neg(RMSE_Temporal),
    TEMP_TSS = norm_pos(TSS_Temporal),
    
    INDEP_MPS = norm_pos(MPS),
    INDEP_FPPS = norm_pos(FPPS),
    INDEP_FNR = norm_neg(FNR),
    INDEP_CBI = norm_pos(CBI),
    
    FIELD_AUC = norm_pos(AUC_Field),
    FIELD_Tjur = norm_pos(Tjur_Field),
    FIELD_RMSE = norm_neg(RMSE_Field),
    FIELD_TSS = norm_pos(TSS_Field),
    #FIELD_SenCV = norm_pos(Sensitivity_FieldCV),
    #FIELD_SpeCV = norm_pos(Specificity_FieldCV),
    FIELD_TSSCV = norm_pos(TSS_FieldCV)
    #FIELD_CliffsD = norm_pos(Cliffs_delta)
  )
# Compute final composite score
rank_table_scaled$FinalScore <- rowMeans(rank_table_scaled %>%
                                           select(CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
                                                  TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
                                                  INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
                                                  FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV),
                                         na.rm = TRUE)
# Rank models based on the final score
final_rank <- rank_table_scaled %>%
  arrange(desc(FinalScore)) %>%
  select(Model, FinalScore)


# ---- Build grouped dataset with spacers ----
rank_table_scaled <- rank_table_scaled %>%
  mutate(Model = recode(Model,
                        "bccm_nospatial" = "GLMM_bccm",
                        "bccm_spatial" = "GLMM_spatial_bccm",
                        "GBM_bccm" = "GBM_bccm",
                        "XGBOOST_bccm" = "XGB_bccm",
                        "nep_nospatial" = "GLMM_nep",
                        "nep_spatial" = "GLMM_spatial_nep",
                        "GBM_nep" = "GBM_nep",
                        "XGBOOST_nep" = "XGB_nep" ))

rank_table_scaled <- rank_table_scaled %>%
  mutate(Model = factor(Model, levels = model_order)) %>%
  arrange(Model)

radar_df <- rank_table_scaled %>%
  select(Model,
         CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
         TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
         INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI,
         FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV) %>%
  mutate(
    spacer1 = NA,
    spacer2 = NA,
    spacer3 = NA,
    spacer4 = NA
  ) %>%
  select(Model,
         CV_AUC, CV_Tjur, CV_RMSE, CV_TSS, spacer1,
         TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS, spacer2,
         INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, spacer3,
         FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV, spacer4)

# Replace spacers with 0
radar_df[is.na(radar_df)] <- 0

# ---- Labels ----
colnames(radar_df) <- c(
  "Model",
  "AUC", "Tjur", "RMSE", "TSS", " ",
  "AUC ", "Tjur ", "RMSE ", "TSS ", "  ",
  "MPS", "FPPS", "FNR", "CBI", "   ",
  "AUC  ", "Tjur  ", "RMSE  ",
  "TSS  ", "    "
)

# ---- Radar format ----
radar_ready <- rbind(
  rep(1, ncol(radar_df) - 1),
  rep(0, ncol(radar_df) - 1),
  radar_df[, -1]
)

rownames(radar_ready) <- c("max", "min", radar_df$Model)
plot_cols <- unname(model_cols[as.character(radar_df$Model)])
# ---- Colors ----
model_cols <- c(
  "GLMM_bccm" = "#66C2A5",
  "GLMM_spatial_bccm" = "#FC8D62",
  "GBM_bccm" = "#8DA0CB",
  "XGB_bccm" = "#E78AC3",
  "GLMM_nep" = "#A6D854",
  "GLMM_spatial_nep" = "#FFD92F",
  "GBM_nep" = "#E5C494",
  "XGB_nep" = "#B3B3B3"
)

# ---- Plot ----
par(mar = c(2, 2, 3, 2))


# Also test if we make a score based on each of the 4 validaitons. 
  
  # Compute final composite scores for each category
  rank_table_scaled_cat <- rank_table_scaled %>%
    mutate(
      Spatial = rowMeans(select(., CV_AUC, CV_Tjur, CV_RMSE, CV_TSS), na.rm = TRUE),
      Temporal = rowMeans(select(., TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS), na.rm = TRUE),
      Independent = rowMeans(select(., INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI), na.rm = TRUE),
      Field = rowMeans(select(., FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV), na.rm = TRUE)
    )
  # Rank models based on the final score for each category
  final_rank_cat <- rank_table_scaled_cat %>%
    arrange(desc(Spatial)) %>%
    select(Model, Spatial, Temporal, Independent, Field)
  
  # Heat map
  heat_df_cat <- rank_table_scaled_cat %>%
    select(Model, Model, Spatial, Temporal, Independent, Field, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  
  ###surfgrass
  rank_table_surf <- combined_metrics_all_surfgrass %>%
    select(-threshold_spatial, -threshold_field, -sensitivity_field, - specificity_field) %>%  # Exclude all columns with 'threshold' in their names
    rename(Model = model, 
           AUC_Spatial = auc_spatial,
           Tjur_Spatial = tjur_spatial,
           RMSE_Spatial = rmse_spatial,
           TSS_Spatial = tss_spatial,
           AUC_Temporal = auc_temporal,
           Tjur_Temporal = tjur_temporal,
           RMSE_Temporal = rmse_temporal,
           TSS_Temporal = tss_temporal,
           MPS = mps_independent,
           FPPS = fpps_independent,
           FNR = fnr_independent,
           CBI = cbi_independent,
           AUC_Field = auc_field,
           Tjur_Field = tjur_r2_field,
           RMSE_Field = rmse_field,
           TSS_Field = TSS_field,
           Sensitivity_FieldCV = sensitivity_field_cv,
           Specificity_FieldCV = specificity_field_cv,
           TSS_FieldCV = TSS_field_cv,
           Cliffs_delta = cliffs_delta_field) 
  # Normalize metrics
  norm_pos <- function(x) (x - min(x)) / (max(x) - min(x))
  norm_neg <- function(x) (max(x) - x) / (max(x) - min(x))
 
   # Scale metrics and create a final composite score
  rank_table_scaled_surf <- rank_table_surf %>%
    mutate(
      CV_AUC = norm_pos(AUC_Spatial),
      CV_Tjur = norm_pos(Tjur_Spatial),
      CV_RMSE = norm_neg(RMSE_Spatial),
      CV_TSS = norm_pos(TSS_Spatial),
      
      TEMP_AUC = norm_pos(AUC_Temporal),
      TEMP_Tjur = norm_pos(Tjur_Temporal),
      TEMP_RMSE = norm_neg(RMSE_Temporal),
      TEMP_TSS = norm_pos(TSS_Temporal),
      
      INDEP_MPS = norm_pos(MPS),
      INDEP_FPPS = norm_pos(FPPS),
      INDEP_FNR = norm_neg(FNR),
      INDEP_CBI = norm_pos(CBI),
      
      FIELD_AUC = norm_pos(AUC_Field),
      FIELD_Tjur = norm_pos(Tjur_Field),
      FIELD_RMSE = norm_neg(RMSE_Field),
      FIELD_TSS = norm_pos(TSS_Field),
      #FIELD_SenCV = norm_pos(Sensitivity_FieldCV),
      #FIELD_SpeCV = norm_pos(Specificity_FieldCV),
      FIELD_TSSCV = norm_pos(TSS_FieldCV)
      #FIELD_CliffsD = norm_pos(Cliffs_delta)
    )
  # Compute final composite score
  rank_table_scaled_surf$FinalScore <- rowMeans(rank_table_scaled_surf %>%
                                             select(CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
                                                    TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
                                                    INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
                                                    FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV),
                                           na.rm = TRUE)
  # Rank models based on the final score
  final_rank <- rank_table_scaled_surf %>%
    arrange(desc(FinalScore)) %>%
    select(Model, FinalScore)
  
  # ---- Build grouped dataset with spacers ----
  rank_table_scaled_surf <- rank_table_scaled_surf %>%
    mutate(Model = recode(Model,
                          "bccm_nospatial" = "GLMM_bccm",
                          "bccm_spatial" = "GLMM_spatial_bccm",
                          "GBM_bccm" = "GBM_bccm",
                          "XGBOOST_bccm" = "XGB_bccm",
                          "nep_nospatial" = "GLMM_nep",
                          "nep_spatial" = "GLMM_spatial_nep",
                          "GBM_nep" = "GBM_nep",
                          "XGBOOST_nep" = "XGB_nep" ))
  
  rank_table_scaled_surf <- rank_table_scaled_surf %>%
    mutate(Model = factor(Model, levels = model_order)) %>%
    arrange(Model)
  
  
  radar_df_surf <- rank_table_scaled_surf %>%
    select(Model,
           CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
           TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
           INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI,
           FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV) %>%
    mutate(
      spacer1 = NA,
      spacer2 = NA,
      spacer3 = NA,
      spacer4 = NA
    ) %>%
    select(Model,
           CV_AUC, CV_Tjur, CV_RMSE, CV_TSS, spacer1,
           TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS, spacer2,
           INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, spacer3,
           FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV, spacer4)
  
  # Replace spacers with 0
  radar_df_surf[is.na(radar_df_surf)] <- 0
  
  # ---- Labels ----
  colnames(radar_df_surf) <- c(
    "Model",
    "AUC", "Tjur", "RMSE", "TSS", " ",
    "AUC ", "Tjur ", "RMSE ", "TSS ", "  ",
    "MPS", "FPPS", "FNR", "CBI", "   ",
    "AUC  ", "Tjur  ", "RMSE  ",
    "TSS  ", "    "
  )
  
  # ---- Radar format ----
  radar_ready_surf <- rbind(
    rep(1, ncol(radar_df_surf) - 1),
    rep(0, ncol(radar_df_surf) - 1),
    radar_df_surf[, -1]
  )
  
  rownames(radar_ready_surf) <- c("max", "min", radar_df_surf$Model)
  plot_cols_surf <- unname(model_cols[as.character(radar_df_surf$Model)])
  # ---- Colors ----
  #cols <- brewer.pal(n = nrow(radar_df_surf), name = "Set2")
  
  # ---- Plot ----
  par(mar = c(2, 2, 3, 2))
  
  # Also test if we make a score based on each of the 4 validaitons. 
  
  # Compute final composite scores for each category
  rank_table_scaled_cat_surf <- rank_table_scaled_surf %>%
    mutate(
      Spatial = rowMeans(select(., CV_AUC, CV_Tjur, CV_RMSE, CV_TSS), na.rm = TRUE),
      Temporal = rowMeans(select(., TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS), na.rm = TRUE),
      Independent = rowMeans(select(., INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI), na.rm = TRUE),
      Field = rowMeans(select(., FIELD_AUC, FIELD_Tjur, FIELD_RMSE, FIELD_TSSCV), na.rm = TRUE)
    )
  # Rank models based on the final score for each category
  final_rank_cat_surf <- rank_table_scaled_cat_surf %>%
    arrange(desc(Spatial)) %>%
    select(Model, Spatial, Temporal, Independent, Field)
  
  # Heat map
  heat_df_cat_surf <- rank_table_scaled_cat_surf %>%
    select(Model, Model, Spatial, Temporal, Independent, Field, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  
  
  
  #updated radar maps
  radar_plot_eelgrass <- as.ggplot(function() {
    par(
      mar = c(2.5, 2.5, 2.5, 2.5),  # much tighter
      oma = c(2, 2, 2, 2),
      xpd = NA
    )
    # par(
    #   mar = c(6, 6, 6, 6),
    #   oma = c(2, 2, 2, 2),
    #   xpd = NA
    # )
    radarchart(
      radar_ready,
      axistype = 0,
      caxislabels = NULL,
      pcol = plot_cols,
      plwd = 2,
      plty = 1,
      cglcol = "grey85",
      cglty = 1,
      cglwd = 0.8,
      vlcex = 0.6
    )

    n_axes <- ncol(radar_ready) - 1
    angles <- seq(pi/2, 2*pi + pi/2, length.out = n_axes + 1)[-(n_axes + 1)]
    
    group_pos <- list(
      Spatial = mean(angles[1:4]),
      Temporal = mean(angles[6:9]),
      Independent = mean(angles[11:14]),
      Field = mean(angles[16:19])
    )
    
    r <- 1.7
    
    text(r * cos(group_pos$Spatial), r * sin(group_pos$Spatial), "Spatial", font = 2, cex = 0.9)
    text(r * cos(group_pos$Temporal), r * sin(group_pos$Temporal), "Temporal", font = 2, cex = 0.9)
    text(r * cos(group_pos$Independent), r * sin(group_pos$Independent), "Independent", font = 2, cex = 0.9)
    text(r * cos(group_pos$Field), r * sin(group_pos$Field), "Field", font = 2, cex = 0.9)
  }) + coord_fixed()
  radar_plot_eelgrass
  
  
  radar_plot_surfgrass <- as.ggplot(function() {
    par(
      mar = c(2.5, 2.5, 2.5, 2.5),  # much tighter
      oma = c(2, 2, 2, 2),
      xpd = NA
    )
    # par(
    #   mar = c(6, 6, 6, 6),
    #   oma = c(2, 2, 2, 2),
    #   xpd = NA
    # )
    radarchart(
      radar_ready_surf,
      axistype = 0,
      caxislabels = NULL,
      pcol = plot_cols_surf,
      plwd = 2,
      plty = 1,
      cglcol = "grey85",
      cglty = 1,
      cglwd = 0.8,
      vlcex = 0.6
    )
    
    n_axes <- ncol(radar_ready_surf) - 1
    angles <- seq(pi/2, 2*pi + pi/2, length.out = n_axes + 1)[-(n_axes + 1)]
    
    group_pos <- list(
      Spatial = mean(angles[1:4]),
      Temporal = mean(angles[6:9]),
      Independent = mean(angles[11:14]),
      Field = mean(angles[16:19])
    )
    
    r <- 1.7
    
    text(r * cos(group_pos$Spatial), r * sin(group_pos$Spatial), "Spatial", font = 2, cex = 0.9)
    text(r * cos(group_pos$Temporal), r * sin(group_pos$Temporal), "Temporal", font = 2, cex = 0.9)
    text(r * cos(group_pos$Independent), r * sin(group_pos$Independent), "Independent", font = 2, cex = 0.9)
    text(r * cos(group_pos$Field), r * sin(group_pos$Field), "Field", font = 2, cex = 0.9)
  }) + coord_fixed()
  radar_plot_surfgrass
  
  
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
      cellwidth = 40,
      cellheight = 24,
      fontsize = 10,
      fontsize_row = 9,
      fontsize_col = 10
    )
  })
  
    heat_plot_surfgrass <- as.ggplot(function() {
      pheatmap(
        heat_df_cat_surf,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        treeheight_row = 0,
        treeheight_col = 0,
        legend = FALSE,
        color = my_cols,
        breaks = breaks,
        cellwidth = 40,
        cellheight = 24,
        fontsize = 10,
        fontsize_row = 9,
        fontsize_col = 10
      )
    })
  
    
    radar_plot_eelgrass  <- radar_plot_eelgrass  + theme(plot.margin = margin(1, 1, 1, 1))
    radar_plot_surfgrass <- radar_plot_surfgrass + theme(plot.margin = margin(1, 1, 1, 1))
    heat_plot_eelgrass  <- heat_plot_eelgrass  + theme(plot.margin = margin(1, 1, 1, 1))
    heat_plot_surfgrass <- heat_plot_surfgrass + theme(plot.margin = margin(1, 1, 1, 1))
    # Add species labels as separate title panels
    eelgrass_title <- ggdraw() +
      draw_label("Eelgrass", fontface = "bold", size = 14, x = 0, hjust = 0)
    surfgrass_title <- ggdraw() +
      draw_label("Surfgrass", fontface = "bold", size = 14, x = 0, hjust = 0)
    # Row 1: eelgrass radar + heatmap
    eelgrass_row <- plot_grid(
      radar_plot_eelgrass,
      heat_plot_eelgrass,
      ncol = 2,
      rel_widths = c(1.25, 1),
      align = "h"
    )
    # Row 2: surfgrass radar + heatmap
    surfgrass_row <- plot_grid(
      radar_plot_surfgrass,
      heat_plot_surfgrass,
      ncol = 2,
      rel_widths = c(1.25, 1),
      align = "h"
    )
    # Rebuild radar legend
    radar_leg_plot <- ggplot(
      data.frame(Model = model_order),
      aes(x = Model, y = 1, color = Model)
    ) +
      geom_point(size = 3) +
      scale_color_manual(
        values = model_cols,
        name = "Model",
        guide = guide_legend(ncol = 4, byrow = TRUE)
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    radar_leg <- cowplot::get_legend(radar_leg_plot)
    # Rebuild heatmap legend
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
        name = "Relative performance",
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          barwidth = unit(5.5, "cm"),
          barheight = unit(0.35, "cm")
        )
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    heat_leg <- cowplot::get_legend(heat_leg_plot)
    legend_row <- plot_grid(
      radar_leg,
      heat_leg,
      ncol = 2,
      rel_widths = c(1.4, 0.75)
    )
    # Final figure
    final_plot <- plot_grid(
      eelgrass_title,
      eelgrass_row,
      surfgrass_title,
      surfgrass_row,
      legend_row,
      ncol = 1,
      rel_heights = c(0.05, 1, 0.05, 1, 0.18)
    )
    final_plot
    
    
    ggsave(
      filename = "figures/model_ranking_radar_heatmap.png",
      plot = final_plot,
      width = 8.6,
      height = 9,
      dpi = 300,
      bg = "white"
    )
    
    
   
  # For eelgrass probably pick XGBOOST NEP, even thought it doesn't do well temporally but should do well for predicting rpeent day eelgrass.  
  #For surfgrass pick bccm spatial, it does good in every validation type. XGBOOST nep does best in field , but bccm_spatial is not horrible
  
  