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
    FIELD_TSS = norm_pos(TSS_Field),
    FIELD_SenCV = norm_pos(Sensitivity_FieldCV),
    FIELD_SpeCV = norm_pos(Specificity_FieldCV),
    FIELD_TSSCV = norm_pos(TSS_FieldCV),
    FIELD_CliffsD = norm_pos(Cliffs_delta)
  )
# Compute final composite score
rank_table_scaled$FinalScore <- rowMeans(rank_table_scaled %>%
                                           select(CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
                                                  TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
                                                  INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
                                                  FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV,
                                                  FIELD_TSSCV, FIELD_CliffsD),
                                         na.rm = TRUE)
# Rank models based on the final score
final_rank <- rank_table_scaled %>%
  arrange(desc(FinalScore)) %>%
  select(Model, FinalScore)



# Radar plots
radar_df <- rank_table_scaled %>%
  select(Model, CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
         TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
         INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
         FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV,
         FIELD_TSSCV, FIELD_CliffsD)
# Prepare for radar chart
radar_ready <- rbind(
  rep(1, ncol(radar_df)-1),
  rep(0, ncol(radar_df)-1),
  radar_df[, -1]
)
rownames(radar_ready) <- c("max", "min", radar_df$Model)

palette <- brewer.pal(n = nrow(radar_df), name = "Set3")  # Using Set3 for high visibility

# Generate radar chart
radarchart(radar_ready,
           plwd = 3,
           pcol = palette,  # Use rainbow colors for different models
           plty = 1)
legend("topright", 
       legend = radar_df$Model,  # Names of the models
       col = palette,  # Same colors as used in the plot
       pch = 20,  # Point type (filled circle)
       bty = "n",  # No box around the legend
       pt.cex = 2,  # Point size
       cex = 0.8)  # Font size


# Heat map
heat_df <- rank_table_scaled %>%
  select(Model, CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
         TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
         INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
         FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV,
         FIELD_TSSCV, FIELD_CliffsD, FinalScore) %>%
  arrange(desc(FinalScore)) %>%
  select(-FinalScore) %>%
  tibble::column_to_rownames("Model")
# Create heatmap
  pheatmap(heat_df, 
           cluster_rows = FALSE,  # Do not cluster the rows (models)
           cluster_cols = FALSE, 
           show_rownames = TRUE,    # Show row names
           treeheight_row = 0,
           treeheight_col = 0)
  
  
  # Also test if we make a score based on each of the 4 validaitons. 
  
  # Compute final composite scores for each category
  rank_table_scaled_cat <- rank_table_scaled %>%
    mutate(
      CV_Score = rowMeans(select(., CV_AUC, CV_Tjur, CV_RMSE, CV_TSS), na.rm = TRUE),
      TEMP_Score = rowMeans(select(., TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS), na.rm = TRUE),
      INDEP_Score = rowMeans(select(., INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI), na.rm = TRUE),
      FIELD_Score = rowMeans(select(., FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV, FIELD_TSSCV, FIELD_CliffsD), na.rm = TRUE)
    )
  # Rank models based on the final score for each category
  final_rank_cat <- rank_table_scaled_cat %>%
    arrange(desc(CV_Score)) %>%
    select(Model, CV_Score, TEMP_Score, INDEP_Score, FIELD_Score)
  # Prepare data for radar plot
  radar_df_cat <- final_rank_cat %>%
    select(Model, CV_Score, TEMP_Score, INDEP_Score, FIELD_Score)
  # Prepare for radar chart
  radar_ready_cat <- rbind(
    rep(1, ncol(radar_df_cat) - 1),  # max range
    rep(0, ncol(radar_df_cat) - 1),  # min range
    radar_df_cat[, -1]               # the models to plot
  )
  rownames(radar_ready_cat) <- c("max", "min", radar_df_cat$Model)
  # Choose a color palette with sufficient distinct colors
  library(RColorBrewer)
  palette <- brewer.pal(n = nrow(radar_df_cat), name = "Set3")
  
  # Generate radar chart
  radarchart(radar_ready_cat,
             plwd = 3,  # Set the line width to 3 (thicker)
             pcol = palette,  # Use the custom color palette
             plty = 1)
  # Adding Legend
  legend("topright", 
         legend = radar_df_cat$Model,  # Names of the models
         col = palette,  # Use the same colors for the legend
         pch = 20,  # Point type (filled circle)
         bty = "n",  # No box around the legend
         pt.cex = 2,  # Point size
         cex = 0.8)  # Font size
  
  # Heat map
  heat_df_cat <- rank_table_scaled_cat %>%
    select(Model, Model, CV_Score, TEMP_Score, INDEP_Score, FIELD_Score, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  # Create heatmap
  pheatmap(heat_df_cat, 
           cluster_rows = FALSE,  # Do not cluster the rows (models)
           cluster_cols = FALSE, 
           show_rownames = TRUE,    # Show row names
           treeheight_row = 0,
           treeheight_col = 0)

  
  # For eelgrass probably pick XGBOOST NEP, even thought it doesn't do well temporally but should do well for predicting rpeent day eelgrass.  
  
  
  
  
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
      FIELD_TSS = norm_pos(TSS_Field),
      FIELD_SenCV = norm_pos(Sensitivity_FieldCV),
      FIELD_SpeCV = norm_pos(Specificity_FieldCV),
      FIELD_TSSCV = norm_pos(TSS_FieldCV),
      FIELD_CliffsD = norm_pos(Cliffs_delta)
    )
  # Compute final composite score
  rank_table_scaled_surf$FinalScore <- rowMeans(rank_table_scaled_surf %>%
                                             select(CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
                                                    TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
                                                    INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
                                                    FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV,
                                                    FIELD_TSSCV, FIELD_CliffsD),
                                           na.rm = TRUE)
  # Rank models based on the final score
  final_rank_surf <- rank_table_scaled_surf %>%
    arrange(desc(FinalScore)) %>%
    select(Model, FinalScore)
  
  
  
  # Radar plots
  radar_df_surf <- rank_table_scaled_surf %>%
    select(Model, CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
           TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
           INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
           FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV,
           FIELD_TSSCV, FIELD_CliffsD)
  # Prepare for radar chart
  radar_ready_surf <- rbind(
    rep(1, ncol(radar_df_surf)-1),
    rep(0, ncol(radar_df_surf)-1),
    radar_df_surf[, -1]
  )
  rownames(radar_ready_surf) <- c("max", "min", radar_df_surf$Model)
  
  palette <- brewer.pal(n = nrow(radar_df_surf), name = "Set3")  # Using Set3 for high visibility
  
  # Generate radar chart
  radarchart(radar_ready_surf,
             plwd = 3,
             pcol = palette,  # Use rainbow colors for different models
             plty = 1)
  legend("topright", 
         legend = radar_df_surf$Model,  # Names of the models
         col = palette,  # Same colors as used in the plot
         pch = 20,  # Point type (filled circle)
         bty = "n",  # No box around the legend
         pt.cex = 2,  # Point size
         cex = 0.8)  # Font size
  
  
  # Heat map
  heat_df_surf <- rank_table_scaled_surf %>%
    select(Model, CV_AUC, CV_Tjur, CV_RMSE, CV_TSS,
           TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS,
           INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI, 
           FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV,
           FIELD_TSSCV, FIELD_CliffsD, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  # Create heatmap
  pheatmap(heat_df_surf, 
           cluster_rows = FALSE,  # Do not cluster the rows (models)
           cluster_cols = FALSE, 
           show_rownames = TRUE,    # Show row names
           treeheight_row = 0,
           treeheight_col = 0)
  
  
  # Also test if we make a score based on each of the 4 validaitons. 
  
  # Compute final composite scores for each category
  rank_table_scaled_surf_cat <- rank_table_scaled_surf %>%
    mutate(
      CV_Score = rowMeans(select(., CV_AUC, CV_Tjur, CV_RMSE, CV_TSS), na.rm = TRUE),
      TEMP_Score = rowMeans(select(., TEMP_AUC, TEMP_Tjur, TEMP_RMSE, TEMP_TSS), na.rm = TRUE),
      INDEP_Score = rowMeans(select(., INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI), na.rm = TRUE),
      FIELD_Score = rowMeans(select(., FIELD_AUC, FIELD_TSS, FIELD_SenCV, FIELD_SpeCV, FIELD_TSSCV, FIELD_CliffsD), na.rm = TRUE)
    )
  # Rank models based on the final score for each category
  final_rank_surf_cat <- rank_table_scaled_surf_cat %>%
    arrange(desc(CV_Score)) %>%
    select(Model, CV_Score, TEMP_Score, INDEP_Score, FIELD_Score)
  # Prepare data for radar plot
  radar_df_surf_cat <- final_rank_surf_cat %>%
    select(Model, CV_Score, TEMP_Score, INDEP_Score, FIELD_Score)
  # Prepare for radar chart
  radar_ready_surf_cat <- rbind(
    rep(1, ncol(radar_df_surf_cat) - 1),  # max range
    rep(0, ncol(radar_df_surf_cat) - 1),  # min range
    radar_df_surf_cat[, -1]               # the models to plot
  )
  rownames(radar_ready_surf_cat) <- c("max", "min", radar_df_surf_cat$Model)
  # Choose a color palette with sufficient distinct colors
  library(RColorBrewer)
  palette <- brewer.pal(n = nrow(radar_df_surf_cat), name = "Set3")
  
  # Generate radar chart
  radarchart(radar_ready_surf_cat,
             plwd = 3,  # Set the line width to 3 (thicker)
             pcol = palette,  # Use the custom color palette
             plty = 1)
  # Adding Legend
  legend("topright", 
         legend = radar_df_surf_cat$Model,  # Names of the models
         col = palette,  # Use the same colors for the legend
         pch = 20,  # Point type (filled circle)
         bty = "n",  # No box around the legend
         pt.cex = 2,  # Point size
         cex = 0.8)  # Font size
  
  # Heat map
  heat_df_surf_cat <- rank_table_scaled_surf_cat %>%
    select(Model, Model, CV_Score, TEMP_Score, INDEP_Score, FIELD_Score, FinalScore) %>%
    arrange(desc(FinalScore)) %>%
    select(-FinalScore) %>%
    tibble::column_to_rownames("Model")
  # Create heatmap
  pheatmap(heat_df_surf_cat, 
           cluster_rows = FALSE,  # Do not cluster the rows (models)
           cluster_cols = FALSE, 
           show_rownames = TRUE,    # Show row names
           treeheight_row = 0,
           treeheight_col = 0)
  
  
  # For surfgrass pick bccm spatial, it does good in every validation type. XGBOOST nep does best in field , but bccm_spatial is not horrible