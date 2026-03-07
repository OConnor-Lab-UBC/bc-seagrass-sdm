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



model_names <- c(
  "bccm_nospatial",
  "bccm_spatial",
  "nep_nospatial",
  "nep_spatial"
)

cv_clean <- cv_summary %>%
  filter(algo %in% model_names) %>%
  select(algo, AUC_mean, Tjur, TSS_mean, RMSE) %>%
  rename(Model = algo)


temporal_clean <- forecast_summary %>%
  filter(model %in% model_names) %>%
  select(model, AUC_forecast, Tjur_forecast, TSS_forecast, RMSE_forecast) %>%
  rename(Model = model)

indep_clean <- independent_results %>%
  select(Model, MPS, FPPS, FNR, CBI)

field_clean <- field_validation %>%
  rename(Model = model)

#normalize metrics
norm_pos <- function(x) (x - min(x)) / (max(x) - min(x))

norm_neg <- function(x) (max(x) - x) / (max(x) - min(x))

#build unified table

# Join everything
rank_table <- cv_clean %>%
  left_join(temporal_clean, by = "Model") %>%
  left_join(indep_clean, by = "Model")

if (!is.null(field_clean)) {
  rank_table <- rank_table %>% left_join(field_clean, by = "Model")
}

# Scale metrics
rank_table_scaled <- rank_table %>%
  mutate(
    CV_AUC    = norm_pos(AUC_mean),
    CV_Tjur   = norm_pos(Tjur),
    CV_TSS    = norm_pos(TSS_mean),
    CV_RMSE   = norm_neg(RMSE),
    
    TEMP_AUC    = norm_pos(AUC_forecast),
    TEMP_Tjur   = norm_pos(Tjur_forecast),
    TEMP_TSS    = norm_pos(TSS_forecast),
    TEMP_RMSE   = norm_neg(RMSE_forecast),
    
    INDEP_MPS  = norm_pos(MPS),
    INDEP_FPPS = norm_pos(FPPS),
    INDEP_FNR  = norm_neg(FNR),
    INDEP_CBI  = norm_pos(CBI)
  )

#compute final composite score

rank_table_scaled$FinalScore <- rowMeans(rank_table_scaled %>%
                                           select(CV_AUC, CV_Tjur, CV_TSS, CV_RMSE,
                                                  TEMP_AUC, TEMP_Tjur, TEMP_TSS, TEMP_RMSE,
                                                  INDEP_MPS, INDEP_FPPS, INDEP_FNR, INDEP_CBI),
                                         na.rm = TRUE
)

#rank models
final_rank <- rank_table_scaled %>%
  arrange(desc(FinalScore)) %>%
  select(Model, FinalScore)


#radar plots

library(fmsb)

radar_df <- rank_table_scaled %>%
  select(Model, CV_AUC, TEMP_AUC, INDEP_MPS, INDEP_CBI)  # choose key metrics

# prepare for radar
radar_ready <- rbind(
  rep(1, ncol(radar_df)-1),
  rep(0, ncol(radar_df)-1),
  radar_df[, -1]
)
rownames(radar_ready) <- c("max", "min", radar_df$Model)

radarchart(radar_ready,
           plwd = 2,
           title = "Model Performance Summary")

#heat map
library(pheatmap)

heat_df <- rank_table_scaled %>%
  select(Model, CV_AUC, TEMP_AUC, INDEP_MPS, INDEP_CBI) %>%
  column_to_rownames("Model")

pheatmap(heat_df, cluster_rows=TRUE, cluster_cols=FALSE)
