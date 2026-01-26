###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# predict sdm models with sdmTMB 


###############################################################################
#### Load modelling functions ####
source("code/modelling-functions.R")

####load packages####
UsePackages(c("sdmTMB", "sdmTMBextra", "tidyverse", "sf", "future", "terra", "future.apply"))


####load prediction data####
load("code/output_data/prediction_model_inputs.RData")
load("code/output_data/seagrass_model_inputs.RData")
load("code/output_data/model_results/final_eelgrass_model.RData")

# make predictions and get standard error

# Parameters to set
survey_type <- c("ABL", "BHM", "Cuk", "GDK", "GSU", "MSE", "Mul", "RSU")
model_name<- "nep_spatial"
species <- "eelgrass"

#Predict for each survey
PredictSDM_bySurvey(
  env = env_20m_all,
  model = fmodel_e_nep_spatial,
  survey_type = survey_type,
  species = species,
  model_name  = model_name
)

#mean predictions across all surveys
mean_pred <- AverageSurveyPredictions(
  species    = species,
  model_name = model_name
)

# change to 0-1 away from log-odds (logit) space
mean_pred <- mean_pred %>%
  mutate(est_p = plogis(est))

#### save rasters ####
outdir <- file.path("./raster", model_name)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- mean_pred %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, est_p)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0(species, "_predictions_hg_", model_name, ".tif")), overwrite = TRUE)

raster_ss <- mean_pred %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, est_p)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0(species, "_predictions_ss_", model_name, ".tif")), overwrite = TRUE)

raster_wcvi <- mean_pred %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, est_p)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0(species, "_predictions_wcvi_", model_name, ".tif")), overwrite = TRUE)

raster_ncc <- mean_pred %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, est_p)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0(species, "_predictions_ncc_", model_name, ".tif")), overwrite = TRUE)

raster_qcs <- mean_pred %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, est_p)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0(species, "_predictions_qcs_", model_name, ".tif")), overwrite = TRUE)

raster_hg_se <- mean_pred %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, SE)
raster_hg_se <- rast(x = raster_hg_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg_se, file.path(outdir, paste0(species, "_predictions_hg_se_", model_name, ".tif")), overwrite = TRUE)

raster_ss_se <- mean_pred %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, SE)
raster_ss_se <- rast(x = raster_ss_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss_se, file.path(outdir, paste0(species, "_predictions_ss_se_", model_name, ".tif")), overwrite = TRUE)

raster_wcvi_se <- mean_pred %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, SE)
raster_wcvi_se <- rast(x = raster_wcvi_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi_se, file.path(outdir, paste0(species, "_predictions_wcvi_se_", model_name, ".tif")), overwrite = TRUE)

raster_ncc_se <- mean_pred %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, SE)
raster_ncc_se <- rast(x = raster_ncc_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc_se, file.path(outdir, paste0(species, "_predictions_ncc_se_", model_name, ".tif")), overwrite = TRUE)

raster_qcs_se <- mean_pred %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, SE)
raster_qcs_se <- rast(x = raster_qcs_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs_se, file.path(outdir, paste0(species, "_predictions_qcs_se_", model_name, ".tif")), overwrite = TRUE)


##### plots ####
#these take a long time to make 
ggplot((mean_pred), aes(X, Y, fill = median)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(trans = "sqrt")

prediction_plot <- ggplot(mean_pred)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=est_p, width=20,height=20))+
  scale_colour_gradient(low = "#f7fcb9", high = "#006837")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_sf(expand = FALSE)+
  ylab("")+
  xlab("") 
prediction_plot
ggsave(paste0("./figures/", species, ".png", height = 6, width = 6))

se_plot <- ggplot(mean_pred)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=SE, width=20,height=20))+
  scale_colour_viridis_b()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_sf(expand = FALSE)+
  ylab("")+
  xlab("") 
se_plot
ggsave(paste0("./figures/", species, "_se.png", height = 6, width = 6))


# refer to https://pbs-assess.github.io/sdmTMB/articles/basic-intro.html to make plots of random spatial fields etc










#NEED TO WORK ON THIS STUFF BELOW!

#Percent cover eelgrass predictions

eelgrass_predictions <- env_20m_all

e_bin <- predict(fmodel_e_bccm_spatial, newdata = env_20m_all)
e_per <- predict(m_e_per_bccm_final, newdata = env_20m_all)
e_bin_prob <- fmodel_e_bccm_spatial$family$linkinv(e_bin$est)
e_per_exp <- m_e_per_bccm_final$family$linkinv(e_per$est)
eelgrass_predictions$est_exp <- e_bin_prob * e_per_exp

set.seed(28239)
p_bin_sim <- predict(fmodel_e_bccm_spatial, newdata = env_20m_all, nsim = 100)
p_per_sim <- predict(m_e_per_bccm_final, newdata = env_20m_all, nsim = 100)
p_bin_prob_sim <- fmodel_e_bccm_spatial$family$linkinv(p_bin_sim)
p_per_exp_sim <- m_e_per_bccm_final$family$linkinv(p_per_sim)
p_combined_sim <- p_bin_prob_sim * p_per_exp_sim

eelgrass_predictions$median <- apply(p_combined_sim, 1, median)
eelgrass_predictions$median_binary <- apply(p_bin_prob_sim, 1, median)
#hold_all$SE <- apply(sims, 1, sd)
eelgrass_predictions$cv <- apply(p_combined_sim, 1, function(x) sd(x) / mean(x))
eelgrass_predictions$cv_binary <- apply(p_bin_prob_sim, 1, function(x) sd(x) / mean(x))

plot(eelgrass_predictions$est_exp, eelgrass_predictions$median)

eelgrass_plot_percent <- ggplot(eelgrass_predictions)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=median, width=20,height=20))+
  scale_colour_gradient(low = "#f7fcb9", high = "#006837")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_sf(expand = FALSE)+
  ylab("")+
  xlab("") 
eelgrass_plot_percent
ggsave("./figures/eelgrass_percent.png", height = 6, width = 6)

eelgrass_se_plot_percent <- ggplot(eelgrass_predictions)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=cv, width=20,height=20))+
  scale_colour_viridis_b()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_sf(expand = FALSE)+
  ylab("")+
  xlab("") 
eelgrass_se_plot_percent
ggsave("./figures/eelgrass_se_percent.png", height = 6, width = 6)