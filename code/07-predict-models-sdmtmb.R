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
#

# STILL TO DO, MAKE PREDICTIONS WITH SURVEY AND YEAR
###############################################################################
#### Load modelling functions ####
source("code/modelling-functions.R")

####load packages####
UsePackages(c("sdmTMB", "sdmTMBextra", "tidyverse", "sf", "future", "terra", "future.apply"))


####load prediction data####
load("code/output_data/prediction_model_inputs.RData")
load("code/output_data/seagrass_model_inputs.RData")
# make predictions and get standard error
# fmodel_e_bccm_nospatial
# testing forst for model with no year, just survey random effect
#env_20m_all_survey <- replicate_df(env_20m_all, "Survey", unique(data$Survey))

survey_type <- c("ABL", "BHM", "Cuk", "GDK", "GSU", "MSE", "Mul", "RSU")



#This below was just for presentation, need to update with random factors 
eelgrass_predictions <- env_20m_all

e_bin <- predict(fmodel_e_bccm_spatial, newdata = env_20m_all)
e_per <- predict(m_e_per_bccm_final, newdata = env_20m_all)
e_bin_prob <- fmodel_e_bccm_spatial$family$linkinv(e_bin$est)
e_per_exp <- m_e_per_bccm_final$family$linkinv(e_per$est)
eelgrass_predictions$est_exp <- e_bin_prob * e_per_exp


#just grab a quarter to test
env_20m_all_fhalf <- env_20m_all %>% filter(ID < 8249999)
#env_20m_all_shalf <- env_20m_all %>% filter(ID > 8249999)

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

ggplot((eelgrass_predictions), aes(X, Y, fill = median)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(trans = "sqrt")

eelgrass_plot <- ggplot(eelgrass_predictions)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=median_binary, width=20,height=20))+
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
eelgrass_plot
ggsave("./figures/eelgrass.png", height = 6, width = 6)


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

eelgrass_se_plot <- ggplot(eelgrass_predictions)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=cv_binary, width=20,height=20))+
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
eelgrass_se_plot
ggsave("./figures/eelgrass_se.png", height = 6, width = 6)

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



eelgrass_raster_qcs <- eelgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, median_binary)
eelgrass_raster_qcs <- rast(x = eelgrass_raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_qcs, file.path("./raster/eelgrass_predictions_qcs_binary.tif"), overwrite=TRUE)

eelgrass_raster_qcs_percent <- eelgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, median)
eelgrass_raster_qcs_percent <- rast(x = eelgrass_raster_qcs_percent %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_qcs_percent, file.path("./raster/eelgrass_predictions_qcs_percent.tif"), overwrite=TRUE)





# Run function
pred <- PredictSDM(
  env = env_20m_all_fhalf,
  model = fmodel_e_bccm_nospatial,
  survey_type = survey_type,
  years = 2013:2023,
  species = "eelgrass"
)


# will need to recombine predictions to full study area


# change to 0-1 away from log-odds (logit) space
eelgrass_predictions <- eelgrass_predictions %>%
  mutate(est_p = plogis(est))

mean_preds <- mean_preds %>%
  mutate(est_p = plogis(est))

save(eelgrass_predictions, file = "code/output_data/eelgrass_predictions.RData")

####Plots####         
eelgrass_plot <- ggplot(eelgrass_predictions)+
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
eelgrass_plot
ggsave("./figures/eelgrass.png", height = 6, width = 6)

eelgrass_se_plot <- ggplot(eelgrass_predictions)+
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
eelgrass_se_plot
ggsave("./figures/eelgrass_se.png", height = 6, width = 6)


#### save rasters ####
#this changes to 100 m resolution though
# eelgrass_raster <- eelgrass_predictions %>%
#   mutate(x_round = round(X_m, -2), y_round = round(Y_m, -2)) %>%
#   select(x_round, y_round, est_p)
# 
# eelgrass_raster <- rast(x = eelgrass_raster %>% as.matrix, type = "xyz", crs = "EPSG:3005")
# writeRaster(eelgrass_raster, file.path("./raster/eelgrass_predictions.tif"), overwrite=TRUE)

eelgrass_raster_hg <- eelgrass_predictions %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, est_p)
eelgrass_raster_hg <- rast(x = eelgrass_raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_hg, file.path("./raster/eelgrass_predictions_hg.tif"), overwrite=TRUE)

eelgrass_raster_ss <- eelgrass_predictions %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, est_p)
eelgrass_raster_ss <- rast(x = eelgrass_raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_ss, file.path("./raster/eelgrass_predictions_ss.tif"), overwrite=TRUE)

eelgrass_raster_wcvi <- eelgrass_predictions %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, est_p)
eelgrass_raster_wcvi <- rast(x = eelgrass_raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_wcvi, file.path("./raster/eelgrass_predictions_wcvi.tif"), overwrite=TRUE)

eelgrass_raster_ncc <- eelgrass_predictions %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, est_p)
eelgrass_raster_ncc <- rast(x = eelgrass_raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_ncc, file.path("./raster/eelgrass_predictions_ncc.tif"), overwrite=TRUE)

eelgrass_raster_qcs <- eelgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, est_p)
eelgrass_raster_qcs <- rast(x = eelgrass_raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_qcs, file.path("./raster/eelgrass_predictions_qcs.tif"), overwrite=TRUE)

eelgrass_raster_hg_se <- eelgrass_predictions %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, SE)
eelgrass_raster_hg_se <- rast(x = eelgrass_raster_hg_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_hg_se, file.path("./raster/eelgrass_predictions_hg_se.tif"), overwrite=TRUE)

eelgrass_raster_ss_se <- eelgrass_predictions %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, SE)
eelgrass_raster_ss_se <- rast(x = eelgrass_raster_ss_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_ss_se, file.path("./raster/eelgrass_predictions_ss_se.tif"), overwrite=TRUE)

eelgrass_raster_wcvi_se <- eelgrass_predictions %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, SE)
eelgrass_raster_wcvi_se <- rast(x = eelgrass_raster_wcvi_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_wcvi_se, file.path("./raster/eelgrass_predictions_wcvi_se.tif"), overwrite=TRUE)

eelgrass_raster_ncc_se <- eelgrass_predictions %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, SE)
eelgrass_raster_ncc_se <- rast(x = eelgrass_raster_ncc_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_ncc_se, file.path("./raster/eelgrass_predictions_ncc_se.tif"), overwrite=TRUE)

eelgrass_raster_qcs_se <- eelgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, SE)
eelgrass_raster_qcs_se <- rast(x = eelgrass_raster_qcs_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_qcs_se, file.path("./raster/eelgrass_predictions_qcs_se.tif"), overwrite=TRUE)

eelgrass_raster_hg_sd <- eelgrass_predictions %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, SD)
eelgrass_raster_hg_sd <- rast(x = eelgrass_raster_hg_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_hg_sd, file.path("./raster/eelgrass_predictions_hg_sd.tif"), overwrite=TRUE)

eelgrass_raster_ss_sd <- eelgrass_predictions %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, SD)
eelgrass_raster_ss_sd <- rast(x = eelgrass_raster_ss_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_ss_sd, file.path("./raster/eelgrass_predictions_ss_sd.tif"), overwrite=TRUE)

eelgrass_raster_wcvi_sd <- eelgrass_predictions %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, SD)
eelgrass_raster_wcvi_sd <- rast(x = eelgrass_raster_wcvi_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_wcvi_sd, file.path("./raster/eelgrass_predictions_wcvi_sd.tif"), overwrite=TRUE)

eelgrass_raster_ncc_sd <- eelgrass_predictions %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, SD)
eelgrass_raster_ncc_sd <- rast(x = eelgrass_raster_ncc_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_ncc_sd, file.path("./raster/eelgrass_predictions_ncc_sd.tif"), overwrite=TRUE)

eelgrass_raster_qcs_sd <- eelgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, SD)
eelgrass_raster_qcs_sd <- rast(x = eelgrass_raster_qcs_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster_qcs_sd, file.path("./raster/eelgrass_predictions_qcs_sd.tif"), overwrite=TRUE)





## surfgrass


#make predictions and get SE
hold <- predict(fmodel_s_bccm, env_20m_all)
sims <- predict(fmodel_s_bccm, newdata = env_20m_all, nsim = 100) #sim needs to be 500? ram is not working for this right now at 20m prediction cells
hold$SE <- apply(sims, 1, sd)
hold$SD <- apply(pclog(sims), 1, sd)
surfgrass_predictions <- env_20m_all
surfgrass_predictions <- bind_cols(surfgrass_predictions, hold %>% select(est:SD))

# hold <- surfgrass_predictions %>%
#   group_by(X_m, Y_m) %>%
#   summarise(CV = sd(exp(est), na.rm = TRUE)/mean(exp(est)), counts_ln = mean(est), SE = mean(SE))

# change to 0-1 away from log-odds (logit) space
surfgrass_predictions <- surfgrass_predictions %>%
  mutate(est_p = plogis(est))

#save outputs####
save(surfgrass_predictions, file = "code/output_data/surfgrass_predictions.RData")


####Plots####         


surfgrass_plot <- ggplot(surfgrass_predictions)+
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
surfgrass_plot
ggsave("./figures/surfgrass.png", height = 6, width = 6)

surfgrass_se_plot <- ggplot(surfgrass_predictions)+
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
surfgrass_se_plot
ggsave("./figures/surfgrass_se.png", height = 6, width = 6)

# refer to https://pbs-assess.github.io/sdmTMB/articles/basic-intro.html to make plots of random spatial fields etc




#save as raster, this changes to 100 m resolution though
# surfgrass_raster <- surfgrass_predictions %>%
#   mutate(x_round = round(X_m, -2), y_round = round(Y_m, -2)) %>%
#   select(x_round, y_round, est_p)
# 
# surfgrass_raster <- rast(x = surfgrass_raster %>% as.matrix, type = "xyz", crs = "EPSG:3005")
# writeRaster(surfgrass_raster, file.path("./raster/surfgrass_predictions.tif"), overwrite=TRUE)

surfgrass_raster_hg <- surfgrass_predictions %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, est_p)
surfgrass_raster_hg <- rast(x = surfgrass_raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_hg, file.path("./raster/surfgrass_predictions_hg.tif"), overwrite=TRUE)

surfgrass_raster_ss <- surfgrass_predictions %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, est_p)
surfgrass_raster_ss <- rast(x = surfgrass_raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_ss, file.path("./raster/surfgrass_predictions_ss.tif"), overwrite=TRUE)

surfgrass_raster_wcvi <- surfgrass_predictions %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, est_p)
surfgrass_raster_wcvi <- rast(x = surfgrass_raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_wcvi, file.path("./raster/surfgrass_predictions_wcvi.tif"), overwrite=TRUE)

surfgrass_raster_ncc <- surfgrass_predictions %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, est_p)
surfgrass_raster_ncc <- rast(x = surfgrass_raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_ncc, file.path("./raster/surfgrass_predictions_ncc.tif"), overwrite=TRUE)

surfgrass_raster_qcs <- surfgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, est_p)
surfgrass_raster_qcs <- rast(x = surfgrass_raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_qcs, file.path("./raster/surfgrass_predictions_qcs.tif"), overwrite=TRUE)

surfgrass_raster_hg_se <- surfgrass_predictions %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, SE)
surfgrass_raster_hg_se <- rast(x = surfgrass_raster_hg_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_hg_se, file.path("./raster/surfgrass_predictions_hg_se.tif"), overwrite=TRUE)

surfgrass_raster_ss_se <- surfgrass_predictions %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, SE)
surfgrass_raster_ss_se <- rast(x = surfgrass_raster_ss_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_ss_se, file.path("./raster/surfgrass_predictions_ss_se.tif"), overwrite=TRUE)

surfgrass_raster_wcvi_se <- surfgrass_predictions %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, SE)
surfgrass_raster_wcvi_se <- rast(x = surfgrass_raster_wcvi_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_wcvi_se, file.path("./raster/surfgrass_predictions_wcvi_se.tif"), overwrite=TRUE)

surfgrass_raster_ncc_se <- surfgrass_predictions %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, SE)
surfgrass_raster_ncc_se <- rast(x = surfgrass_raster_ncc_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_ncc_se, file.path("./raster/surfgrass_predictions_ncc_se.tif"), overwrite=TRUE)

surfgrass_raster_qcs_se <- surfgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, SE)
surfgrass_raster_qcs_se <- rast(x = surfgrass_raster_qcs_se %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_qcs_se, file.path("./raster/surfgrass_predictions_qcs_se.tif"), overwrite=TRUE)

surfgrass_raster_hg_sd <- surfgrass_predictions %>%
  filter(region == "Haida Gwaii") %>%
  select(X_m, Y_m, SD)
surfgrass_raster_hg_sd <- rast(x = surfgrass_raster_hg_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_hg_sd, file.path("./raster/surfgrass_predictions_hg_sd.tif"), overwrite=TRUE)

surfgrass_raster_ss_sd <- surfgrass_predictions %>%
  filter(region == "Salish Sea") %>%
  select(X_m, Y_m, SD)
surfgrass_raster_ss_sd <- rast(x = surfgrass_raster_ss_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_ss_sd, file.path("./raster/surfgrass_predictions_ss_sd.tif"), overwrite=TRUE)

surfgrass_raster_wcvi_sd <- surfgrass_predictions %>%
  filter(region == "West Coast Vancouver Island") %>%
  select(X_m, Y_m, SD)
surfgrass_raster_wcvi_sd <- rast(x = surfgrass_raster_wcvi_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_wcvi_sd, file.path("./raster/surfgrass_predictions_wcvi_sd.tif"), overwrite=TRUE)

surfgrass_raster_ncc_sd <- surfgrass_predictions %>%
  filter(region == "North Central Coast") %>%
  select(X_m, Y_m, SD)
surfgrass_raster_ncc_sd <- rast(x = surfgrass_raster_ncc_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_ncc_sd, file.path("./raster/surfgrass_predictions_ncc_sd.tif"), overwrite=TRUE)

surfgrass_raster_qcs_sd <- surfgrass_predictions %>%
  filter(region == "Queen Charlotte Strait") %>%
  select(X_m, Y_m, SD)
surfgrass_raster_qcs_sd <- rast(x = surfgrass_raster_qcs_sd %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(surfgrass_raster_qcs_sd, file.path("./raster/surfgrass_predictions_qcs_sd.tif"), overwrite=TRUE)


