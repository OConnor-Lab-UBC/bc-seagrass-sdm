###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
#
#
# Objective:
# ---------
# fit sdm models with sdmTMB 
#
###############################################################################


#load packages####
library(sdmTMB)
library(sdmTMBextra)
library(tidyverse)
library(sf)
library(future)
library(terra)


#load seagrass data
load("code/output_data/seagrass_model_inputs.RData")
load("code/output_data/prediction_model_inputs.RData")


#eelgrass model###
eelgrass <- filter(seagrass_data_long, species == "ZO")
print(paste("eelgrass present in ", round((sum(eelgrass$presence)/nrow(eelgrass))*100,2), "% of observations", sep = ""))

library(corrplot)
eelgrass_env<- eelgrass %>%
  select(presence, depth_stnd, rei_sqrt_stnd, tidal_sqrt_stnd, freshwater_sqrt_stnd, slope_sqrt_stnd, NH4_stnd, NO3_stnd, 
         saltmean_sq_stnd, saltmin_sq_stnd, PARmean_stnd, PARmin_stnd, PARmax_stnd, surftempmean_stnd, surftempmin_stnd, 
         surftempmax_stnd, tempmean_stnd, tempmin_stnd, tempmax_stnd, DOmean_stnd, DOmin_stnd)
corrplot(cor(eelgrass_env), method = "number")
# Correlations close to-1 or +1 might indicate the existence of multicollinearity. 
# As a rule of thumb, one might suspect multicollinearity when the correlation between two (predictor) variables is below -0.9 or above +0.9.
# correlated: PAR max and mean, salt min and mean, surftempmean and temp mean, surftempmax and tempmax, tempmean and temp max (but not temp min)
# keep PAR mean and min, salt min, 5 m temps max and min, NH4, NO3

## CREATE A LINEAR REGRESSION MODEL
#my_model <- lm(presence~., data = eelgrass_env)
my_model <- lm(presence~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + freshwater_sqrt_stnd + slope_sqrt_stnd +
                 NH4_stnd + NO3_stnd + saltmin_sq_stnd + PARmean_stnd + tempmin_stnd + tempmax_stnd +
                 DOmin_stnd, data = eelgrass_env)
library(olsrr)
ols_vif_tol(my_model)
# As a general guideline, a Tolerance of <0.1 might indicate multicollinearity.
# As a rule of thumb, a VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. 
# Ideally, the Variance Inflation Factors are below 3.
# there is no multi multicollinearity in variables selevered in my_model.

# Continue if prevalence is greater than 0.5%
substrates_present <- eelgrass %>%
    group_by(substrate) %>%
    summarise(n_present = sum(presence))
substrates_present 
# there are eelgrass presence observations on all substrates


#make mesh
mesh_eelgrass <- make_mesh(data = eelgrass, xy_cols = c("X", "Y"), cutoff = 15)
plot(mesh_eelgrass)

barrier_mesh_eelgrass <- add_barrier_mesh(mesh_eelgrass, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)


m_eelgrass_1 <- sdmTMB_cv(presence ~ s(depth_stnd, k = 3) + freshwater_sqrt_stnd + slope_sqrt_stnd + substrate + s(rei_sqrt_stnd, k=3) + s(tidal_sqrt_stnd, k = 3) + NH4_stnd + NO3_stnd + saltmin_sq_stnd + PARmean_stnd + tempmin_stnd + tempmax_stnd + DOmin_stnd,
                      mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = TRUE, data = eelgrass, fold_ids = "fold")
#  
m_eelgrass_1$sum_loglik
m_eelgrass_1$AIC
#eelgrass$ID<-as.factor(eelgrass$ID)

# m_eelgrass <- sdmTMB(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|ID),
#                           mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = "on", time = "Year", spatiotemporal = "IID", data = eelgrass)
# 
# m_eelgrass <- sdmTMB(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|HKey),
#                      mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = "on", data = eelgrass)

m_eelgrass <- sdmTMB(presence ~ s(depth_stnd, k = 3) + freshwater_sqrt_stnd + slope_sqrt_stnd + substrate + s(rei_sqrt_stnd, k=3) + s(tidal_sqrt_stnd, k = 3) +
                       NH4_stnd + NO3_stnd + saltmin_sq_stnd + PARmean_stnd + tempmin_stnd + tempmax_stnd + DOmin_stnd, mesh = barrier_mesh_eelgrass, 
                     family = binomial(link = "logit"), spatial = "on", data = eelgrass)

sanity(m_eelgrass)
#sdmTMB::run_extra_optimization(m_eelgrass, nlminb_loops = 1L, newton_loops = 1L)

m_eelgrass


pred_fixed <- m_eelgrass$family$linkinv(predict(m_eelgrass)$est_non_rf)
r_pois <- DHARMa::createDHARMa(
  simulatedResponse = simulate(m_eelgrass, nsim = 500),
  observedResponse = eelgrass$presence,
  fittedPredictedResponse = pred_fixed
)
plot(r_pois)

DHARMa::testDispersion(r_pois)


simulate(m_eelgrass, nsim = 500) %>%
  dharma_residuals(m_eelgrass)

predict(m_eelgrass) %>%
  ggplot(aes(x = presence, y = m_eelgrass$family$linkinv(est)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_jitter(width = 0.05, height = 0)

hold <- predict(m_eelgrass, env_20m_all)
sims <- predict(m_eelgrass, newdata = env_20m_all, nsim = 500) #ram is not working for this right now at 20m prediction cells
hold$SE <- apply(sims, 1, sd)
eelgrass_predictions <- env_20m_all
eelgrass_predictions <- bind_cols(eelgrass_predictions, hold %>% select(est:SE))
#eelgrass_predictions <- bind_cols(eelgrass_predictions, hold %>% select(est:omega_s))

eelgrass_predictions <- eelgrass_predictions %>%
  mutate(est_p = plogis(est))
         
         
ggplot(eelgrass_predictions)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=est_p, width=20,height=20))+
  scale_colour_viridis_c()+
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

ggsave("./figures/eelgrass.png", height = 6, width = 20)

# hold <- eelgrass_predictions %>%
#   group_by(X_m, Y_m) %>%
#   summarise(CV = sd(exp(est), na.rm = TRUE)/mean(exp(est)), counts_ln = mean(est), SE = mean(SE))
# 

#save as raster, this changes to 100 m resolution though
eelgrass_raster <- eelgrass_predictions %>%
  mutate(x_round = round(X_m, -2), y_round = round(Y_m, -2)) %>%
  select(x_round, y_round, est_p)

eelgrass_raster <- rast(x = eelgrass_raster %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster, file.path("./raster/eelgrass_predictions.tif"), overwrite=TRUE)

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




surfgrass <- filter(seagrass_data_long, species == "PH")
print(paste("surfgrass present in ", round((sum(surfgrass$presence)/nrow(surfgrass))*100,2), "% of transects", sep = ""))

# Continue if prevalence is greater than 0.5%
substrates_present <- surfgrass %>%
  group_by(substrate) %>%
  summarise(n_present = sum(presence))
substrates_present 
# there are surfgrass presence observations on all substrates


#make mesh
mesh_surfgrass <- make_mesh(data = surfgrass, xy_cols = c("X", "Y"), cutoff = 15)
plot(mesh_surfgrass)
barrier_mesh_surfgrass <- add_barrier_mesh(mesh_surfgrass, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)


# m_surfgrass_1 <- sdmTMB_cv(presence ~1 + s(depth_stnd, k = 3) + slope_stnd + s(temperature_stnd, k = 3) + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3),
#                       mesh = barrier_mesh_surfgrass, family = binomial(link = "logit"), spatial = TRUE, data = surfgrass, fold_ids = "fold")
#  
# m_surfgrass_1$sum_loglik
#surfgrass$ID<-as.factor(surfgrass$ID)

# m_surfgrass <- sdmTMB(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|ID),
#                           mesh = barrier_mesh_surfgrass, family = binomial(link = "logit"), spatial = "on", time = "Year", spatiotemporal = "IID", data = surfgrass)
# 
# m_surfgrass <- sdmTMB(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|HKey),
#                      mesh = barrier_mesh_surfgrass, family = binomial(link = "logit"), spatial = "on", data = surfgrass)

#remove freshwater compared to eelgrass
m_surfgrass <- sdmTMB(presence ~ 1 + s(depth_stnd, k = 3) + slope_sqrt_stnd + substrate + s(rei_sqrt_stnd, k=3) + s(tidal_sqrt_stnd, k = 3) +
                        NH4_stnd + NO3_stnd + saltmin_sq_stnd + PARmean_stnd + tempmin_stnd + tempmax_stnd + DOmin_stnd,
                     mesh = barrier_mesh_surfgrass, family = binomial(link = "logit"), spatial = "on", data = surfgrass)

sanity(m_surfgrass)
#sdmTMB::run_extra_optimization(m_surfgrass, nlminb_loops = 1L, newton_loops = 1L)

m_surfgrass


pred_fixed <- m_surfgrass$family$linkinv(predict(m_surfgrass)$est_non_rf)
r_pois <- DHARMa::createDHARMa(
  simulatedResponse = simulate(m_surfgrass, nsim = 500),
  observedResponse = surfgrass$presence,
  fittedPredictedResponse = pred_fixed
)
plot(r_pois)

DHARMa::testDispersion(r_pois)


simulate(m_surfgrass, nsim = 500) %>%
  dharma_residuals(m_surfgrass)

predict(m_surfgrass) %>%
  ggplot(aes(x = presence, y = m_surfgrass$family$linkinv(est)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_jitter(width = 0.05, height = 0)

hold <- predict(m_surfgrass, env_20m_all)
sims <- predict(m_surfgrass, newdata = env_20m_all, nsim = 500) #ram is not working for this right now at 20m prediction cells
hold$SE <- apply(sims, 1, sd)
surfgrass_predictions <- env_20m_all
surfgrass_predictions <- bind_cols(surfgrass_predictions, hold %>% select(est:SE))
#surfgrass_predictions <- bind_cols(surfgrass_predictions, hold %>% select(est:omega_s))

surfgrass_predictions <- surfgrass_predictions %>%
  mutate(est_p = plogis(est))


ggplot(surfgrass_predictions)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=est_p, width=20,height=20))+
  scale_colour_viridis_c()+
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

ggsave("./figures/surfgrass.png", height = 6, width = 20)

# hold <- surfgrass_predictions %>%
#   group_by(X_m, Y_m) %>%
#   summarise(CV = sd(exp(est), na.rm = TRUE)/mean(exp(est)), counts_ln = mean(est), SE = mean(SE))
# 

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
