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

#set up cv model using spatial cross validation###
eelgrass <- filter(seagrass_data_long, species == "ZO")
print(paste("eelgrass present in ", round((sum(eelgrass$presence)/nrow(eelgrass))*100,2), "% of transects", sep = ""))

library(corrplot)
eelgrass_env<- eelgrass %>%
  select(presence, depth_stnd, slope_stnd, temperature_stnd, salinity_stnd, rei_sqrt_stnd, tidal_sqrt_stnd)
corrplot(cor(eelgrass_env), method = "number")
# Correlations close to-1 or +1 might indicate the existence of multicollinearity. 
# As a rule of thumb, one might suspect multicollinearity when the correlation between two (predictor) variables is below -0.9 or above +0.9.

## CREATE A LINEAR REGRESSION MODEL
my_model <- lm(presence~., data = eelgrass_env)
library(olsrr)
ols_vif_tol(my_model)
# As a general guideline, a Tolerance of <0.1 might indicate multicollinearity.
# As a rule of thumb, a VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. 
# Ideally, the Variance Inflation Factors are below 3.
# there is no multi multicollinearity in this dataset.

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


# m_eelgrass_1 <- sdmTMB_cv(presence ~1 + s(depth_stnd, k = 3) + slope_stnd + s(temperature_stnd, k = 3) + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3),
#                       mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = TRUE, data = eelgrass, fold_ids = "fold")
#  
# m_eelgrass_1$sum_loglik
#eelgrass$ID<-as.factor(eelgrass$ID)

# m_eelgrass <- sdmTMB(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|ID),
#                           mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = "on", time = "Year", spatiotemporal = "IID", data = eelgrass)
# 
# m_eelgrass <- sdmTMB(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|HKey),
#                      mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = "on", data = eelgrass)

m_eelgrass <- sdmTMB(presence ~ 1 + s(depth_stnd, k = 3) + slope_stnd + s(temperature_stnd, k = 3) + substrate + s(rei_sqrt_stnd, k=3) + s(salinity_stnd, k = 3) + s(tidal_sqrt_stnd, k = 3),
                     mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = "on", data = eelgrass)

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
sims <- predict(m_eelgrass, newdata = env_20m_all, nsim = 500)
hold$SE <- apply(sims, 1, sd)
eelgrass_predictions <- env_20m_all
eelgrass_predictions <- bind_cols(eelgrass_predictions, hold %>% select(est:SE))

eelgrass_predictions <- eelgrass_predictions %>%
  mutate(est_p = plogis(est))
         
         
ggplot(eelgrass_predictions)+
  geom_sf(data = coastline, linewidth = 0.1)+
  geom_tile(aes(x = X_m, y = Y_m, colour=est_p, width=20,height=20))+
  scale_colour_gradient(
    low = "#fee5d9",
    high = "#cb181d")+
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

#save as raster
eelgrass_raster <- eelgrass_predictions %>%
  mutate(x_round = round(X_m, -2), y_round = round(Y_m, -2)) %>%
  select(x_round, y_round, est_p)

eelgrass_raster <- rast(x = eelgrass_raster %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(eelgrass_raster, file.path("./raster/eelgrass_predictions.tif"), overwrite=TRUE)




surfgrass <- filter(seagrass_data_long, species == "PH")
print(paste("surfgrass present in ", round((sum(surfgrass$presence)/nrow(surfgrass))*100,2), "% of transects", sep = ""))

