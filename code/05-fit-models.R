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

#load seagrass data
load("code/output_data/seagrass_model_inputs.RData")



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
mesh_eelgrass <- make_mesh(data = eelgrass, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh_eelgrass)
barrier_mesh_eelgrass <- add_barrier_mesh(mesh_eelgrass, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)


# m_eelgrass_1 <- sdmTMB_cv(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|HKey),
#                       mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = TRUE, data = eelgrass, fold_ids = "fold")
# 
# m_eelgrass_1$sum_loglik
eelgrass$ID<-as.factor(eelgrass$ID)

m_eelgrass <- sdmTMB(presence ~ 0 + s(depth_stnd, k = 3) + slope_stnd + temperature_stnd + substrate + s(rei_sqrt_stnd, k=3) + salinity_stnd + s(tidal_sqrt_stnd, k = 3) + (1|ID),
                          mesh = barrier_mesh_eelgrass, family = binomial(link = "logit"), spatial = "on", time = "Year", spatiotemporal = "IID", data = eelgrass)

sanity(m_eelgrass)
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
  geom_jitter(width = 0.05, height = 0)+
  scale_x_log10()+
  scale_y_log10()

hold <- predict(m_eelgrass, prediction_0m.df)
sims <- predict(m_eelgrass, newdata = prediction_0m.df, nsim = 500)
hold$SE <- apply(sims, 1, sd)
eelgrass_predictions <- prediction_0m.df
eelgrass_predictions <- bind_cols(eelgrass_predictions, hold %>% select(est:SE))

ggplot(eelgrass_predictions)+
  geom_sf(data = coastline %>% st_crop(prediction_sf_0m))+
  geom_tile(aes(x = X_m, y = Y_m, fill = exp(est), color = exp(est)), lwd= 0)+
  geom_sf(data = all_data %>% filter(`CHUM SALMON`>0, avg_tow_depth <16) %>% st_crop(prediction_sf_0m), aes(size = `CHUM SALMON`), pch = 1)+
  geom_sf(data = all_data %>% filter(`CHUM SALMON`==0, avg_tow_depth <16) %>% st_crop(prediction_sf_0m), pch = 4)+
  scale_shape_manual(values = c(4,1), guide = "none") +
  scale_size_continuous("catch", breaks = c(1, 100, 1000)) +
  scale_color_viridis_c(trans = "log10", "ind. per\n30 min", guide = "none") +
  scale_fill_viridis_c(trans = "log10", "ind. per\n30 min", breaks = c(1, 10, 100, 1000)) +
  facet_grid(season~year)+
  coord_sf(expand = FALSE)+
  theme_bw()+
  xlab("")+
  ylab("")+
  scale_x_continuous(breaks = seq(-125, -122, by = 1))+
  ggtitle("Chum Salmon")+
  theme(strip.background = element_rect(fill = NA, color = NA))
ggsave("./figures/Chum_counts.png", height = 6, width = 20)

hold <- eelgrass_predictions %>%
  group_by(X_m, Y_m) %>%
  summarise(CV = sd(exp(est), na.rm = TRUE)/mean(exp(est)), counts_ln = mean(est), SE = mean(SE))





surfgrass <- filter(seagrass_data_long, species == "PH")
print(paste("surfgrass present in ", round((sum(surfgrass$presence)/nrow(surfgrass))*100,2), "% of transects", sep = ""))

