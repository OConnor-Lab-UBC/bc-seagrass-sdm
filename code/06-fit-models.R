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

#functions####
tjur <- function(y, pred) {
  categories <- sort(unique(y))
  m1 <- mean(pred[which(y == categories[1])], na.rm = TRUE)
  m2 <- mean(pred[which(y == categories[2])], na.rm = TRUE)
  abs(m2 - m1)
}

ll_binomial <- function(withheld_y, withheld_mu) {
  stats::dbinom(x = withheld_y, size = 1, prob = withheld_mu, log = TRUE)
  }
  
  
#eelgrass model###
eelgrass <- filter(seagrass_data_long, species == "ZO")
print(paste("eelgrass present in ", round((sum(eelgrass$presence)/nrow(eelgrass))*100,2), "% of observations", sep = ""))

library(corrplot)
eelgrass_env<- eelgrass %>%
  select(presence, depth_stnd, rei_sqrt_stnd, tidal_sqrt_stnd, freshwater_sqrt_stnd, slope_sqrt_stnd, NH4_stnd, NO3_stnd, 
         saltmean_sq_stnd, saltmin_sq_stnd, saltcv_stnd, PARmean_stnd, PARmin_stnd, PARmax_stnd, surftempmean_stnd, surftempmin_stnd, 
         surftempmax_stnd, surftempcv_stnd, surftempdiff_stnd, tempmean_stnd, tempmin_stnd, tempmax_stnd, tempcv_stnd, tempdiff_stnd, 
         DOmean_stnd, DOmin_stnd)
corrplot(cor(eelgrass_env), method = "number")
# Correlations close to-1 or +1 might indicate the existence of multicollinearity. 
# As a rule of thumb, one might suspect multicollinearity when the correlation between two (predictor) variables is below -0.9 or above +0.9.
# correlated: PAR max and mean, salt min and mean and cv, surf temp min and surf temp cv, surf temp diff and surf temp max, surftempmean and temp mean, surftempmax and tempmax, temp diff and temp max, tempmean and temp max (but not temp min)
# keep PAR mean and min, salt min, 5 m temps max and min, NH4, NO3

## CREATE A LINEAR REGRESSION MODEL
#my_model <- lm(presence~., data = eelgrass_env)
my_model <- lm(presence~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + freshwater_sqrt_stnd + slope_sqrt_stnd +
                 NH4_stnd + NO3_stnd + saltmin_sq_stnd + saltcv_stnd + PARmean_stnd + surftempcv_stnd  + surftempdiff_stnd +
                 tempmax_stnd + tempcv_stnd + tempdiff_stnd + DOmin_stnd, data = eelgrass_env)
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


### Forward feature selection test eelgrass
k<-n_distinct(eelgrass$fold_eelgrass)
#Set up a list to hold the output for each fold
Fold_Outputs <- list()
for(u in 1:k){
  eelgrass_data <- eelgrass %>% dplyr::filter(fold_eelgrass <=k) %>% dplyr::filter(fold_eelgrass != u)
  Test_Data <- eelgrass %>% dplyr::filter(fold_eelgrass <=k) %>% dplyr::filter(fold_eelgrass == u)
  #create a vector for fold membership
  foldmem <- seq(1:1:length(eelgrass_data$fold_eelgrass))
  folds <- unique(eelgrass_data$fold_eelgrass)
  new_folds <- seq(1:1:(k-1))
  for (i in 1:length(eelgrass_data$fold_eelgrass)){
    for (j in 1:(k-1)){
      if (eelgrass_data$fold_eelgrass[[i]] == folds[j]){
        foldmem[i] <- new_folds[j]
      }
    }
  }
  #Setting up an index list for the folds in the caret model training
  index_list <- list()
  for (i in 1:(k-1)){
    index_list[[i]] <- which(foldmem == i)
  }
    #Setting up parameters for how my model is going to be fitted
  fitControl <- caret::trainControl(method = "cv",
                                    number = (k-1),
                                    index = index_list)
  #setting a seed in case that matters
  set.seed(2024)
  #performing model selection by glmStepAIC I think I may want to use glmboost instead
  caret_model <- CAST::ffs(response = eelgrass_data$presence, 
                           predictors = eelgrass_data[,6:37], 
                           method = "glm", 
                           family = "binomial",
                           trControl = fitControl)
    #Create the final glm model using the above determined formula
  #Create a final formula using the selected variables
  selectedVars <- caret_model$selectedvars
  final_formula <- paste0("presence~",selectedVars[1])
  for (i in 2: length(selectedVars)){
    final_formula <- paste0(final_formula,"+",selectedVars[i])
  }
  #fit the model
  caretmodel <- glm(as.formula(final_formula), data = eelgrass_data, family = binomial(link = "logit"))
  #Calculate final model AUC on the testing fold
  pred.caretModel <- predict(caretmodel, newdata = Test_Data, type = "response")
  roc.caretModel <- pROC::roc(Test_Data$presence, pred.caretModel)
  auc.caretModel <- pROC::auc(roc.caretModel)
  Output <- list(caretmodel, auc.caretModel, caretmodel$formula)
  Fold_Outputs[[u]] <- Output
}
return(Fold_Outputs)
### variables that came out include depth, substrate, slope_stnd,  rei_stnd (not rei sqrt so need to assess)


#SDMtmb model
#make mesh
mesh_eelgrass <- make_mesh(data = eelgrass, xy_cols = c("X", "Y"), cutoff = 12) # going to 10 km makes the model not run, likely will need to reduce fixed effects to make it work
plot(mesh_eelgrass)

barrier_mesh_eelgrass <- add_barrier_mesh(mesh_eelgrass, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit cv model
plan(multisession)
m_e_1 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd,
                          mesh = barrier_mesh_eelgrass, 
                          family = binomial(link = "logit"), 
                          spatial = FALSE, 
                          data = eelgrass, 
                          fold_ids = "fold_eelgrass")
m_e_2 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd,
                   mesh = barrier_mesh_eelgrass, 
                   family = binomial(link = "logit"), 
                   spatial = TRUE, 
                   data = eelgrass, 
                   fold_ids = "fold_eelgrass")
roc <- pROC::roc(m_e_2$data$presence, plogis(m_e_2$data$cv_predicted))
auc <- pROC::auc(roc)
auc 
# auc with no spatial field is 0.9326, with spatial field is 0.9427
# with depth and substrate AUC = 0.9307, adding slope made 0.9337, adding rei made 0.9361
#STILL NEED TO DECIDE BETWEEN REI SQRT AND REI, WITH SQRT AUC IS 0.9431, WITH OUT IS 0.9427, tJUR IS BETTER WITHOUT,  FFC ONLY SELECTED REI WITHOUT SQRT
# addition of depth vastly increases auc from 0.6223 to 0.863, model does better with spline
#slope_stnd produces slightly better auc than with it sqrt so will keep the stnd
# adding Par max made marginal improvement, but PAR mean and PAR min made it go down
# adding salinitymin made model decline, adding salinity mean made it stay same
# adding tidal did not improve the model auc is still 0.9361, don't add spline, relationship is not quadratic

folds <- unique(eelgrass$fold_eelgrass)
CV=cv_list_eelgrass$cv
fitted.df <- data.frame()
for(i in folds){
  mu <-  plogis(predict(m_e_2$models[[i]])$est[eelgrass$fold_eelgrass == i])
  fitted.df <- rbind(fitted.df, data.frame(presence = eelgrass$presence[eelgrass$fold_eelgrass == i], mu = mu, fold = i))    
}
#compare test and train statistics
traintest.df <- data.frame()
for(i in folds){
  # Get train obs
  train <- CV[[i]][["train"]]
  sp_data_cv <- filter(seagrass_data_long, species == "ZO")
  trainobs <- sp_data_cv[ train, 44]
  #trainobs <- obs$presence
  # Get train preds
  trainpred <- plogis(predict(m_e_2$models[[i]])$est[train])
  # Calculate area under the receiver-operator curve (AUC))
  train.AUC <-  pROC::auc(trainobs, trainpred)
  # Calculate tjur R2
  train.tjur <- tjur(trainobs, trainpred )
  # Get test indices
  test <- CV[[i]][["test"]]
  testobs <- sp_data_cv[ test, 44]
  #testobs <- tobs$presence
  # Get test preds
  testpred <- plogis(predict(m_e_2$models[[i]])$est[test])
  # Calculate area under the receiver-operator curve (AUC))
  test.AUC <- pROC::auc(testobs, testpred)
  # Calculate tjur R2
  test.tjur <- tjur(testobs, testpred)
  # sum log likelihood
  ll<- m_e_2$sum_loglik
  traintest.df <- data.frame(ll=ll, train.AUC = train.AUC, test.AUC=test.AUC, train.tjur = train.tjur, test.tjur = test.tjur, species = "ZO", fold = i)
  traintest.df <- traintest.df %>% dplyr::summarise(mean_AUC_train = mean(train.AUC, na.rm = TRUE), mean_AUC_test = mean(test.AUC, na.rm = TRUE), mean_Tjur_train = mean(train.tjur, na.rm = TRUE), mean_Tjur_test = mean(test.tjur, na.rm = TRUE), sum_loglike = mean(ll, na.rm = TRUE))
}






# fit full model

m_eelgrass <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd,
                     mesh = barrier_mesh_eelgrass, 
                     family = binomial(link = "logit"), 
                     spatial = "on", 
                     data = eelgrass)
ggeffects::ggeffect(model = m_eelgrass,  terms = "depth_stnd[-4:0]") %>% plot() 
ggeffects::ggeffect(model = m_eelgrass,  terms = "rei_stnd[-1:8]") %>% plot() 
ggeffects::ggeffect(model = m_eelgrass,  terms = "slope_stnd[-1:8]") %>% plot() 
ggeffects::ggeffect(model = m_eelgrass,  terms = "substrate") %>% plot() 
visreg::visreg(m_eelgrass, "depth_stnd")

tidy(m_eelgrass, conf.int = TRUE)
sanity(m_eelgrass)

eelgrass$resids <- residuals(m_eelgrass) # randomized quantile residuals
qqnorm(eelgrass$resids)
qqline(eelgrass$resids)

#sdmTMB::run_extra_optimization(m_eelgrass, nlminb_loops = 1L, newton_loops = 1L)

m_eelgrass
# refere to https://pbs-assess.github.io/sdmTMB/articles/basic-intro.html to make plots of random spatial fields etc

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
sims <- predict(m_eelgrass, newdata = env_20m_all, nsim = 100) #ram is not working for this right now for 500 sims (or 250) at 20m prediction cells
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


#### Surfgrass model ####

surfgrass <- filter(seagrass_data_long, species == "PH")
print(paste("surfgrass present in ", round((sum(surfgrass$presence)/nrow(surfgrass))*100,2), "% of 20 m cells", sep = ""))

my_model_surfgrass <- lm(presence~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + substrate+ freshwater_sqrt_stnd + slope_sqrt_stnd +
                 NH4_stnd + NO3_stnd + saltmin_stnd + PARmean_stnd + surftempcv_stnd +
                   tempcv_stnd + DOmin_stnd +tempmin_stnd, data = surfgrass)
ols_vif_tol(my_model_surfgrass)
# As a general guideline, a Tolerance of <0.1 might indicate multicollinearity.
# As a rule of thumb, a VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. 
# Ideally, the Variance Inflation Factors are below 3.
# there is no multi multicollinearity in variables selevered in my_model.

### Forward feature selection test eelgrass
k<-n_distinct(surfgrass$fold_seagrass)
#Set up a list to hold the output for each fold
Fold_Outputs_seagrass <- list()
for(u in 1:k){
  surfgrass_data <- surfgrass %>% dplyr::filter(fold_seagrass <=k) %>% dplyr::filter(fold_seagrass != u)
  Test_Data <- surfgrass %>% dplyr::filter(fold_seagrass <=k) %>% dplyr::filter(fold_seagrass == u)
  #create a vector for fold membership
  foldmem <- seq(1:1:length(surfgrass_data$fold_seagrass))
  folds <- unique(surfgrass_data$fold_seagrass)
  new_folds <- seq(1:1:(k-1))
  for (i in 1:length(surfgrass_data$fold_seagrass)){
    for (j in 1:(k-1)){
      if (surfgrass_data$fold_seagrass[[i]] == folds[j]){
        foldmem[i] <- new_folds[j]
      }
    }
  }
  #Setting up an index list for the folds in the caret model training
  index_list <- list()
  for (i in 1:(k-1)){
    index_list[[i]] <- which(foldmem == i)
  }
  #Setting up parameters for how my model is going to be fitted
  fitControl <- caret::trainControl(method = "cv",
                                    number = (k-1),
                                    index = index_list)
  #setting a seed in case that matters
  set.seed(2024)
  #performing model selection by glmStepAIC I think I may want to use glmboost instead
  caret_model <- CAST::ffs(response = surfgrass_data$presence, 
                           predictors = surfgrass_data[,6:37], 
                           method = "glm", 
                           family = "binomial",
                           trControl = fitControl)
  #Create the final glm model using the above determined formula
  #Create a final formula using the selected variables
  selectedVars <- caret_model$selectedvars
  final_formula <- paste0("presence~",selectedVars[1])
  for (i in 2: length(selectedVars)){
    final_formula <- paste0(final_formula,"+",selectedVars[i])
  }
  #fit the model
  caretmodel <- glm(as.formula(final_formula), data = surfgrass_data, family = binomial(link = "logit"))
  #Calculate final model AUC on the testing fold
  pred.caretModel <- predict(caretmodel, newdata = Test_Data, type = "response")
  roc.caretModel <- pROC::roc(Test_Data$presence, pred.caretModel)
  auc.caretModel <- pROC::auc(roc.caretModel)
  Output <- list(caretmodel, auc.caretModel, caretmodel$formula)
  Fold_Outputs_seagrass[[u]] <- Output
}
return(Fold_Outputs_seagrass)
### variables that came out include depth, rei sqrt, subtrate. Surftempmin_stnd (2); DOmin (1); Freshwater sqrt (1); Tempcvstn(1); Tempmin (1); Slope  (1); PAR mean (1)
# so need to include some temp variable

# present on all
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
m_s_1 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd,
                   mesh = barrier_mesh_surfgrass, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = surfgrass, 
                   fold_ids = "fold_seagrass")
m_s_2 <- sdmTMB_cv(formula = presence ~ depth_stnd  + substrate + rei_sqrt_stnd + surftempcv_stnd + tempcv_stnd,
                   mesh = barrier_mesh_surfgrass, 
                   family = binomial(link = "logit"), 
                   spatial = TRUE, 
                   data = surfgrass, 
                   fold_ids = "fold_seagrass")
roc <- pROC::roc(m_s_2$data$presence, plogis(m_s_2$data$cv_predicted))
auc <- pROC::auc(roc)
auc 

# without spatial field AUC is 0.9443, with spatial field is 0.9551
#need depth, substrate, rei sqrt. Test addin these: Surftempmin_stnd (2); DOmin (1); Tempmin (1); 
#don't need spline for depth or rei
# added suftempcv_stnd increased to 0.96
# added tempcv_stnd increased to 0.9607
#adding slope, PAR, DO, freshwater, tempmin_stnd, surftempmin doesn't help or hinders after other factors accounted for
# tested several mess sizes between 20- 10 km and 15 had highest AUC
# cv is calculated as the sd for a year divided by the mean of the year *100. Then the mean across the decade is calculated from each year

#fit full model
m_surfgrass <- sdmTMB(formula = presence ~ depth_stnd  + substrate + rei_sqrt_stnd + surftempcv_stnd + tempcv_stnd,
                     mesh = barrier_mesh_surfgrass, 
                     family = binomial(link = "logit"), 
                     spatial = "on", 
                     data = surfgrass)
ggeffects::ggeffect(model = m_surfgrass,  terms = "depth_stnd[-4:0]") %>% plot() # found highest at shallow
ggeffects::ggeffect(model = m_surfgrass,  terms = "rei_sqrt_stnd[-1:8]") %>% plot()  # increases with exposure
ggeffects::ggeffect(model = m_surfgrass,  terms = "substrate") %>% plot() # found on rock
ggeffects::ggeffect(model = m_surfgrass,  terms = "surftempcv_stnd[-3:8]") %>% plot() # declines with increased variation
ggeffects::ggeffect(model = m_surfgrass,  terms = "tempcv_stnd[-3:6]") %>% plot() # declines with increased variation
visreg::visreg(m_surfgrass, "depth_stnd")
visreg::visreg(m_surfgrass, "surftempcv_stnd")

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
