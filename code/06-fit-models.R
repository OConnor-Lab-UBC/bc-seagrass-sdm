###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# fit sdm models with sdmTMB 
#
###############################################################################

#### Load modelling functions ####
source("code/modelling-functions.R")

####load packages####
UsePackages(c("sdmTMB", "sdmTMBextra", "tidyverse", "sf", "future", "terra", "future.apply"))

####load model input and prediction data####
load("code/output_data/seagrass_model_inputs.RData")
load("code/output_data/prediction_model_inputs.RData")
seagrass_data_long <- seagrass_data_long %>% select(-saltmean_sq_stnd, -slope_sqrt_stnd, -saltmin_sq_stnd)

# No categorical predictors in environmental layers
facVars <- "substrate"

####Eelgrass model####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))

#test variable correlation
data_env<- data %>%  select(6:34)
enames <- names(data_env)
corrplot::corrplot(cor(data_env), method = "number")
# Correlations close to-1 or +1 might indicate the existence of multicollinearity. one might suspect multicollinearity when the correlation between two (predictor) variables is below -0.9 or above +0.9.
# correlated: PAR max and mean, salt min and mean and cv, surf temp min and surf temp cv, surf temp diff and surf temp max, surftempmean and temp mean, surftempmax and tempmax, temp diff and temp max, tempmean and temp max (but not temp min)

## test VIF
my_model <- lm(presence~ depth_stnd + substrate + slope_stnd + rei_stnd + DOmin_stnd + saltmin_stnd + tempcv_stnd +tempmean_stnd, data = data)
olsrr::ols_vif_tol(my_model)
# Tolerance of <0.1 might indicate multicollinearity. VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. Ideally, the Variance Inflation Factors are below 3.
VIFs <- CalcVIFs( dat=data_env[enames], VIFThresh=10 )
# ones recomended to move  "tempmax_stnd"  "surftempmax_stnd" "surftempmin_stnd" "saltmean_stnd" "PARmean_stnd" "tempmean_stnd" "tempcv_stnd" "tidal_sqrt_stnd" "DOmin_stnd" "saltmin_stnd"  

# ensure ZO on all substrates
substrates_present <- data %>% group_by(substrate) %>% summarise(n_present = sum(presence))
substrates_present # there are eelgrass presence observations on all substrates

### Forward feature selection test (testing one method to limit variables )
#tested removing depth, makes horrible models
eelgrass_ffs <- glm_ffs(data)  ### variables that came out include depth, substrate, slope_stnd, rei_stnd

####SDMtmb model####
#make mesh
mesh<- make_mesh(data = data, xy_cols = c("X", "Y"), cutoff = 12) # going to 10 km makes the model not run, likely will need to reduce fixed effects to make it work
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit cv model of spatial blocking 
plan(multisession)
#model indicated by forward feature selection with no spatial field, AUC is 0.93
m_e_1 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd,
                          mesh = barrier_mesh, 
                          family = binomial(link = "logit"), 
                          spatial = FALSE, 
                          data = data, 
                          fold_ids = "fold")
#model indicated by forward feature selection with spatial field, AUC is 0.9427
m_e_2 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")
# m_e_3 auc final is 0.9435. model indicated by looking at variable relative importance and also considering what is important for future change
# spatial field is causing under dispersion (means model is too complex)!  tried fitting model with spatial field and just depth, substrate, slope and rei but still underdispersion
## AUC (0.944 vs 0.932) and TSS drops when don't have spatial random field.
m_e_3 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + DOmin_stnd + saltmin_stnd + tempcv_stnd + tempmean_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")
roc <- pROC::roc(m_e_3$data$presence, plogis(m_e_3$data$cv_predicted))
auc <- pROC::auc(roc)
auc

evaldat <- evalStats( folds=1:numFolds,
                      m=m_e_3,
                      CV=cv_list_eelgrass$cv)


# fit full model
fmodel <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + DOmin_stnd + saltmin_stnd + tempcv_stnd +tempmean_stnd,
                     mesh = barrier_mesh, 
                     family = binomial(link = "logit"), 
                     spatial = FALSE, 
                     data = data)

#have a look at marginal effects
ggeffects::ggeffect(model = fmodel,  terms = "rei_stnd[-1:8]") %>% plot() 
ggeffects::ggeffect(model = fmodel,  terms = "slope_stnd[-1:8]") %>% plot()
ggeffects::ggeffect(model = fmodel,  terms = "saltmin_stnd[-11:2]") %>% plot() # as salinity min increases so does presence
ggeffects::ggeffect(model = fmodel,  terms = "tempmean_stnd[-4:5]") %>% plot() # as mean temp increases presence goes down
ggeffects::ggeffect(model = fmodel,  terms = "tempcv_stnd[-3:6]") %>% plot() # as temp varability increases so does presence
ggeffects::ggeffect(model = fmodel,  terms = "DOmin_stnd[-7:2]") %>% plot() # as min DO increases presence goes down
ggeffects::ggeffect(model = fmodel,  terms = "substrate") %>% plot() 
visreg::visreg(fmodel, "depth_stnd")
visreg::visreg(fmodel, "DOmin_stnd")
visreg::visreg(fmodel, "slope_stnd")
visreg::visreg(fmodel, "rei_stnd")
visreg::visreg(fmodel, "saltmin_stnd")
visreg::visreg(fmodel, "tempmean_stnd")
visreg::visreg(fmodel, "tempcv_stnd")

# Model check
tidy(fmodel, conf.int = TRUE)
sanity(fmodel)

# Add fitted values to data
data$fitted_vals <- predict(fmodel, type="response")$est

# Calculate optimal thresholds
thresh <- calcThresh( x=data ) # TSS threshold is 0.03

pred <- predict(fmodel)
pred$p <- plogis(pred$est)
pred$pred_01 <- ifelse(pred$p < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1) # use TSS threshold
conmat <- table(pred$pred_01, pred$presence)
true_neg <- conmat[1, 1]
false_neg <- conmat[1, 2]
false_pos <- conmat[2, 1]
true_pos <- conmat[2, 2]

# Calculate TSS:
true_pos_rate <- true_pos / (true_pos + false_neg)
true_neg_rate <- true_neg / (true_neg + false_pos)
TSS <- true_pos_rate + true_neg_rate - 1
TSS #0.719
# Values of TSS greater than 0.6 are considered good, between 0.2 and 0.6 moderate, and less than 0.2 poor (Jones et al. 2010; Landis and Koch 1977).

#get relative importance
prednames <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "tempcv_stnd", "tempmean_stnd", "DOmin_stnd", "saltmin_stnd" )
relimp <- varImp( model=fmodel,
                  dat=data,
                  preds=prednames,
                  permute=10 ) # Number of permutations
# depth 80.3, substrate 14.6, slope 2.0, rei 1.9, do min 1.7, salt min 0.3, tempcv 0.1, tempmean 0.1

####check residuals####
# MCMC based randomized quantile residuals (takes a while to compute)
# set.seed(123)
# samps <- sdmTMBextra::predict_mle_mcmc(fmodel, mcmc_iter = 800, mcmc_warmup = 400)
# mcmc_res <- residuals(fmodel, type = "mle-mcmc", mcmc_samples = samps)
# qqnorm(mcmc_res)
# abline(0, 1)


#analytical randomized quantile approach
data$resids <- residuals(fmodel, type = "mle-mvn") # randomized quantile residuals
# check
ggplot(data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + theme_bw()
hist(data$resids)
qqnorm(data$resids);abline(a = 0, b = 1)

# simulation-based randomized quantile residuals, no under dispersion when spatial random field removed
set.seed(123)
ret<- simulate(fmodel, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

predict(fmodel) %>%
  ggplot(aes(x = presence, y = fmodel$family$linkinv(est)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_jitter(width = 0.05, height = 0)

# make predictions and get standard error
hold <- predict(fmodel , env_20m_all)
sims <- predict(fmodel , newdata = env_20m_all, nsim = 100) #ram is not working for this right now for >100 sims at 20m prediction cells
hold$SE <- apply(sims, 1, sd)
hold$SD <- apply(pclog(sims), 1, sd)
eelgrass_predictions <- env_20m_all
eelgrass_predictions <- bind_cols(eelgrass_predictions, hold %>% select(est:SD))

# hold <- eelgrass_predictions %>%
#   group_by(X_m, Y_m) %>%
#   summarise(CV = sd(exp(est), na.rm = TRUE)/mean(exp(est)), counts_ln = mean(est), SE = mean(SE))

# change to 0-1 away from log-odds (logit) space
eelgrass_predictions <- eelgrass_predictions %>%
  mutate(est_p = plogis(est))
 
save(eelgrass_predictions, file = "code/output_data/eelgrass_predictions.RData")
save(data, fmodel, relimp, thresh, r_ret, evaldat, TSS, file = "code/output_data/final_eelgrass_model.RData")


#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009
data_pre2013 <- data %>% filter(Year < 2010)
mesh_pre2013 <- make_mesh(data = data_pre2013, xy_cols = c("X", "Y"), cutoff = 12) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_pre2013)
barrier_mesh_pre2013 <- add_barrier_mesh(mesh_pre2013, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

m_eelgrass_forecast <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + saltmin_stnd + tempmean_stnd + tempcv_stnd + DOmin_stnd,
                               mesh = barrier_mesh_pre2013, 
                               family = binomial(link = "logit"), 
                               spatial = FALSE, 
                               data = data_pre2013) 
data.df <- data %>% select(presence, X, Y, depth_stnd, slope_stnd, rei_stnd, substrate, saltmin_stnd, DOmin_stnd, tempcv_stnd, tempmean_stnd, Year)

forecast <- plogis(predict(m_eelgrass_forecast, newdata = data.df %>% filter(Year > 2012))$est)
pre_2013 <- plogis(predict(m_eelgrass_forecast, newdata = data.df %>% filter(Year < 2010))$est)

forecast_predict_eelgrass <- data.frame(TjurR2 = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast),
                                                    tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013)),
                                         AUC = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast),
                                                 ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013)),
                                         type = factor(c("forecast", "training"), levels = c("training", "forecast"), ordered =  TRUE))
## AUC dropped from 0.932 to 0.92 and Tjur dropped from 0.147 to 0.141 #tjur is quite a bit less when no spatial field





#### Surfgrass model ####
surfgrass <- filter(seagrass_data_long, species == "PH")
print(paste("surfgrass present in ", round((sum(surfgrass$presence)/nrow(surfgrass))*100,2), "% of 20 m cells", sep = ""))

# test correlation
surfgrass_env<- surfgrass %>%
  select(presence, depth_stnd, rei_sqrt_stnd, slope_stnd, saltcv_stnd, tidal_sqrt_stnd,
         surftempmean_stnd, surftempmax_stnd, PARmin_stnd,
         surftempcv_stnd, DOmean_stnd)
corrplot::corrplot(cor(surfgrass_env), method = "number")

my_model_surfgrass <- lm(presence~  depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd  + 
                           saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                           surftempmax_stnd, data = surfgrass_env)
olsrr::ols_vif_tol(my_model_surfgrass)
# As a general guideline, a Tolerance of <0.1 might indicate multicollinearity.
# As a rule of thumb, a VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. 
# Ideally, the Variance Inflation Factors are below 3.


### Forward feature selection test surfgrass
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
mesh_surfgrass <- make_mesh(data = surfgrass, xy_cols = c("X", "Y"), cutoff = 15) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_surfgrass)
barrier_mesh_surfgrass <- add_barrier_mesh(mesh_surfgrass, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)
# model selected by ffs without spatial field AUC = 0.9443
m_s_1 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd,
                   mesh = barrier_mesh_surfgrass, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = surfgrass, 
                   fold_ids = "fold_seagrass")
# model selected by ffs with spatial field AUC = 0.9551
m_s_2 <- sdmTMB_cv(formula = presence ~depth_stnd + substrate + rei_sqrt_stnd,
                   mesh = barrier_mesh_surfgrass, 
                   family = binomial(link = "logit"), 
                   spatial = TRUE, 
                   data = surfgrass, 
                   fold_ids = "fold_seagrass")
# model selected by relimp with spatial field AUC = 0.9659
m_s_3 <- sdmTMB_cv(formula = presence ~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + substrate + 
                     saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                     surftempmax_stnd,
                   mesh = barrier_mesh_surfgrass, 
                   family = binomial(link = "logit"), 
                   spatial = TRUE, 
                   data = surfgrass, 
                   fold_ids = "fold_seagrass")
roc <- pROC::roc(m_s_3$data$presence, plogis(m_s_3$data$cv_predicted))
auc <- pROC::auc(roc)
auc 

#fit full model
m_surfgrass <- sdmTMB(formula = presence ~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + substrate + 
                        saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                        surftempmax_stnd,
                     mesh = barrier_mesh_surfgrass, 
                     family = binomial(link = "logit"), 
                     spatial = "on", 
                     data = surfgrass)
# look ar marginal effects plots, note temperature variables are correlated, so these are not correct 
ggeffects::ggeffect(model = m_surfgrass,  terms = "depth_stnd[-3:4]") %>% plot() # found highest at shallow
ggeffects::ggeffect(model = m_surfgrass,  terms = "rei_sqrt_stnd[-1:9]") %>% plot()  # increases with exposure
ggeffects::ggeffect(model = m_surfgrass,  terms = "tidal_sqrt_stnd[-3:8]") %>% plot() 
ggeffects::ggeffect(model = m_surfgrass,  terms = "substrate") %>% plot() # found on rock
ggeffects::ggeffect(model = m_surfgrass,  terms = "saltcv_stnd[-1:57]") %>% plot()
ggeffects::ggeffect(model = m_surfgrass,  terms = "PARmin_stnd[-3:2]") %>% plot()
ggeffects::ggeffect(model = m_surfgrass,  terms = "DOmean_stnd[-7:3]") %>% plot()
ggeffects::ggeffect(model = m_surfgrass,  terms = "surftempcv_stnd[-3:8]") %>% plot() # declines with increased variation
ggeffects::ggeffect(model = m_surfgrass,  terms = "surftempmax_stnd[-4:5]") %>% plot() # declines with increased variation
ggeffects::ggeffect(model = m_surfgrass,  terms = "surftempmean_stnd[-4:4]") %>% plot()

visreg::visreg(m_surfgrass, "depth_stnd")
visreg::visreg(m_surfgrass, "rei_sqrt_stnd")
visreg::visreg(m_surfgrass, "tidal_sqrt_stnd")
visreg::visreg(m_surfgrass, "saltcv_stnd")
visreg::visreg(m_surfgrass, "PARmin_stnd")
visreg::visreg(m_surfgrass, "DOmean_stnd")
visreg::visreg(m_surfgrass, "surftempcv_stnd")
visreg::visreg(m_surfgrass, "surftempmax_stnd")
visreg::visreg(m_surfgrass, "surftempmean_stnd")

prednames <- c("depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd", "substrate",  
               "saltcv_stnd", "PARmin_stnd", "DOmean_stnd", "surftempcv_stnd",  
               "surftempmax_stnd", "surftempmean_stnd")

# Variable importance (randomization and permutation method)
relimp <- varImp( model=m_surfgrass,
                  dat=surfgrass,
                  preds=prednames,
                  permute=10 ) # Number of permutations

# test model convergence
sanity(m_surfgrass)

set.seed(123)
set<- simulate(m_surfgrass, nsim = 500, type = "mle-mvn") |>
  dharma_residuals(m_surfgrass, return_DHARMa = TRUE)
plot(set)
set
ggsave("./figures/surfgrass_DHARMa_res.png", height = 6, width = 6)

DHARMa::testDispersion(set)

predict(m_surfgrass) %>%
  ggplot(aes(x = presence, y = m_surfgrass$family$linkinv(est)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_jitter(width = 0.05, height = 0)

#NOT ABLE TO MAKE WORK ON COMPUTER, NEED HIGHER COMPUTING
#make predictions and get SE
hold <- predict(m_surfgrass, env_20m_all)
sims <- predict(m_surfgrass, newdata = env_20m_all, nsim = 100) #ram is not working for this right now at 20m prediction cells
hold$SE <- apply(sims, 1, sd)
surfgrass_predictions <- env_20m_all
surfgrass_predictions <- bind_cols(surfgrass_predictions, hold %>% select(est:SE))

# hold <- surfgrass_predictions %>%
#   group_by(X_m, Y_m) %>%
#   summarise(CV = sd(exp(est), na.rm = TRUE)/mean(exp(est)), counts_ln = mean(est), SE = mean(SE))

surfgrass_predictions <- surfgrass_predictions %>%
  mutate(est_p = plogis(est))

#save outputs####
save(surfgrass_predictions, file = "code/output_data/surfgrass_predictions.RData")

#### test forecasting
# left a few years gap 2010-2012
#trained model with 1993-2009
surfgrass_pre2013 <- surfgrass %>% filter(Year < 2010)
mesh_surfgrass_pre2013 <- make_mesh(data = surfgrass_pre2013, xy_cols = c("X", "Y"), cutoff = 15) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_surfgrass_pre2013)
barrier_mesh_surfgrass_pre2013 <- add_barrier_mesh(mesh_surfgrass_pre2013, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

m_surfgrass_forecast <- sdmTMB(formula = presence ~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + substrate + 
                        saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                        surftempmax_stnd,
                      mesh = barrier_mesh_surfgrass_pre2013, 
                      family = binomial(link = "logit"), 
                      spatial = "on", 
                      data = surfgrass_pre2013) 
surfgrass.df <- surfgrass %>% select(presence, X, Y, depth_stnd, rei_sqrt_stnd, tidal_sqrt_stnd, substrate, 
                                                              saltcv_stnd, PARmin_stnd, DOmean_stnd, surftempcv_stnd, surftempmean_stnd,
                                                              surftempmax_stnd, Year)

forecast <- plogis(predict(m_surfgrass_forecast, newdata = surfgrass.df %>% filter(Year > 2012))$est)
pre_2013 <- plogis(predict(m_surfgrass_forecast, newdata = surfgrass.df %>% filter(Year < 2010))$est)

forecast_predict_surfgrass <- data.frame(TjurR2 = c(tjur(y = surfgrass.df$presence[surfgrass.df$Year > 2012], pred = forecast),
                                             tjur(y = surfgrass.df$presence[surfgrass.df$Year < 2010], pred = pre_2013)),
                                  AUC = c(ModelMetrics::auc(surfgrass.df$presence[surfgrass.df$Year > 2012], forecast),
                                          ModelMetrics::auc(surfgrass.df$presence[surfgrass.df$Year < 2010], pre_2013)),
                                  type = factor(c("forecast", "training"), levels = c("training", "forecast"), ordered =  TRUE))
## AUC dropped from 0.983 to 0.89 and Tjur dropped from 0.24 to 0.07







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

# refere to https://pbs-assess.github.io/sdmTMB/articles/basic-intro.html to make plots of random spatial fields etc


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
