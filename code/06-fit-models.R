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
# m_e_3 auc final is 0.932. model indicated by looking at variable relative importance and also considering what is important for future change
# spatial field is causing under dispersion (means model is too complex)!  tried fitting model with spatial field and just depth, substrate, slope and rei but still underdispersion
## AUC (0.944 vs 0.932) and TSS drops when don't have spatial random field.
m_e_3 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + DOmin_stnd + saltmin_stnd + tempcv_stnd + tempmean_stnd + cul_eff_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")
roc <- pROC::roc(m_e_3$data$presence, plogis(m_e_3$data$cv_predicted))
auc <- pROC::auc(roc)
auc

eval_cv <- evalStats( folds=1:numFolds,
                      m=m_e_3,
                      CV=cv_list_eelgrass$cv)

# fit full model
fmodel <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + DOmin_stnd + saltmin_stnd + tempcv_stnd +tempmean_stnd + cul_eff_stnd,
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
ggeffects::ggeffect(model = fmodel,  terms = "cul_eff_stnd[-2:5]") %>% plot() 
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

# Add fitted values (preds) to data
data$fitted_vals <- predict(fmodel, type="response")$est

# Calculate optimal thresholds
thresh <- calcThresh( x=data ) 

#TSS pred thresh add to data
data$pred_TSS_thresh <- ifelse(data$fitted_vals < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1) 

##Kappa pred thresh add to data
data$pred_kappa_thresh <- ifelse(data$fitted_vals < thresh$Predicted[thresh$Method == "MaxKappa"], 0, 1) 

## PCC thresh add to data
data$pred_PCC_thresh <- ifelse(data$fitted_vals < thresh$Predicted[thresh$Method == "MaxPCC"], 0, 1) 

eval_fmod <- evalfmod( x=data, thresh = thresh )
# this model has good TSS, is well calibrated (miller). for calibration (Hosmer & Lemeshow goodness-of-fit) model seems to have issues at mid to higher predicted probabilities

###Notes on evaluation statistics
# Values of TSS greater than 0.6 are considered good, between 0.2 and 0.6 moderate, and less than 0.2 poor (Jones et al. 2010; Landis and Koch 1977).

# miller calibration statistic. If the model is well calibrated, the line should lie along (or at least be nearly parallel to) the reference diagonal, i.e. the slope should ideally equal 1 (i.e., 45 degrees) and the intercept 0
# Miller's calibration statistics are mainly useful when projecting a model outside those training data.
# A slope greater than 1 indicates that predicted values above 0.5 are underestimating, and predicted values below 0.5 are overestimating, the probability of presence. A slope smaller than 1 (while greater than 0) implies that predicted values below 0.5 are underestimating, and values above 0.5 are overestimating, the probability of presence (Pearce & Ferrier 2000). 
#A Miller slope very different from 1 indicates a poorly calibrated model. Baquero et al. (2021) proposed that values between 0.5 and 1.5 can be considered not very different from 1.

# Hosmer & Lemeshow goodness-of-fit. compares predicted probability to observed occurrence frequency at each portion of the probability range
# an important facet of model evaluation is calibration or reliability, i.e., the relationship between predicted probability and observed occurrence frequency (Pearce & Ferrier 2000; Jimenez-Valverde et al. 2013).
# The HLfit function measures model reliability with the Hosmer & Lemeshow goodness-of-fit statistic (Hosmer & Lemeshow 1980).
# and the strong influence of the binning method on the results. try 'HLfit' with different binning methods to see how if the results are robust.
# p-value of the Hosmer-Lemeshow test. Note that this is one of those tests for which higher p-values are better

# The smaller the error measure (eer) returned values, the better the model predictions fit the observations.

#get relative importance
prednames <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "tempcv_stnd", "tempmean_stnd", "DOmin_stnd", "saltmin_stnd", "cul_eff_stnd" )
relimp <- varImp( model=fmodel,
                  dat=data,
                  preds=prednames,
                  permute=10 ) # Number of permutations
# depth 80.3, substrate 14.6, slope 2.0, rei 1.9, do min 1.7, salt min 0.3, tempcv 0.1, tempmean 0.1, cul eff 0.1

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

#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009
data_pre2013 <- data %>% filter(Year < 2010)
mesh_pre2013 <- make_mesh(data = data_pre2013, xy_cols = c("X", "Y"), cutoff = 12) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_pre2013)
barrier_mesh_pre2013 <- add_barrier_mesh(mesh_pre2013, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

m_eelgrass_forecast <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + saltmin_stnd + tempmean_stnd + tempcv_stnd + DOmin_stnd + cul_eff_stnd,
                               mesh = barrier_mesh_pre2013, 
                               family = binomial(link = "logit"), 
                               spatial = FALSE, 
                               data = data_pre2013) 
m_eelgrass_forecast_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + saltmin_stnd + tempmean_stnd + tempcv_stnd + DOmin_stnd + cul_eff_stnd,
                              mesh = barrier_mesh_pre2013, 
                              family = binomial(link = "logit"), 
                              spatial = TRUE, 
                              data = data_pre2013) 
data.df <- data %>% select(presence, X, Y, depth_stnd, slope_stnd, rei_stnd, substrate, saltmin_stnd, DOmin_stnd, tempcv_stnd, tempmean_stnd, cul_eff_stnd, Year)

forecast <- plogis(predict(m_eelgrass_forecast, newdata = data.df %>% filter(Year > 2012))$est)
forecast_spatial <- plogis(predict(m_eelgrass_forecast_spatial, newdata = data.df %>% filter(Year > 2012))$est)

pre_2013 <- plogis(predict(m_eelgrass_forecast, newdata = data.df %>% filter(Year < 2010))$est)
pre_2013_spatial <- plogis(predict(m_eelgrass_forecast_spatial, newdata = data.df %>% filter(Year < 2010))$est)

forecast_predict_eelgrass <- data.frame(TjurR2_no_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast),
                                                    tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013)),
                                        TjurR2_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast_spatial),
                                                              tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013_spatial)),
                                         AUC_no_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast),
                                                 ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013)),
                                        AUC_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast_spatial),
                                                ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013_spatial)),
                                         type = factor(c("forecast", "training"), levels = c("training", "forecast"), ordered =  TRUE))
## AUC dropped from 0.932 to 0.92 and Tjur dropped from 0.147 to 0.141 when forecasted #tjur is quite a bit less when no spatial field
# model did much better forcasting without spatial field
save(data, fmodel, relimp, thresh, r_ret, eval_cv, eval_fmod, forecast_predict_eelgrass, file = "code/output_data/final_eelgrass_model.RData")




#### Surfgrass model ####

sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))

## test VIF
my_model <- lm(presence~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + 
                 saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                 surftempmax_stnd, data = data)
olsrr::ols_vif_tol(my_model)
# Tolerance of <0.1 might indicate multicollinearity. VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. Ideally, the Variance Inflation Factors are below 3.

# ensure PH on all substrates
substrates_present <- data %>% group_by(substrate) %>% summarise(n_present = sum(presence))
substrates_present # there are surgrass presence observations on all substrates

### Forward feature selection test (testing one method to limit variables )
#tested removing depth, makes horrible models
surfgrass_ffs <- glm_ffs(data)  ### variables that came out include depth, substrate, slope_stnd, rei_stnd
### variables that came out include depth, rei sqrt, subtrate. Surftempmin_stnd (2); DOmin (1); Freshwater sqrt (1); Tempcvstn(1); Tempmin (1); Slope  (1); PAR mean (1)
# so need to include some temp variable

#make mesh # tested several mesh sizes between 20- 10 km and 15 had highest AUC
mesh<- make_mesh(data = data, xy_cols = c("X", "Y"), cutoff = 12) 
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)
# model selected by ffs without spatial field AUC = 0.9443
m_s_1 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")
# model selected by ffs with spatial field AUC = 0.9551
m_s_2 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = TRUE, 
                   data = data, 
                   fold_ids = "fold")
# model selected by relimp without spatial field AUC = 0.9613
#having spatial field makes model underdispersed (means too complicated)
# cul_eff not important for surfgrass model so not included (makes auc go down and 0.0 in relimp. 
m_s_3 <- sdmTMB_cv(formula = presence ~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + substrate + 
                     saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                     surftempmax_stnd ,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")

roc <- pROC::roc(m_s_3$data$presence, plogis(m_s_3$data$cv_predicted))
auc <- pROC::auc(roc)
auc

eval_cv <- evalStats( folds=1:numFolds,
                      m=m_s_3,
                      CV=cv_list_seagrass$cv)
# If a model is unbiased bias should be close to zero
# MAE want low values
# AUC want high values above 0.9. According to Pearce and Ferrier (2000) and Jones et al. (2010) values of AUC greater than 0.9 are considered good, between 0.7 and 0.9 moderate, and less than 0.7 poor. values of 0.5 indicate that the model is no better than random.
# TSS balances sensitivity (proportion of presence observations that are correctly classified) and specificity (proportion of absence observations that are correctly classified) and is independent of the prevalence of the observations (Allouche et al. 2006). Values of TSS greater than 0.6 are considered good, between 0.2 and 0.6 moderate, and less than 0.2 poor (Jones et al. 2010; Landis and Koch 1977). 
# Accuracy is the percent of predictions which are correctly classified and varies from values of 0 to 1 where 1 is the highest accuracy.
# Kappa is a measure of agreement between observed and predicted values that accounts for chance agreements and is dependent on prevalence of the observations. Kappa range from -1 to 1 with values less than 0 representing models that are no better than random and values of 1 indicating perfect agreement (Allouche et al. 2006). 

#fit full model
fmodel <- sdmTMB(formula = presence ~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + substrate + 
                     saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                     surftempmax_stnd,
                     mesh = barrier_mesh, 
                     family = binomial(link = "logit"), 
                     spatial = FALSE, 
                     data = data)
# look ar marginal effects plots, note temperature variables are correlated, so these are not correct 
ggeffects::ggeffect(model = fmodel,  terms = "depth_stnd[-3:4]") %>% plot() # found highest at shallow
ggeffects::ggeffect(model = fmodel,  terms = "rei_sqrt_stnd[-1:9]") %>% plot()  # increases with exposure
ggeffects::ggeffect(model = fmodel,  terms = "tidal_sqrt_stnd[-3:8]") %>% plot() # increases with tidal infleunce
ggeffects::ggeffect(model = fmodel,  terms = "substrate") %>% plot() # found on rock
ggeffects::ggeffect(model = fmodel,  terms = "saltcv_stnd[-1:57]") %>% plot() # highest in areas with stable salinity
ggeffects::ggeffect(model = fmodel,  terms = "PARmin_stnd[-3:2]") %>% plot() # increases with min par
ggeffects::ggeffect(model = fmodel,  terms = "DOmean_stnd[-7:3]") %>% plot() #increases with do
ggeffects::ggeffect(model = fmodel,  terms = "surftempcv_stnd[-3:8]") %>% plot() # increases with increased variation
ggeffects::ggeffect(model = fmodel,  terms = "surftempmax_stnd[-4:5]") %>% plot() # declines with increased max temps
ggeffects::ggeffect(model = fmodel,  terms = "surftempmean_stnd[-4:4]") %>% plot() # increases with mean temps


visreg::visreg(fmodel, "depth_stnd")
visreg::visreg(fmodel, "rei_sqrt_stnd")
visreg::visreg(fmodel, "tidal_sqrt_stnd")
visreg::visreg(fmodel, "saltcv_stnd")
visreg::visreg(fmodel, "PARmin_stnd")
visreg::visreg(fmodel, "DOmean_stnd")
visreg::visreg(fmodel, "surftempcv_stnd")
visreg::visreg(fmodel, "surftempmax_stnd")
visreg::visreg(fmodel, "surftempmean_stnd")

# Model check
tidy(fmodel, conf.int = TRUE)
sanity(fmodel)

# Add fitted values (preds) to data
data$fitted_vals <- predict(fmodel, type="response")$est

# Calculate optimal thresholds
thresh <- calcThresh( x=data ) 

#TSS pred thresh add to data
data$pred_TSS_thresh <- ifelse(data$fitted_vals < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1) 

##Kappa pred thresh add to data
data$pred_kappa_thresh <- ifelse(data$fitted_vals < thresh$Predicted[thresh$Method == "MaxKappa"], 0, 1) 

## PCC thresh add to data
data$pred_PCC_thresh <- ifelse(data$fitted_vals < thresh$Predicted[thresh$Method == "MaxPCC"], 0, 1) 

eval_fmod <- evalfmod( x=data, thresh = thresh )
# this model has good TSS, is well calibrated (miller). for calibration (Hosmer & Lemeshow goodness-of-fit) model seems to have issues at higher predicted probabilities

# Variable importance (randomization and permutation method)
prednames <- c("depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd", "substrate",  
               "saltcv_stnd", "PARmin_stnd", "DOmean_stnd", "surftempcv_stnd",  
               "surftempmax_stnd", "surftempmean_stnd")

relimp <- varImp( model=fmodel,
                  dat=data,
                  preds=prednames,
                  permute=10 ) # Number of permutations

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

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

predict(fmodel) %>%
  ggplot(aes(x = presence, y = fmodel$family$linkinv(est)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_jitter(width = 0.05, height = 0)

#make predictions and get SE
hold <- predict(fmodel, env_20m_all)
sims <- predict(fmodel, newdata = env_20m_all, nsim = 100) #sim needs to be 500? ram is not working for this right now at 20m prediction cells
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

#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009
data_pre2013 <- data %>% filter(Year < 2010)
mesh_pre2013 <- make_mesh(data = data_pre2013, xy_cols = c("X", "Y"), cutoff = 15) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_pre2013)
barrier_mesh_pre2013 <- add_barrier_mesh(mesh_pre2013, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

m_surfgrass_forecast <- sdmTMB(formula = presence ~ depth_stnd + rei_sqrt_stnd + tidal_sqrt_stnd + substrate + 
                                 saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd  +
                                 surftempmax_stnd,
                      mesh = barrier_mesh_pre2013, 
                      family = binomial(link = "logit"), 
                      spatial = FALSE, 
                      data = data_pre2013) 
m_surfgrass_forecast_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + saltmin_stnd + tempmean_stnd + tempcv_stnd + DOmin_stnd,
                                      mesh = barrier_mesh_pre2013, 
                                      family = binomial(link = "logit"), 
                                      spatial = TRUE, 
                                      data = data_pre2013) 
data.df <- data %>% select(presence, X, Y, depth_stnd, slope_stnd, rei_stnd, substrate, saltmin_stnd, DOmin_stnd, tempcv_stnd, tempmean_stnd, Year)

forecast <- plogis(predict(m_surfgrass_forecast, newdata = data.df %>% filter(Year > 2012))$est)
forecast_spatial <- plogis(predict(m_surfgrass_forecast_spatial, newdata = data.df %>% filter(Year > 2012))$est)

pre_2013 <- plogis(predict(m_surfgrass_forecast, newdata = data.df %>% filter(Year < 2010))$est)
pre_2013_spatial <- plogis(predict(m_surfgrass_forecast_spatial, newdata = data.df %>% filter(Year < 2010))$est)

forecast_predict_surfgrass <- data.frame(TjurR2_no_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast),
                                                              tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013)),
                                        TjurR2_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast_spatial),
                                                           tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013_spatial)),
                                        AUC_no_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast),
                                                           ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013)),
                                        AUC_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast_spatial),
                                                        ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013_spatial)),
                                        type = factor(c("forecast", "training"), levels = c("training", "forecast"), ordered =  TRUE))
## AUC dropped from 0.96 to 0.95 and Tjur dropped from 0.18 to 0.16
# removing spatial field makes model make  better predictions. So while having spatial random field makes current day predictions better, having it for forecasting makes it worse

save(data, fmodel, relimp, thresh, r_ret, eval_cv, eval_fmod, forecast_predict_surfgrass, file = "code/output_data/final_surfgrass_model.RData")


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

# refer to https://pbs-assess.github.io/sdmTMB/articles/basic-intro.html to make plots of random spatial fields etc


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
