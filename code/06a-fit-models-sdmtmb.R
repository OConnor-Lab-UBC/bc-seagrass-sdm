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
# test spatiotemporal models but they are always underdispersed no matter the mesh size or number of variables
#realised that becasue we don't have repeat sampling and we are not asking quesitons about temporal variability in eelgrass 
#it doesn't make sense to add a spatiotemporal field.

###############################################################################
#### Load modelling functions ####
source("code/modelling-functions.R")

####load packages####
UsePackages(c("sdmTMB", "sdmTMBextra", "ModelMetrics", "tidyverse", "sf", "future", "terra", "future.apply"))

####load model input ####
load("code/output_data/seagrass_model_inputs.RData")
seagrass_data_long <- seagrass_data_long %>% select(-saltmean_bccm_sq_stnd, -saltmean_nep_sq_stnd, -slope_sqrt_stnd, -saltmin_bccm_sq_stnd, -saltmin_nep_sq_stnd)
seagrass_data_long <- seagrass_data_long %>%
  mutate(Survey = as.factor(substr(HKey, 1, 3)),
         HKey = as.factor(HKey),
         Year_factor = as.factor(Year))
#  categorical predictors in environmental layers
facVars <- c("substrate", "Survey")

####Eelgrass model####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
data <- data %>% mutate(Survey = as.factor(Survey)) # was still not recognizing i changed survye to factor so did it again
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))
# eelgrass in 3.71% of 20 m aggregated observations


####test variable correlation####
data_env_bccm<- data[, !grepl("nep", names(data), ignore.case = TRUE)] %>%  dplyr::select(7:49)
#enames <- names(data_env)
corrplot::corrplot(cor(data_env_bccm, use = "pairwise.complete.obs"), method = "color",  col = colorRampPalette(c("red", "orange", "white", "blue", "purple"))(200), is.corr = TRUE, tl.cex = 0.6, tl.col = "black", number.cex = 0.5, order = "hclust", type = "upper")
# Correlations close to-1 or +1 might indicate the existence of multicollinearity. one might suspect multicollinearity when the correlation between two (predictor) variables is below -0.9 or above +0.9.
data_env_nep<- data[, !grepl("bccm", names(data), ignore.case = TRUE)] %>%  dplyr::select(7:49)
#enames <- names(data_env)
corrplot::corrplot(cor(data_env_nep, use = "pairwise.complete.obs"), method = "color",  col = colorRampPalette(c("red", "orange", "white", "blue", "purple"))(200), is.corr = TRUE, tl.cex = 0.6, tl.col = "black", number.cex = 0.5, order = "hclust", type = "upper")
# 
# 
# ## test VIF
# my_model <- lm(presence~ substrate + depth_stnd +rei_sqrt_stnd + tidal_sqrt_stnd + freshwater_sqrt_stnd + 
#                  slope_stnd + NO3_stnd + saltmin_stnd + PARmin_stnd + 
#                  surftempmax_stnd + surftempcv_stnd + 
#                  tempdiff_stnd + cul_eff_stnd, data = data)
# #my_model <- lm(presence~ depth_stnd + substrate + slope_stnd + rei_stnd + DOmin_stnd + saltmean_stnd + tempcv_stnd + tempmean_stnd + freshwater_sqrt_stnd + tidal_sqrt_stnd, data = data)
# olsrr::ols_vif_tol(my_model)
# # Tolerance of <0.1 might indicate multicollinearity. VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. Ideally, the Variance Inflation Factors are below 3.
# VIFs <- CalcVIFs( dat=data_env[enames], VIFThresh=10 )
# # ones recomended to move  "tempmax_stnd"  "surftempmax_stnd" "surftempmin_stnd" "saltmean_stnd" "PARmean_stnd" "tempmean_stnd" "tempcv_stnd" "tidal_sqrt_stnd" "DOmin_stnd" "saltmin_stnd"  
# 
# # ensure ZO on all substrates
# substrates_present <- data %>% group_by(substrate) %>% summarise(n_present = sum(presence))
# substrates_present # there are eelgrass presence observations on all substrates
# 
# ### Forward feature selection test (testing one method to limit variables )
# #tested removing depth, makes horrible models
# eelgrass_ffs <- glm_ffs(data)  
# save(eelgrass_ffs, file = "code/output_data/model_results/eelgrass_ffs_variables.RData")
## most important variables in all 10 folds are substrate, depth, slope. Airtemp min and rei was important in 8 folds, DO_nep_stnd in 1 fold



####SDMtmb cv model####
#make mesh
mesh<- make_mesh(data = data, xy_cols = c("X", "Y"), cutoff = 53) 
#anything 30 and under makes the model underdispersed if no other variables (too complicated).
# anythign under 53 with other variables is underdispersed
#going to 10 km makes the model not run, likely will need to reduce fixed effects to make it work
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit cv model of spatial blocking 
plan(multisession)

# figure out random effects and spatiotemporal structure
#Do not want to add year as a fixed effect as I am not interested in temporal trends and not expecting temporal trend, just want to account that there might be interannual variability
# having spline on depth makes model better, spline on any other variable makes no change
#AUC is 0.825, tjur = 0.047, loglike -12168
m_e_0 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#add spatial
#AUC is 0.853, tjur = 0.121, loglike -13071, improves model
m_e_1 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#add spatiotemporal idd
#AUC is 0.846, tjur = 0.140, loglike -13289. Doesn't improve
m_e_2 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

#add spatiotemporal ar1
#AUC is 0.849, tjur = 0.126, loglike -13186. Doesn't improve, but better than IDD
m_e_2a <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

#spatial with random effect for survey
#but add random effect for survey type (observations from the same survey may be similar, difference in detection, observer training)
#AUC is 0.864, tjur = 0.142, loglike -12441. This improves the model
m_e_3 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + (1|Survey), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#spatial with random effect for survey and Hkey
#AUC is 0.854, tjur = 0.087, loglike -15999.
#makes Tjur and loglike so bad, so not worth including. 
#Likely because a transect goes from deep to shallow so observations on a transect would not be similar to each other in terms of eelgrass as there would be no eelgrass deep
m_e_4 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + (1|Survey) + (1 | HKey), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#spatial with random effect for survey and year
#AUC is 0.862, tjur = 0.147, loglike -12461. 
#This doesn't improve the model but doesn't make it much worse so might be good to include to account for interannual variability
m_e_5 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#spatial and spatiotemporal with random effect for survey and year
#AUC is 0.854, tjur = 0.149, loglike -12562. This doesn't improve the model
m_e_6 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

#add in fixed effects important from forward feature selection with random effects. No spatial field, 
#AUC is 0.931, tjur = 0.198, loglike -9384
m_e_7 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#add in fixed effects important from forward feature selection with random effects. test splines on rei. No spatial field, 
#AUC is 0.931, tjur = 0.198, loglike -9384, model is not improved adding a spline on rei
m_e_8 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + s(rei_stnd, k = 3) + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

# fss include random effects, yes spatial
#AUC is 0.938, tjur = 0.309, loglike -9544
m_e_9 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# fss include random effects, yes spatial, yes spatial temporal
#AUC is 0.930, tjur = 0.319, loglike -9645, this doesn't improve the model
m_e_10 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE,  time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

# fss include survey random effects but not year, yes spatial, yes spatial temporal
#AUC is 0.930, tjur = 0.319, loglike -9645. this doesn't improve the model
#This model is very similar in metrics to m_e_10. But it is also more simple than m_e_10
m_e_11 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + airtempmin_stnd + (1|Survey), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE,  time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

# fss include survey random effects but not year, yes spatial, yes spatial temporal with ar1 (each year correlated with the previous)
#AUC is 0.932, tjur = 0.326, loglike -9620. this model is better than iid
m_e_11a <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + airtempmin_stnd + (1|Survey), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE,  time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

# model indicated by looking at ffs and at variable relative importance and also considering what is important for future change, and also what resulted in highest AUC, Tjur and sum loglikelihood
# no spatial or spatial temporal
# AUC = 0.928, tjur = 0.234, loglike -9298
m_e_12 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                     airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                     saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd +#bccm variables
                     (1|Survey) + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")
# add rsdsmin AUC is 0.931, tjur = 0.206, loglike -9376 , best marginally of rsds
# add prmin AUC is 0.933, tjur = 0.205, loglike -9371, best marginally of pr
# adding another air temp values did not improve model
# adding freshwater, tidal or cul_eff did not improve model
# add salt_cv AUC is 0.934, tjur = 0.224, loglike -9327, best out of salts
# add DO or NO3 or PAR or surftemp does not improve model
# add NH4 AUC is 0.933, tjur = 0.225, loglike -9308
# add tempcv AUC is 0.928, tjur = 0.234, loglike -9298,
# if we remove prmin at this point AUC is 0.933, tjur = 0.231, loglike -9290, so slightly better 

# find nep model that works best 
#AUC = 0.932, tjur = 0.247, loglike -9315 
m_e_13 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + PARmin_nep_stnd + tempcv_nep_stnd + # + #nep variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")
# add salt_cv AUC is 0.935, tjur = 0.222, loglike -9336, best out of salts
# Add DO, NH4, NO3, surftemp does not improve model
# add PARmin AUC is 0.934, tjur = 0.238, loglike -9332
# add tempcv AUC is 0.932, tjur = 0.247, loglike -9315 

#bccm model with all same variables for nep and bccm no year factor
#AUC is 0.929, tjur 0.223, loglike -9231
m_e_14 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                     (1|Survey), # + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#bccm model with all same variables for nep and bccm with year factor
#AUC is 0.927, tjur 0.229, loglike -9316
m_e_15 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

# bccm model with all same variables for nep and bccm with spatial no year factor
# AUC 0.931, tjur 0.316, loglike -9458
m_e_16 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                      (1|Survey), # + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# bccm model with all same variables for nep and bccm with spatial with year factor
# AUC 0.930, tjur 0.312, loglike -9521
m_e_17 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# bccm model with all same variables for nep and bccm with spatial and spatial temporal (AR1)
#AUC 0.926, tjur 0.336, loglike -9626, 
m_e_18 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

#bccm model with all same variables for nep and bccm with spatial and spatial temporal (AR1) with no year random factor
#AUC 0.926, tjur 0.336, loglike -9626, same metrics to include year factor or not 
m_e_19 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                      (1|Survey),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

#bccm model with all same variables for nep and bccm with spatial and spatial temporal (IID) with no year random factor
#AUC 0.924, tjur 0.330, loglike -9646, IDD is worse that AR1
m_e_20 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                      (1|Survey),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

#nep model with all same variables for nep and bccm no year factor
#AUC is 0.930, tjur 0.231, loglike -9219
m_e_21 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                      (1|Survey), # + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#nep model with all same variables for nep and bccm with year factor
#AUC is 0.932, tjur 0.242, loglike -9320
m_e_22 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

# nep model with all same variables for nep and bccm with spatial no year factor
# AUC 0.936, tjur 0.326, loglike -9473
m_e_23 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                      (1|Survey), # + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# nep model with all same variables for nep and bccm with spatial with year factor
# AUC 0.933, tjur 0.322, loglike -9536
m_e_24 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# nep model with all same variables for nep and bccm with spatial and spatial temporal (AR1)
#AUC 0.932, tjur 0.340, loglike -9588, 
m_e_25 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

#nep model with all same variables for nep and bccm with spatial and spatial temporal (AR1) with no year random factor
#AUC 0.932, tjur 0.340, loglike -9588, same metrics to include year factor or not 
m_e_26 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                      (1|Survey),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

#nep model with all same variables for nep and bccm with spatial and spatial temporal (IID) with no year random factor
#AUC 0.930, tjur 0.331, loglike -9609, IDD is worse that AR1
m_e_27 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                      (1|Survey),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

# eval_cv <- evalStats( folds=1:numFolds,m=m_e_27, CV=cv_list_eelgrass$cv)
# eval_cv




# cv stats from all best models and save
eval_cv_bccm_nospatial <- evalStats( folds=1:numFolds,m=m_e_14, CV=cv_list_eelgrass$cv)
eval_cv_bccm_spatial <- evalStats( folds=1:numFolds,m=m_e_16, CV=cv_list_eelgrass$cv)
eval_cv_bccm_spatialtemp <- evalStats( folds=1:numFolds,m=m_e_19, CV=cv_list_eelgrass$cv)
eval_cv_nep_nospatial <- evalStats( folds=1:numFolds,m=m_e_21,CV=cv_list_eelgrass$cv)
eval_cv_nep_spatial <- evalStats( folds=1:numFolds,m=m_e_23,CV=cv_list_eelgrass$cv)
eval_cv_nep_spatialtemp <- evalStats( folds=1:numFolds,m=m_e_26,CV=cv_list_eelgrass$cv)
eval_cv_list <- list(eval_cv_bccm_nospatial, eval_cv_bccm_spatial, eval_cv_bccm_spatialtemp, eval_cv_nep_nospatial, eval_cv_nep_spatial, eval_cv_nep_spatialtemp)
save(eval_cv_list, file = "code/output_data/model_results/eval_cv.RData")

####SDMtmb full model####
# fit full model bccm 
#remove year for testing
fmodel_e_bccm_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                                    airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                    saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                                    (1|Survey), #random effect
                                mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data)

fmodel_e_bccm_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                                  airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                  saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                                  (1|Survey), #random effect
                                mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data)

fmodel_e_nep_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                                   airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                   saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                                   (1|Survey), #random effect
                                 mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data)

fmodel_e_nep_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                                 airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                 saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd  + #bccm variables
                                 (1|Survey), #random effect
                               mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data)

save(data, fmodel_e_bccm_nospatial, fmodel_e_bccm_spatial, fmodel_e_nep_nospatial, fmodel_e_nep_spatial, file = "code/output_data/model_results/final_eelgrass_model.RData")



#have a look at marginal effects
ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "depth_stnd[all]") %>% plot()
ggeffects::ggpredict(model = fmodel_e_bccm_spatial,  terms = "depth_stnd[all]") %>% plot()
ggeffects::ggpredict(model = fmodel_e_bccm_spatialtemporal,  terms = "depth_stnd[all]") %>% plot()

ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "depth_stnd[all]") %>% plot()

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "substrate") %>% plot() # more in sand and mud
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "substrate") %>% plot() # more in sand and mud

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "slope_stnd[all]") %>% plot() # presence declines with slope
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "slope_stnd[all]") %>% plot() # presence declines with slope

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "rei_stnd[all]") %>% plot() #presence declines with exposure
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "rei_stnd[all]") %>% plot() #presence declines with exposure

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "airtempmin_stnd[all]") %>% plot() # presence increases as air temp min increases
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "airtempmin_stnd[all]") %>% plot() # presence increases as air temp min increases

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "rsdsmin_stnd[all]") %>% plot() # presence increases as rsds min increases
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "rsdsmin_stnd[all]") %>% plot() # presence increases as rsds min increases

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "prmin_stnd[all]") %>% plot() # as prmin increases presence goes up
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "prmin_stnd[all]") %>% plot() # as prmin increases presence goes up

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "saltcv_bccm_stnd[all]") %>% plot() # as salinity variability increases presence goes down
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "saltcv_nep_stnd[all]") %>% plot() # as salinity variability increases presence goes down

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "NH4_bccm_stnd[all]") %>% plot() # as NH4 increases presence goes up
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "NH4_nep_stnd[all]") %>% plot() # as NH4 increases presence goes down opposite!!

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "tempcv_bccm_stnd[all]") %>% plot() # as tempcv increases presence goes up
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "tempcv_nep_stnd[all]") %>% plot() # as tempcv increases presence goes up

ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "PARmin_bccm_stnd[all]") %>% plot() # PARmin and presence flat
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "PARmin_nep_stnd[all]") %>% plot() # as PARmin increases presence goes down



visreg::visreg(fmodel_e_bccm, "depth_stnd")
visreg::visreg(fmodel_e_bccm, "DOmin_stnd")
visreg::visreg(fmodel_e_bccm, "slope_stnd")
visreg::visreg(fmodel_e_bccm, "rei_stnd")
visreg::visreg(fmodel_e_bccm, "saltmin_stnd")
visreg::visreg(fmodel_e_bccm, "tempmean_stnd")
visreg::visreg(fmodel_e_bccm, "airtempmin_stnd")
visreg::visreg(fmodel_e_bccm_nospatial, "Survey")

# Model check
tidy(fmodel_e_bccm_nospatial, conf.int = TRUE)
sanity(fmodel_e_bccm_nospatial)
tidy(fmodel_e_bccm_spatial, conf.int = TRUE)
sanity(fmodel_e_bccm_spatial)
tidy(fmodel_e_bccm_spatial, "ran_pars", conf.int = TRUE)
tidy(fmodel_e_nep_nospatial, conf.int = TRUE)
sanity(fmodel_e_nep_nospatial)
tidy(fmodel_e_nep_spatial, conf.int = TRUE)

models <- list(
  bccm_nospatial = fmodel_e_bccm_nospatial,
  bccm_spatial   = fmodel_e_bccm_spatial,
  nep_nospatial  = fmodel_e_nep_nospatial,
  nep_spatial    = fmodel_e_nep_spatial
  )

# Create a list to store evaluation results
eval_results <- list()

# Loop through each model
for (m_name in names(models)) {
  model <- models[[m_name]]
  
  # Calculate fitted values
  data[[paste0("fitted_vals_", m_name)]] <- predict(model, type = "response")$est
  
  # Prepare temporary data for threshold calculation
  data_tmp <- data
  data_tmp$fitted_vals <- data[[paste0("fitted_vals_", m_name)]]
  
  # Calculate optimal thresholds
  thresh <- calcThresh(x = data_tmp)
  
  # Add threshold-based predictions directly to original data
  data[[paste0("pred_TSS_thresh_", m_name)]]   <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1
  )
  data[[paste0("pred_kappa_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxKappa"], 0, 1
  )
  data[[paste0("pred_PCC_thresh_", m_name)]]   <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxPCC"], 0, 1
  )
  
  # Evaluate model
  eval_metrics <- evalfmod(x = data_tmp, thresh = thresh)
  
  # Save both evaluation metrics AND thresholds in the list
  eval_results[[m_name]] <- list(
    metrics = eval_metrics,
    thresholds = thresh
  )
}

save(eval_results, file = "code/output_data/model_results/eelgrass_eval_final_models.RData")

# these models have good TSS (>0.7), well calibrated (miller). for calibration (Hosmer & Lemeshow goodness-of-fit) model seems to have issues at higher predicted probabilities

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
prednames_bccm <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "Survey", "rsdsmin_stnd", "airtempmin_stnd", "prmin_stnd", "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempcv_bccm_stnd", "PARmin_bccm_stnd")
prednames_nep <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "Survey", "rsdsmin_stnd", "airtempmin_stnd", "prmin_stnd", "saltcv_nep_stnd", "NH4_nep_stnd", "tempcv_nep_stnd", "PARmin_nep_stnd")

groups <- list(
  Geomorphic = c("depth_stnd", "slope_stnd", "rei_stnd", "substrate"),
  Atmospheric = c("airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd"),
  Oceanographic = c("saltcv_bccm_stnd", "NH4_bccm_stnd", "tempcv_bccm_stnd", "PARmin_bccm_stnd"),
  Random = c("Survey")
)

plan(sequential) # the spatial models don't run well on multisession

relimp_e_bccm_nospatial <- varImp_sdmTMB(
  model = fmodel_e_bccm_nospatial,
  dat = data,
  preds = prednames_bccm,
  groups = groups,
  permute = 10
)

relimp_e_bccm_nospatial$individual
relimp_e_bccm_nospatial$grouped

# depth 70, substrate 22, slope 3.6, rei 0.5, air temp min 0.5, rsdsmin 0.6, prmin 0.3, salt cv 0.9, NH4 0.2, PARmin 0.0, Survey 1.2, tempcv 0.3
save(relimp_e_bccm_nospatial, file = "code/output_data/model_results/relimp_e_bccm_nospatial.RData")

relimp_e_bccm_spatial <- varImp_sdmTMB( model=fmodel_e_bccm_spatial,
                                   dat=data,
                                   preds=prednames_bccm,
                                   groups = groups,
                                   permute=10 ) # Number of permutations
save(relimp_e_bccm_spatial, file = "code/output_data/model_results/relimp_e_bccm_spatial.RData")
# depth 70.9, substrate 21.3, slope 4.0, rei 0.3, air temp min 0.1, rsdsmin 0.1, prmin 1.0, salt cv 0.3, NH4 0.2, PARmin 0.0, Survey 1.0, tempcv 0.7

groups <- list(
  Geomorphic = c("depth_stnd", "slope_stnd", "rei_stnd", "substrate"),
  Atmospheric = c("airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd"),
  Oceanographic = c("saltcv_nep_stnd", "NH4_nep_stnd", "tempcv_nep_stnd", "PARmin_nep_stnd"),
  Random = c("Survey")
)

relimp_e_nep_nospatial <- varImp_sdmTMB( model=fmodel_e_nep_nospatial,
                         dat=data,
                         preds=prednames_nep,
                         groups = groups,
                         permute=10 ) # Number of permutations
# depth 68.1, substrate 21.4, slope 3.4, rei 0.6, air temp min 0.8, rsdsmin 1.4, prmin 0.2, salt cv 1.0, NH4 0.5, PARmin 1.1, Survey 1.1, tempcv 0.4
save(relimp_e_nep_nospatial, file = "code/output_data/model_results/relimp_e_nep_nospatial.RData")

relimp_e_nep_spatial <- varImp_sdmTMB( model=fmodel_e_nep_spatial,
                                  dat=data,
                                  preds=prednames_nep,
                                  groups = groups,
                                  permute=10 ) # Number of permutations
# depth 69.0, substrate 20.7, slope 3.9, rei 0.5, air temp min 0.4, rsdsmin 0.2, prmin 0.5, salt cv 0.4, NH4 1.1, PARmin 1.9, Survey 1.2, tempcv 0.5
save(relimp_e_nep_spatial, file = "code/output_data/model_results/relimp_e_nep_spatial.RData")

####check residuals####
# MCMC based randomized quantile residuals (takes a while to compute)
set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_e_bccm_nospatial, mcmc_iter = 800, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_e_bccm_nospatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_e_bccm_nospatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_e_bccm_nospatial, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_e_bccm_nospatial.RData")


set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_e_bccm_spatial, mcmc_iter = 800, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_e_bccm_spatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_e_bccm_spatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_e_bccm_spatial, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_e_bccm_spatial.RData")

set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_e_nep_nospatial, mcmc_iter = 800, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_e_nep_nospatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_e_nep_nospatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_e_nep_nospatial, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_e_nep_nospatial.RData")


set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_e_nep_spatial, mcmc_iter = 800, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_e_nep_spatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_e_nep_spatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_e_nep_spatial, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_e_nep_spatial.RData")


#analytical randomized quantile approach
# data$resids_bccm_spatial <- residuals(fmodel_e_bccm_spatial, type = "mle-mvn") # randomized quantile residuals
# # check
# ggplot(data, aes(X, Y, col = resids_bccm_spatial)) + scale_colour_gradient2() +
#   geom_point() + theme_bw()
# hist(data$resids)
# qqnorm(data$resids);abline(a = 0, b = 1)

# predict(fmodel_e_bccm_nospatial) %>%
#   ggplot(aes(x = presence, y = fmodel_e_bccm_nospatial$family$linkinv(est)))+
#   geom_abline(slope = 1, intercept = 0)+
#   geom_jitter(width = 0.05, height = 0)


#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009
# model does better if you don't include oceanographic or atmospheric data. which makes sense, becasue there are no differences between past and present as geomorphic are  static
#cannot account for surveys or years in these models as random factor because only ABL, Cuke, GDK, GSU, and RSU surveys were conducted prior to 2013 and you cant make predictions of years not fit in model

data_pre2013 <- data %>% filter(Year < 2010)
unique(data_pre2013$Survey)
mesh_pre2013 <- make_mesh(data = data_pre2013, xy_cols = c("X", "Y"), cutoff = 53) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_pre2013)
barrier_mesh_pre2013 <- add_barrier_mesh(mesh_pre2013, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)
#BCCM
m_eelgrass_forecast_bccm <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                     airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                     saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd, #bccm variables
                                   mesh = barrier_mesh_pre2013, 
                                   family = binomial(link = "logit"), 
                                   spatial = FALSE, 
                                   data = data_pre2013) 
m_eelgrass_forecast_spatial_bccm <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                             airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                             saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + PARmin_bccm_stnd,
                                           mesh = barrier_mesh_pre2013, 
                                           family = binomial(link = "logit"), 
                                           spatial = TRUE, 
                                           data = data_pre2013) 
data.df <- data %>% select(presence, X, Y, depth_stnd, slope_stnd, rei_stnd, prmin_stnd, substrate, saltcv_bccm_stnd, NH4_bccm_stnd, airtempmin_stnd, rsdsmin_stnd, tempcv_bccm_stnd, PARmin_bccm_stnd, Survey, Year)

forecast <- plogis(predict(m_eelgrass_forecast_bccm, newdata = data.df %>% filter(Year > 2012))$est)
forecast_spatial <- plogis(predict(m_eelgrass_forecast_spatial_bccm, newdata = data.df %>% filter(Year > 2012))$est)

pre_2013 <- plogis(predict(m_eelgrass_forecast_bccm, newdata = data.df %>% filter(Year < 2010))$est)
pre_2013_spatial <- plogis(predict(m_eelgrass_forecast_spatial_bccm, newdata = data.df %>% filter(Year < 2010))$est)

forecast_predict_eelgrass_bccm <- data.frame(TjurR2_no_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast),
                                                    tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013)),
                                        TjurR2_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast_spatial),
                                                              tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013_spatial)),
                                         AUC_no_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast),
                                                 ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013)),
                                        AUC_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast_spatial),
                                                ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013_spatial)),
                                         type = factor(c("forecast", "training"), levels = c("training", "forecast"), ordered =  TRUE))


#NEP36
m_eelgrass_forecast_nep <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                    airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                    saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd,# + #nep variables
                                  mesh = barrier_mesh_pre2013, 
                                  family = binomial(link = "logit"), 
                                  spatial = FALSE, 
                                  data = data_pre2013) 
m_eelgrass_forecast_spatial_nep <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                            airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                            saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + PARmin_nep_stnd,# + #nep variables
                                          mesh = barrier_mesh_pre2013, 
                                          family = binomial(link = "logit"), 
                                          spatial = TRUE, 
                                          data = data_pre2013) 
data.df <- data %>% select(presence, X, Y, depth_stnd, slope_stnd, rei_stnd, prmin_stnd, substrate, saltcv_nep_stnd, NH4_nep_stnd, airtempmin_stnd, rsdsmin_stnd, tempcv_nep_stnd, PARmin_nep_stnd, Year)

forecast <- plogis(predict(m_eelgrass_forecast_nep, newdata = data.df %>% filter(Year > 2012))$est)
forecast_spatial <- plogis(predict(m_eelgrass_forecast_spatial_nep, newdata = data.df %>% filter(Year > 2012))$est)

pre_2013 <- plogis(predict(m_eelgrass_forecast_nep, newdata = data.df %>% filter(Year < 2010))$est)
pre_2013_spatial <- plogis(predict(m_eelgrass_forecast_spatial_nep, newdata = data.df %>% filter(Year < 2010))$est)

forecast_predict_eelgrass_nep <- data.frame(TjurR2_no_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast),
                                                                   tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013)),
                                             TjurR2_spatial = c(tjur(y = data.df$presence[data.df$Year > 2012], pred = forecast_spatial),
                                                                tjur(y = data.df$presence[data.df$Year < 2010], pred = pre_2013_spatial)),
                                             AUC_no_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast),
                                                                ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013)),
                                             AUC_spatial = c(ModelMetrics::auc(data.df$presence[data.df$Year > 2012], forecast_spatial),
                                                             ModelMetrics::auc(data.df$presence[data.df$Year < 2010], pre_2013_spatial)),
                                             type = factor(c("forecast", "training"), levels = c("training", "forecast"), ordered =  TRUE))

save(m_eelgrass_forecast_bccm, m_eelgrass_forecast_spatial_bccm, forecast_predict_eelgrass_bccm, m_eelgrass_forecast_nep, m_eelgrass_forecast_spatial_nep, forecast_predict_eelgrass_nep, file = "code/output_data/model_results/forecast_eelgrass_models.RData")

#bccm model with spatial does best overall at forecasting






####Eelgrass delta model with percent cover ####
dat2 <- subset(data, mean_PerCovZO > 0)
#not all surveys had records of percent cover
dat2$Survey <- factor(dat2$Survey,
                      levels = c("ABL", "BHM", "Cuk", "GSU", "Mul", "RSU"))
mesh2 <- make_mesh(dat2,
                   xy_cols = c("X", "Y"),
                   mesh = mesh$mesh
)
plot(mesh2)
barrier_mesh2 <- add_barrier_mesh(mesh2, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)
#better to have more flexible spline on depth
#better with spatial
m_e_per_1 <- sdmTMB_cv(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + freshwater_sqrt_stnd +
                         surftempmin_bccm_stnd + #bccm variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)
m_e_per_1$sum_loglik # -9048

m_e_per_2 <- sdmTMB_cv(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + freshwater_sqrt_stnd +
                         surftempmin_nep_stnd + PARmean_nep_stnd + #nep variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)
m_e_per_2$sum_loglik # -9043

#So best models that include all best factors from both nep and bccm
m_e_per_bccm_final <- sdmTMB(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + freshwater_sqrt_stnd +
                               surftempmin_bccm_stnd  + PARmean_bccm_stnd + #bccm variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)


m_e_per_nep_final <- sdmTMB(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + freshwater_sqrt_stnd +
                         surftempmin_nep_stnd + PARmean_nep_stnd +  #nep variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)

##ENDED HERE!!!
sanity(m_e_per_bccm_final)

#print(fmodel_e_delta_bccm_nospatial)

tidy(m_e_per_bccm_final)
tidy(m_e_per_bccm_final, "ran_pars", conf.int = TRUE)
tidy(m_e_per_nep_final)
tidy(m_e_per_nep_final, "ran_pars", conf.int = TRUE)

visreg::visreg(m_e_per_bccm_final, "depth_stnd")
visreg::visreg(m_e_per_bccm_final, "substrate")
visreg::visreg(m_e_per_bccm_final, "surftempmin_bccm_stnd")
visreg::visreg(m_e_per_bccm_final, "Survey")

visreg::visreg(m_e_per_nep_final, "depth_stnd")
visreg::visreg(m_e_per_nep_final, "substrate")
visreg::visreg(m_e_per_nep_final, "surftempmin_nep_stnd")
visreg::visreg(m_e_per_nep_final, "Survey")












#### Surfgrass model ####

sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)
data <- data %>% mutate(Survey = as.factor(Survey))
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))
# present in 1.4% of 20 m cells

## test VIF
# my_model <- lm(presence~ depth_stnd + substrate + rei_sqrt_stnd + tempmin_stnd +
#                  saltcv_stnd + PARmin_stnd + surftempmax_stnd + NO3_stnd, data = data)
# olsrr::ols_vif_tol(my_model)
# Tolerance of <0.1 might indicate multicollinearity. VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. Ideally, the Variance Inflation Factors are below 3.

# ensure PH on all substrates
substrates_present <- data %>% group_by(substrate) %>% summarise(n_present = sum(presence))
substrates_present # there are surgrass presence observations on all substrates

### Forward feature selection test (testing one method to limit variables )
#tested removing depth, makes horrible models
# surfgrass_ffs <- glm_ffs(data)  
# save(surfgrass_ffs, file = "code/output_data/surfgrass_ffs_variables.RData")
## variables that came out include most often depth,  airtempcv, substrate, air temp min, rei sqrt, freshwater sqrt, rei, parmin nep, freshwater. A few times surftempcv bccm and tempmean bccm and tidal came out


#make mesh # 30 size mesh means there is no underdispersion
mesh<- make_mesh(data = data, xy_cols = c("X", "Y"), cutoff = 30) 
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)

# bccm model indicated by looking at ffs and at variable relative importance and also considering what is important for future change, and also what resulted in highest AUC, Tjur and sum loglikelihood
# no spatial or spatial temporal
# AUC = 0.966, tjur = 0.276, loglike -3550
m_s_1 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + tidal_sqrt_stnd +
                      airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                      saltcv_bccm_stnd + NO3_bccm_stnd + surftempmax_bccm_stnd + #bccm variables
                      (1|Survey),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")
# best metrics to just have spline on depth with no basis function
# substrate important, rei_sqrt is better than rei (can't include both), tidal makes it better, and also when we add on the spline with basis = 3
# slope is not important, freshwater doesn't improve model, freshwater sqrt only mariginally increases Tjur and AUC but makes logliklihood marignally go down
# rsds does not improve model, NH4 does not improve model, PAR did not improve, DO did not improve
# temp min has the best metrics of the temps
# NO3 improves model, cul eff improves model
# surf temp max makes biggest improvement
# remove cul eff, tidal sqr and temp min based on relimp < 0.3

# nep model indicated by looking at ffs and at variable relative importance and also considering what is important for future change, and also what resulted in highest AUC, Tjur and sum loglikelihood
# no spatial or spatial temporal
# AUC = 0.966, tjur = 0.276, loglike -3550
m_s_2 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_nep_stnd +  surftempdiff_nep_stnd + tempmean_nep_stnd + PARmin_nep_stnd + #  PARmin_nep_stnd +  + #nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")
# will do salt min and suftempmax for both nep and bccm
#PARmin has low relimp so can be excluded. Will have to add in temp mean though 
#NH4 does not improve, DO doesnt improve, NO3 doesn't make worse

#bccm model that has features best for both nep and bccm, no spatial
# AUC = 0.967, tjur = 0.265, loglike -3563
m_s_3 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_bccm_stnd +  NO3_bccm_stnd + surftempmax_bccm_stnd + tempmean_bccm_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#bccm model that has features best for both nep and bccm with spatial
# AUC = 0.968, tjur = 0.297, loglike -3616
m_s_4 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_bccm_stnd +  NO3_bccm_stnd + surftempmax_bccm_stnd + tempmean_bccm_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#bccm model that has features best for both nep and bccm with spatialtemporal idd
# AUC = 0.964, tjur = 0.313, loglike -3535
m_s_5 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_bccm_stnd +  NO3_bccm_stnd + surftempmax_bccm_stnd + tempmean_bccm_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

#bccm model that has features best for both nep and bccm with spatialtemporal ar1
# AUC = 0.964, tjur = 0.312, loglike -3576
m_s_6 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_bccm_stnd +  NO3_bccm_stnd + surftempmax_bccm_stnd + tempmean_bccm_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

#nep model that has features best for both nep and bccm, no spatial
# AUC = 0.965, tjur = 0.259, loglike -3526
m_s_7 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_nep_stnd +  NO3_nep_stnd + surftempmax_nep_stnd + tempmean_nep_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#nep model that has features best for both nep and bccm with spatial
# AUC = 0.964, tjur = 0.302, loglike -3643
m_s_8 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_nep_stnd +  NO3_nep_stnd + surftempmax_nep_stnd + tempmean_nep_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#nep model that has features best for both nep and bccm with spatialtemporal idd
# AUC = 0.960, tjur = 0.331, loglike -3613
m_s_9 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_nep_stnd +  NO3_nep_stnd + surftempmax_nep_stnd + tempmean_nep_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

#nep model that has features best for both nep and bccm with spatialtemporal ar1
# AUC = 0.960, tjur = 0.331, loglike -3631
m_s_10 <- sdmTMB_cv(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                     airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                     saltmin_nep_stnd +  NO3_nep_stnd + surftempmax_nep_stnd + tempmean_nep_stnd +  # nep variables
                     (1|Survey),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "ar1", data = data, fold_ids = "fold")

eval_cv <- evalStats( folds=1:numFolds,m=m_s_10, CV=cv_list_seagrass$cv)
eval_cv

# If a model is unbiased bias should be close to zero
# MAE want low values
# AUC want high values above 0.9. According to Pearce and Ferrier (2000) and Jones et al. (2010) values of AUC greater than 0.9 are considered good, between 0.7 and 0.9 moderate, and less than 0.7 poor. values of 0.5 indicate that the model is no better than random.
# TSS balances sensitivity (proportion of presence observations that are correctly classified) and specificity (proportion of absence observations that are correctly classified) and is independent of the prevalence of the observations (Allouche et al. 2006). Values of TSS greater than 0.6 are considered good, between 0.2 and 0.6 moderate, and less than 0.2 poor (Jones et al. 2010; Landis and Koch 1977). 
# Accuracy is the percent of predictions which are correctly classified and varies from values of 0 to 1 where 1 is the highest accuracy.
# Kappa is a measure of agreement between observed and predicted values that accounts for chance agreements and is dependent on prevalence of the observations. Kappa range from -1 to 1 with values less than 0 representing models that are no better than random and values of 1 indicating perfect agreement (Allouche et al. 2006). 

#fit full model


fmodel_s_bccm_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                                    airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                                    saltmin_bccm_stnd +  NO3_bccm_stnd + surftempmax_bccm_stnd + tempmean_bccm_stnd+  # bccm variables
                                    (1|Survey),  #random effect
                                mesh = barrier_mesh, #random effect
                                family = binomial(link = "logit"), 
                                spatial = FALSE, 
                                data = data)
fmodel_s_bccm_spatial <- sdmTMB(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                                  airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                                  saltmin_bccm_stnd +  NO3_bccm_stnd + surftempmax_bccm_stnd + tempmean_bccm_stnd +  # bccm variables
                                  (1|Survey),  #random effect
                                  mesh = barrier_mesh, 
                                  family = binomial(link = "logit"), 
                                  spatial = TRUE, 
                                  data = data)

fmodel_s_nep_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd) + substrate + rei_sqrt_stnd + 
                                   airtempmean_stnd + airtempcv_stnd + airtempmin_stnd + prmean_stnd + #chelsa variables
                                   saltmin_nep_stnd +  NO3_nep_stnd + surftempmax_nep_stnd + tempmean_nep_stnd +  # nep variables
                                   (1|Survey),  #random effect
                                  mesh = barrier_mesh, 
                                  family = binomial(link = "logit"), 
                                  spatial = FALSE, 
                                  data = data)

# look ar marginal effects plots, note temperature variables are correlated, so these are not correct 
#have a look at marginal effects
ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "depth_stnd[all]") %>% plot()
ggeffects::ggpredict(model = fmodel_s_bccm_spatial,  terms = "depth_stnd[all]") %>% plot()
ggeffects::ggpredict(model = fmodel_s_bccm_spatialtemporal,  terms = "depth_stnd[all]") %>% plot()

ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "depth_stnd[all]") %>% plot()

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "substrate") %>% plot() # more in rock
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "substrate") %>% plot() # more in rock

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "rei_sqrt_stnd[all]") %>% plot() #presence increases with exposure
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "rei_sqrt_stnd[all]") %>% plot() #presence ? with exposure

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "airtempmin_stnd[all]") %>% plot() # presence increases as air temp min decreases
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "airtempmin_stnd[all]") %>% plot() # presence increases as air temp min 

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "airtempmean_stnd[all]") %>% plot() # presence increases as air temp mean increases
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "airtempmean_stnd[all]") %>% plot() # presence increases as air temp mean

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "airtempcv_stnd[all]") %>% plot() # presence increases as air temp cv increases
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "airtempcv_stnd[all]") %>% plot() # presence increases as air temp cv

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "prmean_stnd[all]") %>% plot() # as prmean increases presence goes up
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "prmean_stnd[all]") %>% plot() # as prmeann increases presence goes 

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "saltmin_bccm_stnd[all]") %>% plot() # as salinity min increases presence goes up
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "saltmin_nep_stnd[all]") %>% plot() # as salinity min increases presence goes

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "NO3_bccm_stnd[all]") %>% plot() # as NO3 increases presence goes down
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "NO3_nep_stnd[all]") %>% plot() # as NO3 increases presence goes 

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "tempmean_bccm_stnd[all]") %>% plot() # as tempmean increases presence goes up
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "tempmean_nep_stnd[all]") %>% plot() # as tempmean increases presence goes

ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "surftempmax_bccm_stnd[all]") %>% plot() # as surf temp max increases presences goes down
ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "surftempmax_nep_stnd[all]") %>% plot() # 


# Model check
tidy(fmodel_s_bccm_nospatial, conf.int = TRUE)
sanity(fmodel_s_bccm_nospatial)
tidy(fmodel_s_bccm_spatial, conf.int = TRUE)
sanity(fmodel_s_bccm_spatial)
tidy(fmodel_s_bccm_spatialtemporal, conf.int = TRUE)
sanity(fmodel_s_bccm_spatialtemporal)
tidy(fmodel_s_bccm_spatial, "ran_pars", conf.int = TRUE)
tidy(fmodel_s_nep_nospatial, conf.int = TRUE)
sanity(fmodel_s_nep_nospatial)
tidy(fmodel_s_nep_spatial, conf.int = TRUE)
sanity(fmodel_s_nep_spatialtemporal)
tidy(fmodel_s_nep_spatial, conf.int = TRUE)
sanity(fmodel_s_nep_spatialtemporal)

models <- list(
  bccm_nospatial = fmodel_s_bccm_nospatial,
  bccm_spatial   = fmodel_s_bccm_spatial,
  bccm_spatialtemporal   = fmodel_s_bccm_spatialtemporal,
  nep_nospatial  = fmodel_s_nep_nospatial,
  nep_spatial    = fmodel_s_nep_spatial,
  nep_spatialtemporal    = fmodel_s_nep_spatialtemporal
)

# Create a list to store evaluation results
eval_results <- list()

# Loop through each model
for (m_name in names(models)) {
  model <- models[[m_name]]
  
  # Calculate fitted values
  data[[paste0("fitted_vals_", m_name)]] <- predict(model, type = "response")$est
  
  # Prepare temporary data for threshold calculation
  data_tmp <- data
  data_tmp$fitted_vals <- data[[paste0("fitted_vals_", m_name)]]
  
  # Calculate optimal thresholds
  thresh <- calcThresh(x = data_tmp)
  
  # Add threshold-based predictions directly to original data
  data[[paste0("pred_TSS_thresh_", m_name)]]   <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1
  )
  data[[paste0("pred_kappa_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxKappa"], 0, 1
  )
  data[[paste0("pred_PCC_thresh_", m_name)]]   <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxPCC"], 0, 1
  )
  
  # Evaluate model
  eval_metrics <- evalfmod(x = data_tmp, thresh = thresh)
  
  # Save both evaluation metrics AND thresholds in the list
  eval_results[[m_name]] <- list(
    metrics = eval_metrics,
    thresholds = thresh
  )
}


save(eval_results, file = "code/output_data/model_results/surfgrass_eval_final_models.RData")

# this model has good TSS, is well calibrated (miller). for calibration (Hosmer & Lemeshow goodness-of-fit) model seems to have issues at higher predicted probabilities

plan(sequential) # the spatial models don't run well on multisession
# Variable importance (randomization and permutation method)
prednames_bccm <- c("depth_stnd", "rei_sqrt_stnd", "substrate", "airtempmean_stnd", "airtempcv_stnd",
               "airtempmin_stnd", "prmean_stnd", "saltmin_bccm_stnd", "NO3_bccm_stnd", "tempmean_bccm_stnd", "surftempmax_bccm_stnd", "Survey")
prednames_nep <- c("depth_stnd", "rei_sqrt_stnd", "substrate", "airtempmean_stnd", "airtempcv_stnd",
                    "airtempmin_stnd", "prmean_stnd", "saltmin_nep_stnd", "surftempdiff_nep_stnd", "tempmean_nep_stnd", "PARmin_nep_stnd", "Survey")
relimp_s_bccm_nospatial <- varImp( model=fmodel_s_bccm_nospatial,
                  dat=data,
                  preds=prednames_bccm,
                  permute=10 ) # Number of permutations

# depth 31, airtempcv 19, airtemp mean 6, airtempmin 3, NO3 0.4, prmean 2, rei sqrt 0.7, saltmin 4.3, substrte 2.9, surfmax 6.1, survey 20.9, temp mean 2.6
save(relimp_s_bccm_nospatial, file = "code/output_data/model_results/relimp_s_bccm_nospatial.RData")

relimp_s_bccm_spatial <- varImp( model=fmodel_s_bccm_spatial,
                                 dat=data,
                                 preds=prednames_bccm,
                                 permute=10 ) # Number of permutations
save(relimp_s_bccm_spatial, file = "code/output_data/model_results/relimp_s_bccm_spatial.RData")

relimp_s_nep_nospatial <- varImp( model=fmodel_s_nep_nospatial,
                                 dat=data,
                                 preds=prednames_nep,
                                 permute=10 ) # Number of permutations
save(relimp_s_bccm_spatial, file = "code/output_data/model_results/relimp_s_bccm_spatial.RData")

# depth 55.9, PARmin 10.1, rei_sqrt 1.8, salt cv 1.9, substrate 8.2, surftemp max 10.3, temp min 7.9, NO3 3.8

####check residuals####
# MCMC based randomized quantile residuals (takes a while to compute)
# set.seed(123)
# samps <- sdmTMBextra::predict_mle_mcmc(fmodel_s_bccm, mcmc_iter = 800, mcmc_warmup = 400)
# mcmc_res <- residuals(fmodel_s_bccm, type = "mle-mcmc", mcmc_samples = samps)
# qqnorm(mcmc_res)
# abline(0, 1)

#analytical randomized quantile approach
data$resids <- residuals(fmodel_s_bccm, type = "mle-mvn") # randomized quantile residuals
# check
ggplot(data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + theme_bw()
hist(data$resids)
qqnorm(data$resids);abline(a = 0, b = 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_s_bccm_spatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_s_bccm_spatial, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

predict(fmodel_s_bccm) %>%
  ggplot(aes(x = presence, y = fmodel_s_bccm$family$linkinv(est)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_jitter(width = 0.05, height = 0)



#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009
data_pre2013 <- data %>% filter(Year < 2010)
mesh_pre2013 <- make_mesh(data = data_pre2013, xy_cols = c("X", "Y"), cutoff = 15) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_pre2013)
barrier_mesh_pre2013 <- add_barrier_mesh(mesh_pre2013, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

m_surfgrass_forecast <- sdmTMB(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd + tempmin_stnd +
                                 saltcv_stnd + PARmin_stnd + surftempmax_stnd + NO3_stnd,
                      mesh = barrier_mesh_pre2013, 
                      family = binomial(link = "logit"), 
                      spatial = FALSE, 
                      data = data_pre2013) 
m_surfgrass_forecast_spatial <- sdmTMB(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd + tempmin_stnd +
                                         saltcv_stnd + PARmin_stnd + surftempmax_stnd + NO3_stnd,
                                      mesh = barrier_mesh_pre2013, 
                                      family = binomial(link = "logit"), 
                                      spatial = TRUE, 
                                      data = data_pre2013) 
data.df <- data %>% dplyr::select(presence, X, Y, depth_stnd, rei_sqrt_stnd, substrate, tempmin_stnd, saltcv_stnd, PARmin_stnd, surftempmax_stnd, NO3_stnd,  Year)

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
## AUC doesnt drop, stays at 0.95 and Tjur drops from 0.17 to 0.14
# having spatial random field makes current day predictions better, having it for forecasting makes it worse

save(data, fmodel_s_bccm, relimp, thresh, r_ret, eval_cv, eval_fmod, forecast_predict_surfgrass, file = "code/output_data/final_surfgrass_model.RData")
