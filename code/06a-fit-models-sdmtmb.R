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
# data_env_bccm<- data[, !grepl("nep", names(data), ignore.case = TRUE)] %>%  dplyr::select(7:47)
# #enames <- names(data_env)
# corrplot::corrplot(cor(data_env_bccm, use = "pairwise.complete.obs"), method = "color",  col = colorRampPalette(c("red", "orange", "white", "blue", "purple"))(200), is.corr = TRUE, tl.cex = 0.6, tl.col = "black", number.cex = 0.5, order = "hclust", type = "upper")
# # Correlations close to-1 or +1 might indicate the existence of multicollinearity. one might suspect multicollinearity when the correlation between two (predictor) variables is below -0.9 or above +0.9.
# data_env_nep<- data[, !grepl("bccm", names(data), ignore.case = TRUE)] %>%  dplyr::select(7:47)
# #enames <- names(data_env)
# corrplot::corrplot(cor(data_env_nep, use = "pairwise.complete.obs"), method = "color",  col = colorRampPalette(c("red", "orange", "white", "blue", "purple"))(200), is.corr = TRUE, tl.cex = 0.6, tl.col = "black", number.cex = 0.5, order = "hclust", type = "upper")
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
# eelgrass_ffs <- glm_ffs(data)  ### variables that came out include depth, substrate, slope_stnd, rei_stnd, and twice tidal sqrt, and once each salt mean and freshwater and temp max
# save(eelgrass_ffs, file = "code/output_data/model_results/eelgrass_ffs_variables.RData")
# # most important variables in all 10 folds are substrate, depth, slope. Airtemp min and rei was important in 5 folds, tidal sqrt in 2 folds and surf temp min bccm in 1 fold



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
#AUC is 0.853, tjur = 0.121, loglike -13071
m_e_1 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#add spatiotemporal
#AUC is 0.846, tjur = 0.140, loglike -13289. Doesn't improve
m_e_2 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

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

#add in fixed effects important from forward feature selection. No spatial field, and no random effect
#AUC is 0.926, tjur = 0.190, loglike -9327
m_e_7 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd  + tidal_sqrt_stnd + airtempmin_stnd, mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

# ffs with include random effects, no spatial
#AUC is 0.930, tjur = 0.200, loglike -9379
m_e_8 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + tidal_sqrt_stnd + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

# fss include random effects, yes spatial
#AUC is 0.930, tjur = 0.200, loglike -9379
m_e_9 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + tidal_sqrt_stnd + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# fss include random effects, yes spatial
#AUC is 0.938, tjur = 0.259, loglike -9569
m_e_9 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + tidal_sqrt_stnd + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# fss include random effects, yes spatial, yes spatial temporal
#AUC is 0.929, tjur = 0.313, loglike -9698
m_e_10 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + tidal_sqrt_stnd + airtempmin_stnd + (1|Survey) + (1 | Year_factor), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE,  time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

# fss include survey random effects but not year, yes spatial, yes spatial temporal
#AUC is 0.930, tjur = 0.317, loglike -9668
m_e_11 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + tidal_sqrt_stnd + airtempmin_stnd + (1|Survey), mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE,  time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

# model indicated by looking at ffs and at variable relative importance and also considering what is important for future change, and also what resulted in highest AUC, Tjur and sum loglikelihood
#AUC is 0.932, tjur 0.227, loglike -9321
m_e_12 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                     saltcv_bccm_stnd + NH4_bccm_stnd + #bccm variables
                     (1|Survey) + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")
# m_e_3 bccm model add freshwater AUC is 0.930, tjur 0.196, loglike -9288, so don't keep it doesn't make model better or worse
# m_e_3 add domin AUC is 0.931, tjur 0.197, loglike -9290, # m_e_3 add domean AUC is 0.930, tjur 0.197, loglike -9303, Do make model worse domin has 0 relimp
# m_e_3 add salt min AUC is 0.932, tjur 0.200, loglike -9279 # m_e_3 add salt mean AUC is 0.932, tjur 0.198, loglike -9281 
# m_e_3 add salt cv AUC is 0.932, tjur 0.204, loglike -9264; this makes model better
# m_e_3 add temp diff AUC is 0.926, tjur 0.213, loglike -9285 # m_e_3 add temp cv AUC is 0.924, tjur 0.208, loglike -9276 
# m_e_3 add temp max AUC is 0.922, tjur 0.209, loglike -9288# m_e_3 add temp mean AUC is 0.927, tjur 0.208, loglike -9283
# m_e_3 add temp min AUC is 0.932, tjur 0.203, loglike -9270; no temps make the model better, temp min makes less worse
# m_e_3 add NH4 AUC is 0.931, tjur 0.211, loglike -9241 ,  makes model better add NH4
# m_e_3 add NO3 AUC is 0.929, tjur 0.211, loglike -9274, don't add N03
# m_e_3 add cul_eff AUC is 0.930, tjur 0.212, loglike -9254 , also has 0 relimp
# m_e_3 add PARmean AUC is 0.931, tjur 0.211, loglike -9253# m_e_3 add PARmax AUC is 0.931, tjur 0.209, loglike -9250
# m_e_3 add PARmin AUC is 0.931, tjur 0.214, loglike -9262 # PAR doesn't make the models better
#m_e_3 add surftempmax_bccm_stnd AUC is 0.927, tjur 0.210, loglike -9257#m_e_3 add surftempmean_bccm_stnd AUC is 0.929, tjur 0.210, loglike -9258
#m_e_3 add surftempmin_bccm_stnd AUC is 0.931, tjur 0.212, loglike -9255#m_e_3 addsurftempdiff_bccm_stnd AUC is 0.927, tjur 0.214, loglike -9252
#m_e_3 add surftempcv_bccm_stnd AUC is 0.930, tjur 0.212, loglike -9251 # surf temp doesn't make the models better, and cv has 0 relimp
#m_e_3 prmax_stnd AUC is 0.931, tjur 0.207, loglike -9254 #m_e_3 prmean_stnd AUC is 0.931, tjur 0.208, loglike -9253
#m_e_3 prmin_stnd AUC is 0.931, tjur 0.208, loglike -9248 
#m_e_3 prcv_stnd AUC is 0.931, tjur 0.210, loglike -9245 # pr doesn't make model better but if had to add do prmcv
#m_e_3 rsdsmin_stnd AUC is 0.932, tjur 0.216, loglike -9238, makes model better! add rsdsmin
#m_e_3 rsdsmean_stnd AUC is 0.931, tjur 0.215, loglike -9243 #m_e_3 rsdsmax_stnd AUC is 0.931, tjur 0.215, loglike -9246
#m_e_3 rsdscv_stnd AUC is 0.932, tjur 0.216, loglike -9238 (rsds min and cv highly correlated so pick one)
# tidal ends up not important, if you remove tidal log likelihood increases, also remove prcv, it doesn't improve model

#model with all same variables for nep and bccm
#m_e_3a AUC is 0.933, tjur 0.234, loglike -9319
m_e_13 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempmin_bccm_stnd + PARmin_bccm_stnd  + #bccm variables
                     (1|Survey) + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

# best bccm model with spatial AUC 0.937, tjur 0.322, loglike -9543
m_e_14 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                     saltcv_bccm_stnd + NH4_bccm_stnd + #bccm variables
                     (1|Survey) + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#model with all same variables for nep and bccm with spatial AUC 0.936, tjur 0.311, loglike -9534
m_e_15 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempmin_bccm_stnd + PARmin_bccm_stnd + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

#model with all same variables for nep and bccm with spatial AUC 0.928, tjur 0.324, loglike -9702, with spatiotemporal
m_e_16 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd +  
                      airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                      saltcv_bccm_stnd + NH4_bccm_stnd + tempmin_bccm_stnd + PARmin_bccm_stnd + #bccm variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

# no nep variables auc 0.931, tjur = 0.196, loglike -9293
# m_e_4 nep add saltmin auc 0.935, tjur = 0.216, loglike -9251 # m_e_4 nep add saltmean auc 0.935, tjur = 0.216, loglike -9251
# m_e_4 nep add saltcv auc 0.933, tjur = 0.214, loglike -9247; all salts basically the same, so just add same as bccm
# m_e_4 nep add nh4 auc 0.932, tjur = 0.220, loglike -9271;doesn't make model better
# m_e_4 nep add no3 auc 0.933, tjur = 0.228, loglike -9270 ;doesn't make model better
# m_e_3 add cul_eff AUC is 0.933, tjur 0.214, loglike -9263; deoesn't make model better 
# m_e_3 add temp diff AUC is 0.929, tjur 0.224, loglike -9264 # m_e_3 add temp cv AUC is 0.930, tjur 0.226, loglike -9237
# m_e_3 add temp max AUC is 0.927, tjur 0.213, loglike -9258 # m_e_3 add temp mean AUC is 0.930, tjur 0.208, loglike -9272
# m_e_3 add temp min AUC is 0.935, tjur 0.224, loglike -9232; makes model better
# m_e_3 add PARmean AUC is 0.934, tjur 0.224, loglike -9238 # m_e_3 add PARmax AUC is 0.935, tjur 0.223, loglike -9241
# m_e_3 add PARmin AUC is 0.933, tjur 0.237, loglike -9221 ; makes model better
#m_e_3 add surftempmax AUC is 0.93, tjur 0.232, loglike -9234 #m_e_3 add surftempmean AUC is 0.930, tjur 0.228, loglike -9232
#m_e_3 add surftempmin AUC is 0.932, tjur 0.227, loglike -9224 #m_e_3 add surftempcv AUC is 0.933, tjur 0.237, loglike -9228 
#m_e_3 add surftempdiff AUC is 0.931, tjur 0.236, loglike -9239; doesn't make model better
# m_e_3 add domin AUC is 0.932, tjur 0.231, loglike -9244, # m_e_3 add domean AUC is 0.933, tjur 0.236, loglike -9246 ; doesn't make model better 
# drop prcv and model goes up 
#best nep model auc 0.933 , tjur = 0.245, loglike -9330
m_e_17 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                     airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                     saltcv_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd + #nep36 variables
                     (1|Survey) + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#model with all same variables between bccm and nep auc 0.934 , tjur = 0.240, loglike -9333
m_e_18 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                     airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                     saltcv_nep_stnd + NH4_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd + #nep36 variables
                     (1|Survey) + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

#best nep model with spatial auc 0.934 , tjur = 0.326, loglike -9522
m_e_19 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                     airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                     saltcv_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd + #nep36 variables
                     (1|Survey) + (1|Year_factor),  #random effect
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# model with all same variables between bccm and nep with spatial auc 0.936 , tjur = 0.326, loglike -9514
m_e_20 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                      airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd + #nep36 variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

# model with all same variables between bccm and nep with spatial auc 0.930 , tjur = 0.328, loglike -9634 and spatiotemporal
m_e_21 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                      airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                      saltcv_nep_stnd + NH4_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd + #nep36 variables
                      (1|Survey) + (1|Year_factor),  #random effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), spatial = TRUE, time = "Year", spatiotemporal = "IID", data = data, fold_ids = "fold")

# cv stats from all best models and save
eval_cv_bccm_nospatial <- evalStats( folds=1:numFolds,m=m_e_13, CV=cv_list_eelgrass$cv)
eval_cv_bccm_spatial <- evalStats( folds=1:numFolds,m=m_e_15, CV=cv_list_eelgrass$cv)
eval_cv_nep_nospatial <- evalStats( folds=1:numFolds,m=m_e_18,CV=cv_list_eelgrass$cv)
eval_cv_nep_spatial <- evalStats( folds=1:numFolds,m=m_e_20,CV=cv_list_eelgrass$cv)
eval_cv_list <- list(eval_cv_bccm_nospatial, eval_cv_bccm_spatial, eval_cv_nep_nospatial, eval_cv_nep_spatial)
save(eval_cv_list, file = "code/output_data/model_results/eval_cv.RData")

####SDMtmb full model####
# fit full model bccm 
#remove year for testing
fmodel_e_bccm_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                  airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                  saltcv_bccm_stnd + NH4_bccm_stnd + tempmin_bccm_stnd + PARmin_bccm_stnd + #bccm variables
                                  (1|Survey) + (1|Year_factor),  #random effect
                                mesh = barrier_mesh, 
                                family = binomial(link = "logit"), 
                                spatial = FALSE, 
                                data = data)

fmodel_e_bccm_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                  airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                  saltcv_bccm_stnd + NH4_bccm_stnd + tempmin_bccm_stnd + PARmin_bccm_stnd + #bccm variables
                                  (1|Survey) + (1|Year_factor),  #random effect
                                mesh = barrier_mesh, 
                                family = binomial(link = "logit"), 
                                spatial = TRUE, 
                                data = data)

fmodel_e_nep_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                   airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                   saltcv_nep_stnd + NH4_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd + #nep36 variables
                                   (1|Survey) + (1|Year_factor),  #random effect
                                 mesh = barrier_mesh, 
                                 family = binomial(link = "logit"), 
                                 spatial = FALSE, 
                                 data = data)

fmodel_e_nep_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                   airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                   saltcv_nep_stnd + NH4_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd + #nep36 variables
                                   (1|Survey) + (1|Year_factor),  #random effect
                                 mesh = barrier_mesh, 
                                 family = binomial(link = "logit"), 
                                 spatial = TRUE, 
                                 data = data)

save(data, fmodel_e_bccm_nospatial, fmodel_e_bccm_spatial, fmodel_e_nep_nospatial, fmodel_e_nep_spatial, file = "code/output_data/model_results/final_eelgrass_model.RData")


#have a look at marginal effects
ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "depth_stnd[all]") %>% plot()
ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "depth_stnd[all]") %>% plot()
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "substrate") %>% plot() # more in sand and mud
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "substrate") %>% plot() # more in sand and mud
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "slope_stnd[-2:8]") %>% plot() # presence declines with slope
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "slope_stnd[-2:8]") %>% plot() # presence declines with slope
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "rei_stnd[-1:20]") %>% plot() #presence declines with exposure
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "rei_stnd[-1:20]") %>% plot() #presence declines with exposure
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "airtempmin_stnd[-6:2]") %>% plot() # presence increases as air temp min increases
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "airtempmin_stnd[-6:2]") %>% plot() # presence increases as air temp min increases
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "rsdsmin_stnd[-4:5]") %>% plot() # presence increases as rsds min increases
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "rsdsmin_stnd[-4:5]") %>% plot() # presence increases as rsds min increases
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "saltcv_bccm_stnd[-1:57]") %>% plot() # as salinity variability increases presence goes down
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "saltcv_nep_stnd[-1:56]") %>% plot() # as salinity variability increases presence goes down
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "NH4_bccm_stnd[-2:19]") %>% plot() # as ammonium increases presence goes up
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "NH4_nep_stnd[-1:7]") %>% plot() # as ammonium increases presence goes down DIFFERENCE BETWEEN BCCM AND NEP36
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "tempmin_bccm_stnd[-5:3]") %>% plot() # as tempmin increases presence goes down
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "tempmin_nep_stnd[-8:3]") %>% plot() # as tempmin increases presence goes down
ggeffects::ggeffect(model = fmodel_e_bccm_nospatial,  terms = "PARmin_bccm_stnd[-3:2]") %>% plot() # as PARmin increases presence goes down
ggeffects::ggeffect(model = fmodel_e_nep_nospatial,  terms = "PARmin_nep_stnd[-3:2]") %>% plot() # as PARmin increases presence goes down



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
sanity(fmodel_e_nep_spatial)

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
  data[[paste0("pred_TSS_thresh_", m_name)]]   <- ifelse(data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1)
  data[[paste0("pred_kappa_thresh_", m_name)]] <- ifelse(data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxKappa"], 0, 1)
  data[[paste0("pred_PCC_thresh_", m_name)]]   <- ifelse(data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxPCC"], 0, 1)
  # Evaluate model
  eval_results[[m_name]] <- evalfmod(x = data_tmp, thresh = thresh)
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
prednames_bccm <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "Survey", "Year_factor", "rsdsmin_stnd", "airtempmin_stnd", "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempmin_bccm_stnd", "PARmin_bccm_stnd")
prednames_nep <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "Survey", "Year_factor", "rsdsmin_stnd", "airtempmin_stnd", "saltcv_nep_stnd", "NH4_nep_stnd", "tempmin_nep_stnd", "PARmin_nep_stnd")

relimp_e_bccm_nospatial <- varImp( model=fmodel_e_bccm_nospatial,
                  dat=data,
                  preds=prednames_bccm,
                  permute=10 ) # Number of permutations
# depth 69.9, substrate 22.1, slope 3.4, rei 0.4, air temp min 0.3, rsdsmin 0.6, salt cv 0.8, NH4 0.2, PARmin 0.0, Survey 1.1, tempmin 0.3, year 0.8
save(relimp_e_bccm_nospatial, file = "code/output_data/model_results/relimp_e_bccm_nospatial.RData")
plan(sequential) # the spatial models don't run well on multisession
relimp_e_bccm_spatial <- varImp( model=fmodel_e_bccm_spatial,
                                   dat=data,
                                   preds=prednames_bccm,
                                   permute=10 ) # Number of permutations
save(relimp_e_bccm_spatial, file = "code/output_data/model_results/relimp_e_bccm_spatial.RData")
# depth 71.0, substrate 21.5, slope 3.9, rei 0.3, air temp min 0.0, rsdsmin 0.1, salt cv 0.3, NH4 0.3, PARmin 0.2, Survey 1.2, tempmin 0.5, year 0.8


relimp_e_nep_nospatial <- varImp( model=fmodel_e_nep_nospatial,
                         dat=data,
                         preds=prednames_nep,
                         permute=10 ) # Number of permutations
# depth 68.8, substrate 21.7, slope 3.2, rei 0.6, air temp min 0.8, rsdsmin 1.1, salt cv 0.7, NH4 0.3, PARmin 0.9, Survey 1.0, tempmin 0.4, year 0.6
save(relimp_e_nep_nospatial, file = "code/output_data/model_results/relimp_e_nep_nospatial.RData")

relimp_e_nep_spatial <- varImp( model=fmodel_e_nep_spatial,
                                  dat=data,
                                  preds=prednames_nep,
                                  permute=10 ) # Number of permutations
# depth 69.3, substrate 20.8, slope 3.8, rei 0.5, air temp min 0.0, rsdsmin 0.1, salt cv 0.5, NH4 1.0, PARmin 1.8, Survey 1.1, tempmin 0.6, year 0.5
save(relimp_e_nep_spatial, file = "code/output_data/model_results/relimp_e_nep_spatial.RData")

####check residuals####
# MCMC based randomized quantile residuals (takes a while to compute)
# set.seed(123)
# samps <- sdmTMBextra::predict_mle_mcmc(fmodel_e_bccm, mcmc_iter = 800, mcmc_warmup = 400)
# mcmc_res <- residuals(fmodel_e_bccm, type = "mle-mcmc", mcmc_samples = samps)
# qqnorm(mcmc_res)
# abline(0, 1)

#analytical randomized quantile approach
data$resids_bccm_spatial <- residuals(fmodel_e_bccm_spatial, type = "mle-mvn") # randomized quantile residuals
# check
ggplot(data, aes(X, Y, col = resids_bccm_spatial)) + scale_colour_gradient2() +
  geom_point() + theme_bw()
hist(data$resids)
qqnorm(data$resids);abline(a = 0, b = 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_e_bccm_spatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_e_bccm_spatial, return_DHARMa = TRUE)
plot(r_ret)
DHARMa::testResiduals(r_ret)

predict(fmodel_e_bccm) %>%
  ggplot(aes(x = presence, y = fmodel_e_bccm$family$linkinv(est)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_jitter(width = 0.05, height = 0)
# need to save residuals from all models still!!

#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009

#NEED TO REWRITE THIS TO ACCOUNT FOR RANDOM FACTORS!!

data_pre2013 <- data %>% filter(Year < 2010)
mesh_pre2013 <- make_mesh(data = data_pre2013, xy_cols = c("X", "Y"), cutoff = 53) # tested several mesh sizes between 20- 10 km and 15 had highest AUC
plot(mesh_pre2013)
barrier_mesh_pre2013 <- add_barrier_mesh(mesh_pre2013, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)
#BCCM
m_eelgrass_forecast_bccm <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                saltcv_bccm_stnd + NH4_bccm_stnd + tempmin_bccm_stnd + PARmin_bccm_stnd,# + #bccm variables
                               # (1|Survey) + (1|Year_factor),  # doesn't work with random factor yet
                               mesh = barrier_mesh_pre2013, 
                               family = binomial(link = "logit"), 
                               spatial = FALSE, 
                               data = data_pre2013) 
m_eelgrass_forecast_spatial_bccm <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                        airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                        saltcv_bccm_stnd + NH4_bccm_stnd + tempmin_bccm_stnd + PARmin_bccm_stnd,# + #bccm variables
                                      # (1|Survey) + (1|Year_factor),  # doesn't work with random factor yet
                              mesh = barrier_mesh_pre2013, 
                              family = binomial(link = "logit"), 
                              spatial = TRUE, 
                              data = data_pre2013) 
data.df <- data %>% select(presence, X, Y, depth_stnd, slope_stnd, rei_stnd, substrate, saltcv_bccm_stnd, NH4_bccm_stnd, airtempmin_stnd, rsdsmin_stnd, tempmin_bccm_stnd, PARmin_bccm_stnd, Year)

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
                                     airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                     saltcv_nep_stnd + NH4_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd,# + #nep variables
                                   # (1|Survey) + (1|Year_factor),  # doesn't work with random factor yet
                                   mesh = barrier_mesh_pre2013, 
                                   family = binomial(link = "logit"), 
                                   spatial = FALSE, 
                                   data = data_pre2013) 
m_eelgrass_forecast_spatial_nep <- sdmTMB(formula = presence ~ s(depth_stnd, k = 3) + substrate + slope_stnd + rei_stnd + 
                                             airtempmin_stnd + rsdsmin_stnd + #chelsa variables
                                             saltcv_nep_stnd + NH4_nep_stnd + tempmin_nep_stnd + PARmin_nep_stnd,# + #nep variables
                                           # (1|Survey) + (1|Year_factor),  # doesn't work with random factor yet
                                           mesh = barrier_mesh_pre2013, 
                                           family = binomial(link = "logit"), 
                                           spatial = TRUE, 
                                           data = data_pre2013) 
data.df <- data %>% select(presence, X, Y, depth_stnd, slope_stnd, rei_stnd, substrate, saltcv_nep_stnd, NH4_nep_stnd, airtempmin_stnd, rsdsmin_stnd, tempmin_nep_stnd, PARmin_nep_stnd, Year)

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
#better to have fleximble spline on depth
# don't have random effect of year as not all datasets had all years and dataset is not significantly smaller
#better with spatial
m_e_per_1 <- sdmTMB_cv(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + #slope_stnd + 
                         surftempmin_bccm_stnd + #saltcv_bccm_stnd + #NO3_bccm_stnd + # + #PARmin_bccm_stnd +  #tempcv_bccm_stnd +  #bccm variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)
m_e_per_1$sum_loglik # -9050

m_e_per_2 <- sdmTMB_cv(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + #slope_stnd + 
                         surftempmin_nep_stnd +  #bccm variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)
m_e_per_2$sum_loglik # -9055

#if using spatial can get rid of slope, tempcv, PAR min, salt cv, NO3
m_e_per_bccm_final <- sdmTMB(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + 
                               surftempmin_bccm_stnd  +  #bccm variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)


m_e_per_nep_final <- sdmTMB(formula = mean_PerCovZO ~ s(depth_stnd) + substrate + 
                         surftempmin_nep_stnd +  #bccm variables
                         (1|Survey),  #random effect
                       mesh = barrier_mesh2, 
                       family = Gamma(link = "log"), 
                       spatial = TRUE, 
                       data = dat2)


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
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))
# present in 1.39% of 20 m cells

## test VIF
my_model <- lm(presence~ depth_stnd + substrate + rei_sqrt_stnd + tempmin_stnd +
                 saltcv_stnd + PARmin_stnd + surftempmax_stnd + NO3_stnd, data = data)
olsrr::ols_vif_tol(my_model)
# Tolerance of <0.1 might indicate multicollinearity. VIF exceeding 5 requires further investigation, whereas VIFs above 10 indicate multicollinearity. Ideally, the Variance Inflation Factors are below 3.

# ensure PH on all substrates
substrates_present <- data %>% group_by(substrate) %>% summarise(n_present = sum(presence))
substrates_present # there are surgrass presence observations on all substrates

### Forward feature selection test (testing one method to limit variables )
#tested removing depth, makes horrible models
surfgrass_ffs <- glm_ffs(data)  ### variables that came out include depth, substrate, slope_stnd, rei_stnd
save(surfgrass_ffs, file = "code/output_data/surfgrass_ffs_variables.RData")
## variables that came out include most often depth, rei sqrt, rei, temp min, tidal sqrt, substrate. A few times slope, cul eff, freshwater came out

library(MASS)
# test stepwise model
datasub <- data %>% dplyr::select(substrate:cul_eff_stnd, presence) %>% dplyr::select(-rei_stnd, -tidal_stnd, -saltmean_stnd, -PARmax_stnd)
full_model <- glm(presence ~ ., data = datasub, family = binomial(link = "logit"))

# Perform stepwise selection (both directions)
stepwise_model <- step(full_model, direction = "both")

# View the selected model
summary(stepwise_model)


#make mesh # tested several mesh sizes between 20- 10 km and 15 had highest AUC
mesh<- make_mesh(data = data, xy_cols = c("X", "Y"), cutoff = 12) 
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)
# model selected by ffs without spatial field, none as splines AUC = 0.85, tjur 0.122. 
#Making depth a spline makes it worse (AUC = 0.847; Tjur = 0.122), 
#rei sqrt spline MAKES auc BETTER AND TJUR WORSE  AUC = 0.87, tjur = 0.104), #CHOOSE REI SQRT OVER REI rei instead of rei sqrt 0.83, tjur 0.125 # rei spline auc = 0.84 tjur = 0.07
# make tidal a spline AUC = 0.87, tjur = 0.104 doesnt improve model
m_s_1 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + s(rei_sqrt_stnd, k = 3) + tempmin_stnd + tidal_sqrt_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")
roc <- pROC::roc(m_s_1$data$presence, plogis(m_s_1$data$cv_predicted))
auc <- pROC::auc(roc)
auc

# model selected by ffs with spatial field  AUC = 0.90, tjur 0.20
m_s_2 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + s(rei_sqrt_stnd, k = 3) + tempmin_stnd + tidal_sqrt_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = TRUE, 
                   data = data, 
                   fold_ids = "fold")
roc <- pROC::roc(m_s_2$data$presence, plogis(m_s_2$data$cv_predicted))
auc <- pROC::auc(roc)
auc

# model with all ffs auc = 0.874, tjur = 0.121. This is better than m_s_1
m_s_3 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + s(rei_sqrt_stnd, k = 3) + tempmin_stnd + tidal_sqrt_stnd +
                     slope_stnd + freshwater_sqrt_stnd + cul_eff_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")

# model with all ffs and extra auc = 0.95, tjur 0.268
m_s_4 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + s(rei_sqrt_stnd, k = 3) + tempmin_stnd + tidal_sqrt_stnd +
                     slope_stnd + freshwater_sqrt_stnd + cul_eff_stnd +
                     saltcv_stnd + PARmin_stnd + DOmean_stnd + surftempcv_stnd + surftempmean_stnd + surftempmax_stnd ,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")
# remove freshwater auc=95, tjur = 0.258
# remove cul effect auc = 0.948, tjur 0.225
#remove slope auc =0.95, tjur 0.261, and relim <0.09 so remove
# remove surf temp mean auc =0.95, tjur = 0.23 # if remove temp cv auc = 0.95 tjur 0.23
# salt cv slightly higher over sal min
# do min and mean the same
# removed spline from rei auc 0.949 tjur 0.242
# adding NO3 auc 0.952 tjur 0.242
#remove tidal 0.95 tjur 0.242
# remove domin auc 0.95 tjur 0.246
#remove freshwater auc 0.95 tjur 0.243 (so tjur goes slightly down)
# remove cul eff auc 0.95 tjur 0.233 is final, use this model
m_s_5 <- sdmTMB_cv(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd + tempmin_stnd +
                                         saltcv_stnd + PARmin_stnd + surftempmax_stnd + NO3_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")

# adding in spatial field now makes it go down! auc 0.926 tjur 0.199
m_s_5a <- sdmTMB_cv(formula = presence  ~ depth_stnd + substrate + rei_sqrt_stnd + tempmin_stnd +
                      cul_eff_stnd + 
                      saltcv_stnd + PARmin_stnd + surftempmax_stnd + NO3_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = TRUE, 
                   data = data, 
                   fold_ids = "fold")
eval_cv <- evalStats( folds=1:numFolds,
                      m=m_s_5,
                      CV=cv_list_seagrass$cv)

# model from stepwise, auc = 0.928, tjur =0.23
m_s_6 <- sdmTMB_cv(formula = presence ~ substrate + depth_stnd  + s(rei_sqrt_stnd, k = 3) + tidal_sqrt_stnd + freshwater_sqrt_stnd + 
                     slope_stnd + NO3_stnd + saltmin_stnd + PARmin_stnd + 
                     surftempmax_stnd + surftempcv_stnd + 
                     tempdiff_stnd + cul_eff_stnd,
                   mesh = barrier_mesh, 
                   family = binomial(link = "logit"), 
                   spatial = FALSE, 
                   data = data, 
                   fold_ids = "fold")



roc <- pROC::roc(m_s_3$data$presence, plogis(m_s_3$data$cv_predicted))
auc <- pROC::auc(roc)
auc


# If a model is unbiased bias should be close to zero
# MAE want low values
# AUC want high values above 0.9. According to Pearce and Ferrier (2000) and Jones et al. (2010) values of AUC greater than 0.9 are considered good, between 0.7 and 0.9 moderate, and less than 0.7 poor. values of 0.5 indicate that the model is no better than random.
# TSS balances sensitivity (proportion of presence observations that are correctly classified) and specificity (proportion of absence observations that are correctly classified) and is independent of the prevalence of the observations (Allouche et al. 2006). Values of TSS greater than 0.6 are considered good, between 0.2 and 0.6 moderate, and less than 0.2 poor (Jones et al. 2010; Landis and Koch 1977). 
# Accuracy is the percent of predictions which are correctly classified and varies from values of 0 to 1 where 1 is the highest accuracy.
# Kappa is a measure of agreement between observed and predicted values that accounts for chance agreements and is dependent on prevalence of the observations. Kappa range from -1 to 1 with values less than 0 representing models that are no better than random and values of 1 indicating perfect agreement (Allouche et al. 2006). 

#fit full model
fmodel_s_bccm <- sdmTMB(formula = presence ~ depth_stnd + substrate + rei_sqrt_stnd + tempmin_stnd +
                   saltcv_stnd + PARmin_stnd + surftempmax_stnd + NO3_stnd,
                     mesh = barrier_mesh, 
                     family = binomial(link = "logit"), 
                     spatial = FALSE, 
                     data = data)
# look ar marginal effects plots, note temperature variables are correlated, so these are not correct 
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "depth_stnd[-3:6]") %>% plot() # found highest at shallow
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "rei_sqrt_stnd[-1:9]") %>% plot()  # increases with exposure
#ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "tidal_sqrt_stnd[-3:8]") %>% plot() # increases with tidal infleunce
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "tempmin_stnd[-5:3]") %>% plot() # increases with min temp
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "substrate") %>% plot() # found on rock
#ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "slope_stnd[-2:9]") %>% plot() # decreases with higher slopes
#ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "freshwater_sqrt_stnd[-1:18]") %>% plot() # increases with freshwater
#ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "cul_eff_stnd[-2:5]") %>% plot() # increases with cul eff
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "saltcv_stnd[-1:57]") %>% plot() # highest in areas with stable salinity
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "PARmin_stnd[-3:2]") %>% plot() # increases with min par
#ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "DOmin_stnd[-7:2]") %>% plot() #increases with do
#ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "surftempcv_stnd[-3:8]") %>% plot() # increases with increased variation
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "surftempmax_stnd[-4:5]") %>% plot() # declines with increased max temps
#ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "surftempmean_stnd[-4:4]") %>% plot() # increases with mean temps
ggeffects::ggeffect(model = fmodel_s_bccm,  terms = "NO3_stnd[-3:42]") %>% plot() # as NO3 increases presence decreases

visreg::visreg(fmodel_s_bccm, "depth_stnd")
visreg::visreg(fmodel_s_bccm, "rei_sqrt_stnd")
#visreg::visreg(fmodel_s_bccm, "tidal_sqrt_stnd")
visreg::visreg(fmodel_s_bccm, "saltcv_stnd")
visreg::visreg(fmodel_s_bccm, "PARmin_stnd")
visreg::visreg(fmodel_s_bccm, "tempmin_stnd")
#visreg::visreg(fmodel_s_bccm, "surftempcv_stnd")
visreg::visreg(fmodel_s_bccm, "surftempmax_stnd")
visreg::visreg(fmodel_s_bccm, "NO3_stnd")
# all of these look like linear relationships. There are some outliars of high NO3 and high salt cv (does this matter??)

# Model check
tidy(fmodel_s_bccm, conf.int = TRUE)
sanity(fmodel_s_bccm)

# Add fitted values (preds) to data
data$fitted_vals <- predict(fmodel_s_bccm, type="response")$est

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
prednames <- c("depth_stnd", "rei_sqrt_stnd", "substrate", "tempmin_stnd",
               "saltcv_stnd", "PARmin_stnd", "surftempmax_stnd", "NO3_stnd")


relimp <- varImp( model=fmodel_s_bccm,
                  dat=data,
                  preds=prednames,
                  permute=10 ) # Number of permutations

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
ret<- simulate(fmodel_s_bccm, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_s_bccm, return_DHARMa = TRUE)
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
