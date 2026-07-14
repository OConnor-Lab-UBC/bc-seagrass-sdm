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
UsePackages(c("sdmTMB", "sdmTMBextra", "ModelMetrics", "tidyverse", "sf", "future", "terra", "purrr", "future.apply", "rsample"))

####load model input ####
load("code/output_data/seagrass_model_inputs.RData")
seagrass_data_long <- seagrass_data_long %>% select(-saltmean_bccm_sq_stnd, -saltmean_nep_sq_stnd, -slope_sqrt_stnd, -saltmin_bccm_sq_stnd, -saltmin_nep_sq_stnd)
seagrass_data_long <- seagrass_data_long %>%
  mutate(Survey = as.factor(substr(HKey, 1, 3)),
         HKey = as.factor(HKey),
         Year_factor = as.factor(Year))
#  categorical predictors in environmental layers
facVars <- c("substrate", "Survey")

#having trouble with MSE survey causing b_j standard error warning as it was only absences for sea grass so combined it with BHm because people doing data would be similar and there would be care in iding correctly
#Also GSU and ABL was missing from a number of folds so combined with RSU as these are all set in index sites
seagrass_data_long <- seagrass_data_long %>%
  mutate(Survey = recode(Survey, MSE = "BHM")) %>%
  mutate(Survey = recode(Survey, GSU = "RSU")) %>%
  mutate(Survey = recode(Survey, ABL = "RSU")) %>%
  mutate(Survey = factor(Survey))


####Eelgrass model####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass, inner_fold = innerfold_eelgrass)
data <- data %>% mutate(Survey = as.factor(Survey)) # was still not recognizing i changed survye to factor so did it again
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))
# eelgrass in 3.7% of 20 m aggregated observations


####test variable correlation####
# data_env_bccm<- data[, !grepl("nep", names(data), ignore.case = TRUE)] %>%  dplyr::select(7:49)
# #enames <- names(data_env)
# corrplot::corrplot(cor(data_env_bccm, use = "pairwise.complete.obs"), method = "color",  col = colorRampPalette(c("red", "orange", "white", "blue", "purple"))(200), is.corr = TRUE, tl.cex = 0.6, tl.col = "black", number.cex = 0.5, order = "hclust", type = "upper")
# # Correlations close to-1 or +1 might indicate the existence of multicollinearity. one might suspect multicollinearity when the correlation between two (predictor) variables is below -0.9 or above +0.9.
# data_env_nep<- data[, !grepl("bccm", names(data), ignore.case = TRUE)] %>%  dplyr::select(7:49)
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
# eelgrass_ffs <- glm_ffs(data)  
# save(eelgrass_ffs, file = "code/output_data/model_results/eelgrass_ffs_variables.RData")
## most important variables in all 10 folds are substrate, depth, slope. Airtemp min and rei was important in 8 folds, DO_nep_stnd in 1 fold



####SDMtmb cv model####
#make mesh
mesh<- make_mesh(data = data, xy_cols = c("X", "Y"), cutoff = 55) 
#anything 30 and under makes the model underdispersed if no other variables (too complicated).
# anythign under 53 with other variables is underdispersed
#going to 10 km makes the model not run, likely will need to reduce fixed effects to make it work
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit cv model of spatial blocking 
plan(multisession)

# run nested spatial cross validation on two non spatial models to determine appropriate spline terms
depth_k_values <- list(
  "linear" = NULL,
  "2" = 2,
  "3" = 3,
  "4" = 4
)

airtemp_k_values <- list(
  "linear" = NULL,
  "2" = 2,
  "3" = 3,
  "4" = 4
)


res_bccm_nospatial <- nested_sdmTMB_cv(
  data = data,
  mesh = barrier_mesh,
  spatial_setting = FALSE,
  ocean_model = "bccm",
  outer_fold_col = "fold",
  inner_fold_col = "inner_fold"
)

save(res_bccm_nospatial, file = "code/output_data/model_results/eelgrass_eval_cv_bccm_nospatial.RData")
tss_thresh_bccm_nospatial<- mean(res_bccm_nospatial[["fold_results"]]$tss_threshold)
# threshold 0.04564801, linear on air and k = 4 on depth
# only way to do similar on spatial model is to refit the mesh for every inner fold. More justified to retain the smooth from the non spatial model
#then to add spatial field so you are isolating effect of spatial field

res_nep_nospatial <- nested_sdmTMB_cv(
  data = data,
  mesh = barrier_mesh,
  spatial_setting = FALSE,
  ocean_model = "nep",
  outer_fold_col = "fold",
  inner_fold_col = "inner_fold"
)
save(res_nep_nospatial, file = "code/output_data/model_results/eelgrass_eval_cv_nep_nospatial.RData")
tss_thresh_nep_nospatial<- mean(res_nep_nospatial[["fold_results"]]$tss_threshold)
# threshold 0.0417, linear on air and k = 4 on depth



m_e_1 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + prmin_stnd +#chelsa variables
                     saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + #bccm variables
                     Survey,  #fixed effect
                   mesh = barrier_mesh, priors = sdmTMBpriors(b = normal(0, 1)), family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

eval_cv_bccm_nospatial <- evalStats(folds = 1:numFolds,  m = m_e_1,  CV = cv_list_eelgrass$cv,  response_col = "presence")
eval_cv_bccm_nospatial

m_e_2 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + prmin_stnd +#chelsa variables
                     saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + #bccm variables
                     Survey,  #fixed effect
                   mesh = barrier_mesh, priors = sdmTMBpriors(b = normal(0, 1)), family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

eval_cv_bccm_spatial <- evalStats(folds = 1:numFolds,  m = m_e_2,  CV = cv_list_eelgrass$cv,  response_col = "presence")
eval_cv_bccm_spatial

m_e_3 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + prmin_stnd +#chelsa variables
                     saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + #nep variables
                     Survey,  #fixed effect
                   mesh = barrier_mesh, priors = sdmTMBpriors(b = normal(0, 1)), family = binomial(link = "logit"), spatial = FALSE, data = data, fold_ids = "fold")

eval_cv_nep_nospatial <- evalStats(folds = 1:numFolds,  m = m_e_3,  CV = cv_list_eelgrass$cv,  response_col = "presence")
eval_cv_nep_nospatial

m_e_4 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                     airtempmin_stnd + rsdsmin_stnd + prmin_stnd +#chelsa variables
                     saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + #nep variables
                     Survey,  #fixed effect
                   mesh = barrier_mesh, priors = sdmTMBpriors(b = normal(0, 1)), family = binomial(link = "logit"), spatial = TRUE, data = data, fold_ids = "fold")

eval_cv_nep_spatial <- evalStats(folds = 1:numFolds,  m = m_e_4,  CV = cv_list_eelgrass$cv,  response_col = "presence")
eval_cv_nep_spatial

eval_cv_list <- list(eval_cv_bccm_nospatial, eval_cv_bccm_spatial, eval_cv_nep_nospatial, eval_cv_nep_spatial)
save(eval_cv_list, file = "code/output_data/model_results/eelgrass_eval_cv.RData")


#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009
# model does better if you don't include oceanographic or atmospheric data. which makes sense, becasue there are no differences between past and present as geomorphic are  static
#cannot account for surveys or years in these models as random factor because only ABL, Cuke, GDK, GSU, and RSU surveys were conducted prior to 2013 and you cant make predictions of years not fit in model

data_pre2013 <- data %>%
  filter(Year < 2010) %>%
  mutate(row_id = row_number())
names(data_pre2013)

unique(data_pre2013$Survey)
test_set <- data %>% filter(Year > 2012)
obs_test <- test_set$presence

formula_bccm <- presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
  airtempmin_stnd + rsdsmin_stnd + prmin_stnd + 
  saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd
formula_nep <- presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
  airtempmin_stnd + rsdsmin_stnd + prmin_stnd + 
  saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd
cv_thr_bccm <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_bccm,
  coastline = coastline,
  spatial = FALSE
)
cv_thr_bccm_spatial <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_bccm,
  coastline = coastline,
  spatial = TRUE
)
cv_thr_nep <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_nep,
  coastline = coastline,
  spatial = FALSE
)
cv_thr_nep_spatial <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_nep,
  coastline = coastline,
  spatial = TRUE
)


# Create 10 random folds on pre-2010 data (for model fit)
set.seed(123)
folds <- rsample::vfold_cv(data_pre2013, v = 10, strata = "presence")

results_list <- list()

for (i in seq_along(folds$splits)) {
  split <- folds$splits[[i]]
  train_data <- analysis(split)
  obs_train <- train_data$presence
  
  mesh <- make_mesh(data = train_data, xy_cols = c("X", "Y"), cutoff = 55)
  barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = FALSE)
  
  m_bccm <- sdmTMB(
    formula = formula_bccm,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = FALSE,
    data = train_data
  )
  
  m_bccm_spatial <- sdmTMB(
    formula = formula_bccm,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = TRUE,
    data = train_data
  )
  
  m_nep <- sdmTMB(
    formula = formula_nep,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = FALSE,
    data = train_data
  )
  
  m_nep_spatial <- sdmTMB(
    formula = formula_nep,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = TRUE,
    data = train_data
  )
  
  pred_train_bccm <- plogis(predict(m_bccm, newdata = train_data)$est)
  pred_train_spatial_bccm <- plogis(predict(m_bccm_spatial, newdata = train_data)$est)
  pred_train_nep <- plogis(predict(m_nep, newdata = train_data)$est)
  pred_train_spatial_nep <- plogis(predict(m_nep_spatial, newdata = train_data)$est)
  
  pred_test_bccm <- plogis(predict(m_bccm, newdata = test_set)$est)
  pred_test_spatial_bccm <- plogis(predict(m_bccm_spatial, newdata = test_set)$est)
  pred_test_nep <- plogis(predict(m_nep, newdata = test_set)$est)
  pred_test_spatial_nep <- plogis(predict(m_nep_spatial, newdata = test_set)$est)
  
  eval_nonspatial_bccm <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_bccm,
    obs_test = obs_test,
    pred_test = pred_test_bccm,
    threshold = cv_thr_bccm$threshold
  )
  eval_nonspatial_bccm$model <- "BCCM_no_spatial"
  eval_nonspatial_bccm$fold <- i
  
  eval_spatial_bccm <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_spatial_bccm,
    obs_test = obs_test,
    pred_test = pred_test_spatial_bccm,
    threshold = cv_thr_bccm_spatial$threshold
  )
  eval_spatial_bccm$model <- "BCCM_spatial"
  eval_spatial_bccm$fold <- i
  
  eval_nonspatial_nep <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_nep,
    obs_test = obs_test,
    pred_test = pred_test_nep,
    threshold = cv_thr_nep$threshold
  )
  eval_nonspatial_nep$model <- "NEP_no_spatial"
  eval_nonspatial_nep$fold <- i
  
  eval_spatial_nep <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_spatial_nep,
    obs_test = obs_test,
    pred_test = pred_test_spatial_nep,
    threshold = cv_thr_nep_spatial$threshold
  )
  eval_spatial_nep$model <- "NEP_spatial"
  eval_spatial_nep$fold <- i
  
  results_list[[length(results_list) + 1]] <- eval_nonspatial_bccm
  results_list[[length(results_list) + 1]] <- eval_spatial_bccm
  results_list[[length(results_list) + 1]] <- eval_nonspatial_nep
  results_list[[length(results_list) + 1]] <- eval_spatial_nep
}
forecast_predict_eelgrass <- do.call(rbind, results_list)
forecast_only <- forecast_predict_eelgrass %>%
  filter(dataset == "forecast")
forecast_by_model <- forecast_only %>%
  group_by(model) %>%
  summarise(
    AUC = mean(AUC, na.rm = TRUE),
    TjurR2 = mean(TjurR2, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    Brier = mean(Brier, na.rm = TRUE),
    LogLoss = mean(LogLoss, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    TSS = mean(TSS, na.rm = TRUE),
    Threshold = mean(threshold, na.rm = TRUE),
    .groups = "drop"
  )
forecast_by_model

save(forecast_predict_eelgrass, forecast_by_model, file = "code/output_data/model_results/forecast_eelgrass_models.RData")


####SDMtmb full model####
# fit full model bccm 
fmodel_e_bccm_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                                    airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                    saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + #bccm variables
                                    Survey,  #fixed effect
                                mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = FALSE, data = data)

fmodel_e_bccm_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                                  airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                  saltcv_bccm_stnd + NH4_bccm_stnd + tempcv_bccm_stnd + #bccm variables
                                  Survey,  #fixed effect
                                mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = TRUE, data = data)

fmodel_e_nep_nospatial <- sdmTMB(formula = presence ~  s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                                   airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                   saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + #nep variables
                                   Survey,  #fixed effect
                                 mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = FALSE, data = data)

fmodel_e_nep_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 4) + substrate + slope_stnd + rei_stnd +  
                                 airtempmin_stnd + rsdsmin_stnd + prmin_stnd + #chelsa variables
                                 saltcv_nep_stnd + NH4_nep_stnd + tempcv_nep_stnd + #nep variables
                                 Survey,  #fixed effect
                               mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = TRUE, data = data)

save(data, fmodel_e_bccm_nospatial, fmodel_e_bccm_spatial, fmodel_e_nep_nospatial, fmodel_e_nep_spatial, file = "code/output_data/model_results/final_eelgrass_model.RData")



#have a look at marginal effects
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "depth_stnd[all]") %>% plot() # decreases with depth
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "depth_stnd[all]") %>% plot() # decreases with depth
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "substrate") %>% plot() # more in sand and mud
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "substrate") %>% plot() # more in sand and mud
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "slope_stnd[all]") %>% plot() # presence declines with slope
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "slope_stnd[all]") %>% plot() # presence declines with slope
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "rei_stnd[all]") %>% plot() #presence declines with exposure
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "rei_stnd[all]") %>% plot() #presence declines with exposure
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "airtempmin_stnd[all]") %>% plot() # presence increases as air temp min increases
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "airtempmin_stnd[all]") %>% plot() # presence increases as air temp min increases
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "rsdscv_stnd[all]") %>% plot() # presence decreases as rsds cv increases
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "rsdscv_stnd[all]") %>% plot() # presence decreases as rsds cv increases
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "saltcv_bccm_stnd[all]") %>% plot() # as salinity variability increases presence goes down
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "saltcv_nep_stnd[all]") %>% plot() # as salinity variability increases presence goes down
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "NH4_bccm_stnd[all]") %>% plot() # as NH4 increases presence goes up
# 
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "tempcv_nep_stnd[all]") %>% plot() # as tempcv increases presence goes up
# 
# ggeffects::ggpredict(model = fmodel_e_bccm_nospatial,  terms = "Survey") %>% plot() # # Remote sensed adds survey effect 
# ggeffects::ggpredict(model = fmodel_e_nep_nospatial,  terms = "Survey") %>% plot() # # Remote sensed adds survey effect 

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
  data[[paste0("pred_TSS_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1
  )
  
  data[[paste0("pred_kappa_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxKappa"], 0, 1
  )
  
  data[[paste0("pred_PCC_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxPCC"], 0, 1
  )
  
  # Evaluate model
  eval_metrics <- evalfmod(x = data_tmp, thresh = thresh, sp = m_name)
  
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
# brier score: 0 is perfect, 1 is bad
#log loss: 0 is perfect, higher is bad
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
prednames_bccm <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "Survey", "rsdsmin_stnd", "airtempmin_stnd", "prmin_stnd", "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempcv_bccm_stnd")
prednames_nep <- c("depth_stnd", "substrate", "rei_stnd", "slope_stnd", "Survey", "rsdsmin_stnd", "airtempmin_stnd", "prmin_stnd", "saltcv_nep_stnd", "NH4_nep_stnd", "tempcv_nep_stnd")




plan(sequential) # the spatial models don't run well on multisession

relimp_e_bccm_nospatial <- varImp_sdmTMB(
  model = fmodel_e_bccm_nospatial,
  dat = data,
  preds = prednames_bccm,
  permute = 10
)


save(relimp_e_bccm_nospatial, file = "code/output_data/model_results/relimp_e_bccm_nospatial.RData")

relimp_e_bccm_spatial <- varImp_sdmTMB( model=fmodel_e_bccm_spatial,
                                   dat=data,
                                   preds=prednames_bccm,
                                   permute=10 ) # Number of permutations
save(relimp_e_bccm_spatial, file = "code/output_data/model_results/relimp_e_bccm_spatial.RData")

relimp_e_nep_nospatial <- varImp_sdmTMB( model=fmodel_e_nep_nospatial,
                         dat=data,
                         preds=prednames_nep,
                         permute=10 ) # Number of permutations
save(relimp_e_nep_nospatial, file = "code/output_data/model_results/relimp_e_nep_nospatial.RData")


relimp_e_nep_spatial <- varImp_sdmTMB( model=fmodel_e_nep_spatial,
                                  dat=data,
                                  preds=prednames_nep,
                                  permute=10 ) # Number of permutations
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









# have not updated surfgrass models yet!!!


#### Surfgrass model ####

sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass, inner_fold = innerfold_seagrass)
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
mesh<- make_mesh(data = data, xy_cols = c("X", "Y"), cutoff = 55) 
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)

#fit model
plan(multisession)

# run nested spatial cross validation on two non spatial models to determine appropriate spline terms
depth_k_values <- list(
  "linear" = NULL,
  "2" = 2,
  "3" = 3,
  "4" = 4
)

rei_sqrt_k_values <- list(
  "linear" = NULL,
  "2" = 2,
  "3" = 3,
  "4" = 4
)


res_bccm_nospatial <- nested_sdmTMB_cv_surf(
  data = data,
  mesh = barrier_mesh,
  spatial_setting = FALSE,
  ocean_model = "bccm",
  outer_fold_col = "fold",
  inner_fold_col = "inner_fold"
)

save(res_bccm_nospatial, file = "code/output_data/model_results/surfgrass_eval_cv_bccm_nospatial.RData")
tss_thresh_bccm_nospatial<- mean(res_bccm_nospatial[["fold_results"]]$tss_threshold)
# threshold 0.0232, k=4 on rei (7 were 4 and 3 were 2) and k = 4 on depth (9 were 4 and 1 was 2)

res_nep_nospatial <- nested_sdmTMB_cv_surf(
  data = data,
  mesh = barrier_mesh,
  spatial_setting = FALSE,
  ocean_model = "nep",
  outer_fold_col = "fold",
  inner_fold_col = "inner_fold"
)
save(res_nep_nospatial, file = "code/output_data/model_results/surfgrass_eval_cv_nep_nospatial.RData")
tss_thresh_nep_nospatial<- mean(res_nep_nospatial[["fold_results"]]$tss_threshold)
# threshold 0.019, linear on rei (7 linear, 2 were 2 and 1 was 4) and k = 4 on depth (9 were 4 and 1 was linear)

# having issues with overdispersion on keeping rei k = 4 for the bccm models, 
#to keep it similar between the two we will make rei k=2 for both bccm and nep models to split the difference and account for too high mdoel complexity

# NO3 was significant but after reviewing variable between SSC and BCCM there was no good distribution between borders

# no spatial 
m_s_1 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd + 
                     airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                     saltmean_bccm_stnd + tempmean_bccm_stnd + surftempcv_bccm_stnd + #bccm variables
                     Survey,  #fixed effect
                    mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = FALSE, data = data, fold_ids = "fold")
eval_cv_bccm_nospatial <- evalStats( folds=1:numFolds,m=m_s_1,CV=cv_list_seagrass$cv,  response_col = "presence")
eval_cv_bccm_nospatial

# bccm model with spatial
m_s_2 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd + 
                     airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                     saltmean_bccm_stnd + tempmean_bccm_stnd + surftempcv_bccm_stnd + #bccm variables
                     Survey,  #fixed effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = TRUE, data = data, fold_ids = "fold")
eval_cv_bccm_spatial <- evalStats( folds=1:numFolds,m=m_s_2,CV=cv_list_seagrass$cv,  response_col = "presence")
eval_cv_bccm_spatial

# nep model no spatial
m_s_3 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd +  
                     airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                     saltmean_nep_stnd + tempmean_nep_stnd + surftempcv_nep_stnd + #nep variables
                     Survey,  #fixed effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = FALSE, data = data, fold_ids = "fold")
eval_cv_nep_nospatial <- evalStats( folds=1:numFolds,m=m_s_3,CV=cv_list_seagrass$cv,  response_col = "presence")
eval_cv_nep_nospatial

# nep model spatial 
m_s_4 <- sdmTMB_cv(formula = presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd +  
                     airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                     saltmean_nep_stnd + tempmean_nep_stnd + surftempcv_nep_stnd + #nep variables
                     Survey,  #fixed effect
                   mesh = barrier_mesh, family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)), spatial = TRUE, data = data, fold_ids = "fold")
eval_cv_nep_spatial <- evalStats( folds=1:numFolds,m=m_s_4,CV=cv_list_seagrass$cv,  response_col = "presence")
eval_cv_nep_spatial

# cv stats from all best models and save
eval_cv_list_surfgrass <- list(eval_cv_bccm_nospatial, eval_cv_bccm_spatial, eval_cv_nep_nospatial, eval_cv_nep_spatial)
save(eval_cv_list_surfgrass, file = "code/output_data/model_results/eval_cv_surfgrass.RData")

#HAVING ISSUES FAILING HERE!
#### test forecasting
# left a few years gap 2010-2012 #trained model with 1993-2009
# model does better if you don't include oceanographic or atmospheric data. which makes sense, becasue there are no differences between past and present as geomorphic are  static
#cannot account for surveys or years in these models as random factor because only ABL, Cuke, GDK, GSU, and RSU surveys were conducted prior to 2013 and you cant make predictions of years not fit in model

data_pre2013 <- data %>%
  filter(Year < 2010) %>%
  mutate(row_id = row_number())
names(data_pre2013)

unique(data_pre2013$Survey)
test_set <- data %>% filter(Year > 2012)
obs_test <- test_set$presence

formula_bccm <- presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd + 
  airtempcv_stnd + prmean_stnd + rsdsmin_stnd +
  saltmean_bccm_stnd + tempmean_bccm_stnd + surftempcv_bccm_stnd
formula_nep <- presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd +  
  airtempcv_stnd + prmean_stnd + rsdsmin_stnd +
  saltmean_nep_stnd + tempmean_nep_stnd + surftempcv_nep_stnd

cv_thr_bccm <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_bccm,
  coastline = coastline,
  spatial = FALSE
)
cv_thr_bccm_spatial <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_bccm,
  coastline = coastline,
  spatial = TRUE
)
cv_thr_nep <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_nep,
  coastline = coastline,
  spatial = FALSE
)
cv_thr_nep_spatial <- get_global_cv_threshold_sdmTMB(
  data_train = data_pre2013,
  formula = formula_nep,
  coastline = coastline,
  spatial = TRUE
)


# Create 10 random folds on pre-2010 data (for model fit)
set.seed(123)
folds <- rsample::vfold_cv(data_pre2013, v = 10, strata = "presence")

results_list <- list()

for (i in seq_along(folds$splits)) {
  split <- folds$splits[[i]]
  train_data <- analysis(split)
  obs_train <- train_data$presence
  
  mesh <- make_mesh(data = train_data, xy_cols = c("X", "Y"), cutoff = 55)
  barrier_mesh <- add_barrier_mesh(mesh, barrier_sf = coastline, proj_scaling = 1000, plot = FALSE)
  
  m_bccm <- sdmTMB(
    formula = formula_bccm,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = FALSE,
    data = train_data
  )
  
  m_bccm_spatial <- sdmTMB(
    formula = formula_bccm,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = TRUE,
    data = train_data
  )
  
  m_nep <- sdmTMB(
    formula = formula_nep,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = FALSE,
    data = train_data
  )
  
  m_nep_spatial <- sdmTMB(
    formula = formula_nep,
    mesh = barrier_mesh,
    priors = sdmTMBpriors(b = normal(0, 1)),
    family = binomial(link = "logit"),
    spatial = TRUE,
    data = train_data
  )
  
  pred_train_bccm <- plogis(predict(m_bccm, newdata = train_data)$est)
  pred_train_spatial_bccm <- plogis(predict(m_bccm_spatial, newdata = train_data)$est)
  pred_train_nep <- plogis(predict(m_nep, newdata = train_data)$est)
  pred_train_spatial_nep <- plogis(predict(m_nep_spatial, newdata = train_data)$est)
  
  pred_test_bccm <- plogis(predict(m_bccm, newdata = test_set)$est)
  pred_test_spatial_bccm <- plogis(predict(m_bccm_spatial, newdata = test_set)$est)
  pred_test_nep <- plogis(predict(m_nep, newdata = test_set)$est)
  pred_test_spatial_nep <- plogis(predict(m_nep_spatial, newdata = test_set)$est)
  
  eval_nonspatial_bccm <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_bccm,
    obs_test = obs_test,
    pred_test = pred_test_bccm,
    threshold = cv_thr_bccm$threshold
  )
  eval_nonspatial_bccm$model <- "BCCM_no_spatial"
  eval_nonspatial_bccm$fold <- i
  
  eval_spatial_bccm <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_spatial_bccm,
    obs_test = obs_test,
    pred_test = pred_test_spatial_bccm,
    threshold = cv_thr_bccm_spatial$threshold
  )
  eval_spatial_bccm$model <- "BCCM_spatial"
  eval_spatial_bccm$fold <- i
  
  eval_nonspatial_nep <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_nep,
    obs_test = obs_test,
    pred_test = pred_test_nep,
    threshold = cv_thr_nep$threshold
  )
  eval_nonspatial_nep$model <- "NEP_no_spatial"
  eval_nonspatial_nep$fold <- i
  
  eval_spatial_nep <- evaluate_forecast(
    obs_train = obs_train,
    pred_train = pred_train_spatial_nep,
    obs_test = obs_test,
    pred_test = pred_test_spatial_nep,
    threshold = cv_thr_nep_spatial$threshold
  )
  eval_spatial_nep$model <- "NEP_spatial"
  eval_spatial_nep$fold <- i
  
  results_list[[length(results_list) + 1]] <- eval_nonspatial_bccm
  results_list[[length(results_list) + 1]] <- eval_spatial_bccm
  results_list[[length(results_list) + 1]] <- eval_nonspatial_nep
  results_list[[length(results_list) + 1]] <- eval_spatial_nep
}
forecast_predict_surfgrass <- do.call(rbind, results_list)
forecast_only <- forecast_predict_surfgrass %>%
  filter(dataset == "forecast")
forecast_by_model <- forecast_only %>%
  group_by(model) %>%
  summarise(
    AUC = mean(AUC, na.rm = TRUE),
    TjurR2 = mean(TjurR2, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    Brier = mean(Brier, na.rm = TRUE),
    LogLoss = mean(LogLoss, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    TSS = mean(TSS, na.rm = TRUE),
    Threshold = mean(threshold, na.rm = TRUE),
    .groups = "drop"
  )
forecast_by_model

save(forecast_predict_surfgrass, forecast_by_model, file = "code/output_data/model_results/forecast_surfgrass_models.RData")


#fit full model
fmodel_s_bccm_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd + 
                                    airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                                    saltmean_bccm_stnd + tempmean_bccm_stnd + surftempcv_bccm_stnd + #bccm variables
                                    Survey,  #fixed effect
                                mesh = barrier_mesh, #random effect
                                family = binomial(link = "logit"), priors = sdmTMBpriors(b = normal(0, 1)),
                                spatial = FALSE, 
                                data = data)
fmodel_s_bccm_spatial <- sdmTMB(formula =  presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd +  
                                  airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                                  saltmean_bccm_stnd + tempmean_bccm_stnd + surftempcv_bccm_stnd + #bccm variables
                                  Survey,  #fixed effect
                                mesh = barrier_mesh, priors = sdmTMBpriors(b = normal(0, 1)),
                                family = binomial(link = "logit"), 
                                spatial = TRUE, 
                                data = data)

fmodel_s_nep_nospatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd +  
                                   airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                                   saltmean_nep_stnd + tempmean_nep_stnd + surftempcv_nep_stnd + #nep variables
                                   Survey,  #fixed effect
                                  mesh = barrier_mesh, priors = sdmTMBpriors(b = normal(0, 1)),
                                  family = binomial(link = "logit"), 
                                  spatial = FALSE, 
                                  data = data)

fmodel_s_nep_spatial <- sdmTMB(formula = presence ~ s(depth_stnd, k = 4) + s(rei_sqrt_stnd, k = 2) + substrate + tidal_sqrt_stnd +  
                                 airtempcv_stnd + prmean_stnd + rsdsmin_stnd + #chelsa variables
                                 saltmean_nep_stnd + tempmean_nep_stnd + surftempcv_nep_stnd + #nep variables
                                 Survey,  #fixed effect
                                 mesh = barrier_mesh, priors = sdmTMBpriors(b = normal(0, 1)),
                                 family = binomial(link = "logit"), 
                                 spatial = TRUE, 
                                 data = data)

save(data, fmodel_s_bccm_nospatial, fmodel_s_bccm_spatial, fmodel_s_nep_nospatial, fmodel_s_nep_spatial, file = "code/output_data/model_results/final_surfgrass_model.RData")


#have a look at marginal effects
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "depth_stnd[all]") %>% plot() # unimodal depth distribution
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "depth_stnd[all]") %>% plot() # unimodal depth distribution
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "substrate") %>% plot() # more rock
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "substrate") %>% plot() # more rock
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "rei_sqrt_stnd[all]") %>% plot() #unimodal distribution
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "rei_sqrt_stnd[all]") %>% plot() #unimodal distribution
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "tidal_sqrt_stnd[all]") %>% plot() # presence increases with tidal current
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "tidal_sqrt_stnd[all]") %>% plot() 
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "cul_eff_stnd[all]") %>% plot() #presence higher in areas with impacts??
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "cul_eff_stnd[all]") %>% plot() #presence higher in areas with impacts??
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "airtempcv_stnd[all]") %>% plot() # presence decreases as air temp variability increases
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "airtempcv_stnd[all]") %>% plot() # presence decreases as air temp variability increases
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "rsdsmin_stnd[all]") %>% plot() # presence increases as rsds min increases
#  
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "prmean_stnd[all]") %>% plot() # presence increases as precipitation increases
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "prmean_stnd[all]") %>% plot() 
#  
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "saltcv_bccm_stnd[all]") %>% plot() # preence decreases as salt variability increases
# 
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "saltmean_nep_stnd[all]") %>% plot() # presence increases as salt mean icnrease
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "NO3_bccm_stnd[all]") %>% plot() # presence decreases as NO3 decreases
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "tempmin_bccm_stnd[all]") %>% plot() #presence increases as temp min increase
# 
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "tempmean_nep_stnd[all]") %>% plot() # presence increases as temp mean increases
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "surftempcv_bccm_stnd[all]") %>% plot() # presence decreases as surf temp variability increases
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "surftempcv_nep_stnd[all]") %>% plot() 
# 
# ggeffects::ggpredict(model = fmodel_s_bccm_nospatial,  terms = "Survey") %>% plot() # doesn't look like difference between surveys even though it has a lot of influence
# ggeffects::ggpredict(model = fmodel_s_nep_nospatial,  terms = "Survey") %>% plot() 



# Model check
tidy(fmodel_s_bccm_nospatial, conf.int = TRUE)
sanity(fmodel_s_bccm_nospatial)
tidy(fmodel_s_bccm_spatial, conf.int = TRUE)
sanity(fmodel_s_bccm_spatial)
tidy(fmodel_s_bccm_spatial, "ran_pars", conf.int = TRUE)
tidy(fmodel_s_nep_nospatial, conf.int = TRUE)
sanity(fmodel_s_nep_nospatial)
tidy(fmodel_s_nep_spatial, conf.int = TRUE)
sanity(fmodel_s_nep_spatial)

models <- list(
  bccm_nospatial = fmodel_s_bccm_nospatial,
  bccm_spatial   = fmodel_s_bccm_spatial,
  nep_nospatial  = fmodel_s_nep_nospatial,
  nep_spatial    = fmodel_s_nep_spatial
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
  data[[paste0("pred_TSS_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxSens+Spec"], 0, 1
  )
  
  data[[paste0("pred_kappa_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxKappa"], 0, 1
  )
  
  data[[paste0("pred_PCC_thresh_", m_name)]] <- ifelse(
    data_tmp$fitted_vals < thresh$Predicted[thresh$Method == "MaxPCC"], 0, 1
  )
  
  # Evaluate model
  eval_metrics <- evalfmod(x = data_tmp, thresh = thresh, sp = m_name)
  
  # Save both evaluation metrics AND thresholds in the list
  eval_results[[m_name]] <- list(
    metrics = eval_metrics,
    thresholds = thresh
  )
}

save(eval_results, file = "code/output_data/model_results/surfgrass_eval_final_models.RData")


# this model has good TSS, is well calibrated (miller). for calibration (Hosmer & Lemeshow goodness-of-fit) model seems to have issues at higher predicted probabilities

#get relative importance
plan(sequential) # the spatial models don't run well on multisession

prednames_bccm <- c("depth_stnd", "substrate", "rei_sqrt_stnd", "tidal_sqrt_stnd", "Survey", "airtempcv_stnd", 
                    "prmean_stnd", "rsdsmin_stnd", "saltmean_bccm_stnd", "tempmean_bccm_stnd", "surftempcv_bccm_stnd")
prednames_nep <- c("depth_stnd", "substrate", "rei_sqrt_stnd", "tidal_sqrt_stnd", "Survey", "airtempcv_stnd", 
                   "prmean_stnd", "rsdsmin_stnd", "saltmean_nep_stnd", "tempmean_nep_stnd", "surftempcv_nep_stnd")

plan(sequential) # the spatial models don't run well on multisession

relimp_s_bccm_nospatial <- varImp_sdmTMB(model = fmodel_s_bccm_nospatial,
                                         dat = data,
                                         preds = prednames_bccm,
                                         permute = 10)
save(relimp_s_bccm_nospatial, file = "code/output_data/model_results/relimp_s_bccm_nospatial.RData")

relimp_s_bccm_spatial <- varImp_sdmTMB( model=fmodel_s_bccm_spatial,
                                        dat=data,
                                        preds=prednames_bccm,
                                        permute=10 ) # Number of permutations
save(relimp_s_bccm_spatial, file = "code/output_data/model_results/relimp_s_bccm_spatial.RData")

relimp_s_nep_nospatial <- varImp_sdmTMB( model=fmodel_s_nep_nospatial,
                                         dat=data,
                                         preds=prednames_nep,
                                         permute=10 ) # Number of permutations
save(relimp_s_nep_nospatial, file = "code/output_data/model_results/relimp_s_nep_nospatial.RData")


relimp_s_nep_spatial <- varImp_sdmTMB( model=fmodel_s_nep_spatial,
                                       dat=data,
                                       preds=prednames_nep,
                                       permute=10 ) # Number of permutations
save(relimp_s_nep_spatial, file = "code/output_data/model_results/relimp_s_nep_spatial.RData")


####check residuals####
# MCMC based randomized quantile residuals (takes a while to compute)
set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_s_bccm_nospatial, mcmc_iter = 1600, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_s_bccm_nospatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_s_bccm_nospatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_s_bccm_nospatial, return_DHARMa = TRUE)
#plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_s_bccm_nospatial.RData")


set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_s_bccm_spatial, mcmc_iter = 1600, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_s_bccm_spatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_s_bccm_spatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_s_bccm_spatial, return_DHARMa = TRUE)
#plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_s_bccm_spatial.RData")

set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_s_nep_nospatial, mcmc_iter = 1600, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_s_nep_nospatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_s_nep_nospatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_s_nep_nospatial, return_DHARMa = TRUE)
#plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_s_nep_nospatial.RData")


set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(fmodel_s_nep_spatial, mcmc_iter = 1600, mcmc_warmup = 400)
mcmc_res <- residuals(fmodel_s_nep_spatial, type = "mle-mcmc", mcmc_samples = samps)
qqnorm(mcmc_res)
abline(0, 1)

# simulation-based randomized quantile residuals
set.seed(123)
ret<- simulate(fmodel_s_nep_spatial, nsim = 500, type = "mle-mvn") 
r_ret <-  dharma_residuals(ret, fmodel_s_nep_spatial, return_DHARMa = TRUE)
#plot(r_ret)
DHARMa::testResiduals(r_ret)

save(samps, mcmc_res, ret, r_ret, file = "code/output_data/model_results/residuals_s_nep_spatial.RData")





















# 
# ####Eelgrass delta model with percent cover ####
# dat2 <- subset(data, mean_PerCovZO > 0)
# #not all surveys had records of percent cover
# 
# # recomend to use beta family which has to have percent cover rescaled to 0 to 1 with no 100% percent cover 
# 
# dat2$cover_beta <- (dat2$mean_PerCovZO - 0.5) / 100
# 
# range(dat2$cover_beta)
# 
# dat2$Survey <- factor(dat2$Survey,
#                       levels = c("ABL", "BHM", "Cuk", "GSU", "Mul", "RSU"))
# mesh2 <- make_mesh(dat2,
#                    xy_cols = c("X", "Y"),
#                    mesh = mesh$mesh
# )
# plot(mesh2)
# barrier_mesh2 <- add_barrier_mesh(mesh2, barrier_sf = coastline, proj_scaling = 1000, plot = TRUE)
# #better to have more flexible spline on depth
# #better with spatial
# m_e_per_1 <- sdmTMB_cv(formula = cover_beta ~ s(depth_stnd, k=5) + substrate + freshwater_sqrt_stnd +
#                          surftempmin_bccm_stnd + PARmean_bccm_stnd+ #bccm variables
#                          (1|Survey),  #random effect
#                        mesh = barrier_mesh2, 
#                        family = Beta(link = "logit"), 
#                        spatial = TRUE, 
#                        fold_ids = "fold",
#                        data = dat2)
# m_e_per_1$sum_loglik # -9038
# 
# m_e_per_2 <- sdmTMB_cv(formula = cover_beta ~ s(depth_stnd, k=5) + substrate + freshwater_sqrt_stnd +
#                          surftempmin_nep_stnd + PARmean_nep_stnd + #nep variables
#                          (1|Survey),  #random effect
#                        mesh = barrier_mesh2, 
#                        family = Beta(link = "logit"),
#                        spatial = TRUE, 
#                        fold_ids = "fold",
#                        data = dat2)
# m_e_per_2$sum_loglik # -9048
# 
# #So best models that include all best factors from both nep and bccm
# # best model has surftemp min, parmean, freshwater sqrt, substrate and depth but model was underdispersed!
# m_e_per_bccm_final <- sdmTMB(formula = cover_beta ~ s(depth_stnd, k=5), #+ substrate + #freshwater_sqrt_stnd +
#                                #surftempmin_bccm_stnd + PARmean_bccm_stnd, # + #+ #bccm variables
#                                #(1|Survey),  #random effect
#                        mesh = barrier_mesh2, 
#                        family = Beta(link = "logit"), 
#                        spatial = FALSE, 
#                        data = dat2)
# # Gamma and lognormal family don't work well
# 
# 
# m_e_per_nep_final <- sdmTMB(formula = mean_PerCovZO ~ s(depth_stnd, k=5) + substrate + freshwater_sqrt_stnd +
#                          surftempmin_nep_stnd + PARmean_nep_stnd +  #nep variables
#                          (1|Survey),  #random effect
#                        mesh = barrier_mesh2, 
#                        family = Gamma(link = "log"), 
#                        spatial = TRUE, 
#                        data = dat2)
# 
# sanity(m_e_per_bccm_final)
# sanity(m_e_per_nep_final)
# #print(fmodel_e_delta_bccm_nospatial)
# 
# tidy(m_e_per_bccm_final)
# tidy(m_e_per_bccm_final, "ran_pars", conf.int = TRUE)
# tidy(m_e_per_nep_final)
# tidy(m_e_per_nep_final, "ran_pars", conf.int = TRUE)
# 
# visreg::visreg(m_e_per_bccm_final, "depth_stnd")
# visreg::visreg(m_e_per_nep_final, "depth_stnd")
# 
# visreg::visreg(m_e_per_bccm_final, "substrate")
# visreg::visreg(m_e_per_nep_final, "substrate")
# 
# visreg::visreg(m_e_per_bccm_final, "freshwater_sqrt_stnd")
# visreg::visreg(m_e_per_nep_final, "freshwater_sqrt_stnd")
# 
# visreg::visreg(m_e_per_bccm_final, "surftempmin_bccm_stnd")
# visreg::visreg(m_e_per_nep_final, "surftempmin_nep_stnd")
# 
# visreg::visreg(m_e_per_bccm_final, "PARmean_bccm_stnd")
# visreg::visreg(m_e_per_nep_final, "PARmean_nep_stnd")
# 
# visreg::visreg(m_e_per_bccm_final, "Survey")
# visreg::visreg(m_e_per_nep_final, "Survey")
# 
# library(DHARMa)
# sim <- simulateResiduals(
#   fittedModel = m_e_per_bccm_final,
#   n = 250,
#   refit = FALSE
# )
# fitted_vals <- fitted(m_e_per_bccm_final)
# plot(sim)
# summary(sim$scaledResiduals)
# plotResiduals(sim, fitted_vals)
# plotResiduals(sim, dat2$depth)
# 
# set.seed(123)
# rq_res <- residuals(m_e_per_bccm_final, type = "mle-mvn")
# rq_res <- rq_res[is.finite(rq_res)] # some Inf
# qqnorm(rq_res);abline(0, 1)
# 
# set.seed(123)
# ret<- simulate(m_e_per_bccm_final, nsim = 500, type = "mle-mvn") 
# r_ret <-  dharma_residuals(ret, m_e_per_bccm_final, return_DHARMa = TRUE)
# plot(r_ret)
# DHARMa::testResiduals(r_ret)
# 
# DHARMa::testSpatialAutocorrelation(r_ret, x = dat2$X, y = dat2$Y)
# 
# 
# # residuals are not plotting right, need to figure out right family!!!
# 
