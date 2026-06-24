###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# predict ML models 


###############################################################################
library(dplyr)
library(terra)
library(gbm)
library(pROC)
library(xgboost)
library(ModelMetrics)  


# Load functions
source("code/modelling-functions.R")
# Load data
load("code/output_data/seagrass_model_inputs.RData")
load("code/output_data/prediction_model_inputs.RData")

# Preprocess
seagrass_data_long <- seagrass_data_long %>%
  dplyr::select(
    -saltmean_bccm_sq_stnd, -saltmean_nep_sq_stnd,
    -slope_sqrt_stnd, -saltmin_bccm_sq_stnd, -saltmin_nep_sq_stnd
  ) %>%
  mutate(
    Survey = as.factor(substr(HKey, 1, 3)),
    HKey = as.factor(HKey),
    Year_factor = as.factor(Year)
  )

seagrass_data_long <- seagrass_data_long %>%
  mutate(Survey = recode(Survey, MSE = "BHM")) %>%
  mutate(Survey = recode(Survey, GSU = "RSU")) %>%
  mutate(Survey = recode(Survey, ABL = "RSU")) %>%
  mutate(Survey = factor(Survey))



#eelgrass models

sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
dat2 <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)

pred_vars = c("depth_stnd", "slope_stnd", "rei_stnd", "substrate",
  "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
  "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempcv_bccm_stnd", "Survey")

#BCCM model
env<- env_20m_all %>% dplyr::select(depth_stnd, slope_stnd, rei_stnd, substrate, airtempmin_stnd, rsdsmin_stnd, prmin_stnd, 
                                    saltcv_bccm_stnd, NH4_bccm_stnd, tempcv_bccm_stnd)

env$substrate <- factor(
  env$substrate,
  levels = levels(dat2$substrate)
)

env$Survey <- factor(
  "BHM",
  levels = levels(dat2$Survey)
)


set.seed(123)
# make a full gbm model without cv
gbm_fit <- gbm(
  presence ~ .,
  data = dat2[, c("presence", pred_vars)],
  distribution = "bernoulli",
  n.trees = 5000,
  interaction.depth = 2, #capture main effects, less risk of overfitting
  shrinkage = 0.005, # Low values (0.001–0.01) → slower learning, usually better generalization, but requires more trees.
  bag.fraction = 0.7,
  n.minobsinnode = 10,
  train.fraction = 1,
  verbose = FALSE
)

gbm_varimp_bccm_eelgrass <- summary(
  gbm_fit,
  n.trees = 5000,
  method = permutation.test.gbm,
  plotit = FALSE
)

save(gbm_varimp_bccm_eelgrass, file = "code/output_data/model_results/relimp_e_bccm_gbm.RData")

gbm_prob <- predict(
  gbm_fit,
  newdata = env,
  n.trees = 5000,
  type = "response"   # returns probabilities for bernoulli 
)

range(gbm_prob)


env_20m_all$gbm_bccm_prob <- gbm_prob

#standard error with bootstrapping
n_boot <- 50
pred_boot <- matrix(NA, nrow = nrow(env), ncol = n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  dat_boot <- dat2[boot_idx, ]
  
  gbm_boot <- gbm(
    presence ~ ., data = dat_boot[, c("presence", pred_vars)],
    distribution = "bernoulli", n.trees = 5000,
    interaction.depth = 2, shrinkage = 0.005,
    bag.fraction = 0.7, n.minobsinnode = 10,
    train.fraction = 1, verbose = FALSE
  )
  
  pred_boot[, i] <- predict(gbm_boot, newdata = env[, pred_vars], n.trees = 5000, type = "response")
}

# Mean and standard error
#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$gbm_bccm_se <- apply(pred_boot, 1, sd)


outdir <- file.path("./raster/eelgrass/gbm_bccm")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_gbm_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_gbm_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_gbm_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_gbm_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_gbm_bccm.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_se_gbm_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_se_gbm_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_se_gbm_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_se_gbm_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_se_gbm_bccm.tif")), overwrite = TRUE)


#make a full xgb model with no cv
# 1. Training matrix
X <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y <- dat2$presence

dtrain <- xgb.DMatrix(data = X, label = y)

# 2. Train XGBoost
xgb_fit <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 5,
    eta = 0.05,
    gamma = 1,
    min_child_weight = 5,
    subsample = 0.7,
    colsample_bytree = 0.7,
    nthread = parallel::detectCores() - 1
  ),
  data = dtrain,
  nrounds = 1500,
  verbose = 0
)

# 3. Variable importance using gain method
xgb_varimp <- xgb.importance(model = xgb_fit)

# 3. Variable importance using permutation method
xgb_perm_imp <- perm_importance_xgb(model = xgb_fit, X = X, y = y, nrep = 10)

# identify substrate-related rows
substrate_rows <- grepl("substrate", xgb_perm_imp$Feature)

# aggregate substrate gain
substrate_gain <- sum(xgb_perm_imp$RelImportance[substrate_rows])

# keep non-substrate variables
xgb_varimp_collapsed <- xgb_perm_imp[!substrate_rows, ]

# add aggregated substrate row
xgb_varimp_collapsed <- rbind(
  xgb_varimp_collapsed,
  data.frame(
    Feature = "substrate",
    RelImportance = substrate_gain,
    Importance = NA
  )
)

xgb_varimp_bccm_eelgrass <- xgb_varimp_collapsed[
  order(-xgb_varimp_collapsed$RelImportance), ]

save(xgb_varimp_bccm_eelgrass,
     file = "code/output_data/model_results/relimp_e_bccm_xgb.RData")



# 4. Predict onto your environmental dataframe
X_pred <- model.matrix(~ . - 1, data = env)
dtest <- xgb.DMatrix(X_pred)
env_20m_all$xgb_bccm_prob <- predict(xgb_fit, dtest)

# calculate xgb se
n_boot <- 50
X_train <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y_train <- dat2$presence

X_env <- model.matrix(~ . - 1, data = env)

#-------------------------------
# Bootstrapping loop
#-------------------------------
pred_boot <- matrix(NA, nrow = nrow(X_env), ncol = n_boot)

for (i in 1:n_boot) {
  
  # Resample training data with replacement
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  X_boot <- X_train[boot_idx, ]
  y_boot <- y_train[boot_idx]
  
  dtrain <- xgb.DMatrix(data = X_boot, label = y_boot)
  
  # Fit XGBoost
  xgb_model <- xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 5,
      eta = 0.05,
      gamma = 1,
      min_child_weight = 5,
      subsample = 0.7,
      colsample_bytree = 0.7
    ),
    data = dtrain,
    nrounds = 1500,
    verbose = 0
  )
  
  # Predict on full env dataset
  dtest <- xgb.DMatrix(data = X_env)
  pred_boot[, i] <- predict(xgb_model, dtest)
  
  message("Bootstrap iteration ", i, " complete")
}

#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$xgb_bccm_se   <- apply(pred_boot, 1, sd)



outdir <- file.path("./raster/eelgrass/xgb_bccm")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_xgb_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_xgb_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_xgb_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_xgb_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_xgb_bccm.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_se_xgb_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_se_xgb_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_se_xgb_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_se_xgb_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_se_xgb_bccm.tif")), overwrite = TRUE)




#NEP model
pred_vars = c("depth_stnd", "slope_stnd", "rei_stnd", "substrate",
              "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd",
              "saltcv_nep_stnd", "NH4_nep_stnd", "tempcv_nep_stnd", "Survey")

env<- env_20m_all %>% dplyr::select(depth_stnd, slope_stnd, rei_stnd, substrate, airtempmin_stnd, rsdsmin_stnd, prmin_stnd, 
                                    saltcv_nep_stnd, NH4_nep_stnd, tempcv_nep_stnd)

env$substrate <- factor(
  env$substrate,
  levels = levels(dat2$substrate)
)

env$Survey <- factor(
  "BHM",
  levels = levels(dat2$Survey)
)


set.seed(123)
# make a full gbm model without cv
gbm_fit <- gbm(
  presence ~ .,
  data = dat2[, c("presence", pred_vars)],
  distribution = "bernoulli",
  n.trees = 5000,
  interaction.depth = 2, #capture main effects, less risk of overfitting
  shrinkage = 0.005, # Low values (0.001–0.01) → slower learning, usually better generalization, but requires more trees.
  bag.fraction = 0.7,
  n.minobsinnode = 10,
  train.fraction = 1,
  verbose = FALSE
)

gbm_varimp_nep_eelgrass <- summary(
  gbm_fit,
  n.trees = 5000,
  method = permutation.test.gbm,
  plotit = FALSE
)

save(gbm_varimp_nep_eelgrass, file = "code/output_data/model_results/relimp_e_nep_gbm.RData")

gbm_prob <- predict(
  gbm_fit,
  newdata = env,
  n.trees = 5000,
  type = "response"   # returns probabilities for bernoulli 
)

range(gbm_prob)


env_20m_all$gbm_nep_prob <- gbm_prob

#standard error with bootstrapping
n_boot <- 50
pred_boot <- matrix(NA, nrow = nrow(env), ncol = n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  dat_boot <- dat2[boot_idx, ]
  
  gbm_boot <- gbm(
    presence ~ ., data = dat_boot[, c("presence", pred_vars)],
    distribution = "bernoulli", n.trees = 5000,
    interaction.depth = 2, shrinkage = 0.005,
    bag.fraction = 0.7, n.minobsinnode = 10,
    train.fraction = 1, verbose = FALSE
  )
  
  pred_boot[, i] <- predict(gbm_boot, newdata = env[, pred_vars], n.trees = 5000, type = "response")
}

# Mean and standard error
#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$gbm_nep_se <- apply(pred_boot, 1, sd)


outdir <- file.path("./raster/eelgrass/gbm_nep")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_gbm_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_gbm_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_gbm_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_gbm_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_gbm_nep.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_se_gbm_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_se_gbm_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_se_gbm_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_se_gbm_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_se_gbm_nep.tif")), overwrite = TRUE)


#make a full xgb model with no cv
# 1. Training matrix
X <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y <- dat2$presence

dtrain <- xgb.DMatrix(data = X, label = y)

# 2. Train XGBoost
xgb_fit <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 5,
    eta = 0.05,
    gamma = 1,
    min_child_weight = 5,
    subsample = 0.7,
    colsample_bytree = 0.7,
    nthread = parallel::detectCores() - 1
  ),
  data = dtrain,
  nrounds = 1500,
  verbose = 0
)

# 3. Variable importance using gain method
xgb_varimp <- xgb.importance(model = xgb_fit)

# 3. Variable importance using permutation method
xgb_perm_imp <- perm_importance_xgb(model = xgb_fit, X = X, y = y, nrep = 10)

# identify substrate-related rows
substrate_rows <- grepl("substrate", xgb_perm_imp$Feature)

# aggregate substrate gain
substrate_gain <- sum(xgb_perm_imp$RelImportance[substrate_rows])

# keep non-substrate variables
xgb_varimp_collapsed <- xgb_perm_imp[!substrate_rows, ]

# add aggregated substrate row
xgb_varimp_collapsed <- rbind(
  xgb_varimp_collapsed,
  data.frame(
    Feature = "substrate",
    RelImportance = substrate_gain,
    Importance = NA
  )
)

xgb_varimp_nep_eelgrass <- xgb_varimp_collapsed[
  order(-xgb_varimp_collapsed$RelImportance), ]

save(xgb_varimp_nep_eelgrass, file = "code/output_data/model_results/relimp_e_nep_xgb.RData")

# 4. Predict onto your environmental dataframe
X_pred <- model.matrix(~ . - 1, data = env)
dtest <- xgb.DMatrix(X_pred)
env_20m_all$xgb_nep_prob <- predict(xgb_fit, dtest)

# calculate xgb se
n_boot <- 50
X_train <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y_train <- dat2$presence

X_env <- model.matrix(~ . - 1, data = env)

#-------------------------------
# Bootstrapping loop
#-------------------------------
pred_boot <- matrix(NA, nrow = nrow(X_env), ncol = n_boot)

for (i in 1:n_boot) {
  
  # Resample training data with replacement
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  X_boot <- X_train[boot_idx, ]
  y_boot <- y_train[boot_idx]
  
  dtrain <- xgb.DMatrix(data = X_boot, label = y_boot)
  
  # Fit XGBoost
  xgb_model <- xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 5,
      eta = 0.05,
      gamma = 1,
      min_child_weight = 5,
      subsample = 0.7,
      colsample_bytree = 0.7
    ),
    data = dtrain,
    nrounds = 1500,
    verbose = 0
  )
  
  # Predict on full env dataset
  dtest <- xgb.DMatrix(data = X_env)
  pred_boot[, i] <- predict(xgb_model, dtest)
  
  message("Bootstrap iteration ", i, " complete")
}

#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$xgb_nep_se   <- apply(pred_boot, 1, sd)



outdir <- file.path("./raster/eelgrass/xgb_nep")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_xgb_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_xgb_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_xgb_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_xgb_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_xgb_nep.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("eelgrass_predictions_hg_se_xgb_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("eelgrass_predictions_ss_se_xgb_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("eelgrass_predictions_wcvi_xgb_gbm_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("eelgrass_predictions_ncc_se_xgb_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("eelgrass_predictions_qcs_se_xgb_nep.tif")), overwrite = TRUE)




#HAVE NO MODIFIED ANYTHING FOR SURFGRASS YET!!


#surfgrass models

sp = "PH"
numFolds <- length(unique(seagrass_data$fold_seagrass))
dat2 <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_seagrass)


pred_vars = c("depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd",
  "substrate", "cul_eff_stnd", "airtempcv_stnd",
  "rsdsmin_stnd", "prmean_stnd", "saltcv_bccm_stnd",
  "tempmin_bccm_stnd", "surftempcv_bccm_stnd")

#BCCM model
env<- env_20m_all %>% dplyr::select(depth_stnd,  rei_sqrt_stnd, tidal_sqrt_stnd, substrate, cul_eff_stnd, 
                                    airtempcv_stnd, rsdsmin_stnd, prmean_stnd, saltcv_bccm_stnd,
                                    tempmin_bccm_stnd, surftempcv_bccm_stnd)
env$substrate <- as.factor(env$substrate)


set.seed(123)
# make a full gbm model without cv
gbm_fit <- gbm(
  presence ~ .,
  data = dat2[, c("presence", pred_vars)],
  distribution = "bernoulli",
  n.trees = 5000,
  interaction.depth = 2, #capture main effects, less risk of overfitting
  shrinkage = 0.005, # Low values (0.001–0.01) → slower learning, usually better generalization, but requires more trees.
  bag.fraction = 0.7,
  n.minobsinnode = 10,
  train.fraction = 1,
  verbose = FALSE
)

gbm_varimp_bccm_surfgrass <- summary(
  gbm_fit,
  n.trees = 5000,
  method = permutation.test.gbm,
  plotit = FALSE
)

save(gbm_varimp_bccm_surfgrass, file = "code/output_data/model_results/relimp_s_bccm_gbm.RData")



gbm_prob <- predict(
  gbm_fit,
  newdata = env,
  n.trees = 5000,
  type = "response"   # returns probabilities for bernoulli 
)

range(gbm_prob)


env_20m_all$gbm_bccm_prob <- gbm_prob

#standard error with bootstrapping
n_boot <- 50
pred_boot <- matrix(NA, nrow = nrow(env), ncol = n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  dat_boot <- dat2[boot_idx, ]
  
  gbm_boot <- gbm(
    presence ~ ., data = dat_boot[, c("presence", pred_vars)],
    distribution = "bernoulli", n.trees = 5000,
    interaction.depth = 2, shrinkage = 0.005,
    bag.fraction = 0.7, n.minobsinnode = 10,
    train.fraction = 1, verbose = FALSE
  )
  
  pred_boot[, i] <- predict(gbm_boot, newdata = env[, pred_vars], n.trees = 5000, type = "response")
}

# Mean and standard error
#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$gbm_bccm_se <- apply(pred_boot, 1, sd)


outdir <- file.path("./raster/surfgrass/gbm_bccm")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_gbm_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_gbm_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_gbm_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_gbm_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_gbm_bccm.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_se_gbm_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_se_gbm_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_se_gbm_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_se_gbm_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_bccm_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_se_gbm_bccm.tif")), overwrite = TRUE)


#make a full xgb model with no cv
# 1. Training matrix
X <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y <- dat2$presence

dtrain <- xgb.DMatrix(data = X, label = y)

# 2. Train XGBoost
xgb_fit <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 5,
    eta = 0.05,
    gamma = 1,
    min_child_weight = 5,
    subsample = 0.7,
    colsample_bytree = 0.7,
    nthread = parallel::detectCores() - 1
  ),
  data = dtrain,
  nrounds = 1500,
  verbose = 0
)

# 3. Variable importance using gain method
xgb_varimp <- xgb.importance(model = xgb_fit)

# 3. Variable importance using permutation method
xgb_perm_imp <- perm_importance_xgb(model = xgb_fit, X = X, y = y, nrep = 10)

# identify substrate-related rows
substrate_rows <- grepl("substrate", xgb_perm_imp$Feature)

# aggregate substrate gain
substrate_gain <- sum(xgb_perm_imp$RelImportance[substrate_rows])

# keep non-substrate variables
xgb_varimp_collapsed <- xgb_perm_imp[!substrate_rows, ]

# add aggregated substrate row
xgb_varimp_collapsed <- rbind(
  xgb_varimp_collapsed,
  data.frame(
    Feature = "substrate",
    RelImportance = substrate_gain,
    Importance = NA
  )
)

xgb_varimp_bccm_surfgrass <- xgb_varimp_collapsed[
  order(-xgb_varimp_collapsed$RelImportance), ]

save(xgb_varimp_bccm_surfgrass, file = "code/output_data/model_results/relimp_s_bccm_xgb.RData")

# 4. Predict onto your environmental dataframe
X_pred <- model.matrix(~ . - 1, data = env)
dtest <- xgb.DMatrix(X_pred)
env_20m_all$xgb_bccm_prob <- predict(xgb_fit, dtest)

# calculate xgb se
n_boot <- 50
X_train <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y_train <- dat2$presence

X_env <- model.matrix(~ . - 1, data = env)

#-------------------------------
# Bootstrapping loop
#-------------------------------
pred_boot <- matrix(NA, nrow = nrow(X_env), ncol = n_boot)

for (i in 1:n_boot) {
  
  # Resample training data with replacement
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  X_boot <- X_train[boot_idx, ]
  y_boot <- y_train[boot_idx]
  
  dtrain <- xgb.DMatrix(data = X_boot, label = y_boot)
  
  # Fit XGBoost
  xgb_model <- xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 5,
      eta = 0.05,
      gamma = 1,
      min_child_weight = 5,
      subsample = 0.7,
      colsample_bytree = 0.7
    ),
    data = dtrain,
    nrounds = 1500,
    verbose = 0
  )
  
  # Predict on full env dataset
  dtest <- xgb.DMatrix(data = X_env)
  pred_boot[, i] <- predict(xgb_model, dtest)
  
  message("Bootstrap iteration ", i, " complete")
}

#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$xgb_bccm_se   <- apply(pred_boot, 1, sd)



outdir <- file.path("./raster/surfgrass/xgb_bccm")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_xgb_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_xgb_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_xgb_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_xgb_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_xgb_bccm.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_se_xgb_bccm.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_se_xgb_bccm.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_xgb_gbm_bccm.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_se_xgb_bccm.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_bccm_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_se_xgb_bccm.tif")), overwrite = TRUE)




#NEP model
pred_vars = c( "depth_stnd", "rei_sqrt_stnd", "tidal_sqrt_stnd",
               "substrate", "cul_eff_stnd", "airtempcv_stnd",
               "rsdsmin_stnd", "prmean_stnd", "saltmean_nep_stnd",
               "tempmean_nep_stnd", "surftempcv_nep_stnd")

env<- env_20m_all %>% dplyr::select( depth_stnd, rei_sqrt_stnd, tidal_sqrt_stnd,
                                     substrate, cul_eff_stnd, airtempcv_stnd,
                                     rsdsmin_stnd, prmean_stnd, saltmean_nep_stnd,
                                     tempmean_nep_stnd, surftempcv_nep_stnd)
env$substrate <- as.factor(env$substrate)

set.seed(123)
# make a full gbm model without cv
gbm_fit <- gbm(
  presence ~ .,
  data = dat2[, c("presence", pred_vars)],
  distribution = "bernoulli",
  n.trees = 5000,
  interaction.depth = 2, #capture main effects, less risk of overfitting
  shrinkage = 0.005, # Low values (0.001–0.01) → slower learning, usually better generalization, but requires more trees.
  bag.fraction = 0.7,
  n.minobsinnode = 10,
  train.fraction = 1,
  verbose = FALSE
)

gbm_varimp_nep_surfgrass <- summary(
  gbm_fit,
  n.trees = 5000,
  method = permutation.test.gbm,
  plotit = FALSE
)

save(gbm_varimp_nep_surfgrass, file = "code/output_data/model_results/relimp_s_nep_gbm.RData")

gbm_prob <- predict(
  gbm_fit,
  newdata = env,
  n.trees = 5000,
  type = "response"   # returns probabilities for bernoulli 
)

range(gbm_prob)


env_20m_all$gbm_nep_prob <- gbm_prob

#standard error with bootstrapping
n_boot <- 50
pred_boot <- matrix(NA, nrow = nrow(env), ncol = n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  dat_boot <- dat2[boot_idx, ]
  
  gbm_boot <- gbm(
    presence ~ ., data = dat_boot[, c("presence", pred_vars)],
    distribution = "bernoulli", n.trees = 5000,
    interaction.depth = 2, shrinkage = 0.005,
    bag.fraction = 0.7, n.minobsinnode = 10,
    train.fraction = 1, verbose = FALSE
  )
  
  pred_boot[, i] <- predict(gbm_boot, newdata = env[, pred_vars], n.trees = 5000, type = "response")
}

# Mean and standard error
#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$gbm_nep_se <- apply(pred_boot, 1, sd)


outdir <- file.path("./raster/surfgrass/gbm_nep")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_gbm_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_gbm_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_gbm_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_gbm_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_nep_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_gbm_nep.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_se_gbm_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_se_gbm_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_se_gbm_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_se_gbm_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, gbm_nep_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_se_gbm_nep.tif")), overwrite = TRUE)


#make a full xgb model with no cv
# 1. Training matrix
X <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y <- dat2$presence

dtrain <- xgb.DMatrix(data = X, label = y)

# 2. Train XGBoost
xgb_fit <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 5,
    eta = 0.05,
    gamma = 1,
    min_child_weight = 5,
    subsample = 0.7,
    colsample_bytree = 0.7,
    nthread = parallel::detectCores() - 1
  ),
  data = dtrain,
  nrounds = 1500,
  verbose = 0
)

# 3. Variable importance using gain method
xgb_varimp <- xgb.importance(model = xgb_fit)

# 3. Variable importance using permutation method
xgb_perm_imp <- perm_importance_xgb(model = xgb_fit, X = X, y = y, nrep = 10)

# identify substrate-related rows
substrate_rows <- grepl("substrate", xgb_perm_imp$Feature)

# aggregate substrate gain
substrate_gain <- sum(xgb_perm_imp$RelImportance[substrate_rows])

# keep non-substrate variables
xgb_varimp_collapsed <- xgb_perm_imp[!substrate_rows, ]

# add aggregated substrate row
xgb_varimp_collapsed <- rbind(
  xgb_varimp_collapsed,
  data.frame(
    Feature = "substrate",
    RelImportance = substrate_gain,
    Importance = NA
  )
)

xgb_varimp_nep_surfgrass <- xgb_varimp_collapsed[
  order(-xgb_varimp_collapsed$RelImportance), ]

save(xgb_varimp_nep_surfgrass, file = "code/output_data/model_results/relimp_s_nep_xgb.RData")

# 4. Predict onto your environmental dataframe
X_pred <- model.matrix(~ . - 1, data = env)
dtest <- xgb.DMatrix(X_pred)
env_20m_all$xgb_nep_prob <- predict(xgb_fit, dtest)

# calculate xgb se
n_boot <- 50
X_train <- model.matrix(~ . - 1, data = dat2[, pred_vars])
y_train <- dat2$presence

X_env <- model.matrix(~ . - 1, data = env)

#-------------------------------
# Bootstrapping loop
#-------------------------------
pred_boot <- matrix(NA, nrow = nrow(X_env), ncol = n_boot)

for (i in 1:n_boot) {
  
  # Resample training data with replacement
  boot_idx <- sample(1:nrow(dat2), replace = TRUE)
  X_boot <- X_train[boot_idx, ]
  y_boot <- y_train[boot_idx]
  
  dtrain <- xgb.DMatrix(data = X_boot, label = y_boot)
  
  # Fit XGBoost
  xgb_model <- xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 5,
      eta = 0.05,
      gamma = 1,
      min_child_weight = 5,
      subsample = 0.7,
      colsample_bytree = 0.7
    ),
    data = dtrain,
    nrounds = 1500,
    verbose = 0
  )
  
  # Predict on full env dataset
  dtest <- xgb.DMatrix(data = X_env)
  pred_boot[, i] <- predict(xgb_model, dtest)
  
  message("Bootstrap iteration ", i, " complete")
}

#env_20m_all$pred_mean <- rowMeans(pred_boot)
env_20m_all$xgb_nep_se   <- apply(pred_boot, 1, sd)



outdir <- file.path("./raster/surfgrass/xgb_nep")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_xgb_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_xgb_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_xgb_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_xgb_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_nep_prob)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_xgb_nep.tif")), overwrite = TRUE)

raster_hg <- env_20m_all %>%
  filter(region == "Haida Gwaii") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_hg <- rast(x = raster_hg %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_hg, file.path(outdir, paste0("surfgrass_predictions_hg_se_xgb_nep.tif")), overwrite = TRUE)

raster_ss <- env_20m_all %>%
  filter(region == "Salish Sea") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_ss <- rast(x = raster_ss %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ss, file.path(outdir, paste0("surfgrass_predictions_ss_se_xgb_nep.tif")), overwrite = TRUE)

raster_wcvi <- env_20m_all %>%
  filter(region == "West Coast Vancouver Island") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_wcvi <- rast(x = raster_wcvi %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_wcvi, file.path(outdir, paste0("surfgrass_predictions_wcvi_xgb_gbm_nep.tif")), overwrite = TRUE)

raster_ncc <- env_20m_all %>%
  filter(region == "North Central Coast") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_ncc <- rast(x = raster_ncc %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_ncc, file.path(outdir, paste0("surfgrass_predictions_ncc_se_xgb_nep.tif")), overwrite = TRUE)

raster_qcs <- env_20m_all %>%
  filter(region == "Queen Charlotte Strait") %>%
  dplyr::select(X_m, Y_m, xgb_nep_se)
raster_qcs <- rast(x = raster_qcs %>% as.matrix, type = "xyz", crs = "EPSG:3005")
writeRaster(raster_qcs, file.path(outdir, paste0("surfgrass_predictions_qcs_se_xgb_nep.tif")), overwrite = TRUE)



