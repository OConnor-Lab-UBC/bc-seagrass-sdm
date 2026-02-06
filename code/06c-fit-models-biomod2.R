###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# fit sdm models with biomod2 
# # none of the ML models do better than the sdmTMB models  
# but this was a test if we wanted to continue to do further validation beyond cross validation and temporal forecasting
# and we do not
###############################################################################


library(biomod2)
library(dplyr)
library(terra)

load("code/output_data/seagrass_model_inputs.RData")
source("code/modelling-functions.R")

seagrass_data_long <- seagrass_data_long %>% select(-saltmean_bccm_sq_stnd, -saltmean_nep_sq_stnd, -slope_sqrt_stnd, -saltmin_bccm_sq_stnd, -saltmin_nep_sq_stnd)
seagrass_data_long <- seagrass_data_long %>%
  mutate(Survey = as.factor(substr(HKey, 1, 3)),
         HKey = as.factor(HKey),
         Year_factor = as.factor(Year))

####Eelgrass model####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
data <- data %>% mutate(Survey = as.factor(Survey)) # was still not recognizing i changed survye to factor so did it again

pred_vars <- c("depth_stnd", "slope_stnd", "rei_stnd", "substrate", "airtempmin_stnd", "rsdsmin_stnd", "prmin_stnd", "saltcv_bccm_stnd", "NH4_bccm_stnd", "tempmin_bccm_stnd")

dat2 <- data %>%
  dplyr::select(presence, X_m, Y_m, fold, Year, all_of(pred_vars)) %>%
  na.omit()

dat2$presence <- as.numeric(dat2$presence)

myBiomodData <- BIOMOD_FormatingData(
  resp.var  = dat2$presence,
  expl.var  = dat2 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = dat2[, c("X_m", "Y_m")],
  resp.name = "eelgrass",
  PA.nb.rep = 0   )


folds <- sort(unique(dat2$fold))
K <- length(folds)

cv_table <- matrix(TRUE, nrow = nrow(dat2), ncol = K + 1)  # extra column for full model
for (k in seq_along(folds)) {
  test_fold <- folds[k]
  cv_table[dat2$fold == test_fold, k] <- FALSE  # use your fold column
}
colnames(cv_table) <- c(paste0("_allData_RUN", seq_len(K)), "_allData_allRun")

ml_models <- c("GBM", "RF", "XGBOOST", "ANN", "CTA")

myOptions <- bm_ModelingOptions(
  data.type = "binary",
  models = ml_models,
  strategy = "bigboss"  # stronger defaults than the "default" set
)

# need to change var import to 10-20 for final runs
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = ml_models,
  CV.strategy = "user.defined",
  CV.user.table = cv_table,
  CV.do.full.models = TRUE,
  metric.eval = c("AUCroc", "TSS", "KAPPA"),
  var.import = 3,
  seed.val = 123
)

evals <- get_evaluations(myBiomodModelOut)

plot_output <- bm_PlotEvalMean(bm.out = myBiomodModelOut, metric.eval = c("AUCroc", "TSS"), dataset = 'validation') 
ggdat <- plot_output$tab
names(ggdat) <- c("algo", "AUC_mean","TSS_mean","AUC_sd","TSS_sd")
ggdat$algo <- as.factor(ggdat$algo)



get_variables_importance(myBiomodModelOut)



bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# random forest here can get cut out, TSS is horrible

# none of the ML models do better than the sdmTMB models 

# test forecasting
data_pre2013 <- dat2 %>% filter(Year < 2010)
data_post2012 <- dat2 %>% filter(Year > 2012)

# Format data for BIOMOD2
myBiomodData_pre2013 <- BIOMOD_FormatingData(
  resp.var  = data_pre2013$presence,
  expl.var  = data_pre2013 %>% dplyr::select(all_of(pred_vars)),
  resp.xy   = data_pre2013[, c("X_m", "Y_m")],
  resp.name = "eelgrass",
  PA.nb.rep = 0
)

n <- nrow(data_pre2013)
cv_table <- matrix(TRUE, nrow = n, ncol = 1)  # 1 column = 1 run
colnames(cv_table) <- "_allData_allRun"       # required by BIOMOD2


myBiomodModelOut_pre2013 <- BIOMOD_Modeling(
  bm.format       = myBiomodData_pre2013,
  models          = ml_models,
  CV.strategy     = "user.defined",
  CV.user.table   = cv_table,
  CV.do.full.models = TRUE,  # ensures BIOMOD2 builds the “full model”
  metric.eval     = c("AUCroc","TSS","KAPPA"),
  var.import      = 3,
  seed.val        = 123,
  modeling.id     = "pre2013_ML"
)


library(randomForest)
library(xgboost)
library(gbm)
library(ModelMetrics)  
library(nnet)       
library(MASS) 
library(rpart) 


# Random Forest
rf_mod <- randomForest(
  x = data_pre2013[, pred_vars],
  y = as.factor(data_pre2013$presence),
  ntree = 500
)

# GBM
gbm_mod <- gbm(
  formula = presence ~ .,
  distribution = "bernoulli",
  data = data_pre2013[, c("presence", pred_vars)],
  n.trees = 500,
  interaction.depth = 3,
  shrinkage = 0.01,
  verbose = FALSE
)

# XGBoost
xgb_matrix_pre <- model.matrix(~ . -1, data = data_pre2013[, pred_vars])
xgb_matrix_post <- model.matrix(~ . -1, data = data_post2012[, pred_vars])

dtrain <- xgb.DMatrix(
  data = xgb_matrix_pre,
  label = data_pre2013$presence
)

xgb_params <- list(
  objective = "binary:logistic",  # presence/absence
  eval_metric = "auc",
  max_depth = 6,
  eta = 0.1,
  nthread = 4
)

xgb_mod <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 200,
  verbose = 1
)

# --- ANN ---

data_pre2013_num <- data_pre2013 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.)))) %>%
  replace(is.na(.), 0)  
data_post2012_num <- data_post2012 %>%
  dplyr::select(all_of(pred_vars)) %>%
  mutate(across(everything(), ~ if(is.factor(.)) as.numeric(.) else as.numeric(as.character(.))))%>%
  replace(is.na(.), 0)  
#data_pre2013_num$substrate <- as.factor(data_pre2013$substrate)

y_pre_1hot <- class.ind(as.factor(data_pre2013$presence))  # one-hot encoding: 0 -> [1,0], 1 -> [0,1]

# -------------------------------
# 3️⃣ Train ANN
# -------------------------------
set.seed(123)
ann_mod <- nnet(
  x = as.matrix(data_pre2013_num),
  y = y_pre_1hot,
  size = 5,       # hidden neurons
  linout = FALSE,
  entropy = TRUE,
  maxit = 500,
  trace = FALSE
)



# --- CTA (Classification Tree Analysis) ---
cta_mod <- rpart(
  formula = presence ~ .,
  data = cbind(presence = as.factor(data_pre2013$presence), data_pre2013_num),
  method = "class",
  control = rpart.control(minsplit = 5, cp = 0.001)
)




predictions <- list(
  RF = list(
    training = predict(rf_mod, data_pre2013[, pred_vars], type = "prob")[,2],
    forecast = predict(rf_mod, data_post2012[, pred_vars], type = "prob")[,2]
  ),
  GBM = list(
    training = predict(gbm_mod, data_pre2013[, pred_vars], n.trees = 2000, type = "response"),
    forecast = predict(gbm_mod, data_post2012[, pred_vars], n.trees = 2000, type = "response")
  ),
  XGBOOST = list(
    training = predict(xgb_mod, newdata = xgb_matrix_pre),
    forecast = predict(xgb_mod, newdata = xgb_matrix_post)
  ),
  ANN = list(
    training = predict(ann_mod, newdata = as.matrix(data_pre2013_num))[,2],
    forecast = predict(ann_mod, newdata = as.matrix(data_post2012_num))[,2]
  ),
  CTA = list(
    training = predict(cta_mod, newdata = data_pre2013_num, type = "prob")[,2],
    forecast = predict(cta_mod, newdata = data_post2012_num, type = "prob")[,2]
  )
)


eval_results <- lapply(names(predictions), function(mod) {
  data.frame(
    model = mod,
    type = c("training", "forecast"),
    AUC = c(
      ModelMetrics::auc(data_pre2013$presence, predictions[[mod]]$training),
      ModelMetrics::auc(data_post2012$presence, predictions[[mod]]$forecast)
    ),
    TjurR2 = c(
      tjur(data_pre2013$presence, predictions[[mod]]$training),
      tjur(data_post2012$presence, predictions[[mod]]$forecast)
    )
  )
}) %>% bind_rows()

print(eval_results)


# none of the ML models do better than the sdmTMB models 