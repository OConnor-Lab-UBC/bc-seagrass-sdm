###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# fit sdm models with flexsdm 
#
###############################################################################

# devtools::install_github('sjevelazco/flexsdm')
library(flexsdm)
library(dplyr)
library(terra)

load("code/output_data/seagrass_model_inputs.RData")
load("code/output_data/prediction_model_inputs.RData")
seagrass_data_long <- seagrass_data_long %>% select(-saltmean_sq_stnd, -slope_sqrt_stnd, -saltmin_sq_stnd)

rei_df <- env_20m_all %>%
  select(X_m, Y_m, rei_stnd)
names(rei_df) <- c("x", "y", "rei_stnd")  # Rename to expected column names

# Create SpatVector of points
points <- vect(rei_df, geom = c("x", "y"), crs = "EPSG:3005")

# Create an empty raster template
# Define extent and resolution explicitly — here assuming 20m grid
r_template <- rast(ext(points), resolution = 20, crs = "EPSG:3005")

# Rasterize using the 'z' values
rei_raster <- rasterize(points, r_template, field = "rei_stnd")
names(rei_raster) <- "rei_stnd"

depth_df <- env_20m_all %>%
  select(X_m, Y_m, depth_stnd)
names(depth_df) <- c("x", "y", "depth_stnd")  # Rename to expected column names

# Create SpatVector of points
points <- vect(depth_df, geom = c("x", "y"), crs = "EPSG:3005")

# Rasterize using the 'z' values
depth_raster <- rasterize(points, r_template, field = "depth_stnd")
names(depth_raster) <- "depth_stnd"

slope_df <- env_20m_all %>%
  select(X_m, Y_m, slope_stnd)
names(slope_df) <- c("x", "y", "slope_stnd")  # Rename to expected column names

# Create SpatVector of points
points <- vect(slope_df, geom = c("x", "y"), crs = "EPSG:3005")

# Rasterize using the 'z' values
slope_raster <- rasterize(points, r_template, field = "slope_stnd")
names(slope_raster) <- "slope_stnd"

substrate_df <- env_20m_all %>%
  select(X_m, Y_m, substrate)
substrate_levels <- c("Rock", "Mixed", "Sand", "Mud")

# Reassign numeric codes based on position in substrate_levels
substrate_df$substrate <- match(substrate_df$substrate, substrate_levels)

names(substrate_df) <- c("x", "y", "substrate")  # Rename to expected column names

# Create SpatVector of points
points <- vect(substrate_df, geom = c("x", "y"), crs = "EPSG:3005")

# Rasterize using the 'z' values
substrate_raster <- rasterize(points, r_template, field = "substrate")
names(substrate_raster) <- "substrate"

somevar <- terra::rast(c(rei_raster, slope_raster, depth_raster, substrate_raster))


####Eelgrass model####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
data$substrate <- match(data$substrate, substrate_levels)
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))

# predictors
predictors_continuous = c("depth_stnd", "slope_stnd", "rei_stnd")
predictors_cat = c("substrate")

# random forest model
tune_grid <-
  expand.grid(
    mtry = seq(1, 7, 1),
    ntree = c(300, 500, 700, 900)
  )

rf_t <-
  fit_raf(
    data = data,
    response = "presence",
    predictors = predictors_continuous,
    predictors_f = predictors_cat,
    partition = "fold",
    #grid = tune_grid,
    ntree=500,
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

#glm
glm1 <- 
  fit_glm(
    data = data,
    response = "presence",
    predictors = predictors_continuous,
    predictors_f = predictors_cat,
    select_pred = TRUE,
    partition = "fold",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = NULL,
    poly = 2,
    inter_order = 0
  )


# generalized boosted regression methods
gbm_t1 <- 
  fit_gbm(
    data = data,
    response = "presence",
    predictors = predictors_continuous,
    predictors_f = predictors_cat,
    partition = "fold",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    n_trees = 100
  )

# gam
gam1 <-
  fit_gam(
    data = data,
    response = "presence",
    predictors = predictors_continuous,
    predictors_f = predictors_cat,
    partition = "fold",
    select_pred = FALSE, #not working when TRUE
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = NULL,
    k = -1
  )

# neural networks
net_t1<- fit_net(
  data = data,
  response = "presence",
  predictors = predictors_continuous,
  predictors_f = predictors_cat,
  partition = "fold",
  fit_formula = NULL,
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  size = 2,
  decay = 0.1
)

# maxent model, this one takes a long time so maybe notworking?
max_t1 <- 
  fit_max(
    data = data,
    response = "presence",
    predictors = predictors_continuous,
    predictors_f = predictors_cat,
    partition = "fold",
    background = NULL,
    fit_formula = NULL,
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    clamp = TRUE,
    classes = "default",
    pred_type = "cloglog",
    regmult = 1
  )

# support vector machine, this one takes a long time so maybe notworking?
svm_t1 <- 
  fit_svm(
    data = data,
    response = "presence",
    predictors = predictors_continuous,
    predictors_f = predictors_cat,
    partition = "fold",
    fit_formula = NULL,
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    sigma = "automatic",
    C = 1
  )

merge_df <- sdm_summarize(models = list(rf_t, glm1, gbm_t1, gam1, net_t1, svm_t1, max_t1))
knitr::kable(
  merge_df %>% dplyr::select(
    model,
    AUC = AUC_mean,
    TSS = TSS_mean,
    JACCARD = JACCARD_mean,
    BOYCE = BOYCE_mean,
    IMAE = IMAE_mean
  )
)
#random forest can be removed


mensemble <- 
  fit_ensemble(
    models = list(rf_t,  gbm_t1,  net_t1, svm_t1, max_t1),
    ens_method = "meanw",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    thr_model = "max_sens_spec",
    metric = "TSS"
  )

model_perf <- sdm_summarize(list(rf_t, gbm_t1, net_t1,svm_t1, max_t1, mensemble))
model_perf

pr_1 <- sdm_predict(
  models = glm1,
  pred = somevar,
  thr = "max_sens_spec",
  con_thr = FALSE,
  predict_area = NULL,
  nchunk = 10
)
#> Predicting ensembles

unconstrained <- pr_1$meanw[[1]]
names(unconstrained) <- "unconstrained"

cl <- c("#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D", "#3D4988", "#462777", "#440154")
plot(unconstrained, col = cl, legend = FALSE, axes = FALSE)



