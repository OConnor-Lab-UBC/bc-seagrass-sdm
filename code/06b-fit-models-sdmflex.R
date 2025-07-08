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
#library(terra)

load("code/output_data/seagrass_model_inputs.RData")
load("code/output_data/prediction_model_inputs.RData")
seagrass_data_long <- seagrass_data_long %>% select(-saltmean_sq_stnd, -slope_sqrt_stnd, -saltmin_sq_stnd)

####Eelgrass model####
sp = "ZO"
numFolds <- length(unique(seagrass_data$fold_eelgrass))
data <- filter(seagrass_data_long, species == sp) %>% rename(fold = fold_eelgrass)
print(paste(sp, " present in ", round((sum(data$presence)/nrow(data))*100,2), "% of observations", sep = ""))

# random forest model
tune_grid <-
  expand.grid(
    mtry = seq(1, 7, 1),
    ntree = c(300, 500, 700, 900)
  )

rf_t <-
  tune_raf(
    data = data,
    response = "presence",
    predictors = c("depth", "slope_stnd", "rei_stnd"),
    predictors_f = c("substrate"),
    partition = "fold",
    grid = tune_grid,
    thr = "max_sens_spec",
    metric = "TSS",
  )


# generalized boosted regression methods
gbm_t1 <- fit_gbm(
  data = data,
  response = "presence",
  predictors = c("depth", "slope_stnd", "rei_stnd"),
  predictors_f = c("substrate"),
  partition = "fold",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
)
gbm_t1$model
gbm_t1$predictors
gbm_t1$performance
gbm_t1$data_ens


