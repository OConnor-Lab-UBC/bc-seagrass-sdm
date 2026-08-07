###############################################################################
#
# Authors:      Ashley Park
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
##
# Objective:
# ---------
# figures for paper
###############################################################################

library(ggplot2)
library(tidyverse)
library(forcats)
library(patchwork)
library(sf)
library(terra)
library(tidyterra)
library(GGally)
library(reproducible)
library(factoextra)
library(viridis)
library(gridExtra)
library(grid)
library(purrr)
library(tibble)
library(stringr)
library(tidyterra)


#functions
source("code/modelling-functions.R")

coastline_full <- st_read("raw_data/CHS_HWL2015_Coastline.gdb", layer = "Polygon_CHS_Pacific_HWL_2015_5028437_simple")
coastline <- coastline_full %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_crop(st_bbox(coastline_full) + c(0, 100000, -10000, -80000)) %>%
  st_transform(crs = "EPSG:3005")



#spatial and temporal blocks
load("code/output_data/seagrass_sf_file.RData")

jewel_cb_10 <- c(
  "#0F766E", "#3B5BA5", "#7A5195", "#C0567D", "#DD6B20",
  "#2F855A", "#4C6A92", "#9C4221", "#6B46C1", "#B7791F")

okabe_ito_10 <- c(
  "#E69F00", "#56B4E9", "#009E73","#F0E442","#0072B2", 
  "#D55E00", "#CC79A7", "#000000", "#999999", "#882255" )

base_map_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

map_extent <- coord_sf(
  xlim = c(530000, 1228000),
  ylim = c(366000, 1101000),
  expand = FALSE
)

plot_data <- bind_rows(
  spatialised_sf %>%
    mutate(species = "Eelgrass dataset (1993–2023)\n10 spatial folds",
           fold = factor(fold_eelgrass)),
  spatialised_sf %>%
    mutate(species = "Surfgrass dataset (1993–2023)\n10 spatial folds",
           fold = factor(fold_seagrass))
)

spatial_fold <- ggplot() +
  geom_sf(data = coastline,
          fill = "grey96",
          colour = "grey60",
          linewidth = 0.3) +
  geom_sf(data = plot_data,
          aes(colour = fold),
          size = 1.2) +
  facet_wrap(~species) +
  scale_colour_manual(values = jewel_cb_10, name = "Spatial fold") +
  coord_sf(
    xlim = c(530000, 1228000),
    ylim = c(366000, 1101000),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  filename = "figures/spatial_folds.png",
  plot = spatial_fold,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)



# Training data
data_pre2013 <- spatialised_sf %>%
  filter(Year < 2010) %>%
  mutate(row_id = row_number())

# Testing data
data_post2013 <- spatialised_sf %>%
  filter(Year > 2012)

# Create folds
set.seed(123)
folds <- rsample::vfold_cv(data_pre2013, v = 10, strata = "ZO")

# Assign fold numbers
data_pre2013$tempfold <- NA_integer_

for (i in seq_along(folds$splits)) {
  test_ids <- rsample::assessment(folds$splits[[i]])$row_id
  data_pre2013$tempfold[
    match(test_ids, data_pre2013$row_id)
  ] <- i
}

# Add panel labels
data_pre2013 <- data_pre2013 %>%
  mutate(
    panel = "Eelgrass training dataset (1993–2009)\n10 random folds",
    plot_group = factor(tempfold)
  )

data_post2013 <- data_post2013 %>%
  mutate(
    panel = "Eelgrass and surfgrass\ntesting dataset (2013–2023)",
    plot_group = "Testing"
  )

# Combine datasets
plot_data <- bind_rows(data_pre2013, data_post2013)

# Plot
temp_fold <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    colour = "grey60",
    linewidth = 0.3
  ) +
  geom_sf(data = data_pre2013,
          aes(colour = factor(tempfold)),
          size = 1.2) +
  geom_sf(data = data_post2013,
          colour = "#F58518",
          size = 1.2) +
  facet_wrap(~panel) +
  scale_colour_manual(values = okabe_ito_10, name = "Random fold")+
  coord_sf(
    xlim = c(530000, 1228000),
    ylim = c(366000, 1101000),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  filename = "figures/temp_folds.png",
  plot = temp_fold,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)




# make plots of eelgrass independent dataset 
load("code/output_data/model_results/eelgrass_independent_eval.RData")

# Separate into the two panels
obs_plot <- eelgrass_sf %>%
  filter(obs == 1) %>%
  mutate(
    panel = factor(
      "Eelgrass areal observations\n(2013–2023)",
      levels = c(
        "Eelgrass areal observations\n(2013–2023)",
        "Eelgrass\n10 pseudo-absence datasets")))

pa_plot <- eelgrass_sf %>%
  filter(obs == 0) %>%
  mutate(
    panel = factor(
      "Eelgrass\n10 pseudo-absence datasets",
      levels = c(
        "Eelgrass areal observations\n(2013–2023)",
        "Eelgrass\n10 pseudo-absence datasets")))

pa_plot <- pa_plot %>%
  slice_sample(prop = 1)

pa_plot_jitter <- pa_plot %>%
  mutate(
    geometry = st_jitter(geometry, amount = 40)  # 40 m
  )

indep_map_eel <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    colour = "grey60",
    linewidth = 0.3
  ) +
  geom_sf(
    data = obs_plot,
    colour = "#009E73",
    size = 1
  ) +
  geom_sf(
    data = pa_plot_jitter,
    aes(colour = factor(pa_set)),
    size = 1,
    alpha = 0.6
  ) +
  facet_wrap(~panel) +
  scale_colour_manual(values = okabe_ito_10,
    name = "Pseudoabsence set"
  ) +
  coord_sf(
    xlim = c(530000, 1228000),
    ylim = c(366000, 1101000),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )


ggsave(
  filename = "figures/independent_eel.png",
  indep_map_eel,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


# make plots of surfgrass independent dataset 
load("code/output_data/model_results/surfgrass_independent_eval.RData")

# Separate into the two panels
obs_plot_surf <- surfgrass_sf %>%
  filter(obs == 1) %>%
  mutate(
    panel = factor(
      "Surfgrass areal observations\n(2013–2024)",
      levels = c(
        "Surfgrass areal observations\n(2013–2024)",
        "Surfgrass\n10 pseudo-absence datasets")))

pa_plot_surf <- surfgrass_sf %>%
  filter(obs == 0) %>%
  mutate(
    panel = factor(
      "Surfgrass\n10 pseudo-absence datasets",
      levels = c(
        "Surfgrass areal observations\n(2013–2023)",
        "Surfgrass\n10 pseudo-absence datasets")))

pa_plot_surf <- pa_plot_surf %>%
  slice_sample(prop = 1)

pa_plot_jitter_surf <- pa_plot_surf %>%
  mutate(
    geometry = st_jitter(geometry, amount = 40)  # 40 m
  )

indep_map_surf <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    colour = "grey60",
    linewidth = 0.3
  ) +
  geom_sf(
    data = obs_plot_surf,
    colour = "#B22222",
    size = 1
  ) +
  geom_sf(
    data = pa_plot_jitter_surf,
    aes(colour = factor(pa_set)),
    size = 1,
    alpha = 0.6
  ) +
  facet_wrap(~panel) +
  scale_colour_manual(values = okabe_ito_10,
                      name = "Pseudoabsence set"
  ) +
  coord_sf(
    xlim = c(530000, 1228000),
    ylim = c(366000, 1101000),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )


ggsave(
  filename = "figures/independent_surf.png",
  indep_map_surf,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

surfgrass_sf %>%
  summarise(n_presences = sum(obs == 1, na.rm = TRUE))
eelgrass_sf %>%
  summarise(n_presences = sum(obs == 1, na.rm = TRUE))

surfgrass_sf %>%
  count(pa_set, .drop = FALSE)

eelgrass_sf %>%
  count(pa_set, .drop = FALSE)


# make plot of field validation
load("code/output_data/field_validation/validation_dataset.RData")
validation_sf <- validation_sf %>%
  relocate(PH_freq, PH, ZM_freq, ZM, .after = PC_PH)

# need to also drop sites that after reviewing the database should not have been included as they were aborted or vis was too bad to see species
validation_sf <- validation_sf %>%
  filter(!HKey %in% c(74, 84, 85, 161))

training_plot <- spatialised_sf %>%
  mutate(
    panel = "SDM training dataset\nwith targeted validation sites",
    plot_group = "Training"
  )

validation_overlay <- validation_sf %>%
  mutate(
    panel = "SDM training dataset\nwith targeted validation sites",
    plot_group = "Validation"
  )


validation_class <- validation_sf %>%
  mutate(
    habitat = case_when(
      ZM == 1 & PH == 1 ~ "Both",
      ZM == 1 ~ "Eelgrass",
      PH == 1 ~ "Surfgrass",
      TRUE ~ "Neither"
    ),
    panel = "Targeted validation sites\nby observed taxa"
  ) %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        "Neither",
        "Eelgrass",
        "Surfgrass",
        "Both"
      )
    )
  ) %>%
  arrange(habitat)

p1 <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    colour = "grey60",
    linewidth = 0.3
  ) +
  geom_sf(
    data = training_plot,
    aes(colour = "Training data"),
    size = 1
  ) +
  geom_sf(
    data = validation_overlay,
    aes(colour = "Targeted sites"),
    size = 1
  ) +
  scale_colour_manual(
    values = c(
      "Training data" = "grey75",
      "Targeted sites" = "#0072B2"
    ),
    name = NULL
  ) +
  facet_wrap(~panel) +
  coord_sf(
    xlim = c(530000, 1228000),
    ylim = c(366000, 1101000),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = c(0.82, 0.82),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )


# -------------------------------
# Plot 2: Validation classifications
# -------------------------------

p2 <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    colour = "grey60",
    linewidth = 0.3
  ) +
  geom_sf(
    data = validation_class,
    aes(colour = habitat),
    size = 1
  ) +
  facet_wrap(~panel) +
  scale_colour_manual(
    values = c(
      "Eelgrass" = "#009E73",
      "Surfgrass" = "#D55E00",
      "Both" = "#882255",
      "Neither" = "grey40"
    ),
    name = NULL
  ) +
  coord_sf(
    xlim = c(530000, 1228000),
    ylim = c(366000, 1101000),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = c(0.82, 0.82),   
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )

p2 <- p2 +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

# Combine panels
validation_map <- p1 + p2 +
  plot_layout(widths = c(1, 1))


validation_map

ggsave(
  filename = "figures/targeted.png",
  validation_map,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)



#relative variable importance
#eelgrass
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_bccm_gbm.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_bccm_nospatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_bccm_spatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_bccm_xgb.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_nep_gbm.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_nep_nospatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_nep_spatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_e_nep_xgb.RData")

glmm_varimp_bccm_eelgrass <- relimp_e_bccm_nospatial[[2]]
glmm_varimp_bccm_eelgrass <- glmm_varimp_bccm_eelgrass %>%
  rename(variable = term)%>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLM_bccm")

glmmspatial_varimp_bccm_eelgrass <- relimp_e_bccm_spatial[[2]]
glmmspatial_varimp_bccm_eelgrass <- glmmspatial_varimp_bccm_eelgrass %>%
  rename(variable = term) %>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLMM_bccm")

glmm_varimp_nep_eelgrass <- relimp_e_nep_nospatial[[2]]
glmm_varimp_nep_eelgrass <- glmm_varimp_nep_eelgrass %>%
  rename(variable = term) %>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLM_nep")

glmmspatial_varimp_nep_eelgrass <- relimp_e_nep_spatial[[2]]
glmmspatial_varimp_nep_eelgrass <- glmmspatial_varimp_nep_eelgrass %>%
  rename(variable = term) %>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLMM_nep")

gbm_varimp_bccm_eelgrass <- gbm_varimp_bccm_eelgrass %>%
  rename(variable = var,
         relimp = rel.inf) %>%
  mutate(model = "GBM_bccm")

gbm_varimp_nep_eelgrass <- gbm_varimp_nep_eelgrass %>%
  rename(variable = var,
         relimp = rel.inf) %>%
  mutate(model = "GBM_nep")

xgb_varimp_bccm_eelgrass <- xgb_varimp_bccm_eelgrass %>%
  select(Feature, RelImportance) %>%
  rename(variable = Feature, 
         relimp = RelImportance) %>%
  mutate(model = "XGBoost_bccm")

xgb_varimp_nep_eelgrass <- xgb_varimp_nep_eelgrass %>%
  select(Feature, RelImportance) %>%
  rename(variable = Feature, 
         relimp = RelImportance) %>%
  mutate(model = "XGBoost_nep")

eelgrass_relimp <- rbind(glmm_varimp_bccm_eelgrass, glmmspatial_varimp_bccm_eelgrass, glmm_varimp_nep_eelgrass, glmmspatial_varimp_nep_eelgrass,
                         gbm_varimp_bccm_eelgrass, gbm_varimp_nep_eelgrass, xgb_varimp_bccm_eelgrass, xgb_varimp_nep_eelgrass)


library(viridis)
eelgrass_relimp$model <- factor(
  eelgrass_relimp$model,
  levels = c("GLM_bccm", "GLMM_bccm", "GBM_bccm", "XGBoost_bccm", "GLM_nep",  "GLMM_nep", "GBM_nep", "XGBoost_nep")   
)

eelgrass_relimp <- eelgrass_relimp %>%
  mutate(variable = recode(variable,
                           "depth_stnd" = "Depth",
                           "substrate" = "Substrate",
                           "slope_stnd" = "Slope",
                           "prmin_stnd" = "Min precipitation",
                           "saltcv_nep_stnd" = "Subsurface salinity variability",
                           "airtempmin_stnd" = "Min air temperature",
                           "rsdsmin_stnd" = "Min surface downwelling shortwave flux",
                           "saltcv_bccm_stnd" = "Subsurface salinity variability",
                           "rei_stnd" = "Exposure",
                           "tempcv_nep_stnd" = "Subsurface temperature variability",
                           "tempcv_bccm_stnd" = "Subsurface temperature variability",
                           "NH4_bccm_stnd" = "Ammonium",
                           "NH4_nep_stnd" = "Ammonium",
                           "tempmin_bccm_stnd" = "Min subsurface temperature",
                           "Survey" = "Survey"
  ))



# ---- order variables by mean importance ----
eelgrass_relimp <- eelgrass_relimp %>%
  mutate(variable = reorder(variable, relimp, FUN = mean))

eelrelimp<-ggplot(eelgrass_relimp, aes(x = model, y = variable, fill = relimp)) +
  geom_tile(color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "viridis",
    trans = "sqrt",
    name = "Relative importance"
  ) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_blank(),
    axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
    axis.text.x.bottom = element_blank(),
    axis.ticks = element_blank()
  )
eelrelimp


#surfgrass
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_bccm_gbm.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_bccm_nospatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_bccm_spatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_bccm_xgb.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_nep_gbm.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_nep_nospatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_nep_spatial.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/relimp_s_nep_xgb.RData")

glmm_varimp_bccm_surfgrass <- relimp_s_bccm_nospatial[[2]]
glmm_varimp_bccm_surfgrass <- glmm_varimp_bccm_surfgrass %>%
  rename(variable = term)%>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLM_bccm")

glmmspatial_varimp_bccm_surfgrass <- relimp_s_bccm_spatial[[2]]
glmmspatial_varimp_bccm_surfgrass <- glmmspatial_varimp_bccm_surfgrass %>%
  rename(variable = term) %>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLMM_bccm")

glmm_varimp_nep_surfgrass <- relimp_s_nep_nospatial[[2]]
glmm_varimp_nep_surfgrass <- glmm_varimp_nep_surfgrass %>%
  rename(variable = term) %>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLM_nep")

glmmspatial_varimp_nep_surfgrass <- relimp_s_nep_spatial[[2]]
glmmspatial_varimp_nep_surfgrass <- glmmspatial_varimp_nep_surfgrass %>%
  rename(variable = term) %>%
  select(-mean_permuted_auc, -importance) %>%
  mutate(model = "GLMM_nep")

gbm_varimp_bccm_surfgrass <- gbm_varimp_bccm_surfgrass %>%
  rename(variable = var,
         relimp = rel.inf) %>%
  mutate(model = "GBM_bccm")

gbm_varimp_nep_surfgrass <- gbm_varimp_nep_surfgrass %>%
  rename(variable = var,
         relimp = rel.inf) %>%
  mutate(model = "GBM_nep")

xgb_varimp_bccm_surfgrass <- xgb_varimp_bccm_surfgrass %>%
  select(Feature, RelImportance) %>%
  rename(variable = Feature, 
         relimp = RelImportance) %>%
  mutate(model = "XGBoost_bccm")

xgb_varimp_nep_surfgrass <- xgb_varimp_nep_surfgrass %>%
  select(Feature, RelImportance) %>%
  rename(variable = Feature, 
         relimp = RelImportance) %>%
  mutate(model = "XGBoost_nep")


surfgrass_relimp <- rbind(glmm_varimp_bccm_surfgrass, glmmspatial_varimp_bccm_surfgrass, glmm_varimp_nep_surfgrass, glmmspatial_varimp_nep_surfgrass,
                         gbm_varimp_bccm_surfgrass, gbm_varimp_nep_surfgrass, xgb_varimp_bccm_surfgrass, xgb_varimp_nep_surfgrass)


library(viridis)
surfgrass_relimp$model <- factor(
  surfgrass_relimp$model,
  levels = c("GLM_bccm", "GLMM_bccm", "GBM_bccm", "XGBoost_bccm", "GLM_nep",  "GLMM_nep", "GBM_nep", "XGBoost_nep")    
)

surfgrass_relimp <- surfgrass_relimp %>%
  mutate(variable = recode(variable,
                           "depth_stnd" = "Depth",
                           "substrate" = "Substrate",
                           "tidal_sqrt_stnd" = "Tidal current",
                           "prmean_stnd" = "Precipitation",
                           "saltmean_nep_stnd" = "Subsurface salinity",
                           "saltmean_bccm_stnd" = "Subsurface salinity",
                           "airtempcv_stnd" = "Air temperature variability",
                           "rsdsmin_stnd" = "Min surface downwelling shortwave flux",
                           "rei_sqrt_stnd" = "Exposure",
                           "tempmean_nep_stnd" = "Subsurface temperature",
                           "tempmean_bccm_stnd" = "Subsurface temperature",
                           "surftempcv_bccm_stnd" = "Surface temperature variability",
                           "surftempcv_nep_stnd" = "Surface temperature variability",
                           "Survey" = "Survey"
  ))



# ---- order variables by mean importance ----
surfgrass_relimp <- surfgrass_relimp %>%
  mutate(variable = reorder(variable, relimp, FUN = mean))

surfrelimp<-ggplot(surfgrass_relimp, aes(x = model, y = variable, fill = relimp)) +
  geom_tile(color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "viridis",
    trans = "log1p",
    name = "Relative importance"
  ) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_blank(),
    axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
    axis.text.x.bottom = element_blank(),
    axis.ticks = element_blank()
  )

eelrelimp <- eelrelimp + theme(plot.margin = margin(5, 5, 5, 5))
surfrelimp <- surfrelimp + theme(plot.margin = margin(5, 5, 5, 5))

lims <- range(c(eelgrass_relimp$relimp, surfgrass_relimp$relimp), na.rm = TRUE)

eelrelimp <- eelrelimp +
  scale_fill_viridis_c(option = "viridis", trans = "log1p",
                       name = "Relative importance", limits = lims)

surfrelimp <- surfrelimp +
  scale_fill_viridis_c(option = "viridis", trans = "log1p",
                       name = "Relative importance", limits = lims)

surfrelimp <- surfrelimp +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks.x = element_blank()
  )

relimp<- (eelrelimp / surfrelimp) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = list(c("Eelgrass", "Surfgrass"))
  ) &
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 12, face = "bold")
  )
relimp
ggsave("figures/figure_relimp.png", plot = relimp, width = 8, height = 6, dpi = 300)





# probability of occurence compared to test sites, independent data, and field validation data with all thresholds
# just do best and worst models for figure
# eelgrass best model is XGBoost_nep and worst is BGM_nep
# surfgrass best model is and worse is
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/final_metrics_eelgrass.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/final_metrics_surfgrass.RData")


#need to extract predictions from the 2013-2023 test data
load("code/output_data/seagrass_model_inputs.RData")
seagrass_data %>%
  count(Month)

seagrass_data %>%
  summarise(
    pct_may_sep = 100 * mean(Month %in% 5:9)
  )

seagrass_data <- seagrass_data %>% filter (Year >= 2013) 
seagrass_data_sf <- sf::st_as_sf(seagrass_data, coords = c("X_m", "Y_m"), crs = 3005) 

eelgrass_xgb_nep <- terra::vrt(c("raster/eelgrass/xgb_nep/eelgrass_predictions_hg_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ncc_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_qcs_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_ss_xgb_nep.tif", "raster/eelgrass/xgb_nep/eelgrass_predictions_wcvi_xgb_nep.tif"), "eelgrass_xgb_nep.vrt", overwrite=T)   # values 0–1
eelgrass_gbm_nep <- terra::vrt(c("raster/eelgrass/gbm_nep/eelgrass_predictions_hg_gbm_nep.tif", "raster/eelgrass/gbm_nep/eelgrass_predictions_ncc_gbm_nep.tif", "raster/eelgrass/gbm_nep/eelgrass_predictions_qcs_gbm_nep.tif", "raster/eelgrass/gbm_nep/eelgrass_predictions_ss_gbm_nep.tif", "raster/eelgrass/gbm_nep/eelgrass_predictions_wcvi_gbm_nep.tif"), "eelgrass_gbm_nep.vrt", overwrite=T)   # values 0–1
surfgrass_xgb_nep <- terra::vrt(c("raster/surfgrass/xgb_nep/surfgrass_predictions_hg_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_ncc_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_qcs_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_ss_xgb_nep.tif", "raster/surfgrass/xgb_nep/surfgrass_predictions_wcvi_xgb_nep.tif"), "surfgrass_xgb_nep.vrt", overwrite=T)   # values 0–1
surfgrass_gbm_bccm <- terra::vrt(c("raster/surfgrass/gbm_bccm/surfgrass_predictions_hg_gbm_bccm.tif", "raster/surfgrass/gbm_bccm/surfgrass_predictions_ncc_gbm_bccm.tif", "raster/surfgrass/gbm_bccm/surfgrass_predictions_qcs_gbm_bccm.tif", "raster/surfgrass/gbm_bccm/surfgrass_predictions_ss_gbm_bccm.tif", "raster/surfgrass/gbm_bccm/surfgrass_predictions_wcvi_gbm_bccm.tif"), "surfgrass_gbm_bccm.vrt", overwrite=T)   # values 0–1

pred_stack <- c(
  eelgrass_xgb_nep,
  eelgrass_gbm_nep,
  surfgrass_xgb_nep,
  surfgrass_gbm_bccm
)

names(pred_stack) <- c(
  "XGBoost_nep_eel",
  "GBM_nep_eel",
  "XGBoost_nep_surf",
  "GBM_bccm_surf")

prediction_extract <- terra::extract(
  pred_stack,
  terra::vect(seagrass_data_sf)
) %>%
  select(-ID)

seagrass_data_sf <- seagrass_data_sf %>%
  bind_cols(prediction_extract) %>%
  filter(if_all(all_of(names(pred_stack)), ~ !is.na(.)))
summary(seagrass_data_sf)

seagrass_data_sf <- seagrass_data_sf %>%
  mutate(
    Tidal_zone = if_else(depth > -1, "Subtidal", "Intertidal"),
    Tidal_zone = factor(Tidal_zone, levels = c("Subtidal", "Intertidal"))
  )

test_eel_pa_long <- seagrass_data_sf %>%
  sf::st_drop_geometry() %>%
  mutate(
    Presence = factor(ifelse(ZO == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    XGBoost_nep_eel,
    GBM_nep_eel
  ) %>%
  pivot_longer(
    cols = -c(Presence, Tidal_zone),
    names_to = "model",
    values_to = "predicted_suitability"
  ) %>%
  mutate(
    model = recode(
      model,
      "XGBoost_nep_eel"    = "XGBoost_nep",
      "GBM_nep_eel"        = "GBM_nep"  )
  )

test_surf_pa_long <- seagrass_data_sf %>%
  sf::st_drop_geometry() %>%
  mutate(
    Presence = factor(ifelse(PH == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    XGBoost_nep_surf,
    GBM_bccm_surf
    
  ) %>%
  pivot_longer(
    cols = -c(Presence, Tidal_zone),
    names_to = "model",
    values_to = "predicted_suitability"
  ) %>%
  mutate(
    model = recode(
      model,
      "XGBoost_nep_surf"    = "XGBoost_nep",
      "GBM_bccm_surf"        = "GBM_bccm"   )
  )


#independent
load("code/output_data/model_results/surfgrass_independent_eval.RData")
load("code/output_data/model_results/eelgrass_independent_eval.RData")
#need to add bathy data to independent dataset to be able to classify it by tide height

indep_eelgrass_sf <- eelgrass_sf %>%
  mutate(
    Tidal_zone = if_else(bathy > -1, "Subtidal", "Intertidal"),
    Tidal_zone = factor(Tidal_zone, levels = c("Subtidal", "Intertidal"))
  ) %>%
  distinct(obs, geometry, .keep_all = TRUE)

indep_eel_pa_long <- indep_eelgrass_sf %>%  
  sf::st_drop_geometry() %>%
  mutate(
    Presence = factor(ifelse(obs == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    XGBoost_nep,
    GBM_nep
  ) %>%
  pivot_longer(
    cols = -c(Presence, Tidal_zone),
    names_to = "model",
    values_to = "predicted_suitability"
  ) 


indep_surfgrass_sf <- surfgrass_sf %>%
  mutate(
    Tidal_zone = if_else(bathy > -1, "Subtidal", "Intertidal"),
    Tidal_zone = factor(Tidal_zone, levels = c("Subtidal", "Intertidal"))
  ) %>%
  distinct(obs, geometry, .keep_all = TRUE)

indep_surf_pa_long <- indep_surfgrass_sf %>%  
  sf::st_drop_geometry() %>%
  mutate(
    Presence = factor(ifelse(obs == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    XGBoost_nep,
    GBM_bccm
    
  ) %>%
  pivot_longer(
    cols = -c(Presence, Tidal_zone),
    names_to = "model",
    values_to = "predicted_suitability"
  ) 


#field
load("code/output_data/field_validation/validation_dataset.RData")

# need to also drop sites that after reviewing the database should not have been included as they were aborted or vis was too bad to see species
validation_sf <- validation_sf %>% sf::st_drop_geometry() %>%
  filter(!HKey %in% c(74, 84, 85, 161))  

validation_sf <- validation_sf %>% select (HKey, Region, Survey, avgCorDepth_obs, bathy_mod,
                                           ZM,  zo_XGBOOST_nep_pred, zo_GBM_nep_pred,
                                           PH, ph_XGBOOST_nep_pred, ph_GBM_bccm_pred) 
validation_sf <- validation_sf %>%
  mutate(
    Tidal_zone = if_else(bathy_mod > -1, "Subtidal", "Intertidal"),
    Tidal_zone = factor(Tidal_zone, levels = c("Subtidal", "Intertidal"))
  )

# Prepare eelgrass validation data
field_eelgrass_pa_long <- validation_sf %>%
  mutate(
    Presence = factor(ifelse(ZM == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    zo_XGBOOST_nep_pred,
    zo_GBM_nep_pred
  ) %>%
  pivot_longer(
    cols = -c(Presence, Tidal_zone),
    names_to = "model",
    values_to = "predicted_suitability"
  ) %>%
  mutate(
    model = recode(
      model,
      "zo_XGBOOST_nep_pred"    = "XGBoost_nep",
      "zo_GBM_nep_pred"        = "GBM_nep"  )
    )


field_surfgrass_pa_long <- validation_sf %>%
  mutate(
    Presence = factor(ifelse(PH == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    ph_XGBOOST_nep_pred,
    ph_GBM_bccm_pred
  ) %>%
  pivot_longer(
    cols = -c(Presence, Tidal_zone),
    names_to = "model",
    values_to = "predicted_suitability"
  ) %>%
  mutate(
    model = recode(
      model,
      "ph_XGBOOST_nep_pred"    = "XGBoost_nep",
      "ph_GBM_bccm_pred"       = "GBM_bccm" )
  )



best_model <- "XGBoost_nep"
worst_model <- "GBM_nep"

training <- test_eel_pa_long %>%
  filter(model == best_model) %>%
  mutate(
    Dataset = "SDM building",
    Algorithm = "Highest performing"
  )

training_worst <- test_eel_pa_long %>%
  filter(model == worst_model) %>%
  mutate(
    Dataset = "SDM building",
    Algorithm = "Lowest performing"
  )

independent <- indep_eel_pa_long %>%
  filter(model == best_model) %>%
  mutate(
    Dataset = "Independent",
    Algorithm = "Highest performing"
  )

independent_worst <- indep_eel_pa_long %>%
  filter(model == worst_model) %>%
  mutate(
    Dataset = "Independent",
    Algorithm = "Lowest performing"
  )

field <- field_eelgrass_pa_long %>%
  filter(model == best_model) %>%
  mutate(
    Dataset = "Targeted",
    Algorithm = "Highest performing"
  )

field_worst <- field_eelgrass_pa_long %>%
  filter(model == worst_model) %>%
  mutate(
    Dataset = "Targeted",
    Algorithm = "Lowest performing"
  )

thresholds <- final_metrics_eelgrass %>%
  filter(model %in% c("XGBoost_nep", "GBM_nep")) %>%
  select(
    model,
    threshold_spatial,
    threshold_temporal,
    threshold_independent,
    threshold_field
  ) %>%
  pivot_longer(
    cols = starts_with("threshold"),
    names_to = "Threshold_type",
    values_to = "threshold"
  ) %>%
  mutate(
    Algorithm = ifelse(model == "XGBoost_nep", "Highest performing", "Lowest performing"),
    Dataset = case_when(
      Threshold_type %in% c("threshold_spatial",
                            "threshold_temporal") ~ "SDM building",
      Threshold_type == "threshold_independent" ~ "Independent",
      Threshold_type == "threshold_field" ~ "Targeted"
    ),
    Threshold_type = case_when(
      Threshold_type == "threshold_spatial" ~ "Spatial",
      Threshold_type == "threshold_temporal" ~ "Temporal",
      Threshold_type == "threshold_independent" ~ "Independent",
      Threshold_type == "threshold_field" ~ "Targeted"
    )
  )

thresholds

thresholds <- thresholds %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = c("SDM building", "Independent", "Targeted")
    )
  )

plot_dat <-
  bind_rows(
    training,
    independent,
    field,
    training_worst,
    independent_worst,
    field_worst
  )

plot_dat <- plot_dat %>%
  mutate(
    Group = interaction(Tidal_zone, Presence, sep = "\n"),
    Dataset = factor(
      Dataset,
      levels = c("SDM building","Independent","Targeted")
    ),
    Algorithm = factor(
      Algorithm,
      levels = c("Highest performing","Lowest performing")
    )
  )

n_labels <- bind_rows(
  test_eel_pa_long %>%
    filter(model == best_model) %>%
    mutate(Dataset = "SDM building"),
  
  indep_eel_pa_long %>%
    filter(model == best_model) %>%
    mutate(Dataset = "Independent"),
  
  field_eelgrass_pa_long %>%
    filter(model == best_model) %>%
    mutate(Dataset = "Targeted")
) %>%
  count(Dataset, Presence, Tidal_zone, name = "n") %>%
  mutate(
    Algorithm = "Highest performing",
    label = paste0("n = ", n),
    y = 0.74
  )

n_labels <- n_labels %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = c("SDM building", "Independent", "Targeted")
    ),
    Algorithm = factor(
      Algorithm,
      levels = c("Highest performing", "Lowest performing")
    )
  )

pd <- position_dodge(width = 0.9)

panel_labels <- data.frame(
  Dataset = factor(
    c("SDM building", "Independent", "Targeted"),
    levels = levels(plot_dat$Dataset)
  ),
  Algorithm = levels(plot_dat$Algorithm)[1],  # top row
  Presence = factor("Absent", levels = levels(plot_dat$Presence)),
  y = 0.995,
  label = c("A", "B", "C")
)

eelgrass <- ggplot(plot_dat,
       aes(x = Presence,
           y = predicted_suitability,
           fill = Tidal_zone)) +
  geom_text(
    data = panel_labels,
    aes(
      x = Presence,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 6,
    hjust = 2.8,
    vjust = 1.1
  ) +
  geom_boxplot(
    aes(group = interaction(Presence, Tidal_zone)),
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_hline(
    data = thresholds,
    aes(
      yintercept = threshold,
      linetype = Threshold_type
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.5
  ) +
  facet_grid(
    Algorithm ~ Dataset
  ) +
  geom_text(
    data = n_labels,
    aes(
      x = Presence,
      y = y,
      group = Tidal_zone,
      label = label
    ),
    position = position_dodge2(width = 0.75),
    hjust = -0.08,
    angle = 30,
    size = 2.5,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      "Subtidal" = "#0B5D5E",   # dark teal
      "Intertidal" = "#6EC6C4"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Spatial" = "solid",
      "Temporal" = "dashed",
      "Independent" = "dotdash",
      "Targeted" = "twodash"
    ),
    name = "Validation threshold"
  ) +
  labs(
    x = NULL,
    y = "Predicted suitability",
    fill = "Tidal zone"
  ) +
  
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  
  labs(
    title = "Eelgrass",
    x = NULL,
    y = "Relative probability of occurrence"
  ) + 
  
  theme_bw() +
  
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(), 
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    panel.spacing.y = unit(1, "lines"),
    legend.position = "none"
  )

  eelgrass


  
  
  
  
  best_model_surf <- "XGBoost_nep"
  worst_model_surf <- "GBM_bccm"
  
  training_surf <- test_surf_pa_long %>%
    filter(model == best_model_surf) %>%
    mutate(
      Dataset = "SDM building",
      Algorithm = "Highest performing"
    )
  
  training_worst_surf <- test_surf_pa_long %>%
    filter(model == worst_model_surf) %>%
    mutate(
      Dataset = "SDM building",
      Algorithm = "Lowest performing"
    )
  
  independent_surf <- indep_surf_pa_long %>%
    filter(model == best_model_surf) %>%
    mutate(
      Dataset = "Independent",
      Algorithm = "Highest performing"
    )
  
  independent_worst_surf <- indep_surf_pa_long %>%
    filter(model == worst_model_surf) %>%
    mutate(
      Dataset = "Independent",
      Algorithm = "Lowest performing"
    )
  
  field_surf <- field_surfgrass_pa_long %>%
    filter(model == best_model_surf) %>%
    mutate(
      Dataset = "Targeted",
      Algorithm = "Highest performing"
    )
  
  field_worst_surf <- field_surfgrass_pa_long %>%
    filter(model == worst_model_surf) %>%
    mutate(
      Dataset = "Targeted",
      Algorithm = "Lowest performing"
    )
  
  thresholds_surf <- final_metrics_surfgrass %>%
    filter(model %in% c("XGBoost_nep", "GBM_bccm")) %>%
    select(
      model,
      threshold_spatial,
      threshold_temporal,
      threshold_independent,
      threshold_field
    ) %>%
    pivot_longer(
      cols = starts_with("threshold"),
      names_to = "Threshold_type",
      values_to = "threshold"
    ) %>%
    mutate(
      Algorithm = ifelse(model == "XGBoost_nep", "Highest performing", "Lowest performing"),
      Dataset = case_when(
        Threshold_type %in% c("threshold_spatial",
                              "threshold_temporal") ~ "SDM building",
        Threshold_type == "threshold_independent" ~ "Independent",
        Threshold_type == "threshold_field" ~ "Targeted"
      ),
      Threshold_type = case_when(
        Threshold_type == "threshold_spatial" ~ "Spatial",
        Threshold_type == "threshold_temporal" ~ "Temporal",
        Threshold_type == "threshold_independent" ~ "Independent",
        Threshold_type == "threshold_field" ~ "Targeted"
      )
    )
  
  thresholds_surf
  
  thresholds_surf <- thresholds_surf %>%
    mutate(
      Dataset = factor(
        Dataset,
        levels = c("SDM building", "Independent", "Targeted")
      )
    )
  
  plot_dat_surf <-
    bind_rows(
      training_surf,
      independent_surf,
      field_surf,
      training_worst_surf,
      independent_worst_surf,
      field_worst_surf
    )
  
  plot_dat_surf <- plot_dat_surf %>%
    mutate(
      Group = interaction(Tidal_zone, Presence, sep = "\n"),
      Dataset = factor(
        Dataset,
        levels = c("SDM building","Independent","Targeted")
      ),
      Algorithm = factor(
        Algorithm,
        levels = c("Highest performing","Lowest performing")
      )
    )
  
  n_labels_surf <- bind_rows(
    test_surf_pa_long %>%
      filter(model == best_model_surf) %>%
      mutate(Dataset = "SDM building"),
    
    indep_surf_pa_long %>%
      filter(model == best_model_surf) %>%
      mutate(Dataset = "Independent"),
    
    field_surfgrass_pa_long %>%
      filter(model == best_model_surf) %>%
      mutate(Dataset = "Targeted")
  ) %>%
    count(Dataset, Presence, Tidal_zone, name = "n") %>%
    mutate(
      Algorithm = "Highest performing",
      label = paste0("n = ", n),
      y = 0.74
    )
  
  n_labels_surf <- n_labels_surf %>%
    mutate(
      Dataset = factor(
        Dataset,
        levels = c("SDM building", "Independent", "Targeted")
      ),
      Algorithm = factor(
        Algorithm,
        levels = c("Highest performing", "Lowest performing")
      )
    )
  
  panel_labels_surf <- data.frame(
    Dataset = factor(
      c("SDM building", "Independent", "Targeted"),
      levels = levels(plot_dat_surf$Dataset)
    ),
    Algorithm = levels(plot_dat_surf$Algorithm)[1],  # top row
    Presence = factor("Absent", levels = levels(plot_dat_surf$Presence)),
    y = 0.995,
    label = c("D", "E", "F")
  )
  
      
    
  
  surfgrass <- ggplot(plot_dat_surf,
                     aes(x = Presence,
                         y = predicted_suitability,
                         fill = Tidal_zone)) +
    geom_text(
      data = panel_labels_surf,
      aes(
        x = Presence,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 6,
      hjust = 2.8,
      vjust = 1.1
    ) +
    geom_boxplot(
      aes(group = interaction(Presence, Tidal_zone)),
      position = position_dodge(width = 0.75),
      outlier.shape = NA,
      width = 0.6
    ) +
    
    geom_hline(
      data = thresholds_surf,
      aes(
        yintercept = threshold,
        linetype = Threshold_type
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.5
    ) +
    facet_grid(
      Algorithm ~ Dataset
    ) +
    geom_text(
      data = n_labels_surf,
      aes(
        x = Presence,
        y = y,
        label = label,
        group = Tidal_zone
      ),
      position = position_dodge2(width = 0.75),
      hjust = -0.08,
      angle = 30,
      size = 2.5,
      inherit.aes = FALSE
    ) +
    
    scale_fill_manual(
      values = c(
        "Subtidal" = "#0B5D5E",   # dark teal
        "Intertidal" = "#6EC6C4"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Spatial" = "solid",
        "Temporal" = "dashed",
        "Independent" = "dotdash",
        "Targeted" = "twodash"
      ),
      name = "Validation threshold"
    ) +
    labs(
      x = NULL,
      y = "Predicted suitability",
      fill = "Tidal zone"
    ) +
    
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    
    labs(
      title = "Surfgrass",
      x = NULL,
      y = "Relative probability of occurrence"
    ) + 
    
    theme_bw() +
    
    theme(
      strip.background = element_blank(),
      panel.grid = element_blank(), 
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 10),
      panel.spacing.y = unit(1, "lines"),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(
      linetype = guide_legend(
        title = "Validation\nthreshold",
        nrow = 2,
        byrow = TRUE
      )
    )
  
  surfgrass
  
  
  
 
  
  occurrence_plot <- eelgrass / surfgrass +
    plot_layout(
      heights = c(3, 3.5))
  
  occurrence_plot   
 
  ggsave(
    filename = "figures/occurrence_plots.png",
    plot = occurrence_plot,
    width = 180,
    height = 210,
    units = "mm",
    dpi = 600,
    bg = "white" )
  
  
  
  

 
  
  
  #Binary Prediction maps for eelgrass and surfgrass
 
  eelgrass<- rast("raster/eelgrass/eelgrass_predictions.tif")
  eelgrass_plot <- ifel(eelgrass == 1, 1, NA)
  
  
  surfgrass<- rast("raster/surfgrass/surfgrass_predictions.tif")
  surfgrass_plot <- ifel(surfgrass == 1, 1, NA)
  
  eelgrass_poly <- as.polygons(eelgrass_plot, dissolve = TRUE)
  surfgrass_poly <- as.polygons(surfgrass_plot, dissolve = TRUE)
  
  plot <- ggplot() +
    geom_sf(
      data = coastline,
      fill = "grey65",
      colour = "grey35",
      linewidth = 0.1
    ) +
    geom_spatvector(
      data = eelgrass_poly,
      aes(fill = "Eelgrass"),
      colour = "#00A676",
      linewidth = 0.05
    ) +
    geom_spatvector(
      data = surfgrass_poly,
      aes(fill = "Surfgrass"),
      colour = "#E66101",
      linewidth = 0.05
    ) +
    scale_fill_manual(
      name = NULL,
      values = c(
        "Eelgrass" = "#00A676",
        "Surfgrass" = "#E66101"
      )
    ) +
    coord_sf(
      xlim = c(530000, 1228000),
      ylim = c(366000, 1101000),
      expand = FALSE
    ) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "grey88", colour = NA),
      panel.grid = element_blank(),
      legend.position = c(0.97, 0.97),   # top right
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha("white", 0.8), colour = "grey50"),
      legend.text = element_text(size = 9),
      legend.key = element_blank()
    )
   
plot 
  
ggsave(
  filename = "figures/distribution.png",
  plot = plot,
  width = 180,
  height = 210,
  units = "mm",
  dpi = 600,
  bg = "white" )


  
 #make figures of environmental characteristics

 
 #load seagrass data
 load("code/output_data/seagrass_model_inputs.RData")
 
 bathy_hg <- rast("raw_data/envlayers-20m-hg//bathymetry.tif")
 bathy_ncc <- rast("raw_data/envlayers-20m-ncc/bathymetry.tif")
 bathy_qcs <- rast("raw_data/envlayers-20m-qcs/bathymetry.tif")
 bathy_wcvi <- rast("raw_data/envlayers-20m-wcvi/bathymetry.tif")
 bathy_ss <- rast("raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif")
 
 tidal_all <- vrt(c("raw_data/current_20m/Nearshore_CurrentSpeedIndex.tif"))
 names(tidal_all)<-"tidal"
 #change to index 0-1 scale
 tidal_index_all <- tidal_all/(maxFn(tidal_all))
 crs(tidal_index_all) <- "EPSG:3005"
 
 rei_hg <- rast("raw_data/REI/rei_hg.tif")
 rei_ncc <- rast("raw_data/REI/rei_ncc.tif")
 rei_qcs <- rast("raw_data/REI/rei_qcs.tif")
 rei_wcvi <- rast("raw_data/REI/rei_wcvi.tif")
 rei_ss <- rast("raw_data/REI/rei_sog.tif")
 
 folder_2013_2023 <- "code/output_data/processed_ocean_variables/years_2013-2023"
 files_2013_2023 <- list.files(folder_2013_2023, pattern = "\\.tif$", full.names = TRUE)
 selected_files_2013_2023 <- files_2013_2023[c(5, 6, 17, 18, 22, 23, 24, 25, 26, 29, 30, 35, 36, 39, 44, 45, 46)]
 hindcast2013_2023 <- terra::rast(selected_files_2013_2023)
 
 folder_1993_2002 <- "code/output_data/processed_ocean_variables/years_1993-2002"
 files_1993_2002 <- list.files(folder_1993_2002, pattern = "\\.tif$", full.names = TRUE)
 selected_files_1993_2002 <- files_1993_2002[c(5, 6, 17, 18, 22, 23, 24, 25, 26, 29, 30, 35, 36, 39, 44, 45, 46)]
 hindcast1993_2002 <- terra::rast(selected_files_1993_2002)
 
 folder_2003_2012 <- "code/output_data/processed_ocean_variables/years_2003-2012"
 files_2003_2012 <- list.files(folder_2003_2012, pattern = "\\.tif$", full.names = TRUE)
 selected_files_2003_2012 <- files_2003_2012[c(5, 6, 17, 18, 22, 23, 24, 25, 26, 29, 30, 35, 36, 39, 44, 45, 46)]
 hindcast2003_2012 <- terra::rast(selected_files_2003_2012)
 
 
max_depth <- quantile(seagrass_data$depth, probs = c(0.99))
 min_depth <- quantile(seagrass_data$depth, probs = c(0.001))
 max_rei <- quantile(seagrass_data$rei, probs = c(0.9999))
 max_tidal <- quantile(seagrass_data$tidal, probs = c(0.9999))

 #haida gwaii
 env_20m_hg <- as.data.frame(bathy_hg, xy=TRUE)
 names(env_20m_hg) <- c("X_m", "Y_m", "depth")
 env_20m_hg <- env_20m_hg %>% filter(depth <= max_depth, depth >= min_depth)
 env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = rei_hg, y = env_20m_hg_sf)
 env_20m_hg$rei <- hold$rei 
 env_20m_hg <- env_20m_hg %>% filter(rei <= max_rei)
 env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m","Y_m"), crs = "EPSG:3005")
hold <- terra::extract(x = tidal_index_all, y = env_20m_hg_sf)
 env_20m_hg$tidal <- hold$tidal 
 env_20m_hg <- env_20m_hg %>% filter(tidal <= max_tidal)
 env_20m_hg_sf <- st_as_sf(env_20m_hg, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = hindcast2013_2023, y = env_20m_hg_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_hg$precip_mean2013_2023 <- hold$precip_mean
 env_20m_hg$precip_min2013_2023 <- hold$precip_min
 env_20m_hg$salt_5m_cv_bccm2013_2023 <- hold$salt_5m_cv_bccm
 env_20m_hg$salt_5m_cv_nep362013_2023 <- hold$salt_5m_cv_nep36
 env_20m_hg$salt_5m_mean_bccm2013_2023 <- hold$salt_5m_mean_bccm
 env_20m_hg$salt_5m_mean_nep362013_2023 <- hold$salt_5m_mean_nep36
 env_20m_hg$temp_air_min2013_2023 <- hold$temp_air_min
 env_20m_hg$rsds_min2013_2023 <- hold$rsds_min
 env_20m_hg$temp_5m_cv_bccm2013_2023 <- hold$temp_5m_cv_bccm
 env_20m_hg$temp_5m_cv_nep362013_2023 <- hold$temp_5m_cv_nep36
 env_20m_hg$temp_5m_mean_bccm2013_2023 <- hold$temp_5m_mean_bccm
 env_20m_hg$temp_5m_mean_nep362013_2023 <- hold$temp_5m_mean_nep36
 env_20m_hg$NH4_5m_mean_bccm2013_2023 <- hold$NH4_5m_mean_bccm
 env_20m_hg$NH4_5m_mean_nep362013_2023 <- hold$NH4_5m_mean_nep36
 env_20m_hg$temp_air_cv2013_2023 <- hold$temp_air_cv
 env_20m_hg$temp_s_cv_bccm2013_2023 <- hold$temp_s_cv_bccm
 env_20m_hg$temp_s_cv_nep362013_2023 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast1993_2002, y = env_20m_hg_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_hg$precip_mean1993_2002 <- hold$precip_mean
 env_20m_hg$precip_min1993_2002 <- hold$precip_min
 env_20m_hg$salt_5m_cv_bccm1993_2002 <- hold$salt_5m_cv_bccm
 env_20m_hg$salt_5m_cv_nep361993_2002 <- hold$salt_5m_cv_nep36
 env_20m_hg$salt_5m_mean_bccm1993_2002 <- hold$salt_5m_mean_bccm
 env_20m_hg$salt_5m_mean_nep361993_2002 <- hold$salt_5m_mean_nep36
 env_20m_hg$temp_air_min1993_2002 <- hold$temp_air_min
 env_20m_hg$rsds_min1993_2002 <- hold$rsds_min
 env_20m_hg$temp_5m_cv_bccm1993_2002 <- hold$temp_5m_cv_bccm
 env_20m_hg$temp_5m_cv_nep361993_2002 <- hold$temp_5m_cv_nep36
 env_20m_hg$temp_5m_mean_bccm1993_2002 <- hold$temp_5m_mean_bccm
 env_20m_hg$temp_5m_mean_nep361993_2002 <- hold$temp_5m_mean_nep36
 env_20m_hg$NH4_5m_mean_bccm1993_2002 <- hold$NH4_5m_mean_bccm
 env_20m_hg$NH4_5m_mean_nep361993_2002 <- hold$NH4_5m_mean_nep36
 env_20m_hg$temp_air_cv1993_2002 <- hold$temp_air_cv
 env_20m_hg$temp_s_cv_bccm1993_2002 <- hold$temp_s_cv_bccm
 env_20m_hg$temp_s_cv_nep361993_2002 <- hold$temp_s_cv_nep36

 hold <- terra::extract(x = hindcast2003_2012, y = env_20m_hg_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_hg$precip_mean2003_2012 <- hold$precip_mean
 env_20m_hg$precip_min2003_2012 <- hold$precip_min
 env_20m_hg$salt_5m_cv_bccm2003_2012 <- hold$salt_5m_cv_bccm
 env_20m_hg$salt_5m_cv_nep362003_2012 <- hold$salt_5m_cv_nep36
 env_20m_hg$salt_5m_mean_bccm2003_2012 <- hold$salt_5m_mean_bccm
 env_20m_hg$salt_5m_mean_nep362003_2012 <- hold$salt_5m_mean_nep36
 env_20m_hg$temp_air_min2003_2012 <- hold$temp_air_min
 env_20m_hg$rsds_min2003_2012 <- hold$rsds_min
 env_20m_hg$temp_5m_cv_bccm2003_2012 <- hold$temp_5m_cv_bccm
 env_20m_hg$temp_5m_cv_nep362003_2012 <- hold$temp_5m_cv_nep36
 env_20m_hg$temp_5m_mean_bccm2003_2012 <- hold$temp_5m_mean_bccm
 env_20m_hg$temp_5m_mean_nep362003_2012 <- hold$temp_5m_mean_nep36
 env_20m_hg$NH4_5m_mean_bccm2003_2012 <- hold$NH4_5m_mean_bccm
 env_20m_hg$NH4_5m_mean_nep362003_2012 <- hold$NH4_5m_mean_nep36
 env_20m_hg$temp_air_cv2003_2012 <- hold$temp_air_cv
 env_20m_hg$temp_s_cv_bccm2003_2012 <- hold$temp_s_cv_bccm
 env_20m_hg$temp_s_cv_nep362003_2012 <- hold$temp_s_cv_nep36
 

 
 #north central coast
 env_20m_ncc <- as.data.frame(bathy_ncc, xy=TRUE)
 names(env_20m_ncc) <- c("X_m", "Y_m", "depth")
 env_20m_ncc <- env_20m_ncc %>% filter(depth <= max_depth, depth >= min_depth)
 env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = rei_ncc, y = env_20m_ncc_sf)
 env_20m_ncc$rei <- hold$rei 
 env_20m_ncc <- env_20m_ncc %>% filter(rei <= max_rei)
 env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = tidal_index_all, y = env_20m_ncc_sf)
 env_20m_ncc$tidal <- hold$tidal 
 env_20m_ncc <- env_20m_ncc %>% filter(tidal <= max_tidal)
 env_20m_ncc_sf <- st_as_sf(env_20m_ncc, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = hindcast2013_2023, y = env_20m_ncc_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_ncc$precip_mean2013_2023 <- hold$precip_mean
 env_20m_ncc$precip_min2013_2023 <- hold$precip_min
 env_20m_ncc$salt_5m_cv_bccm2013_2023 <- hold$salt_5m_cv_bccm
 env_20m_ncc$salt_5m_cv_nep362013_2023 <- hold$salt_5m_cv_nep36
 env_20m_ncc$salt_5m_mean_bccm2013_2023 <- hold$salt_5m_mean_bccm
 env_20m_ncc$salt_5m_mean_nep362013_2023 <- hold$salt_5m_mean_nep36
 env_20m_ncc$temp_air_min2013_2023 <- hold$temp_air_min
 env_20m_ncc$rsds_min2013_2023 <- hold$rsds_min
 env_20m_ncc$temp_5m_cv_bccm2013_2023 <- hold$temp_5m_cv_bccm
 env_20m_ncc$temp_5m_cv_nep362013_2023 <- hold$temp_5m_cv_nep36
 env_20m_ncc$temp_5m_mean_bccm2013_2023 <- hold$temp_5m_mean_bccm
 env_20m_ncc$temp_5m_mean_nep362013_2023 <- hold$temp_5m_mean_nep36
 env_20m_ncc$NH4_5m_mean_bccm2013_2023 <- hold$NH4_5m_mean_bccm
 env_20m_ncc$NH4_5m_mean_nep362013_2023 <- hold$NH4_5m_mean_nep36
 env_20m_ncc$temp_air_cv2013_2023 <- hold$temp_air_cv
 env_20m_ncc$temp_s_cv_bccm2013_2023 <- hold$temp_s_cv_bccm
 env_20m_ncc$temp_s_cv_nep362013_2023 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast1993_2002, y = env_20m_ncc_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_ncc$precip_mean1993_2002 <- hold$precip_mean
 env_20m_ncc$precip_min1993_2002 <- hold$precip_min
 env_20m_ncc$salt_5m_cv_bccm1993_2002 <- hold$salt_5m_cv_bccm
 env_20m_ncc$salt_5m_cv_nep361993_2002 <- hold$salt_5m_cv_nep36
 env_20m_ncc$salt_5m_mean_bccm1993_2002 <- hold$salt_5m_mean_bccm
 env_20m_ncc$salt_5m_mean_nep361993_2002 <- hold$salt_5m_mean_nep36
 env_20m_ncc$temp_air_min1993_2002 <- hold$temp_air_min
 env_20m_ncc$rsds_min1993_2002 <- hold$rsds_min
 env_20m_ncc$temp_5m_cv_bccm1993_2002 <- hold$temp_5m_cv_bccm
 env_20m_ncc$temp_5m_cv_nep361993_2002 <- hold$temp_5m_cv_nep36
 env_20m_ncc$temp_5m_mean_bccm1993_2002 <- hold$temp_5m_mean_bccm
 env_20m_ncc$temp_5m_mean_nep361993_2002 <- hold$temp_5m_mean_nep36
 env_20m_ncc$NH4_5m_mean_bccm1993_2002 <- hold$NH4_5m_mean_bccm
 env_20m_ncc$NH4_5m_mean_nep361993_2002 <- hold$NH4_5m_mean_nep36
 env_20m_ncc$temp_air_cv1993_2002 <- hold$temp_air_cv
 env_20m_ncc$temp_s_cv_bccm1993_2002 <- hold$temp_s_cv_bccm
 env_20m_ncc$temp_s_cv_nep361993_2002 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast2003_2012, y = env_20m_ncc_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_ncc$precip_mean2003_2012 <- hold$precip_mean
 env_20m_ncc$precip_min2003_2012 <- hold$precip_min
 env_20m_ncc$salt_5m_cv_bccm2003_2012 <- hold$salt_5m_cv_bccm
 env_20m_ncc$salt_5m_cv_nep362003_2012 <- hold$salt_5m_cv_nep36
 env_20m_ncc$salt_5m_mean_bccm2003_2012 <- hold$salt_5m_mean_bccm
 env_20m_ncc$salt_5m_mean_nep362003_2012 <- hold$salt_5m_mean_nep36
 env_20m_ncc$temp_air_min2003_2012 <- hold$temp_air_min
 env_20m_ncc$rsds_min2003_2012 <- hold$rsds_min
 env_20m_ncc$temp_5m_cv_bccm2003_2012 <- hold$temp_5m_cv_bccm
 env_20m_ncc$temp_5m_cv_nep362003_2012 <- hold$temp_5m_cv_nep36
 env_20m_ncc$temp_5m_mean_bccm2003_2012 <- hold$temp_5m_mean_bccm
 env_20m_ncc$temp_5m_mean_nep362003_2012 <- hold$temp_5m_mean_nep36
 env_20m_ncc$NH4_5m_mean_bccm2003_2012 <- hold$NH4_5m_mean_bccm
 env_20m_ncc$NH4_5m_mean_nep362003_2012 <- hold$NH4_5m_mean_nep36
 env_20m_ncc$temp_air_cv2003_2012 <- hold$temp_air_cv
 env_20m_ncc$temp_s_cv_bccm2003_2012 <- hold$temp_s_cv_bccm
 env_20m_ncc$temp_s_cv_nep362003_2012 <- hold$temp_s_cv_nep36
 

 #queen charlotte
 env_20m_qcs <- as.data.frame(bathy_qcs, xy=TRUE)
 names(env_20m_qcs) <- c("X_m", "Y_m", "depth")
 env_20m_qcs <- env_20m_qcs %>% filter(depth <= max_depth, depth >= min_depth)
 env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = rei_qcs, y = env_20m_qcs_sf)
 env_20m_qcs$rei <- hold$rei 
 env_20m_qcs <- env_20m_qcs %>% filter(rei <= max_rei)
 env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = tidal_index_all, y = env_20m_qcs_sf)
 env_20m_qcs$tidal <- hold$tidal 
 env_20m_qcs <- env_20m_qcs %>% filter(tidal <= max_tidal)
 env_20m_qcs_sf <- st_as_sf(env_20m_qcs, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = hindcast2013_2023, y = env_20m_qcs_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_qcs$precip_mean2013_2023 <- hold$precip_mean
 env_20m_qcs$precip_min2013_2023 <- hold$precip_min
 env_20m_qcs$salt_5m_cv_bccm2013_2023 <- hold$salt_5m_cv_bccm
 env_20m_qcs$salt_5m_cv_nep362013_2023 <- hold$salt_5m_cv_nep36
 env_20m_qcs$salt_5m_mean_bccm2013_2023 <- hold$salt_5m_mean_bccm
 env_20m_qcs$salt_5m_mean_nep362013_2023 <- hold$salt_5m_mean_nep36
 env_20m_qcs$temp_air_min2013_2023 <- hold$temp_air_min
 env_20m_qcs$rsds_min2013_2023 <- hold$rsds_min
 env_20m_qcs$temp_5m_cv_bccm2013_2023 <- hold$temp_5m_cv_bccm
 env_20m_qcs$temp_5m_cv_nep362013_2023 <- hold$temp_5m_cv_nep36
 env_20m_qcs$temp_5m_mean_bccm2013_2023 <- hold$temp_5m_mean_bccm
 env_20m_qcs$temp_5m_mean_nep362013_2023 <- hold$temp_5m_mean_nep36
 env_20m_qcs$NH4_5m_mean_bccm2013_2023 <- hold$NH4_5m_mean_bccm
 env_20m_qcs$NH4_5m_mean_nep362013_2023 <- hold$NH4_5m_mean_nep36
 env_20m_qcs$temp_air_cv2013_2023 <- hold$temp_air_cv
 env_20m_qcs$temp_s_cv_bccm2013_2023 <- hold$temp_s_cv_bccm
 env_20m_qcs$temp_s_cv_nep362013_2023 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast1993_2002, y = env_20m_qcs_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_qcs$precip_mean1993_2002 <- hold$precip_mean
 env_20m_qcs$precip_min1993_2002 <- hold$precip_min
 env_20m_qcs$salt_5m_cv_bccm1993_2002 <- hold$salt_5m_cv_bccm
 env_20m_qcs$salt_5m_cv_nep361993_2002 <- hold$salt_5m_cv_nep36
 env_20m_qcs$salt_5m_mean_bccm1993_2002 <- hold$salt_5m_mean_bccm
 env_20m_qcs$salt_5m_mean_nep361993_2002 <- hold$salt_5m_mean_nep36
 env_20m_qcs$temp_air_min1993_2002 <- hold$temp_air_min
 env_20m_qcs$rsds_min1993_2002 <- hold$rsds_min
 env_20m_qcs$temp_5m_cv_bccm1993_2002 <- hold$temp_5m_cv_bccm
 env_20m_qcs$temp_5m_cv_nep361993_2002 <- hold$temp_5m_cv_nep36
 env_20m_qcs$temp_5m_mean_bccm1993_2002 <- hold$temp_5m_mean_bccm
 env_20m_qcs$temp_5m_mean_nep361993_2002 <- hold$temp_5m_mean_nep36
 env_20m_qcs$NH4_5m_mean_bccm1993_2002 <- hold$NH4_5m_mean_bccm
 env_20m_qcs$NH4_5m_mean_nep361993_2002 <- hold$NH4_5m_mean_nep36
 env_20m_qcs$temp_air_cv1993_2002 <- hold$temp_air_cv
 env_20m_qcs$temp_s_cv_bccm1993_2002 <- hold$temp_s_cv_bccm
 env_20m_qcs$temp_s_cv_nep361993_2002 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast2003_2012, y = env_20m_qcs_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_qcs$precip_mean2003_2012 <- hold$precip_mean
 env_20m_qcs$precip_min2003_2012 <- hold$precip_min
 env_20m_qcs$salt_5m_cv_bccm2003_2012 <- hold$salt_5m_cv_bccm
 env_20m_qcs$salt_5m_cv_nep362003_2012 <- hold$salt_5m_cv_nep36
 env_20m_qcs$salt_5m_mean_bccm2003_2012 <- hold$salt_5m_mean_bccm
 env_20m_qcs$salt_5m_mean_nep362003_2012 <- hold$salt_5m_mean_nep36
 env_20m_qcs$temp_air_min2003_2012 <- hold$temp_air_min
 env_20m_qcs$rsds_min2003_2012 <- hold$rsds_min
 env_20m_qcs$temp_5m_cv_bccm2003_2012 <- hold$temp_5m_cv_bccm
 env_20m_qcs$temp_5m_cv_nep362003_2012 <- hold$temp_5m_cv_nep36
 env_20m_qcs$temp_5m_mean_bccm2003_2012 <- hold$temp_5m_mean_bccm
 env_20m_qcs$temp_5m_mean_nep362003_2012 <- hold$temp_5m_mean_nep36
 env_20m_qcs$NH4_5m_mean_bccm2003_2012 <- hold$NH4_5m_mean_bccm
 env_20m_qcs$NH4_5m_mean_nep362003_2012 <- hold$NH4_5m_mean_nep36
 env_20m_qcs$temp_air_cv2003_2012 <- hold$temp_air_cv
 env_20m_qcs$temp_s_cv_bccm2003_2012 <- hold$temp_s_cv_bccm
 env_20m_qcs$temp_s_cv_nep362003_2012 <- hold$temp_s_cv_nep36
 

 #salish sea
 env_20m_ss <- as.data.frame(bathy_ss, xy=TRUE)
 names(env_20m_ss) <- c("X_m", "Y_m", "depth")
 env_20m_ss <- env_20m_ss %>% filter(depth <= max_depth, depth >= min_depth)
 env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = rei_ss, y = env_20m_ss_sf)
 env_20m_ss$rei <- hold$rei 
 env_20m_ss <- env_20m_ss %>% filter(rei <= max_rei)
 env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = tidal_index_all, y = env_20m_ss_sf)
 env_20m_ss$tidal <- hold$tidal 
 env_20m_ss <- env_20m_ss %>% filter(tidal <= max_tidal)
 env_20m_ss_sf <- st_as_sf(env_20m_ss, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = hindcast2013_2023, y = env_20m_ss_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_ss$precip_mean2013_2023 <- hold$precip_mean
 env_20m_ss$precip_min2013_2023 <- hold$precip_min
 env_20m_ss$salt_5m_cv_bccm2013_2023 <- hold$salt_5m_cv_bccm
 env_20m_ss$salt_5m_cv_nep362013_2023 <- hold$salt_5m_cv_nep36
 env_20m_ss$salt_5m_mean_bccm2013_2023 <- hold$salt_5m_mean_bccm
 env_20m_ss$salt_5m_mean_nep362013_2023 <- hold$salt_5m_mean_nep36
 env_20m_ss$temp_air_min2013_2023 <- hold$temp_air_min
 env_20m_ss$rsds_min2013_2023 <- hold$rsds_min
 env_20m_ss$temp_5m_cv_bccm2013_2023 <- hold$temp_5m_cv_bccm
 env_20m_ss$temp_5m_cv_nep362013_2023 <- hold$temp_5m_cv_nep36
 env_20m_ss$temp_5m_mean_bccm2013_2023 <- hold$temp_5m_mean_bccm
 env_20m_ss$temp_5m_mean_nep362013_2023 <- hold$temp_5m_mean_nep36
 env_20m_ss$NH4_5m_mean_bccm2013_2023 <- hold$NH4_5m_mean_bccm
 env_20m_ss$NH4_5m_mean_nep362013_2023 <- hold$NH4_5m_mean_nep36
 env_20m_ss$temp_air_cv2013_2023 <- hold$temp_air_cv
 env_20m_ss$temp_s_cv_bccm2013_2023 <- hold$temp_s_cv_bccm
 env_20m_ss$temp_s_cv_nep362013_2023 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast1993_2002, y = env_20m_ss_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_ss$precip_mean1993_2002 <- hold$precip_mean
 env_20m_ss$precip_min1993_2002 <- hold$precip_min
 env_20m_ss$salt_5m_cv_bccm1993_2002 <- hold$salt_5m_cv_bccm
 env_20m_ss$salt_5m_cv_nep361993_2002 <- hold$salt_5m_cv_nep36
 env_20m_ss$salt_5m_mean_bccm1993_2002 <- hold$salt_5m_mean_bccm
 env_20m_ss$salt_5m_mean_nep361993_2002 <- hold$salt_5m_mean_nep36
 env_20m_ss$temp_air_min1993_2002 <- hold$temp_air_min
 env_20m_ss$rsds_min1993_2002 <- hold$rsds_min
 env_20m_ss$temp_5m_cv_bccm1993_2002 <- hold$temp_5m_cv_bccm
 env_20m_ss$temp_5m_cv_nep361993_2002 <- hold$temp_5m_cv_nep36
 env_20m_ss$temp_5m_mean_bccm1993_2002 <- hold$temp_5m_mean_bccm
 env_20m_ss$temp_5m_mean_nep361993_2002 <- hold$temp_5m_mean_nep36
 env_20m_ss$NH4_5m_mean_bccm1993_2002 <- hold$NH4_5m_mean_bccm
 env_20m_ss$NH4_5m_mean_nep361993_2002 <- hold$NH4_5m_mean_nep36
 env_20m_ss$temp_air_cv1993_2002 <- hold$temp_air_cv
 env_20m_ss$temp_s_cv_bccm1993_2002 <- hold$temp_s_cv_bccm
 env_20m_ss$temp_s_cv_nep361993_2002 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast2003_2012, y = env_20m_ss_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_ss$precip_mean2003_2012 <- hold$precip_mean
 env_20m_ss$precip_min2003_2012 <- hold$precip_min
 env_20m_ss$salt_5m_cv_bccm2003_2012 <- hold$salt_5m_cv_bccm
 env_20m_ss$salt_5m_cv_nep362003_2012 <- hold$salt_5m_cv_nep36
 env_20m_ss$salt_5m_mean_bccm2003_2012 <- hold$salt_5m_mean_bccm
 env_20m_ss$salt_5m_mean_nep362003_2012 <- hold$salt_5m_mean_nep36
 env_20m_ss$temp_air_min2003_2012 <- hold$temp_air_min
 env_20m_ss$rsds_min2003_2012 <- hold$rsds_min
 env_20m_ss$temp_5m_cv_bccm2003_2012 <- hold$temp_5m_cv_bccm
 env_20m_ss$temp_5m_cv_nep362003_2012 <- hold$temp_5m_cv_nep36
 env_20m_ss$temp_5m_mean_bccm2003_2012 <- hold$temp_5m_mean_bccm
 env_20m_ss$temp_5m_mean_nep362003_2012 <- hold$temp_5m_mean_nep36
 env_20m_ss$NH4_5m_mean_bccm2003_2012 <- hold$NH4_5m_mean_bccm
 env_20m_ss$NH4_5m_mean_nep362003_2012 <- hold$NH4_5m_mean_nep36
 env_20m_ss$temp_air_cv2003_2012 <- hold$temp_air_cv
 env_20m_ss$temp_s_cv_bccm2003_2012 <- hold$temp_s_cv_bccm
 env_20m_ss$temp_s_cv_nep362003_2012 <- hold$temp_s_cv_nep36
 


 #west coast vancouver island
 env_20m_wcvi <- as.data.frame(bathy_wcvi, xy=TRUE)
 names(env_20m_wcvi) <- c("X_m", "Y_m", "depth")
 env_20m_wcvi <- env_20m_wcvi %>% filter(depth <= max_depth, depth >= min_depth)
 env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = rei_wcvi, y = env_20m_wcvi_sf)
 env_20m_wcvi$rei <- hold$rei 
 env_20m_wcvi <- env_20m_wcvi %>% filter(rei <= max_rei)
 env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = tidal_index_all, y = env_20m_wcvi_sf)
 env_20m_wcvi$tidal <- hold$tidal 
 env_20m_wcvi <- env_20m_wcvi %>% filter(tidal <= max_tidal)
 env_20m_wcvi_sf <- st_as_sf(env_20m_wcvi, coords = c("X_m","Y_m"), crs = "EPSG:3005")
 hold <- terra::extract(x = hindcast2013_2023, y = env_20m_wcvi_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_wcvi$precip_mean2013_2023 <- hold$precip_mean
 env_20m_wcvi$precip_min2013_2023 <- hold$precip_min
 env_20m_wcvi$salt_5m_cv_bccm2013_2023 <- hold$salt_5m_cv_bccm
 env_20m_wcvi$salt_5m_cv_nep362013_2023 <- hold$salt_5m_cv_nep36
 env_20m_wcvi$salt_5m_mean_bccm2013_2023 <- hold$salt_5m_mean_bccm
 env_20m_wcvi$salt_5m_mean_nep362013_2023 <- hold$salt_5m_mean_nep36
 env_20m_wcvi$temp_air_min2013_2023 <- hold$temp_air_min
 env_20m_wcvi$rsds_min2013_2023 <- hold$rsds_min
 env_20m_wcvi$temp_5m_cv_bccm2013_2023 <- hold$temp_5m_cv_bccm
 env_20m_wcvi$temp_5m_cv_nep362013_2023 <- hold$temp_5m_cv_nep36
 env_20m_wcvi$temp_5m_mean_bccm2013_2023 <- hold$temp_5m_mean_bccm
 env_20m_wcvi$temp_5m_mean_nep362013_2023 <- hold$temp_5m_mean_nep36
 env_20m_wcvi$NH4_5m_mean_bccm2013_2023 <- hold$NH4_5m_mean_bccm
 env_20m_wcvi$NH4_5m_mean_nep362013_2023 <- hold$NH4_5m_mean_nep36
 env_20m_wcvi$temp_air_cv2013_2023 <- hold$temp_air_cv
 env_20m_wcvi$temp_s_cv_bccm2013_2023 <- hold$temp_s_cv_bccm
 env_20m_wcvi$temp_s_cv_nep362013_2023 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast1993_2002, y = env_20m_wcvi_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_wcvi$precip_mean1993_2002 <- hold$precip_mean
 env_20m_wcvi$precip_min1993_2002 <- hold$precip_min
 env_20m_wcvi$salt_5m_cv_bccm1993_2002 <- hold$salt_5m_cv_bccm
 env_20m_wcvi$salt_5m_cv_nep361993_2002 <- hold$salt_5m_cv_nep36
 env_20m_wcvi$salt_5m_mean_bccm1993_2002 <- hold$salt_5m_mean_bccm
 env_20m_wcvi$salt_5m_mean_nep361993_2002 <- hold$salt_5m_mean_nep36
 env_20m_wcvi$temp_air_min1993_2002 <- hold$temp_air_min
 env_20m_wcvi$rsds_min1993_2002 <- hold$rsds_min
 env_20m_wcvi$temp_5m_cv_bccm1993_2002 <- hold$temp_5m_cv_bccm
 env_20m_wcvi$temp_5m_cv_nep361993_2002 <- hold$temp_5m_cv_nep36
 env_20m_wcvi$temp_5m_mean_bccm1993_2002 <- hold$temp_5m_mean_bccm
 env_20m_wcvi$temp_5m_mean_nep361993_2002 <- hold$temp_5m_mean_nep36
 env_20m_wcvi$NH4_5m_mean_bccm1993_2002 <- hold$NH4_5m_mean_bccm
 env_20m_wcvi$NH4_5m_mean_nep361993_2002 <- hold$NH4_5m_mean_nep36
 env_20m_wcvi$temp_air_cv1993_2002 <- hold$temp_air_cv
 env_20m_wcvi$temp_s_cv_bccm1993_2002 <- hold$temp_s_cv_bccm
 env_20m_wcvi$temp_s_cv_nep361993_2002 <- hold$temp_s_cv_nep36
 
 hold <- terra::extract(x = hindcast2003_2012, y = env_20m_wcvi_sf, fun = "mean", touches = TRUE, bind = TRUE) %>% terra::as.data.frame() 
 env_20m_wcvi$precip_mean2003_2012 <- hold$precip_mean
 env_20m_wcvi$precip_min2003_2012 <- hold$precip_min
 env_20m_wcvi$salt_5m_cv_bccm2003_2012 <- hold$salt_5m_cv_bccm
 env_20m_wcvi$salt_5m_cv_nep362003_2012 <- hold$salt_5m_cv_nep36
 env_20m_wcvi$salt_5m_mean_bccm2003_2012 <- hold$salt_5m_mean_bccm
 env_20m_wcvi$salt_5m_mean_nep362003_2012 <- hold$salt_5m_mean_nep36
 env_20m_wcvi$temp_air_min2003_2012 <- hold$temp_air_min
 env_20m_wcvi$rsds_min2003_2012 <- hold$rsds_min
 env_20m_wcvi$temp_5m_cv_bccm2003_2012 <- hold$temp_5m_cv_bccm
 env_20m_wcvi$temp_5m_cv_nep362003_2012 <- hold$temp_5m_cv_nep36
 env_20m_wcvi$temp_5m_mean_bccm2003_2012 <- hold$temp_5m_mean_bccm
 env_20m_wcvi$temp_5m_mean_nep362003_2012 <- hold$temp_5m_mean_nep36
 env_20m_wcvi$NH4_5m_mean_bccm2003_2012 <- hold$NH4_5m_mean_bccm
 env_20m_wcvi$NH4_5m_mean_nep362003_2012 <- hold$NH4_5m_mean_nep36
 env_20m_wcvi$temp_air_cv2003_2012 <- hold$temp_air_cv
 env_20m_wcvi$temp_s_cv_bccm2003_2012 <- hold$temp_s_cv_bccm
 env_20m_wcvi$temp_s_cv_nep362003_2012 <- hold$temp_s_cv_nep36
 
 env_20m_all <- bind_rows(env_20m_hg, env_20m_ncc, env_20m_qcs, env_20m_wcvi,env_20m_ss)
 summary(env_20m_all)
 env_20m_all <- env_20m_all %>%
   filter(salt_5m_mean_bccm2013_2023 > quantile(seagrass_data$saltmean_bccm, probs = 0.001))
 summary(env_20m_all)   
 
 env_20m_all <- na.omit(env_20m_all)
 
 env_long <- env_20m_all %>%
   pivot_longer(
     cols = matches("(1993_2002|2003_2012|2013_2023)$"),
     names_to = c("variable", "decade"),
     names_pattern = "^(.*?)(1993_2002|2003_2012|2013_2023)$",
     values_to = "value"
   ) %>%
   mutate(
     decade = recode(
       decade,
       "1993_2002" = "1993–2002",
       "2003_2012" = "2003–2012",
       "2013_2023" = "2013–2023"
     )
   )
 
 save(env_long, file = "code/output_data/environmental_data_all_decades.RData")
 save(env_20m_all, file = "code/output_data/environmental_data_all_decades_short.RData")
 
 
 obs_decade <- seagrass_data %>%
   mutate(
     decade = case_when(
       Year >= 1993 & Year <= 2002 ~ "1993–2002",
       Year >= 2003 & Year <= 2012 ~ "2003–2012",
       Year >= 2013 & Year <= 2023 ~ "2013–2023",
       TRUE ~ NA_character_
     )
   ) %>%
   filter(!is.na(decade))
 
 save(obs_decade, file = "code/output_data/survey_data_all_decades.RData")
 
 load("code/output_data/environmental_data_all_decades.RData")
 load("code/output_data/survey_data_all_decades.RData")
 load("code/output_data/environmental_data_all_decades_short.RData")
 
 
 
 
 obs_decade_subset <- obs_decade %>%
   select (X_m, Y_m, NH4_bccm, NH4_nep, saltmean_bccm, saltmean_nep, saltcv_bccm, saltcv_nep,
           tempmean_bccm, tempmean_nep, tempcv_bccm, tempcv_nep, surftempcv_bccm, surftempcv_nep,
           prmean, prmin, rsdsmin, airtempmin, airtempcv, decade)

 obs_decade_subset <- obs_decade_subset %>%
   rename(NH4_5m_mean_bccm = NH4_bccm, 
          NH4_5m_mean_nep36 = NH4_nep, 
          salt_5m_mean_bccm = saltmean_bccm, 
          salt_5m_mean_nep36 = saltmean_nep, 
          salt_5m_cv_bccm = saltcv_bccm, 
          salt_5m_cv_nep36 = saltcv_nep,
          temp_5m_mean_bccm= tempmean_bccm, 
          temp_5m_mean_nep36 = tempmean_nep, 
          temp_5m_cv_bccm = tempcv_bccm, 
          temp_5m_cv_nep36 = tempcv_nep, 
          temp_s_cv_bccm= surftempcv_bccm, 
          temp_s_cv_nep36 = surftempcv_nep,
          precip_mean = prmean, 
          precip_min = prmin, 
          rsds_min = rsdsmin, 
          temp_air_min = airtempmin, 
          temp_air_cv = airtempcv)
 
 obs_long <- obs_decade_subset %>%
   pivot_longer(
     cols = -c(X_m, Y_m, decade),
     names_to = "variable",
     values_to = "value"
   )
 
 env_long <- env_long %>% select(-depth, -rei, -tidal)
 
 
 env_stats <- env_long %>%
   filter(!is.na(value)) %>%
   group_by(variable, decade) %>%
   summarise(
     env_mean   = mean(value),
     env_median = median(value),
     env_sd     = sd(value),
     env_p10    = quantile(value, 0.10),
     env_p90    = quantile(value, 0.90),
     .groups = "drop"
   )
 
 
 obs_stats <- obs_long %>%
   filter(!is.na(value)) %>%
   group_by(variable, decade) %>%
   summarise(
     obs_mean   = mean(value),
     obs_median = median(value),
     obs_sd     = sd(value),
     obs_p10    = quantile(value, 0.10),
     obs_p90    = quantile(value, 0.90),
     .groups = "drop"
   )
 
 
 overlap_stats <- obs_long %>%
   filter(!is.na(value)) %>%
   left_join(
     env_stats %>%
       select(variable, decade, env_p10, env_p90),
     by = c("variable", "decade")
   ) %>%
   mutate(
     within_env = value >= env_p10 & value <= env_p90
   ) %>%
   group_by(variable, decade) %>%
   summarise(
     n_obs = n(),
     n_within = sum(within_env),
     prop_overlap = n_within / n_obs,
     .groups = "drop"
   )
 
 summary_table <- env_stats %>%
   left_join(obs_stats,
             by = c("variable", "decade")) %>%
   left_join(overlap_stats,
             by = c("variable", "decade"))
 
 summary_table <- summary_table %>%
   mutate(
     env_range80 = paste0(round(env_p10, 2), "–", round(env_p90, 2)),
     obs_range80 = paste0(round(obs_p10, 2), "–", round(obs_p90, 2))
   )
 
 
 summary_table %>%
   arrange(variable, decade)
 
 change_stats <- summary_table %>%
   group_by(variable) %>%
   arrange(decade) %>%
   mutate(
     mean_change = env_mean - lag(env_mean),
     sd_pooled = sqrt((env_sd^2 + lag(env_sd)^2)/2),
     cohens_d = mean_change / sd_pooled
   )
 
 #Positive d = variable increased relative to the previous decade.
 #Negative d = variable decreased relative to the previous decade.
 #Large absolute values = stronger shifts.
 # 0.2 small, 0.5 moderate, 0.8 large
 

 #plots
 precip_mean <- make_atm_plot(
   "precip_mean",
   expression("Precipitation (kg m"^-2 ~ "month"^-1 * ")")
 )
 
 precip_min <- make_atm_plot(
   "precip_min",
   expression("Min precipitation (kg m"^-2 ~ "month"^-1 * ")")
 )
 
 rsds_min <- make_atm_plot(
   "rsds_min",
   expression("Min Surface Downwelling\nShortwave Flux (MJ m"^-2 * ")")
 )
 
 airtemp_min <- make_atm_plot(
   "temp_air_min",
   expression("Min air temperature (" * degree * C * ")"),
   convert_kelvin = TRUE
 )
 
 airtemp_cv <- make_atm_plot(
   "temp_air_cv",
   "Air temperature CV",
   show_x = TRUE
 )
 
 salt_mean <- make_bccm_nep_plot(
   "salt_5m_mean_bccm",
   "salt_5m_mean_nep36",
   expression("Salinity (psu)"),
   ylim = c(20, NA)
 )
 
 salt_cv <- make_bccm_nep_plot(
   "salt_5m_cv_bccm",
   "salt_5m_cv_nep36",
   expression("Salinity CV"),
   ylim = c(0, 20)
 )
 
 temp_mean <- make_bccm_nep_plot(
   "temp_5m_mean_bccm",
   "temp_5m_mean_nep36",
   expression("Ocean subsurface temperature (" * degree * C * ")")
 )
 
 temp_cv <- make_bccm_nep_plot(
   "temp_5m_cv_bccm",
   "temp_5m_cv_nep36",
   expression("Ocean subsurface temperature CV")
 )
 
 NH4_mean <- make_bccm_nep_plot(
   "NH4_5m_mean_bccm",
   "NH4_5m_mean_nep36",
   expression("Ammonium (mmol m"^-3 * ")"),
   ylim = c(0, 2.5),
   show_x = TRUE
 )
 
 surfacetemp_cv <- make_bccm_nep_plot(
   "temp_s_cv_bccm",
   "temp_s_cv_nep36",
   expression("Ocean surface temperature CV")
 )
 
 
 atmos <- precip_mean / precip_min / rsds_min / airtemp_min / airtemp_cv + plot_layout(guides = "collect") & theme(legend.position = "right")
 
 ggsave(
   filename = "figures/atmos_var.png",
   plot = atmos,
   width = 12,
   height = 20,
   dpi = 300,
   bg = "white"
 )
 
 ocean<- salt_mean / salt_cv / temp_mean / temp_cv / surfacetemp_cv / NH4_mean + plot_layout(guides = "collect") & theme(legend.position = "right")
 
 ggsave(
   filename = "figures/ocean_var.png",
   plot = ocean,
   width = 14,
   height = 20,
   dpi = 300,
   bg = "white"
 )
 


 
 
 # MESS figures
 
 hindcast_predictor_data_nep <- env_20m_all %>% select(temp_s_cv_nep362013_2023, NH4_5m_mean_nep362013_2023, temp_5m_mean_nep362013_2023, 
                                                       temp_5m_cv_nep362013_2023, salt_5m_mean_nep362013_2023, salt_5m_cv_nep362013_2023, 
                                                       precip_mean2013_2023, precip_min2013_2023, temp_air_min2013_2023, temp_air_cv2013_2023,
                                                       rsds_min2013_2023, tidal, rei)
 hindcast_predictor_data_nep <- hindcast_predictor_data_nep %>% rename (surftempcv_nep = temp_s_cv_nep362013_2023, 
                                                                        NH4_nep = NH4_5m_mean_nep362013_2023, 
                                                                        tempmean_nep = temp_5m_mean_nep362013_2023, 
                                                                        tempcv_nep = temp_5m_cv_nep362013_2023, 
                                                                        saltmean_nep = salt_5m_mean_nep362013_2023, 
                                                                        saltcv_nep = salt_5m_cv_nep362013_2023, 
                                                                        prmean = precip_mean2013_2023, 
                                                                        prmin = precip_min2013_2023, 
                                                                        airtempmin = temp_air_min2013_2023, 
                                                                        airtempcv = temp_air_cv2013_2023,
                                                                        rsdsmin = rsds_min2013_2023)
 
 
 transect_predictor_data <- obs_decade %>% select(surftempcv_nep, NH4_nep, tempmean_nep, tempcv_nep, saltmean_nep, saltcv_nep, 
                                                   prmean, prmin, airtempmin, airtempcv, rsdsmin, tidal, rei)
 names(hindcast_predictor_data_nep)
 names(transect_predictor_data)
 mess_all<- predicts::mess(x = hindcast_predictor_data_nep, v= transect_predictor_data, full = FALSE)
 
 transect_predictor_data_1993_2009 <- obs_decade %>% 
   filter(Year< 2010) %>%
   select(surftempcv_nep, NH4_nep, tempmean_nep, tempcv_nep, saltmean_nep, saltcv_nep, 
          prmean, prmin, airtempmin, airtempcv, rsdsmin, tidal, rei)
 
 # MESS from only pre-2010 survey data
  mess_pre2010 <- predicts::mess(x = hindcast_predictor_data_nep, v = transect_predictor_data_1993_2009, full = FALSE)
  
  mess_df <- env_20m_all %>%
    select(X_m, Y_m) %>%
    mutate(
      mess_all = mess_all,
      mess_pre2010 = mess_pre2010,
      extrap_all = mess_all < 0,
      extrap_pre2010 = mess_pre2010 < 0
    )

  mess_df <- mess_df %>%
    mutate(
      extrap_type = case_when(
        extrap_all & extrap_pre2010 ~ "Both",
        !extrap_all & extrap_pre2010 ~ "Pre-2010 only",
        extrap_all & !extrap_pre2010 ~ "All-data only",
        TRUE ~ "Neither"
      )
    )  

  mess_df <- mess_df %>%
    mutate(
      category = case_when(
        mess_all < 0 ~ "Outside survey range",
        mess_all >= 0 & mess_pre2010 < 0 ~ "Outside pre-2010 range only",
        TRUE ~ "Within range"
      )
    )
  
  mess_summary <- mess_df %>%
    count(category) %>%
    mutate(
      percent = 100 * n / sum(n)
    )
  
  mess_summary
  
  save(mess_df, file = "code/output_data/mess.RData")
  
  
  mess_plot <- mess_df |>
    dplyr::filter(category %in% c(
      "Outside survey range",
      "Outside pre-2010 range only"
    ))
  
  mess_fig <-  ggplot() +
    geom_sf(
      data = coastline,
      fill = "grey95",
      colour = "grey55",
      linewidth = 0.3
    ) +
    geom_point(
      data = mess_plot,
      aes(X_m, Y_m, colour = category),
      size = 1
    ) +
    scale_colour_manual(
      values = c(
        "Outside survey range" = "#1B9E77",
        "Outside pre-2010 range only" = "#7570B3"
      ),
      name = NULL
    ) +
    coord_sf(
      xlim = c(530000, 1228000),
      ylim = c(366000, 1101000),
      expand = FALSE
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      legend.position = "bottom"
    )
  

  
  mess_fig
  
  ggsave(
    filename = "figures/mess.png",
    plot = mess_fig,
    width = 180,
    height = 210,
    units = "mm",
    dpi = 600,
    bg = "white" )
  
  
  
  load("code/output_data/mess.RData")
  
  
  mess_extrap <- cbind(
    env_20m_all,
    mess_df %>% dplyr::select(-X_m, -Y_m)
  )
  
  mess_extrap <- mess_extrap %>%
    dplyr::filter(category %in% c(
      "Outside pre-2010 range only"
    ))
  
  mess_extrap <- mess_extrap %>% select(temp_s_cv_nep362013_2023, NH4_5m_mean_nep362013_2023, temp_5m_mean_nep362013_2023, 
                                                        temp_5m_cv_nep362013_2023, salt_5m_mean_nep362013_2023, salt_5m_cv_nep362013_2023, 
                                                        precip_mean2013_2023, precip_min2013_2023, temp_air_min2013_2023, temp_air_cv2013_2023,
                                                        rsds_min2013_2023, tidal, rei)
  
  
  mess_extrap <- mess_extrap %>% 
    rename(NH4_5m_mean_nep36 = NH4_5m_mean_nep362013_2023, 
           salt_5m_mean_nep36 = salt_5m_mean_nep362013_2023, 
           salt_5m_cv_nep36 = salt_5m_cv_nep362013_2023, 
           temp_5m_mean_nep36 = temp_5m_mean_nep362013_2023, 
           temp_5m_cv_nep36 = temp_5m_cv_nep362013_2023, 
           temp_s_cv_nep36 = temp_s_cv_nep362013_2023,
           precip_mean = precip_mean2013_2023, 
           precip_min = precip_min2013_2023, 
           rsds_min = rsds_min2013_2023, 
           temp_air_min = temp_air_min2013_2023, 
           temp_air_cv = temp_air_cv2013_2023)
    
  obs_mess <- obs_decade %>%
    select (Year, X_m, Y_m, NH4_nep, saltmean_nep, saltcv_nep,
            tempmean_nep, tempcv_nep, surftempcv_nep,
            prmean, prmin, rsdsmin, airtempmin, airtempcv, rei, tidal, decade) %>%
    filter(Year < 2010)
  
  obs_mess <- obs_mess %>%
    rename(NH4_5m_mean_nep36 = NH4_nep,  
           salt_5m_mean_nep36 = saltmean_nep,
           salt_5m_cv_nep36 = saltcv_nep, 
           temp_5m_mean_nep36 = tempmean_nep,
           temp_5m_cv_nep36 = tempcv_nep,  
           temp_s_cv_nep36 = surftempcv_nep,
           precip_mean = prmean, 
           precip_min = prmin, 
           rsds_min = rsdsmin, 
           temp_air_min = airtempmin, 
           temp_air_cv = airtempcv)

  obs_mess <- obs_mess %>% select(-Year, -X_m, -Y_m, -decade)
  
  vars <- names(obs_mess)
  
  range_comparison <- purrr::map_dfr(vars, function(v){
    
    obs_min <- min(obs_mess[[v]], na.rm = TRUE)
    obs_max <- max(obs_mess[[v]], na.rm = TRUE)
    
    mess_vals <- mess_extrap[[v]]
    
    tibble(
      variable = v,
      prop_below = mean(mess_vals < obs_min, na.rm = TRUE),
      prop_above = mean(mess_vals > obs_max, na.rm = TRUE),
      prop_outside = mean(
        mess_vals < obs_min | mess_vals > obs_max,
        na.rm = TRUE
      )
    )
    
  })
  
  range_comparison %>%
    arrange(desc(prop_outside))  
  
  
  cohens_d <- purrr::map_dfr(vars, function(v){
    
    x <- obs_mess[[v]]
    y <- mess_extrap[[v]]
    
    m1 <- mean(x, na.rm = TRUE)
    m2 <- mean(y, na.rm = TRUE)
    
    sd_pool <- sqrt(
      (sd(x, na.rm = TRUE)^2 +
         sd(y, na.rm = TRUE)^2) / 2
    )
    
    tibble(
      variable = v,
      cohens_d = (m2 - m1) / sd_pool
    )
    
  })
  
  cohens_d %>%
    arrange(desc(abs(cohens_d)))
  
  
  
  
  
  ##### Independent validation figures
  # load("code/output_data/model_results/eelgrass_independent_dataframe.RData")
  # load("code/output_data/model_results/surfgrass_independent_dataframe.RData")
  # load("code/output_data/model_results/eelgrass_independent_eval.RData")
  # 
  # 
  # # Model columns in the order you want shown
  # model_labels <- c(
  #   "bccm_nospatial" = "GLMM_bccm",
  #   "bccm_spatial"   = "GLMM_spatial_bccm",
  #   "GBM_bccm"       = "GBM_bccm",
  #   "XGBOOST_bccm"   = "XGB_bccm",
  #   "nep_nospatial"  = "GLMM_nep",
  #   "nep_spatial"    = "GLMM_spatial_nep",
  #   "GBM_nep"        = "GBM_nep",
  #   "XGBOOST_nep"    = "XGB_nep"
  # )
  # 
  # 
  # eelgrass_long <- eelgrass_df %>%
  #   #select(all_of(model_names)) %>%
  #   pivot_longer(
  #     cols = -c(obs, ID, substrate, rock_group),
  #     names_to = "model",
  #     values_to = "predicted_suitability"
  #   ) %>%
  #   mutate(species = "Eelgrass")
  # 
  # eelgrass_long <- eelgrass_long %>%
  #   mutate(model = recode(as.character(model), !!!model_labels),
  #          model = factor(model, levels = model_labels))
  # 
  # surfgrass_long <- surfgrass_df %>%
  #   select(-bathy) %>%
  #   pivot_longer(
  #     cols = -c(obs, ID, substrate, rock_group),
  #     names_to = "model",
  #     values_to = "predicted_suitability"
  #   ) %>%
  #   mutate(species = "Surfgrass")
  # 
  # surfgrass_long <- surfgrass_long %>%
  #   mutate(model = recode(as.character(model), !!!model_labels),
  #          model = factor(model, levels = model_labels))
  # 
  # eelgrass_subset <- eelgrass_long %>% filter (model == "XGB_nep" | model == "GLMM_spatial_nep")
  # 
  # model_levels <- c("Lowest-performing model", "Highest-performing model")
  # 
  # thresholds <- tibble(
  #   model = c("GLMM_spatial_nep", "XGB_nep"),
  #   threshold = c(0.037, 0.031)   
  # )
  # 
  # summary_df <- eelgrass_subset %>%
  #   left_join(thresholds, by = "model") %>%
  #   group_by(model) %>%
  #   summarise(
  #     FPPS = mean(predicted_suitability >= threshold, na.rm = TRUE),
  #     FNR  = mean(predicted_suitability < threshold, na.rm = TRUE),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #     model_label = recode(
  #       model,
  #       "GLMM_spatial_nep" = "Lowest-performing model",
  #       "XGB_nep" = "Highest-performing model"
  #     ),
  #     model_label = factor(model_label, levels = model_levels)
  #   ) %>%
  #   pivot_longer(
  #     cols = c(FPPS, FNR),
  #     names_to = "metric",
  #     values_to = "value"
  #   )
  # 
  # mps_df <- eelgrass_subset %>%
  #   mutate(
  #     model_label = recode(
  #       model,
  #       "GLMM_spatial_nep" = "Lowest-performing model",
  #       "XGB_nep" = "Highest-performing model"
  #     ),
  #     model_label = factor(model_label, levels = model_levels),
  #     value = predicted_suitability
  #   ) %>%
  #   select(model_label, value)
  # 
  # 
  # 
  # 
  # p1 <- ggplot(
  #   mps_df,
  #   aes(x = model_label, y = value, fill = model_label)
  # ) +
  #   geom_violin(
  #     alpha = 0.5,
  #     colour = NA,
  #     width = 0.9
  #   ) +
  #   geom_boxplot(
  #     width = 0.15,
  #     outlier.alpha = 0.1,
  #     colour = "black"
  #   ) +
  #   stat_summary(
  #     fun = median,
  #     geom = "point",
  #     size = 2,
  #     colour = "black"
  #   ) +
  #   scale_fill_manual(values = c(
  #     "Highest-performing model" = "#8C6BB1",
  #     "Lowest-performing model" = "#4D4D4D"
  #   )) +
  #   scale_x_discrete(limits = model_levels) +
  #   labs(
  #     y = "Relative probability of occurrence",
  #     x = NULL
  #   ) +
  #   scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  #   theme_classic(base_size = 12) +
  #   theme(
  #     legend.position = "none",
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     plot.margin = margin(3, 3, 3, 3),
  #     axis.ticks.length = unit(2, "pt"),
  #     plot.background = element_rect(fill = "white", colour = NA)
  #   )
  # 
  # 
  # p2 <- ggplot(
  #   summary_df,
  #   aes(
  #     x = model_label,
  #     y = value,
  #     shape = metric
  #   )
  # ) +
  #   geom_point(
  #     aes(colour = model_label),
  #     size = 3.5,
  #     position = position_dodge(width = 0.1)
  #   ) +
  #   
  #   scale_colour_manual(values = c(
  #     "Highest-performing model" = "#8C6BB1",
  #     "Lowest-performing model" = "#4D4D4D"
  #   )) +
  #   
  #   scale_shape_manual(values = c(
  #     "FPPS" = 17,
  #     "FNR" = 15
  #   )) +
  #   
  #   scale_y_continuous(
  #     limits = c(0, 1),
  #     expand = expansion(mult = c(0.02, 0.02))
  #   ) +
  #   
  #   labs(
  #     y = "Proportion",
  #     x = NULL,
  #     shape = "Metric"
  #   ) +
  #   
  #   guides(
  #     colour = "none",
  #     shape = guide_legend(
  #       order = 1,
  #       nrow = 1
  #     )
  #   ) +
  #   
  #   theme_classic(base_size = 12) +
  #   theme(
  #     axis.text.x = element_text(),
  #     
  #     legend.position = "bottom",
  #     legend.direction = "horizontal",
  #     
  #     legend.background = element_rect(
  #       fill = scales::alpha("white", 0.85),
  #       colour = NA
  #     ),
  #     
  #     plot.margin = margin(3, 3, 3, 3),
  #     axis.ticks.length = unit(2, "pt"),
  #     plot.background = element_rect(fill = "white", colour = NA)
  #   )
  # 
  # 
  # 
  # 
  # surfgrass_subset <- surfgrass_long %>% filter (model == "GLMM_spatial_bccm" | model == "GLMM_nep")
  # 
  # model_levels <- c("Lowest-performing model", "Highest-performing model")
  # 
  # thresholds <- tibble(
  #   model = c("GLMM_nep", "GLMM_spatial_bccm"),
  #   threshold = c(0.014, 0.014)   
  # )
  # 
  # summary_df_surf <- surfgrass_subset %>%
  #   left_join(thresholds, by = "model") %>%
  #   group_by(model) %>%
  #   summarise(
  #     FPPS = mean(predicted_suitability >= threshold, na.rm = TRUE),
  #     FNR  = mean(predicted_suitability < threshold, na.rm = TRUE),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #     model_label = recode(
  #       model,
  #       "GLMM_nep" = "Lowest-performing model",
  #       "GLMM_spatial_bccm" = "Highest-performing model"
  #     ),
  #     model_label = factor(model_label, levels = model_levels)
  #   ) %>%
  #   pivot_longer(
  #     cols = c(FPPS, FNR),
  #     names_to = "metric",
  #     values_to = "value"
  #   )
  # 
  # mps_df_surf <- surfgrass_subset %>%
  #   mutate(
  #     model_label = recode(
  #       model,
  #       "GLMM_nep" = "Lowest-performing model",
  #       "GLMM_spatial_bccm" = "Highest-performing model"
  #     ),
  #     model_label = factor(model_label, levels = model_levels),
  #     value = predicted_suitability
  #   ) %>%
  #   select(model_label, value)
  # 
  # 
  # p4 <- ggplot(
  #   mps_df_surf,
  #   aes(x = model_label, y = value, fill = model_label)
  # ) +
  #   geom_violin(
  #     alpha = 0.5,
  #     colour = NA,
  #     width = 0.9
  #   ) +
  #   geom_boxplot(
  #     width = 0.15,
  #     outlier.alpha = 0.1,
  #     colour = "black"
  #   ) +
  #   stat_summary(
  #     fun = median,
  #     geom = "point",
  #     size = 2,
  #     colour = "black"
  #   ) +
  #   scale_fill_manual(values = c(
  #     "Highest-performing model" = "#8C6BB1",
  #     "Lowest-performing model" = "#4D4D4D"
  #   )) +
  #   scale_x_discrete(limits = model_levels) +
  #   labs(
  #     y = "Relative probability of occurrence",
  #     x = NULL
  #   ) +
  #   scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  #   theme_classic(base_size = 12) +
  #   theme(
  #     legend.position = "none",
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     plot.margin = margin(3, 3, 3, 3),
  #     axis.ticks.length = unit(2, "pt"),
  #     plot.background = element_rect(fill = "white", colour = NA)
  #   )
  # 
  # 
  # p5 <- ggplot(
  #   summary_df_surf,
  #   aes(
  #     x = model_label,
  #     y = value,
  #     shape = metric
  #   )
  # ) +
  #   geom_point(
  #     aes(colour = model_label),
  #     size = 3.5,
  #     position = position_dodge(width = 0.1)
  #   ) +
  #   
  #   scale_colour_manual(values = c(
  #     "Highest-performing model" = "#8C6BB1",
  #     "Lowest-performing model" = "#4D4D4D"
  #   )) +
  #   
  #   scale_shape_manual(values = c(
  #     "FPPS" = 17,
  #     "FNR" = 15
  #   )) +
  #   
  #   scale_y_continuous(
  #     limits = c(0, 1),
  #     expand = expansion(mult = c(0.02, 0.02))
  #   ) +
  #   
  #   labs(
  #     y = "Proportion",
  #     x = NULL,
  #     shape = "Metric"
  #   ) +
  #   
  #   guides(
  #     colour = "none",
  #     shape = "none"
  #   ) +
  #   
  #   theme_classic(base_size = 12) +
  #   theme(
  #     axis.text.x = element_text(),
  #     
  #     legend.position = "none",
  #     
  #     plot.margin = margin(3, 3, 3, 3),
  #     axis.ticks.length = unit(2, "pt"),
  #     plot.background = element_rect(fill = "white", colour = NA)
  #   )
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # #make map of individuals
  # #eelgrass
  # eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))
  # eel_bin_xgb_nep <-rast("raster/eelgrass/eelgrass_predictions_xgb_nep_binary_notmasked.tif")
  # eel_bin_nepspatial<-rast("raster/eelgrass/eelgrass_predictions_nepspatial_binary_notmasked.tif")
  # values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1
  # 
  # eelgrass_indep <- resample(eelgrass_indep, eel_bin_xgb_nep, method = "near")
  # 
  # # make sure rasters align
  # stopifnot(compareGeom(eelgrass_indep, eel_bin_xgb_nep, eel_bin_nepspatial))
  # 
  # # create output raster
  # out <- ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 1 & eel_bin_nepspatial == 1, 4,
  #             ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 1 & eel_bin_nepspatial == 0, 3,
  #                  ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 0 & eel_bin_nepspatial == 1, 2,
  #                       ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 0 & eel_bin_nepspatial == 0, 1, NA))))
  # 
  # # write result
  # writeRaster(out, "raster/eelgrass/eelgrass_good&badmodel_overlayindependent.tif", overwrite = TRUE)
  # # looked in arc for best example areas to showcase and clipped it
  # 
  # 
  # eelgrass_indep_val_clipped <- rast(c("raster/independent_figue.tif"))
  # 
  # eelgrass_cat <- as.factor(eelgrass_indep_val_clipped)
  # levels(eelgrass_cat) <- data.frame(
  #   value = c(1, 2, 3, 4),
  #   category = c(
  #     "Observed only",
  #     "Observed + lowest-performing \nmodel",
  #     "Observed + highest-performing \nmodel",
  #     "Observed + both models"
  #   )
  # )
  # 
  # r_ext <- terra::ext(eelgrass_cat)
  # 
  # p3 <- ggplot() +
  #   geom_spatraster(data = eelgrass_cat) +
  #   geom_sf(data = coastline, fill = "grey85", color = "black", linewidth = 0.2) +
  #   
  #   scale_fill_manual(
  #     name = "Agreement category",
  #     values = c(
  #       "Observed only" = "#D73027",
  #       "Observed + lowest-performing \nmodel" = "#FC8D59",
  #       "Observed + highest-performing \nmodel" = "#91BFDB",
  #       "Observed + both models" = "#1A9850"
  #     ),
  #     na.value = "white",
  #     na.translate = FALSE,
  #     drop = FALSE,
  #     guide = guide_legend(ncol = 2)
  #   ) +
  #   
  #   coord_sf(
  #     xlim = c(r_ext$xmin, r_ext$xmax),
  #     ylim = c(r_ext$ymin, r_ext$ymax),
  #     expand = FALSE
  #   ) +
  #   
  #   theme_minimal() +
  #   
  #   theme(
  #     panel.background = element_rect(fill = "white", color = NA),
  #     plot.background  = element_rect(fill = "white", color = NA),
  #     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
  #     
  #     legend.position = "bottom",
  #     legend.direction = "horizontal",
  #     legend.box = "horizontal",
  #     
  #     legend.background = element_rect(
  #       fill = scales::alpha("white", 0.8),
  #       color = NA
  #     ),
  #     
  #     legend.key = element_rect(fill = NA, color = NA),
  #     
  #     plot.margin = margin(3, 3, 3, 3)
  #   )
  # p3
  # 
  # surfgrass_indep <- rast(c("code/output_data/independent_validation/surfgrass_validation_raster_2013_2024.tif"))
  # surf_bin_bccmspatial <-rast("raster/surfgrass/surfgrass_predictions_bccmspatial_binary_notmasked.tif")
  # surf_bin_nepnospatial<-rast("raster/surfgrass/surfgrass_predictions_nepnospatial_binary_notmasked.tif")
  # values(surfgrass_indep)[values(surfgrass_indep) >= 1] <- 1
  # 
  # surfgrass_indep <- resample(surfgrass_indep, surf_bin_bccmspatial, method = "near")
  # 
  # # make sure rasters align
  # stopifnot(compareGeom(surfgrass_indep, surf_bin_bccmspatial, surf_bin_nepnospatial))
  # 
  # # create output raster
  # out <- ifel(surfgrass_indep == 1 & surf_bin_bccmspatial == 1 & surf_bin_nepnospatial == 1, 4,
  #             ifel(surfgrass_indep == 1 & surf_bin_bccmspatial == 1 & surf_bin_nepnospatial == 0, 3,
  #                  ifel(surfgrass_indep == 1 & surf_bin_bccmspatial == 0 & surf_bin_nepnospatial == 1, 2,
  #                       ifel(surfgrass_indep == 1 & surf_bin_bccmspatial == 0 & surf_bin_nepnospatial == 0, 1, NA))))
  # 
  # # write result
  # writeRaster(out, "raster/surfgrass/surfgrass_good&badmodel_overlayindependent.tif", overwrite = TRUE)
  # # looked in arc for best example areas to showcase and clipped it
  # 
  # surfgrass_indep_val_clipped <- rast(c("raster/independent_figue_surfgrass.tif"))
  # 
  # surfgrass_cat <- as.factor(surfgrass_indep_val_clipped)
  # levels(surfgrass_cat) <- data.frame(
  #   value = c(1, 2, 3, 4),
  #   category = c(
  #     "Observed only",
  #     "Observed + lowest-performing model",
  #     "Observed + highest-performing model",
  #     "Observed + both models"
  #   )
  # )
  # 
  # r_ext_surf <- terra::ext(surfgrass_cat)
  # 
  # levels(surfgrass_cat) <- data.frame(
  #   value = 1:4,
  #   category = c(
  #     "Observed only",
  #     "Observed + lowest-performing model",
  #     "Observed + highest-performing model",
  #     "Observed + both models"
  #   )
  # )
  # 
  # #surfgrass_cat <- as.factor(surfgrass_cat)
  # 
  # legend_levels <- c(
  #   "Observed only",
  #   "Observed + lowest-performing model",
  #   "Observed + highest-performing model",
  #   "Observed + both models"
  # )
  # 
  # p6 <-ggplot() +
  #   geom_spatraster(data = surfgrass_cat) +
  #   geom_sf(data = coastline, fill = "grey85", color = "black", linewidth = 0.2) +
  #   scale_fill_manual(
  #     name = "Agreement category",
  #     values = c(
  #       "Observed only" = "#D73027",
  #       "Observed + lowest-performing model" = "#FC8D59",
  #       "Observed + highest-performing model" = "#91BFDB",
  #       "Observed + both models" = "#1A9850"
  #     ),
  #     limits = legend_levels,  
  #     na.value = "white",
  #     na.translate = FALSE,
  #     drop = FALSE
  #   ) +
  #   coord_sf(
  #     xlim = c(r_ext_surf$xmin, r_ext_surf$xmax),
  #     ylim = c(r_ext_surf$ymin, r_ext_surf$ymax),
  #     expand = FALSE
  #   ) +
  #   theme_minimal() +
  #   theme(
  #     panel.background = element_rect(fill = "white", color = NA),
  #     plot.background  = element_rect(fill = "white", color = NA),
  #     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
  #     legend.position = "none",
  #     plot.margin = margin(3, 3, 3, 3)
  #   )  
  # 
  # p6
  # 
  # 
  # 
  # 
  # 
  # 
  # # COMBINE ------------------------------------------------------------
  # # p1 <- p1 + theme(plot.margin = margin(5,5,5,5))
  # # p2 <- p2 + theme(plot.margin = margin(5,5,5,5))
  # # p3 <- p3 + theme(plot.margin = margin(5,5,5,5))
  # 
  # p1 <- p1 +
  #   annotate(
  #     "text",
  #     x = -Inf, y = Inf,
  #     label = "Eelgrass",
  #     hjust = -0.1, vjust = 1.2,
  #     fontface = "bold"
  #   )
  # 
  # p4 <- p4 +
  #   annotate(
  #     "text",
  #     x = -Inf, y = Inf,
  #     label = "Surfgrass",
  #     hjust = -0.1, vjust = 1.2,
  #     fontface = "bold"
  #   )
  # 
  # 
  # indep_plot_eel <- (p1 / p2) | p3
  # indep_plot_surf <- (p4 / p5) | p6
  # 
  # indep_plot = indep_plot_eel/indep_plot_surf
  #  
  # ggsave(
  #   filename = "figures/independent_validation.png",
  #   plot = indep_plot,
  #   width = 12,
  #   height = 11,
  #   dpi = 300,
  #   bg = "white"
  # )