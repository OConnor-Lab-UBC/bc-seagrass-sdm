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
    panel = "Targeted validation sites\nby observed genus"
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





# probability of occurence compared to test sites, independent data, and field validation data with all thresholds
# just do best and worst models for figure
# eelgrass best model is XGBoost_nep and worst is BGM_nep
# surfgrass best model is and worse is
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/final_metrics_eelgrass.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/final_metrics_surfgrass.RData")

seagrass_data %>%
  count(Month)

seagrass_data %>%
  summarise(
    pct_may_sep = 100 * mean(Month %in% 5:9)
  )
#need to extract predictions from the 2013-2023 test data
load("code/output_data/seagrass_model_inputs.RData")
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
    Algorithm = "Best"
  )

training_worst <- test_eel_pa_long %>%
  filter(model == worst_model) %>%
  mutate(
    Dataset = "SDM building",
    Algorithm = "Worst"
  )

independent <- indep_eel_pa_long %>%
  filter(model == best_model) %>%
  mutate(
    Dataset = "Independent",
    Algorithm = "Best"
  )

independent_worst <- indep_eel_pa_long %>%
  filter(model == worst_model) %>%
  mutate(
    Dataset = "Independent",
    Algorithm = "Worst"
  )

field <- field_eelgrass_pa_long %>%
  filter(model == best_model) %>%
  mutate(
    Dataset = "Targeted",
    Algorithm = "Best"
  )

field_worst <- field_eelgrass_pa_long %>%
  filter(model == worst_model) %>%
  mutate(
    Dataset = "Targeted",
    Algorithm = "Worst"
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
    Algorithm = ifelse(model == "XGBoost_nep", "Best", "Worst"),
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
      levels = c("Best","Worst")
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
    Algorithm = "Best",
    label = paste0("n = ", n),
    y = 0.75
  )

n_labels <- n_labels %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = c("SDM building", "Independent", "Targeted")
    ),
    Algorithm = factor(
      Algorithm,
      levels = c("Best", "Worst")
    )
  )

pd <- position_dodge(width = 0.9)

eelgrass <- ggplot(plot_dat,
       aes(x = Presence,
           y = predicted_suitability,
           fill = Tidal_zone)) +
  
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
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    angle = 30,
    hjust = 0,
    vjust = 0,
    size = 2.5
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
      Algorithm = "Best"
    )
  
  training_worst_surf <- test_surf_pa_long %>%
    filter(model == worst_model_surf) %>%
    mutate(
      Dataset = "SDM building",
      Algorithm = "Worst"
    )
  
  independent_surf <- indep_surf_pa_long %>%
    filter(model == best_model_surf) %>%
    mutate(
      Dataset = "Independent",
      Algorithm = "Best"
    )
  
  independent_worst_surf <- indep_surf_pa_long %>%
    filter(model == worst_model_surf) %>%
    mutate(
      Dataset = "Independent",
      Algorithm = "Worst"
    )
  
  field_surf <- field_surfgrass_pa_long %>%
    filter(model == best_model_surf) %>%
    mutate(
      Dataset = "Targeted",
      Algorithm = "Best"
    )
  
  field_worst_surf <- field_surfgrass_pa_long %>%
    filter(model == worst_model_surf) %>%
    mutate(
      Dataset = "Targeted",
      Algorithm = "Worst"
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
      Algorithm = ifelse(model == "XGBoost_nep", "Best", "Worst"),
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
        levels = c("Best","Worst")
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
      Algorithm = "Best",
      label = paste0("n = ", n),
      y = 0.75
    )
  
  n_labels_surf <- n_labels_surf %>%
    mutate(
      Dataset = factor(
        Dataset,
        levels = c("SDM building", "Independent", "Targeted")
      ),
      Algorithm = factor(
        Algorithm,
        levels = c("Best", "Worst")
      )
    )
  
  
  surfgrass <- ggplot(plot_dat_surf,
                     aes(x = Presence,
                         y = predicted_suitability,
                         fill = Tidal_zone)) +
    
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
        group = Tidal_zone,
        label = label
      ),
      position = position_dodge(width = 0.75),
      inherit.aes = FALSE,
      angle = 30,
      hjust = 0,
      vjust = 0,
      size = 2.5
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
  
  
  
 
  
  occurence_plot <- eelgrass / surfgrass +
    plot_layout(
      heights = c(3, 3.5))
  
      occurence_plot   
 
  ggsave(
    filename = "figures/occurence_plots.png",
    plot = occurence_plot,
    width = 180,
    height = 210,
    units = "mm",
    dpi = 600,
    bg = "white" )
  
  
  
  
  






# Plot
# eelgrass_field <- ggplot(eelgrass_pa_long,
#        aes(x = Presence,
#            y = predicted_suitability,
#            fill = Tidal_zone)) +
#   
#   geom_boxplot(
#     aes(group = interaction(Presence, Tidal_zone)),
#     position = position_dodge(width = 0.75),
#     outlier.alpha = 0.3,
#     width = 0.6
#   ) +
#   
#   
#   facet_wrap(~ model, ncol = 4) +
#   
#   geom_hline(
#     data = threshold_df_eelgrass,
#     aes(yintercept = threshold_spatial),
#     linetype = "dashed",
#     colour = "black",
#     linewidth = 0.8,
#     inherit.aes = FALSE
#   ) +
#   scale_fill_manual(
#     name = "Tidal zone",
#     values = c(
#       "Subtidal" = "#0B5D5E",   # dark teal
#       "Intertidal" = "#6EC6C4"  # light teal
#     )
#   ) +
#   
#   scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
#   #coord_cartesian(ylim = c(0, 1))
#   
#   labs(
#     title = "Eelgrass",
#     x = NULL,
#     y = "Relative probability of occurrence"
#   ) +
#   
#   theme_bw() +
#   theme(
#     panel.grid = element_blank(),   # remove grid lines
#     panel.background = element_rect(fill = "white"),
#     plot.background = element_rect(fill = "white"),
#     legend.position = "none",
#     axis.text.x = element_text(size = 10),
#     strip.text = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0, face = "bold"),
#     panel.spacing = unit(0, "lines")  
#   )
# 
# 
# surfgrass_field <- ggplot(surfgrass_pa_long,
#                           aes(x = Presence,
#                               y = predicted_suitability,
#                               fill = Tidal_zone)) +
#   
#   geom_boxplot(
#     aes(group = interaction(Presence, Tidal_zone)),
#     position = position_dodge(width = 0.75),
#     outlier.alpha = 0.3,
#     width = 0.6
#   ) +
#   
#   facet_wrap(~ model, ncol = 4) +
#   
#   geom_hline(
#     data = threshold_df_surfgrass,
#     aes(yintercept = threshold_spatial),
#     linetype = "dashed",
#     colour = "black",
#     linewidth = 0.8,
#     inherit.aes = FALSE
#   ) +
#   
#   scale_fill_manual(
#     name = "Tidal zone",
#     values = c(
#       "Subtidal" = "#0B5D5E",   # dark teal
#       "Intertidal" = "#6EC6C4"  # light teal
#     )
#   ) +
#   
#   scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
#   #coord_cartesian(ylim = c(0, 1))
#   
#   labs(
#     title = "Surfgrass",
#     x = NULL,
#     y = NULL
#   ) +
#   
#   theme_bw() +
#   theme(
#     panel.grid = element_blank(),   # remove grid lines
#     panel.background = element_rect(fill = "white"),
#     plot.background = element_rect(fill = "white"),
#     legend.position = "right",
#     axis.text.x = element_text(size = 10),
#     strip.text = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0, face = "bold"),
#     panel.spacing = unit(0, "lines")  
#   )
# 
#  field_val<-eelgrass_field + surfgrass_field
# 
#  ggsave(
#    filename = "figures/field_validation.png",
#    plot = field_val,
#    width = 13,
#    height = 12,
#    dpi = 300,
#    bg = "white"
#  )
# 
#  
#  
 
 
 
 
 
 
 
 
 
 
 
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
 selected_files_2013_2023 <- files_2013_2023[c(5, 6, 17, 18, 22, 23, 24, 25, 26, 29, 30, 35, 36, 37, 38, 39, 44, 45, 46)]
 hindcast2013_2023 <- terra::rast(selected_files_2013_2023)
 
 folder_1993_2002 <- "code/output_data/processed_ocean_variables/years_1993-2002"
 files_1993_2002 <- list.files(folder_1993_2002, pattern = "\\.tif$", full.names = TRUE)
 selected_files_1993_2002 <- files_1993_2002[c(5, 6, 17, 18, 22, 23, 24, 25, 26, 29, 30, 35, 36, 37, 38, 39, 44, 45, 46)]
 hindcast1993_2002 <- terra::rast(selected_files_1993_2002)
 
 folder_2003_2012 <- "code/output_data/processed_ocean_variables/years_2003-2012"
 files_2003_2012 <- list.files(folder_2003_2012, pattern = "\\.tif$", full.names = TRUE)
 selected_files_2003_2012 <- files_2003_2012[c(5, 6, 17, 18, 22, 23, 24, 25, 26, 29, 30, 35, 36, 37, 38, 39, 44, 45, 46)]
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
 env_20m_hg$temp_5m_min_bccm2013_2023 <- hold$temp_5m_min_bccm
 env_20m_hg$temp_5m_min_nep362013_2023 <- hold$temp_5m_min_nep36
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
 env_20m_hg$temp_5m_min_bccm1993_2002 <- hold$temp_5m_min_bccm
 env_20m_hg$temp_5m_min_nep361993_2002 <- hold$temp_5m_min_nep36
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
 env_20m_hg$temp_5m_min_bccm2003_2012 <- hold$temp_5m_min_bccm
 env_20m_hg$temp_5m_min_nep362003_2012 <- hold$temp_5m_min_nep36
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
 env_20m_ncc$temp_5m_min_bccm2013_2023 <- hold$temp_5m_min_bccm
 env_20m_ncc$temp_5m_min_nep362013_2023 <- hold$temp_5m_min_nep36
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
 env_20m_ncc$temp_5m_min_bccm1993_2002 <- hold$temp_5m_min_bccm
 env_20m_ncc$temp_5m_min_nep361993_2002 <- hold$temp_5m_min_nep36
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
 env_20m_ncc$temp_5m_min_bccm2003_2012 <- hold$temp_5m_min_bccm
 env_20m_ncc$temp_5m_min_nep362003_2012 <- hold$temp_5m_min_nep36
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
 env_20m_qcs$temp_5m_min_bccm2013_2023 <- hold$temp_5m_min_bccm
 env_20m_qcs$temp_5m_min_nep362013_2023 <- hold$temp_5m_min_nep36
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
 env_20m_qcs$temp_5m_min_bccm1993_2002 <- hold$temp_5m_min_bccm
 env_20m_qcs$temp_5m_min_nep361993_2002 <- hold$temp_5m_min_nep36
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
 env_20m_qcs$temp_5m_min_bccm2003_2012 <- hold$temp_5m_min_bccm
 env_20m_qcs$temp_5m_min_nep362003_2012 <- hold$temp_5m_min_nep36
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
 env_20m_ss$temp_5m_min_bccm2013_2023 <- hold$temp_5m_min_bccm
 env_20m_ss$temp_5m_min_nep362013_2023 <- hold$temp_5m_min_nep36
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
 env_20m_ss$temp_5m_min_bccm1993_2002 <- hold$temp_5m_min_bccm
 env_20m_ss$temp_5m_min_nep361993_2002 <- hold$temp_5m_min_nep36
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
 env_20m_ss$temp_5m_min_bccm2003_2012 <- hold$temp_5m_min_bccm
 env_20m_ss$temp_5m_min_nep362003_2012 <- hold$temp_5m_min_nep36
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
 env_20m_wcvi$temp_5m_min_bccm2013_2023 <- hold$temp_5m_min_bccm
 env_20m_wcvi$temp_5m_min_nep362013_2023 <- hold$temp_5m_min_nep36
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
 env_20m_wcvi$temp_5m_min_bccm1993_2002 <- hold$temp_5m_min_bccm
 env_20m_wcvi$temp_5m_min_nep361993_2002 <- hold$temp_5m_min_nep36
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
 env_20m_wcvi$temp_5m_min_bccm2003_2012 <- hold$temp_5m_min_bccm
 env_20m_wcvi$temp_5m_min_nep362003_2012 <- hold$temp_5m_min_nep36
 env_20m_wcvi$temp_5m_mean_bccm2003_2012 <- hold$temp_5m_mean_bccm
 env_20m_wcvi$temp_5m_mean_nep362003_2012 <- hold$temp_5m_mean_nep36
 env_20m_wcvi$NH4_5m_mean_bccm2003_2012 <- hold$NH4_5m_mean_bccm
 env_20m_wcvi$NH4_5m_mean_nep362003_2012 <- hold$NH4_5m_mean_nep36
 env_20m_wcvi$temp_air_cv2003_2012 <- hold$temp_air_cv
 env_20m_wcvi$temp_s_cv_bccm2003_2012 <- hold$temp_s_cv_bccm
 env_20m_wcvi$temp_s_cv_nep362003_2012 <- hold$temp_s_cv_nep36
 
 env_20m_all <- bind_rows(env_20m_hg, env_20m_ncc, env_20m_qcs, env_20m_wcvi,env_20m_ss)
 
 env_20m_all <- env_20m_all %>%
   filter(salt_5m_mean_bccm2013_2023 > quantile(seagrass_data$saltmean_bccm, probs = 0.001))
     
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
 
 
 var_crosswalk <- tribble(
   ~plot_group, ~variable_label, ~model_source, ~env_variable,              ~obs_variable,        ~units,                  ~convert,
   "Atmospheric", "Mean precipitation",        "single", "precip_mean",           "prmean",             "kg m-2 month-1",       "none",
   "Atmospheric", "Minimum precipitation",     "single", "precip_min",            "prmin",              "kg m-2 month-1",       "none",
   "Atmospheric", "Minimum shortwave flux",    "single", "rsds_min",              "rsdsmin",            "MJ m-2",               "none",
   "Atmospheric", "Minimum air temperature",   "single", "temp_air_min",          "airtempmin",         "degrees C",            "kelvin_to_c",
   "Atmospheric", "Air temperature CV",        "single", "temp_air_cv",           "airtempcv",          "CV",                   "none",
   
   "Ocean", "Mean salinity",                   "BCCM",   "salt_5m_mean_bccm",     "saltmean_bccm",      "psu",                  "none",
   "Ocean", "Mean salinity",                   "NEP36",  "salt_5m_mean_nep36",    "saltmean_nep",       "psu",                  "none",
   "Ocean", "Salinity CV",                     "BCCM",   "salt_5m_cv_bccm",       "saltcv_bccm",        "CV",                   "none",
   "Ocean", "Salinity CV",                     "NEP36",  "salt_5m_cv_nep36",      "saltcv_nep",         "CV",                   "none",
   
   "Ocean", "Mean subsurface temperature",     "BCCM",   "temp_5m_mean_bccm",     "tempmean_bccm",      "degrees C",            "none",
   "Ocean", "Mean subsurface temperature",     "NEP36",  "temp_5m_mean_nep36",    "tempmean_nep",       "degrees C",            "none",
   "Ocean", "Minimum subsurface temperature",  "BCCM",   "temp_5m_min_bccm",      "tempmin_bccm",       "degrees C",            "none",
   "Ocean", "Minimum subsurface temperature",  "NEP36",  "temp_5m_min_nep36",     "tempmin_nep",        "degrees C",            "none",
   "Ocean", "Subsurface temperature CV",       "BCCM",   "temp_5m_cv_bccm",       "tempcv_bccm",        "CV",                   "none",
   "Ocean", "Subsurface temperature CV",       "NEP36",  "temp_5m_cv_nep36",      "tempcv_nep",         "CV",                   "none",
   
   "Ocean", "Mean ammonium",                   "BCCM",   "NH4_5m_mean_bccm",      "NH4_bccm",           "mmol m-3",             "none",
   "Ocean", "Mean ammonium",                   "NEP36",  "NH4_5m_mean_nep36",     "NH4_nep",            "mmol m-3",             "none",
   
   "Ocean", "Surface temperature CV",          "BCCM",   "temp_s_cv_bccm",        "surftempcv_bccm",    "CV",                   "none",
   "Ocean", "Surface temperature CV",          "NEP36",  "temp_s_cv_nep36",       "surftempcv_nep",     "CV",                   "none"
 )
 
 
 # Study-area / hindcast values
 hindcast_compare <- var_crosswalk %>%
   select(
     plot_group, variable_label, model_source, env_variable, units, convert
   ) %>%
   pmap_dfr(function(plot_group, variable_label, model_source,
                     env_variable, units, convert) {
     
     env_long %>%
       filter(variable == env_variable) %>%
       transmute(
         plot_group = plot_group,
         variable_label = variable_label,
         model_source = model_source,
         decade = decade,
         source_type = "Hindcast",
         source = if_else(
           model_source == "single",
           "Hindcast",
           paste("Hindcast", model_source)
         ),
         value = value,
         units = units,
         convert = convert
       )
   })
 
 
 # SDM observation / survey values
 survey_compare <- var_crosswalk %>%
   select(
     plot_group, variable_label, model_source, obs_variable, units, convert
   ) %>%
   pmap_dfr(function(plot_group, variable_label, model_source,
                     obs_variable, units, convert) {
     
     obs_decade %>%
       transmute(
         plot_group = plot_group,
         variable_label = variable_label,
         model_source = model_source,
         decade = decade,
         source_type = "Survey observations",
         source = if_else(
           model_source == "single",
           "Survey observations",
           paste("Survey observations", model_source)
         ),
         value = .data[[obs_variable]],
         units = units,
         convert = convert
       ) %>%
       filter(!is.na(value))
   })
 
 
 env_obs_compare <- bind_rows(
   hindcast_compare,
   survey_compare
 ) %>%
   mutate(
     value_plot = case_when(
       convert == "kelvin_to_c" ~ value - 273.15,
       TRUE ~ value
     ),
     decade = factor(
       decade,
       levels = c("1993–2002", "2003–2012", "2013–2023")
     ),
     source_type = factor(
       source_type,
       levels = c("Hindcast", "Survey observations")
     )
   )
 
 
 env_obs_compare %>%
   count(variable_label, model_source, decade, source_type)
 
 
 
 
 
 hindcast_summary <- env_long %>%
   group_by(variable, decade) %>%
   summarise(
     source = "Hindcast",
     n = n(),
     mean = mean(value, na.rm = TRUE),
     sd = sd(value, na.rm = TRUE),
     median = median(value, na.rm = TRUE),
     q10 = quantile(value, 0.10, na.rm = TRUE),
     q25 = quantile(value, 0.25, na.rm = TRUE),
     q75 = quantile(value, 0.75, na.rm = TRUE),
     q90 = quantile(value, 0.90, na.rm = TRUE),
     min = min(value, na.rm = TRUE),
     max = max(value, na.rm = TRUE),
     .groups = "drop"
   )
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 precip_mean_plot <- env_long %>%
   filter(variable == "precip_mean") %>%
   mutate(source = "Hindcast") %>%
   select(decade, value, source)
 
 precip_mean_obs_plot <- obs_decade %>%
   mutate(
     source = "Survey",
     value = prmean
   ) %>%
   filter(!is.na(value)) %>%
   select(decade, value, source)
 

 precip_mean_plot_df <- bind_rows(precip_mean_plot, precip_mean_obs_plot)
 
 precip_mean <- ggplot(precip_mean_plot_df, aes(x = decade, y = value, fill = source)) +
   geom_violin(
     position = position_dodge(width = 0.8),
     trim = TRUE,
     alpha = 0.3
   ) +
   geom_boxplot(
     width = 0.12,
     position = position_dodge(width = 0.8),
     outlier.shape = NA,
     alpha = 0.8
   ) +
   scale_fill_manual(
     values = c(
       "Survey" = "#E69F00",  # orange
       "Hindcast"  = "#56B4E9"  # sky blue
     )
   ) +
   labs(
     x = "Decade",
     y = expression("Precipitation (kg m"^{-2}~"month"^{-1}*")"),
     fill = NULL
   ) +
   theme_bw() +
   theme(
     panel.grid = element_blank(),   # remove grid lines
     panel.background = element_rect(fill = "white"),
     plot.background = element_rect(fill = "white"),
     #legend.position = "right",
     #axis.text.x = element_text(size = 10),
     strip.text = element_text(face = "bold"),
     plot.title = element_text(hjust = 0, face = "bold"),
     axis.title.x = element_blank(),
     axis.text.x  = element_blank(),
     axis.ticks.x = element_blank(),
     panel.spacing = unit(0, "lines"))  
 
 
 
 precip_min_plot <- env_long %>%
   filter(variable == "precip_min") %>%
   mutate(source = "Hindcast") %>%
   select(decade, value, source)
 
 precip_min_obs_plot <- obs_decade %>%
   mutate(
     source = "Survey",
     value = prmin
   ) %>%
   filter(!is.na(value)) %>%
   select(decade, value, source)
 
 
 precip_min_plot_df <- bind_rows(precip_min_plot, precip_min_obs_plot)
 
 precip_min <- ggplot(precip_min_plot_df, aes(x = decade, y = value, fill = source)) +
   
   geom_violin(
     position = position_dodge(width = 0.8),
     trim = TRUE,
     alpha = 0.3
   ) +
   
   geom_boxplot(
     width = 0.12,
     position = position_dodge(width = 0.8),
     outlier.shape = NA,
     alpha = 0.8
   ) +
   scale_fill_manual(
     values = c(
       "Survey" = "#E69F00",  # orange
       "Hindcast"  = "#56B4E9"  # sky blue
     )
   ) +
   labs(
     x = "Decade",
     y = expression("Min precipitation (kg m"^{-2}~"month"^{-1}*")"),
     fill = NULL
   ) +
   
   theme_bw() +
   theme(
     panel.grid = element_blank(),   # remove grid lines
     panel.background = element_rect(fill = "white"),
     plot.background = element_rect(fill = "white"),
     #legend.position = "right",
     #axis.text.x = element_text(size = 10),
     strip.text = element_text(face = "bold"),
     plot.title = element_text(hjust = 0, face = "bold"),
     axis.title.x = element_blank(),
     axis.text.x  = element_blank(),
     axis.ticks.x = element_blank(),
     panel.spacing = unit(0, "lines"))  
 
 

 
 rsds_min_plot <- env_long %>%
   filter(variable == "rsds_min") %>%
   mutate(source = "Hindcast") %>%
   select(decade, value, source)
 
 rsds_min_obs_plot <- obs_decade %>%
   mutate(
     source = "Survey",
     value = rsdsmin
   ) %>%
   filter(!is.na(value)) %>%
   select(decade, value, source)
 
  rsds_min_plot_df <- bind_rows(rsds_min_plot, rsds_min_obs_plot)
  
  rsds_min <- ggplot(rsds_min_plot_df, aes(x = decade, y = value, fill = source)) +
   geom_violin(
     position = position_dodge(width = 0.8),
     trim = TRUE,
     alpha = 0.3
   ) +
   geom_boxplot(
     width = 0.12,
     position = position_dodge(width = 0.8),
     outlier.shape = NA,
     alpha = 0.8
   ) +
    scale_fill_manual(
      values = c(
        "Survey" = "#E69F00",  # orange
        "Hindcast"  = "#56B4E9"  # sky blue
      )
    ) +
   labs(
     x = "Decade",
     y = expression("Min Surface Downwelling \nShortwave Flux (MJ m"^{-2}~")"),
     fill = NULL
   ) +
   theme_bw() +
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      #legend.position = "right",
      #axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"))  
  
  #temp_air_min 
  airtemp_min_plot <- env_long %>%
    filter(variable == "temp_air_min") %>%
    mutate(source = "Hindcast") %>%
    select(decade, value, source)
  
  airtemp_min_obs_plot <- obs_decade %>%
    mutate(
      source = "Survey",
      value = airtempmin
    ) %>%
    filter(!is.na(value)) %>%
    select(decade, value, source)
  
  airtemp_min_plot_df <- bind_rows(airtemp_min_plot, airtemp_min_obs_plot)
  
  airtemp_min <- ggplot(airtemp_min_plot_df, aes(x = decade, y = value - 273.15, fill = source)) + # convert K to degrees C
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Survey" = "#E69F00",  # orange
        "Hindcast"  = "#56B4E9"  # sky blue
      )
    ) +
    labs(
      x = "Decade",
      y = expression("Min air temperature ("*degree*C*")"),
      fill = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      #legend.position = "right",
      #axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"))  
 
  
  #temp_air_cv
  airtemp_cv_plot <- env_long %>%
    filter(variable == "temp_air_cv") %>%
    mutate(source = "Hindcast") %>%
    select(decade, value, source)
  
  airtemp_cv_obs_plot <- obs_decade %>%
    mutate(
      source = "Survey",
      value = airtempcv
    ) %>%
    filter(!is.na(value)) %>%
    select(decade, value, source)
  
  airtemp_cv_plot_df <- bind_rows(airtemp_cv_plot, airtemp_cv_obs_plot)
  
  airtemp_cv <- ggplot(airtemp_cv_plot_df, aes(x = decade, y = value, fill = source)) +
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Survey" = "#E69F00",  # orange
        "Hindcast"  = "#56B4E9"  # sky blue
      )
    ) +
    labs(
      x = "Decade",
      y = expression("Air temperature CV"),
      fill = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.position = "right",
      axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      panel.spacing = unit(0, "lines"))
 
  
  # salt_5m_mean_bccm  salt_5m_mean_nep36
  salt_mean_plot <- env_long %>%
    filter(variable %in% c("salt_5m_mean_bccm", "salt_5m_mean_nep36")) %>%
    mutate(
      source = case_when(
        variable == "salt_5m_mean_bccm" ~ "Hindcast BCCM",
        variable == "salt_5m_mean_nep36" ~ "Hindcast NEP36"
      )
    ) %>%
    select(decade, value, source)
  
  salt_mean_obs_plot <- obs_decade %>%
    select(
      decade,
      saltmean_bccm,
      saltmean_nep
    ) %>%
    pivot_longer(
      cols = c(saltmean_bccm, saltmean_nep),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = case_when(
        source == "saltmean_bccm" ~ "Survey BCCM",
        source == "saltmean_nep" ~ "Survey NEP36"
      )
    ) %>%
    filter(!is.na(value))
  
  
  salt_mean_plot_df <- bind_rows(
    salt_mean_plot,
    salt_mean_obs_plot
  )

  salt_mean <- ggplot(
    salt_mean_plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Hindcast BCCM" = "#56B4E9",  
        "Hindcast NEP36"  = "#009E73",  
        "Survey BCCM"   = "#E69F00",  
        "Survey NEP36"    = "#CC79A7"   
      )
    ) +
    coord_cartesian(ylim = c(20, NA)) +
    labs(
      x = "Decade",
      y = expression("Salinity (psu)"),
      fill = NULL
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      #legend.position = "right",
      #axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"))  
  
   #salt_5m_cv_bccm  salt_5m_cv_nep36
  salt_cv_plot <- env_long %>%
    filter(variable %in% c("salt_5m_cv_bccm", "salt_5m_cv_nep36")) %>%
    mutate(
      source = case_when(
        variable == "salt_5m_cv_bccm" ~ "Hindcast BCCM",
        variable == "salt_5m_cv_nep36" ~ "Hindcast NEP36"
      )
    ) %>%
    select(decade, value, source)
  
  salt_cv_obs_plot <- obs_decade %>%
    select(
      decade,
      saltcv_bccm,
      saltcv_nep
    ) %>%
    pivot_longer(
      cols = c(saltcv_bccm, saltcv_nep),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = case_when(
        source == "saltcv_bccm" ~ "Survey BCCM",
        source == "saltcv_nep" ~ "Survey NEP36"
      )
    ) %>%
    filter(!is.na(value))
  
  
  salt_cv_plot_df <- bind_rows(
    salt_cv_plot,
    salt_cv_obs_plot
  )
  
  salt_cv <- ggplot(
    salt_cv_plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Hindcast BCCM" = "#56B4E9",  
        "Hindcast NEP36"  = "#009E73",  
        "Survey BCCM"   = "#E69F00",  
        "Survey NEP36"    = "#CC79A7"   
      )
    ) +
    coord_cartesian(ylim = c(0, 20)) +
    labs(
      x = "Decade",
      y = expression("Salinity CV"),
      fill = NULL
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      #legend.position = "right",
      #axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"))  
  
  #temp_5m_cv_bccm temp_5m_cv_nep36
  
  temp_cv_plot <- env_long %>%
    filter(variable %in% c("temp_5m_cv_bccm", "temp_5m_cv_nep36")) %>%
    mutate(
      source = case_when(
        variable == "temp_5m_cv_bccm" ~ "Hindcast BCCM",
        variable == "temp_5m_cv_nep36" ~ "Hindcast NEP36"
      )
    ) %>%
    select(decade, value, source)
  
  temp_cv_obs_plot <- obs_decade %>%
    select(
      decade,
      tempcv_bccm,
      tempcv_nep
    ) %>%
    pivot_longer(
      cols = c(tempcv_bccm, tempcv_nep),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = case_when(
        source == "tempcv_bccm" ~ "Survey BCCM",
        source == "tempcv_nep" ~ "Survey NEP36"
      )
    ) %>%
    filter(!is.na(value))
  
  
  temp_cv_plot_df <- bind_rows(
    temp_cv_plot,
    temp_cv_obs_plot
  )
  
  temp_cv <- ggplot(
    temp_cv_plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Hindcast BCCM" = "#56B4E9",  
        "Hindcast NEP36"  = "#009E73",  
        "Survey BCCM"   = "#E69F00",  
        "Survey NEP36"    = "#CC79A7"   
      )
    ) +
    labs(
      x = "Decade",
      y = expression("Ocean subsurface temperature CV"),
      fill = NULL
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      #legend.position = "right",
      #axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"))  
  
  #temp_5m_min_bccm  temp_5m_min_nep36
  
  temp_min_plot <- env_long %>%
    filter(variable %in% c("temp_5m_min_bccm", "temp_5m_min_nep36")) %>%
    mutate(
      source = case_when(
        variable == "temp_5m_min_bccm" ~ "Hindcast BCCM",
        variable == "temp_5m_min_nep36" ~ "Hindcast NEP36"
      )
    ) %>%
    select(decade, value, source)
  
  temp_min_obs_plot <- obs_decade %>%
    select(
      decade,
      tempmin_bccm,
      tempmin_nep
    ) %>%
    pivot_longer(
      cols = c(tempmin_bccm, tempmin_nep),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = case_when(
        source == "tempmin_bccm" ~ "Survey BCCM",
        source == "tempmin_nep" ~ "Survey NEP36"
      )
    ) %>%
    filter(!is.na(value))
  
  
  temp_min_plot_df <- bind_rows(
    temp_min_plot,
    temp_min_obs_plot
  )
  
  temp_min <- ggplot(
    temp_min_plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Hindcast BCCM" = "#56B4E9",  
        "Hindcast NEP36"  = "#009E73",  
        "Survey BCCM"   = "#E69F00",  
        "Survey NEP36"    = "#CC79A7"   
      )
    ) +
    labs(
      x = "Decade",
      y = expression("Ocean subsurface temperature min ("*degree*C*")"),
      fill = NULL
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      #legend.position = "right",
      #axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"))  
  
  #temp_5m_mean_bccm  temp_5m_mean_nep36
  temp_mean_plot <- env_long %>%
    filter(variable %in% c("temp_5m_mean_bccm", "temp_5m_mean_nep36")) %>%
    mutate(
      source = case_when(
        variable == "temp_5m_mean_bccm" ~ "Hindcast BCCM",
        variable == "temp_5m_mean_nep36" ~ "Hindcast NEP36"
      )
    ) %>%
    select(decade, value, source)
  
  temp_mean_obs_plot <- obs_decade %>%
    select(
      decade,
      tempmean_bccm,
      tempmean_nep
    ) %>%
    pivot_longer(
      cols = c(tempmean_bccm, tempmean_nep),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = case_when(
        source == "tempmean_bccm" ~ "Survey BCCM",
        source == "tempmean_nep" ~ "Survey NEP36"
      )
    ) %>%
    filter(!is.na(value))
  
  
  temp_mean_plot_df <- bind_rows(
    temp_mean_plot,
    temp_mean_obs_plot
  )
  
  temp_mean <- ggplot(
    temp_mean_plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Hindcast BCCM" = "#56B4E9",  
        "Hindcast NEP36"  = "#009E73",  
        "Survey BCCM"   = "#E69F00",  
        "Survey NEP36"    = "#CC79A7"   
      )
    ) +
    labs(
      x = "Decade",
      y = expression("Ocean subsurface temperature ("*degree*C*")"),
      fill = NULL
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),   # remove grid lines
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      #legend.position = "right",
      #axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0, "lines"))  
  
  #NH4_5m_mean_bccm  NH4_5m_mean_nep36
  NH4_mean_plot <- env_long %>%
    filter(variable %in% c("NH4_5m_mean_bccm", "NH4_5m_mean_nep36")) %>%
    mutate(
      source = case_when(
        variable == "NH4_5m_mean_bccm" ~ "Hindcast BCCM",
        variable == "NH4_5m_mean_nep36" ~ "Hindcast NEP36"
      )
    ) %>%
    select(decade, value, source)
  
  NH4_mean_obs_plot <- obs_decade %>%
    select(
      decade,
      NH4_bccm,
      NH4_nep
    ) %>%
    pivot_longer(
      cols = c(NH4_bccm, NH4_nep),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = case_when(
        source == "NH4_bccm" ~ "Survey BCCM",
        source == "NH4_nep" ~ "Survey NEP36"
      )
    ) %>%
    filter(!is.na(value))
  
  
  NH4_mean_plot_df <- bind_rows(
    NH4_mean_plot,
    NH4_mean_obs_plot
  )
  
   NH4_mean <- ggplot(
    NH4_mean_plot_df,
    aes(x = decade, y = value, fill = source)
  ) +
    
    geom_violin(
      position = position_dodge(width = 0.8),
      trim = TRUE,
      alpha = 0.3
    ) +
    
    geom_boxplot(
      width = 0.12,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      alpha = 0.8
    ) +
     scale_fill_manual(
       values = c(
         "Hindcast BCCM" = "#56B4E9",  
         "Hindcast NEP36"  = "#009E73",  
         "Survey BCCM"   = "#E69F00",  
         "Survey NEP36"    = "#CC79A7"   
       )
     ) +
     coord_cartesian(ylim = c(0, 2.5)) +
    labs(
      x = "Decade",
      y = expression("Ammonium (mmol m"^{-3}~")"),
      fill = NULL
    ) +
     
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.position = "right",
      axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      panel.spacing = unit(0, "lines")
    )
  
  
    #temp_s_cv_bccm temp_s_cv_nep36
   surfacetemp_cv_plot <- env_long %>%
     filter(variable %in% c("temp_s_cv_bccm", "temp_s_cv_nep36")) %>%
     mutate(
       source = case_when(
         variable == "temp_s_cv_bccm" ~ "Hindcast BCCM",
         variable == "temp_s_cv_nep36" ~ "Hindcast NEP36"
       )
     ) %>%
     select(decade, value, source)
   
   surfacetemp_cv_obs_plot <- obs_decade %>%
     select(
       decade,
       surftempcv_bccm,
       surftempcv_nep
     ) %>%
     pivot_longer(
       cols = c(surftempcv_bccm, surftempcv_nep),
       names_to = "source",
       values_to = "value"
     ) %>%
     mutate(
       source = case_when(
         source == "surftempcv_bccm" ~ "Survey BCCM",
         source == "surftempcv_nep" ~ "Survey NEP36"
       )
     ) %>%
     filter(!is.na(value))
   
   
   surfacetemp_cv_plot_df <- bind_rows(
     surfacetemp_cv_plot,
     surfacetemp_cv_obs_plot
   )
   surfacetemp_cv <- ggplot(
     surfacetemp_cv_plot_df,
     aes(x = decade, y = value, fill = source)
   ) +
     
     geom_violin(
       position = position_dodge(width = 0.8),
       trim = TRUE,
       alpha = 0.3
     ) +
     
     geom_boxplot(
       width = 0.12,
       position = position_dodge(width = 0.8),
       outlier.shape = NA,
       alpha = 0.8
     ) +
     scale_fill_manual(
       values = c(
         "Hindcast BCCM" = "#56B4E9",  
         "Hindcast NEP36"  = "#009E73",  
         "Survey BCCM"   = "#E69F00",  
         "Survey NEP36"    = "#CC79A7"   
       )
     ) +
     labs(
       x = "Decade",
       y = expression("Ocean surface temperature CV"),
       fill = NULL
     ) +
     
     theme_bw() +
     
     theme(
       panel.grid = element_blank(),   # remove grid lines
       panel.background = element_rect(fill = "white"),
       plot.background = element_rect(fill = "white"),
       #legend.position = "right",
       #axis.text.x = element_text(size = 10),
       strip.text = element_text(face = "bold"),
       plot.title = element_text(hjust = 0, face = "bold"),
       axis.title.x = element_blank(),
       axis.text.x  = element_blank(),
       axis.ticks.x = element_blank(),
       panel.spacing = unit(0, "lines"))  
   
   
   
 atmos <- precip_mean / precip_min / rsds_min / airtemp_min / airtemp_cv + plot_layout(guides = "collect") & theme(legend.position = "right")
 
 ggsave(
   filename = "figures/atmos_var.png",
   plot = atmos,
   width = 12,
   height = 20,
   dpi = 300,
   bg = "white"
 )
 
 ocean<- salt_mean / salt_cv / temp_mean / temp_min / temp_cv / surfacetemp_cv / NH4_mean + plot_layout(guides = "collect") & theme(legend.position = "right")
 
 ggsave(
   filename = "figures/ocean_var.png",
   plot = ocean,
   width = 14,
   height = 20,
   dpi = 300,
   bg = "white"
 )
 