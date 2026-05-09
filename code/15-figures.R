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

coastline_full <- st_read("raw_data/CHS_HWL2015_Coastline.gdb", layer = "Polygon_CHS_Pacific_HWL_2015_5028437_simple")
coastline <- coastline_full %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_crop(st_bbox(coastline_full) + c(0, 100000, -10000, -80000)) %>%
  st_transform(crs = "EPSG:3005")



#spatial and temporal blocks
load("code/output_data/seagrass_sf_file.RData")

jewel_cb_10 <- c(
  "#0F766E", # deep teal
  "#3B5BA5", # sapphire
  "#7A5195", # amethyst
  "#C0567D", # rose
  "#DD6B20", # amber
  "#2F855A", # emerald
  "#4C6A92", # steel blue
  "#9C4221", # copper
  "#6B46C1", # violet
  "#B7791F"  # ochre gold
)

base_map_theme <- theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2)
  )
p_eelgrass <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    color = "grey50",
    linewidth = 0.35
  ) +
  geom_sf(
    data = spatialised_sf,
    aes(color = factor(fold_eelgrass)),
    size = 0.8
  ) +
  scale_color_manual(
    values = jewel_cb_10,
    name = "Fold"
  ) +
  coord_sf(expand = FALSE) +
  labs(title = "Eelgrass") +
  base_map_theme

p_surfgrass <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    color = "grey50",
    linewidth = 0.35
  ) +
  geom_sf(
    data = spatialised_sf,
    aes(color = factor(fold_seagrass)),
    size = 0.8
  ) +
  scale_color_manual(
    values = jewel_cb_10,
    name = "Fold"
  ) +
  coord_sf(expand = FALSE) +
  labs(title = "Surfgrass") +
  base_map_theme
combined_plot <- p_eelgrass + p_surfgrass +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
combined_plot


spatialised_sf <- spatialised_sf %>%
  mutate(
    year_group = case_when(
      Year >= 1993 & Year <= 2009 ~ "1993–2009",
      Year >= 2012 & Year <= 2023 ~ "2012–2023",
      TRUE ~ NA_character_
    ),
    year_group = factor(year_group, levels = c("1993–2009", "2012–2023"))
  )
p_year_groups <- ggplot() +
  geom_sf(
    data = coastline,
    fill = "grey96",
    color = "grey50",
    linewidth = 0.35
  ) +
  geom_sf(
    data = filter(spatialised_sf, year_group == "1993–2009"),
    aes(color = year_group),
    size = 1.2
  ) +
  geom_sf(
    data = filter(spatialised_sf, year_group == "2012–2023"),
    aes(color = year_group),
    size = 0.7,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "1993–2009" = "#F58518",
      "2012–2023" = "#4C78A8"
    ),
    name = "Year"
  ) +
  coord_sf(expand = FALSE) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2)
  )
p_year_groups

fold_schemes <- combined_plot/p_year_groups

ggsave(
  filename = "figures/fold_schemes.png",
  plot = fold_schemes,
  width = 6,
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
  mutate(model = "GLMM_bccm")

glmmspatial_varimp_bccm_eelgrass <- relimp_e_bccm_spatial[[2]]
glmmspatial_varimp_bccm_eelgrass <- glmmspatial_varimp_bccm_eelgrass %>%
  rename(variable = term) %>%
  mutate(model = "GLMM_spatial_bccm")

glmm_varimp_nep_eelgrass <- relimp_e_nep_nospatial[[2]]
glmm_varimp_nep_eelgrass <- glmm_varimp_nep_eelgrass %>%
  rename(variable = term) %>%
  mutate(model = "GLMM_nep")

glmmspatial_varimp_nep_eelgrass <- relimp_e_nep_spatial[[2]]
glmmspatial_varimp_nep_eelgrass <- glmmspatial_varimp_nep_eelgrass %>%
  rename(variable = term) %>%
  mutate(model = "GLMM_spatial_nep")

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
  mutate(model = "XGB_bccm")

xgb_varimp_nep_eelgrass <- xgb_varimp_nep_eelgrass %>%
  select(Feature, RelImportance) %>%
  rename(variable = Feature, 
         relimp = RelImportance) %>%
  mutate(model = "XGB_nep")

eelgrass_relimp <- rbind(glmm_varimp_bccm_eelgrass, glmmspatial_varimp_bccm_eelgrass, glmm_varimp_nep_eelgrass, glmmspatial_varimp_nep_eelgrass,
                         gbm_varimp_bccm_eelgrass, gbm_varimp_nep_eelgrass, xgb_varimp_bccm_eelgrass, xgb_varimp_nep_eelgrass)


library(viridis)
eelgrass_relimp$model <- factor(
  eelgrass_relimp$model,
  levels = c("GLMM_bccm", "GLMM_spatial_bccm", "GBM_bccm", "XGB_bccm", "GLMM_nep",  "GLMM_spatial_nep", "GBM_nep", "XGB_nep")   
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
                           "NH4_bccm_stnd" = "Ammonium",
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
  mutate(model = "GLMM_bccm")

glmmspatial_varimp_bccm_surfgrass <- relimp_s_bccm_spatial[[2]]
glmmspatial_varimp_bccm_surfgrass <- glmmspatial_varimp_bccm_surfgrass %>%
  rename(variable = term) %>%
  mutate(model = "GLMM_spatial_bccm")

glmm_varimp_nep_surfgrass <- relimp_s_nep_nospatial[[2]]
glmm_varimp_nep_surfgrass <- glmm_varimp_nep_surfgrass %>%
  rename(variable = term) %>%
  mutate(model = "GLMM_nep")

glmmspatial_varimp_nep_surfgrass <- relimp_s_nep_spatial[[2]]
glmmspatial_varimp_nep_surfgrass <- glmmspatial_varimp_nep_surfgrass %>%
  rename(variable = term) %>%
  mutate(model = "GLMM_spatial_nep")

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
  mutate(model = "XGB_bccm")

xgb_varimp_nep_surfgrass <- xgb_varimp_nep_surfgrass %>%
  select(Feature, RelImportance) %>%
  rename(variable = Feature, 
         relimp = RelImportance) %>%
  mutate(model = "XGB_nep")

surfgrass_relimp <- rbind(glmm_varimp_bccm_surfgrass, glmmspatial_varimp_bccm_surfgrass, glmm_varimp_nep_surfgrass, glmmspatial_varimp_nep_surfgrass,
                         gbm_varimp_bccm_surfgrass, gbm_varimp_nep_surfgrass, xgb_varimp_bccm_surfgrass, xgb_varimp_nep_surfgrass)


library(viridis)
surfgrass_relimp$model <- factor(
  surfgrass_relimp$model,
  levels = c("GLMM_bccm", "GLMM_spatial_bccm", "GBM_bccm", "XGB_bccm", "GLMM_nep",  "GLMM_spatial_nep", "GBM_nep", "XGB_nep")   
)

surfgrass_relimp <- surfgrass_relimp %>%
  mutate(variable = recode(variable,
                           "depth_stnd" = "Depth",
                           "substrate" = "Substrate",
                           "tidal_sqrt_stnd" = "Tidal current",
                           "cul_eff_stnd" = "Cumulative impact",
                           "prmean_stnd" = "Precipitation",
                           "saltcv_nep_stnd" = "Subsurface salinity variability",
                           "saltmean_nep_stnd" = "Subsurface salinity",
                           "airtempcv_stnd" = "Air temperature variability",
                           "rsdsmin_stnd" = "Min surface downwelling shortwave flux",
                           "saltcv_bccm_stnd" = "Subsurface salinity variability",
                           "rei_sqrt_stnd" = "Exposure",
                           "tempcv_nep_stnd" = "Subsurface temperature variability",
                           "tempmean_nep_stnd" = "Subsurface temperature",
                           "tempmin_bccm_stnd" = "Min subsurface temperature",
                           "surftempcv_bccm_stnd" = "Surface temperature CV",
                           "surftempcv_nep_stnd" = "Surface temperature CV",
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

ggsave("figures/figure_relimp.png", plot = relimp, width = 8, height = 6, dpi = 300)



##### Independent validation figures
load("code/output_data/model_results/eelgrass_independent_dataframe.RData")
load("code/output_data/model_results/surfgrass_independent_dataframe.RData")
load("code/output_data/model_results/eelgrass_independent_eval.RData")


# Model columns in the order you want shown
model_labels <- c(
  "bccm_nospatial" = "GLMM_bccm",
  "bccm_spatial"   = "GLMM_spatial_bccm",
  "GBM_bccm"       = "GBM_bccm",
  "XGBOOST_bccm"   = "XGB_bccm",
  "nep_nospatial"  = "GLMM_nep",
  "nep_spatial"    = "GLMM_spatial_nep",
  "GBM_nep"        = "GBM_nep",
  "XGBOOST_nep"    = "XGB_nep"
)


eelgrass_long <- eelgrass_df %>%
  #select(all_of(model_names)) %>%
  pivot_longer(
    cols = -c(obs, ID, substrate, rock_group),
    names_to = "model",
    values_to = "predicted_suitability"
  ) %>%
  mutate(species = "Eelgrass")

eelgrass_long <- eelgrass_long %>%
  mutate(model = recode(as.character(model), !!!model_labels),
         model = factor(model, levels = model_labels))

surfgrass_long <- surfgrass_df %>%
  select(-bathy) %>%
  pivot_longer(
    cols = -c(obs, ID, substrate, rock_group),
    names_to = "model",
    values_to = "predicted_suitability"
  ) %>%
  mutate(species = "Surfgrass")

surfgrass_long <- surfgrass_long %>%
  mutate(model = recode(as.character(model), !!!model_labels),
         model = factor(model, levels = model_labels))

plot_df <- bind_rows(eelgrass_long, surfgrass_long) %>%
  mutate(
    model = factor(model, levels = model_labels)
  )

independent <- ggplot(plot_df, aes(x = model, y = predicted_suitability, fill = species)) +
  geom_boxplot(
    outlier.alpha = 0.15,
    width = 0.7,
    colour = "black",
    linewidth = 0.3
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.8,
    colour = "black"
  ) +
  facet_wrap(~ species, ncol = 1, scales = "fixed") +
  
  scale_fill_manual(values = c(
    "Eelgrass" = "grey",
    "Surfgrass" = "grey"
  )) +
  
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  
  labs(
    x = NULL,
    y = "Relative probability of occurrence"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  theme_classic(base_size = 12) +
  theme(
    # --- facet labels (move left, remove box) ---
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(
      face = "bold",
      hjust = 0   # left-align text
    ),
    
    # --- axes styling ---
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(margin = margin(r = 10)),
    
    # --- remove legend ---
    legend.position = "none",
    
    # --- spacing ---
    panel.spacing = unit(0.8, "lines")
  )

ggsave(
  filename = "figures/independent_validation.png",
  plot = independent,
  width = 6,
  height = 6,
  dpi = 300,
  bg = "white"
)


eelgrass_subset <- eelgrass_long %>% filter (model == "XGB_nep" | model == "GLMM_spatial_nep")

model_levels <- c("Lowest-performing model", "Highest-performing model")


summary_df <- eelgrass_subset %>%
  left_join(thresholds, by = "model") %>%
  group_by(model) %>%
  summarise(
    FPPS = mean(predicted_suitability >= threshold, na.rm = TRUE),
    FNR  = mean(predicted_suitability < threshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model_label = recode(
      model,
      "GLMM_spatial_nep" = "Lowest-performing model",
      "XGB_nep" = "Highest-performing model"
    ),
    model_label = factor(model_label, levels = model_levels)
  ) %>%
  pivot_longer(
    cols = c(FPPS, FNR),
    names_to = "metric",
    values_to = "value"
  )

mps_df <- eelgrass_subset %>%
  mutate(
    model_label = recode(
      model,
      "GLMM_spatial_nep" = "Lowest-performing model",
      "XGB_nep" = "Highest-performing model"
    ),
    model_label = factor(model_label, levels = model_levels),
    value = predicted_suitability
  ) %>%
  select(model_label, value)

thresholds <- tibble(
  model = c("GLMM_spatial_nep", "XGB_nep"),
  threshold = c(0.037, 0.031)   
)


p1 <- ggplot(
  mps_df,
  aes(x = model_label, y = value, fill = model_label)
) +
  geom_violin(
    alpha = 0.5,
    colour = NA,
    width = 0.9
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.alpha = 0.1,
    colour = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2,
    colour = "black"
  ) +
  scale_fill_manual(values = c(
    "Highest-performing model" = "#8C6BB1",
    "Lowest-performing model" = "#4D4D4D"
  )) +
  scale_x_discrete(limits = model_levels) +
  labs(
    y = "Relative probability of occurrence",
    x = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(3, 3, 3, 3),
    axis.ticks.length = unit(2, "pt"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# PANEL B ------------------------------------------------------------

p2 <- ggplot(
  summary_df,
  aes(
    x = model_label,
    y = value,
    shape = metric
  )
) +
  geom_point(
    aes(colour = model_label),
    size = 3.5,
    position = position_dodge(width = 0.1)
  ) +
  
  scale_colour_manual(values = c(
    "Highest-performing model" = "#8C6BB1",
    "Lowest-performing model" = "#4D4D4D"
  )) +
  
  scale_shape_manual(values = c(
    "FPPS" = 17,
    "FNR" = 15
  )) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    y = "Proportion",
    x = NULL,
    shape = "Metric"
  ) +
  
  guides(
    colour = "none",   
    shape = guide_legend(order = 1)
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(),
    legend.position = c(0.15, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      colour = NA),
    plot.margin = margin(3, 3, 3, 3),
    axis.ticks.length = unit(2, "pt"),
    plot.background = element_rect(fill = "white", colour = NA)
  )


#make map of individuals
eelgrass_indep <- rast(c("code/output_data/independent_validation/BCeelgrass_netforce_2013_2023.tif"))
eel_bin_xgb_nep <-rast("raster/eelgrass/eelgrass_predictions_xgb_nep_binary_notmasked.tif")
eel_bin_nepspatial<-rast("raster/eelgrass/eelgrass_predictions_nepspatial_binary_notmasked.tif")
values(eelgrass_indep)[values(eelgrass_indep) >= 1] <- 1

eelgrass_indep <- resample(eelgrass_indep, eel_bin_xgb_nep, method = "near")

# make sure rasters align
stopifnot(compareGeom(eelgrass_indep, eel_bin_xgb_nep, eel_bin_nepspatial))

# create output raster
out <- ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 1 & eel_bin_nepspatial == 1, 4,
            ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 1 & eel_bin_nepspatial == 0, 3,
                 ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 0 & eel_bin_nepspatial == 1, 2,
                      ifel(eelgrass_indep == 1 & eel_bin_xgb_nep == 0 & eel_bin_nepspatial == 0, 1, NA))))

# write result
writeRaster(out, "raster/eelgrass/eelgrass_good&badmodel_overlayindependent.tif", overwrite = TRUE)
# looked in arc for best example areas to showcase and clipped it


eelgrass_indep_val_clipped <- rast(c("raster/independent_figue.tif"))

eelgrass_cat <- as.factor(eelgrass_indep_val_clipped)
levels(eelgrass_cat) <- data.frame(
  value = c(1, 2, 3, 4),
  category = c(
    "Observed only",
    "Observed + lowest-performing model",
    "Observed + highest-performing model",
    "Observed + both models"
  )
)

r_ext <- terra::ext(eelgrass_cat)

p3 <-ggplot() +
  geom_spatraster(data = eelgrass_cat) +
  geom_sf(data = coastline, fill = "grey85", color = "black", linewidth = 0.2) +
  scale_fill_manual(
    name = "Agreement category",
    values = c(
      "Observed only" = "#D73027",
      "Observed + lowest-performing model" = "#FC8D59",
      "Observed + highest-performing model" = "#91BFDB",
      "Observed + both models" = "#1A9850"
    ),
    na.value = "white",
    na.translate = FALSE,
    drop = FALSE
  ) +
  coord_sf(
    xlim = c(r_ext$xmin, r_ext$xmax),
    ylim = c(r_ext$ymin, r_ext$ymax),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.8),
      color = NA
    ),
    legend.key = element_rect(fill = NA, color = NA),
    plot.margin = margin(3, 3, 3, 3)
  )


# COMBINE ------------------------------------------------------------
# p1 <- p1 + theme(plot.margin = margin(5,5,5,5))
# p2 <- p2 + theme(plot.margin = margin(5,5,5,5))
# p3 <- p3 + theme(plot.margin = margin(5,5,5,5))


indep_plot <- (p1 / p2) | p3

ggsave(
  filename = "figures/independent_validation.png",
  plot = indep_plot,
  width = 13,
  height = 7,
  dpi = 300,
  bg = "white"
)





# field validation figures
load("code/output_data/field_validation/validation_dataset.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/combined_metrics_eelgrass.RData")
load("~/PhD/Chapters/Chapter-1-BC-Seagrass-Distribution/bc-seagrass-sdm/code/output_data/model_results/combined_metrics_surfgrass.RData")
threshold_df_eelgrass <- combined_metrics_eelgrass %>%
  select(model, threshold_spatial) %>%
  mutate(
    model = recode(
      model,
      "bccm_nospatial" = "GLMM_bccm",
      "bccm_spatial"   = "GLMM_spatial_bccm",
      "nep_nospatial"  = "GLMM_nep",
      "nep_spatial"    = "GLMM_spatial_nep",
      "XGBOOST_bccm"   = "XGB_bccm",
      "XGBOOST_nep"    = "XGB_nep",
      "GBM_bccm"       = "GBM_bccm",
      "GBM_nep"        = "GBM_nep"  ),
    model = factor(
        model,
        levels = c(
          "GLMM_bccm",
          "GLMM_spatial_bccm",
          "GBM_bccm",
          "XGB_bccm",
          "GLMM_nep",
          "GLMM_spatial_nep",
          "GBM_nep",
          "XGB_nep")))

threshold_df_surfgrass <- combined_metrics_surfgrass %>%
  select(model, threshold_spatial) %>%
  mutate(
    model = recode(
      model,
      "bccm_nospatial" = "GLMM_bccm",
      "bccm_spatial"   = "GLMM_spatial_bccm",
      "nep_nospatial"  = "GLMM_nep",
      "nep_spatial"    = "GLMM_spatial_nep",
      "XGBOOST_bccm"   = "XGB_bccm",
      "XGBOOST_nep"    = "XGB_nep",
      "GBM_bccm"       = "GBM_bccm",
      "GBM_nep"        = "GBM_nep"  ),
    model = factor(
      model,
      levels = c(
        "GLMM_bccm",
        "GLMM_spatial_bccm",
        "GBM_bccm",
        "XGB_bccm",
        "GLMM_nep",
        "GLMM_spatial_nep",
        "GBM_nep",
        "XGB_nep")))

validation_sf <- validation_sf %>%
  relocate(PH_freq, PH, ZM_freq, ZM, .after = PC_PH)

# need to also drop sites that after reviewing the database should not have been included as they were aborted or vis was too bad to see species
validation_sf <- validation_sf %>% sf::st_drop_geometry() %>%
  filter(!HKey %in% c(74, 84, 85, 161))  

# determine difference between observed and modelled environmental layers
# if either observed or modelled layers were mixed these were always retained as we are just looking at predominant substrate and most sites are actually mixed in reality
validation_sf <- validation_sf %>%
  mutate(Depth_diff = abs(bathy_mod - avgCorDepth_obs),
         Slope_diff = abs(slope_mod - slope_obs),
         Substrate_diff =
           if_else(
             substrate_mod == "Mixed" | substrate_mode_obs == "Mixed", TRUE,
             substrate_mod == substrate_mode_obs))

validation_sf <- validation_sf %>% select (HKey, Region, Survey, Visibility, High.vegetation, avgCorDepth_obs, bathy_mod,
                                           Depth_diff, slope_obs, slope_mod, Slope_diff, substrate_mode_obs, 
                                           substrate_mod, Substrate_diff, 
                                           ZM, ZM_freq, PC_ZM, zo_bccm_nospatial_pred, zo_bccm_spatial_pred, zo_nep_nospatial_pred, zo_nep_spatial_pred, zo_XGBOOST_bccm_pred, zo_XGBOOST_nep_pred, zo_GBM_bccm_pred, zo_GBM_nep_pred,
                                           PH, PH_freq, PC_PH, ph_bccm_nospatial_pred, ph_bccm_spatial_pred, ph_nep_nospatial_pred, ph_nep_spatial_pred, ph_XGBOOST_bccm_pred, ph_XGBOOST_nep_pred, ph_GBM_bccm_pred, ph_GBM_nep_pred) # need to add in surfgrass predictions
validation_sf <- validation_sf %>%
  mutate(
    Tidal_zone = if_else(Survey == "Intertidal", "Intertidal", "Subtidal"),
    Tidal_zone = factor(Tidal_zone, levels = c("Subtidal", "Intertidal"))
  )

# Prepare eelgrass validation data
eelgrass_pa_long <- validation_sf %>%
  mutate(
    Presence = factor(ifelse(ZM == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    zo_bccm_nospatial_pred,
    zo_bccm_spatial_pred,
    zo_nep_nospatial_pred,
    zo_nep_spatial_pred,
    zo_XGBOOST_bccm_pred,
    zo_XGBOOST_nep_pred,
    zo_GBM_bccm_pred,
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
      "zo_bccm_nospatial_pred" = "GLMM_bccm",
      "zo_bccm_spatial_pred"   = "GLMM_spatial_bccm",
      "zo_nep_nospatial_pred"  = "GLMM_nep",
      "zo_nep_spatial_pred"    = "GLMM_spatial_nep",
      "zo_XGBOOST_bccm_pred"   = "XGB_bccm",
      "zo_XGBOOST_nep_pred"    = "XGB_nep",
      "zo_GBM_bccm_pred"       = "GBM_bccm",
      "zo_GBM_nep_pred"        = "GBM_nep"  ),
    model = factor(
      model,
      levels = c(
        "GLMM_bccm",
        "GLMM_spatial_bccm",
        "GBM_bccm",
        "XGB_bccm",
        "GLMM_nep",
        "GLMM_spatial_nep",
        "GBM_nep",
        "XGB_nep"
      )
    )
  )

surfgrass_pa_long <- validation_sf %>%
  mutate(
    Presence = factor(ifelse(PH == 0, "Absent", "Present"),
                      levels = c("Absent", "Present"))
  ) %>%
  select(
    Presence,
    Tidal_zone,
    ph_bccm_nospatial_pred,
    ph_bccm_spatial_pred,
    ph_nep_nospatial_pred,
    ph_nep_spatial_pred,
    ph_XGBOOST_bccm_pred,
    ph_XGBOOST_nep_pred,
    ph_GBM_bccm_pred,
    ph_GBM_nep_pred
  ) %>%
  sf::st_drop_geometry() %>%
  pivot_longer(
    cols = -c(Presence, Tidal_zone),
    names_to = "model",
    values_to = "predicted_suitability"
  ) %>%
  mutate(
    model = recode(
      model,
      "ph_bccm_nospatial_pred" = "GLMM_bccm",
      "ph_bccm_spatial_pred"   = "GLMM_spatial_bccm",
      "ph_nep_nospatial_pred"  = "GLMM_nep",
      "ph_nep_spatial_pred"    = "GLMM_spatial_nep",
      "ph_XGBOOST_bccm_pred"   = "XGB_bccm",
      "ph_XGBOOST_nep_pred"    = "XGB_nep",
      "ph_GBM_bccm_pred"       = "GBM_bccm",
      "ph_GBM_nep_pred"        = "GBM_nep"  ),
    model = factor(
      model,
      levels = c(
        "GLMM_bccm",
        "GLMM_spatial_bccm",
        "GBM_bccm",
        "XGB_bccm",
        "GLMM_nep",
        "GLMM_spatial_nep",
        "GBM_nep",
        "XGB_nep"
      )
    )
  )


# Plot
eelgrass_field <- ggplot(eelgrass_pa_long,
       aes(x = Presence,
           y = predicted_suitability,
           fill = Tidal_zone)) +
  
  geom_boxplot(
    aes(group = interaction(Presence, Tidal_zone)),
    position = position_dodge(width = 0.75),
    outlier.alpha = 0.3,
    width = 0.6
  ) +
  
  
  facet_wrap(~ model, ncol = 4) +
  
  geom_hline(
    data = threshold_df_eelgrass,
    aes(yintercept = threshold_spatial),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Tidal zone",
    values = c(
      "Subtidal" = "#0B5D5E",   # dark teal
      "Intertidal" = "#6EC6C4"  # light teal
    )
  ) +
  
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  #coord_cartesian(ylim = c(0, 1))
  
  labs(
    title = "Eelgrass",
    x = NULL,
    y = "Relative probability of occurrence"
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),   # remove grid lines
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0, face = "bold"),
    panel.spacing = unit(0, "lines")  
  )


surfgrass_field <- ggplot(surfgrass_pa_long,
                          aes(x = Presence,
                              y = predicted_suitability,
                              fill = Tidal_zone)) +
  
  geom_boxplot(
    aes(group = interaction(Presence, Tidal_zone)),
    position = position_dodge(width = 0.75),
    outlier.alpha = 0.3,
    width = 0.6
  ) +
  
  facet_wrap(~ model, ncol = 4) +
  
  geom_hline(
    data = threshold_df_surfgrass,
    aes(yintercept = threshold_spatial),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    name = "Tidal zone",
    values = c(
      "Subtidal" = "#0B5D5E",   # dark teal
      "Intertidal" = "#6EC6C4"  # light teal
    )
  ) +
  
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  #coord_cartesian(ylim = c(0, 1))
  
  labs(
    title = "Surfgrass",
    x = NULL,
    y = NULL
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
    panel.spacing = unit(0, "lines")  
  )

 field_val<-eelgrass_field + surfgrass_field

 ggsave(
   filename = "figures/field_validation.png",
   plot = field_val,
   width = 13,
   height = 12,
   dpi = 300,
   bg = "white"
 )
 