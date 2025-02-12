###############################################################################
#
# Authors:      Ashley Park (with some code modified from Matt Csordas and Graham Epstien)
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
#
#
# Objective:
# ---------
# prepare oceanography data from BCCM, NEP36 and Salish Sea Cast at 20 m resolution 
#
#STILL TO DO, FIGURE OUT HOW TO APPLY SSH. dOES THIS NEED TO BE APPLIED TO HINDCAST OR JUST PROJECTIONS?? Or do you just subtract the hindcast from the projection and that is how much you change the depth by?
###############################################################################

#load packages####
library(sf)
library(tidyverse)
library(terra)
library(tidync)

#read in 20m bathymetry layer####
bathy_hg <- rast("raw_data/envlayers-20m-hg//bathymetry.tif")
bathy_ncc <- rast("raw_data/envlayers-20m-ncc/bathymetry.tif")
bathy_qcs <- rast("raw_data/envlayers-20m-qcs/bathymetry.tif")
bathy_wcvi <- rast("raw_data/envlayers-20m-wcvi/bathymetry.tif")
bathy_ss <- rast("raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif")

bathy20m <- mosaic(bathy_hg, bathy_ncc, bathy_qcs, bathy_wcvi, bathy_ss)



#read in BCCM data

#NH4
#read in the nc file
BCCM_h_NH4_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_mon1993to2021_NH4.nc")
#Load in the data and grid
BCCM_h_NH4_data <- tidync::hyper_tibble(BCCM_h_NH4_nc)
BCCM_h_NH4_data$months <- as.numeric(BCCM_h_NH4_data$months)
BCCM_h_NH4_data <- BCCM_h_NH4_data %>%
  mutate(date = as.Date('1992-12-01') + months(months),
         year = year(date),
         month = month(date))

BCCM_h_NH4_grid <- BCCM_h_NH4_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

#NO3
#read in the nc file
BCCM_h_NO3_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_mon1993to2021_NO3.nc")
#Load in the data and grid
BCCM_h_NO3_data <- tidync::hyper_tibble(BCCM_h_NO3_nc)
BCCM_h_NO3_data$months <- as.numeric(BCCM_h_NO3_data$months)
BCCM_h_NO3_data <- BCCM_h_NO3_data %>%
  mutate(date = as.Date('1992-12-01') + months(months),
         year = year(date),
         month = month(date))

BCCM_h_NO3_grid <- BCCM_h_NO3_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

#Salinity
#read in the nc file
BCCM_h_salt_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_mon1993to2021_salt.nc")
#Load in the data and grid
BCCM_h_salt_data <- tidync::hyper_tibble(BCCM_h_salt_nc)
BCCM_h_salt_data$months <- as.numeric(BCCM_h_salt_data$months)
BCCM_h_salt_data <- BCCM_h_salt_data %>%
  mutate(date = as.Date('1992-12-01') + months(months),
         year = year(date),
         month = month(date))
BCCM_h_salt_grid <- BCCM_h_salt_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

#PAR
#read in the nc file
BCCM_h_urad_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_1993to2019_Urad.nc")
#Load in the data and grid
BCCM_h_urad_data <- tidync::hyper_tibble(BCCM_h_urad_nc)
BCCM_h_urad_data$years <- as.numeric(BCCM_h_urad_data$years)
BCCM_h_urad_data <- BCCM_h_urad_data %>%
  mutate(date = as.Date('1992-01-01') + years(years),
         year = year(date))
BCCM_h_urad_grid <- BCCM_h_urad_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

#DO
#read in the nc file
BCCM_h_do_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_mon1993to2021_Oxy.nc")
#Load in the data and grid
BCCM_h_do_data <- tidync::hyper_tibble(BCCM_h_do_nc)
BCCM_h_do_data$months <- as.numeric(BCCM_h_do_data$months)
BCCM_h_do_data <- BCCM_h_do_data %>%
  mutate(date = as.Date('1992-12-01') + months(months),
         year = year(date),
         month = month(date))
BCCM_h_do_grid <- BCCM_h_do_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

#Temperature
#read in the nc file
BCCM_h_temp_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_mon1993to2021_temp.nc")
#Load in the data and grid
BCCM_h_temp_data <- tidync::hyper_tibble(BCCM_h_temp_nc)
BCCM_h_temp_data$months <- as.numeric(BCCM_h_temp_data$months)
BCCM_h_temp_data <- BCCM_h_temp_data %>%
  mutate(date = as.Date('1992-12-01') + months(months),
         year = year(date),
         month = month(date))
BCCM_h_temp_grid <- BCCM_h_temp_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

gc()


## Create climatologies of hindcast variables for prediction to occur onto for present day
#Set the start and end year for the climatology to be calculated based off of
# for seagrass climatologies were created in decadal slices as wanted to match data based on climatologies, not annual data as eelgrass is a perrenial species
#that spreads by colonal growth and beds are usually fairly static even if they expand and contract
# 1993-2002, 2003-2012, 2013-2021. Also did 1993-2021 for prediction grid
start <- 1993
end <- 2002

########
#NH4
########
#subset the NH4 data by these start and end year
NH4_sub <- BCCM_h_NH4_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
NH4_summary_year <- NH4_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% 
  dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE))
#Summarise the NH4 data into the desired climatologies for the entire time period
NH4_summary <- NH4_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(NH4_5_mean = mean(mean_5_year, na.rm = TRUE))
#Combine the above dataframe with the NH4 grid and filter out erroneous entries
NH4_summary_grid <- dplyr::full_join(NH4_summary, BCCM_h_NH4_grid, by=c("xi_rho","eta_rho")) %>% na.omit() 
#Convert to spatial points and transform to grid crs
NH4_points <- sf::st_as_sf(NH4_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_NH4_5_mean_rast <- terra::rasterize(NH4_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(NH4_points), res=3000),field = "NH4_5_mean", fun=mean)

########
#NO3
########
#subset the NO3 data by these start and end year
NO3_sub <- BCCM_h_NO3_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
NO3_summary_year <- NO3_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% 
  dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE))
#Summarise the NO3 data into the desired climatologies for the entire time period
NO3_summary <- NO3_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(NO3_5_mean = mean(mean_5_year, na.rm = TRUE))
#Combine the above dataframe with the NO3 grid and filter out erroneous entries
NO3_summary_grid <- dplyr::full_join(NO3_summary, BCCM_h_NO3_grid, by=c("xi_rho","eta_rho")) %>% na.omit() 
#Convert to spatial points and transform to grid crs
NO3_points <- sf::st_as_sf(NO3_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_NO3_5_mean_rast <- terra::rasterize(NO3_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(NO3_points), res=3000),field = "NO3_5_mean", fun=mean)

########
#salt
########
#subset the salt data by these start and end year
salt_sub <- BCCM_h_salt_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
salt_summary_year <- salt_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE), min_5_year = min(min_5m, na.rm = TRUE))
#Summarise the salt data into the desired climatologies for the entire time period
salt_summary <- salt_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(salt_5_mean = mean(mean_5_year, na.rm = TRUE),  salt_5_min = mean(min_5_year, na.rm = TRUE))
#Combine the above dataframe with the salt grid and filter out erroneous entries
salt_summary_grid <- dplyr::full_join(salt_summary, BCCM_h_salt_grid, by=c("xi_rho","eta_rho")) %>% na.omit()
#Convert to spatial points and transform to grid crs
salt_points <- sf::st_as_sf(salt_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_salt_5_mean_rast <- terra::rasterize(salt_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(salt_points), res=3000),field = "salt_5_mean", fun=mean)
bccm_salt_5_min_rast <- terra::rasterize(salt_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(salt_points), res=3000),field = "salt_5_min", fun=mean)

########
#PAR
########
#subset the sshPAR data by these start and end year
urad_sub <- BCCM_h_urad_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
PAR_summary_year <- urad_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE), max_5_year = max(max_5m, na.rm = TRUE), min_5_year = min(min_5m, na.rm = TRUE))
#Summarise the sshPAR data into the desired climatologies for the entire time period
PAR_summary <- PAR_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(PAR_5_mean = mean(mean_5_year, na.rm = TRUE), PAR_5_max = mean(max_5_year, na.rm = TRUE), PAR_5_min = mean(min_5_year, na.rm = TRUE))
#Combine the above dataframe with the sshPAR grid and filter out erroneous entries
PAR_summary_grid <- dplyr::full_join(PAR_summary, BCCM_h_urad_grid, by=c("xi_rho","eta_rho")) %>% na.omit()
#Convert to spatial points and transform to grid crs
PAR_points <- sf::st_as_sf(PAR_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of interest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_PAR_5_mean_rast <- terra::rasterize(PAR_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(PAR_points), res=3000),field = "PAR_5_mean", fun=mean)
bccm_PAR_5_max_rast <- terra::rasterize(PAR_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(PAR_points), res=3000),field = "PAR_5_max", fun=mean)
bccm_PAR_5_min_rast <- terra::rasterize(PAR_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(PAR_points), res=3000),field = "PAR_5_min", fun=mean)

########
#temp
########
#subset the temp data by these start and end year
temp_sub <- BCCM_h_temp_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
temp_summary_year <- temp_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_sur_year = mean(mean_sur, na.rm = TRUE), max_sur_year = max(max_sur, na.rm = TRUE), min_sur_year = min(min_sur, na.rm = TRUE), mean_5_year = mean(mean_5m, na.rm = TRUE), max_5_year = max(max_5m, na.rm = TRUE), min_5_year = min(min_5m, na.rm = TRUE))
#Summarise the temp data into the desired climatologies for the entire time period
temp_summary <- temp_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(temp_s_mean = mean(mean_sur_year, na.rm = TRUE), temp_s_max = mean(max_sur_year, na.rm = TRUE), temp_s_min = mean(min_sur_year, na.rm = TRUE), temp_5_mean = mean(mean_5_year, na.rm = TRUE), temp_5_max = mean(max_5_year, na.rm = TRUE), temp_5_min = mean(min_5_year, na.rm = TRUE))
#Combine the above dataframe with the temp grid and filter out erroneous entries
temp_summary_grid <- dplyr::full_join(temp_summary, BCCM_h_temp_grid, by=c("xi_rho","eta_rho")) %>% na.omit() # %>% dplyr::filter(temp_s_min >= 0, temp_10_min >= 0)
#Convert to spatial points and transform to grid crs
temp_points <- sf::st_as_sf(temp_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of interest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_temp_s_mean_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_s_mean", fun=mean)
bccm_temp_s_max_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_s_max", fun=max)
bccm_temp_s_min_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_s_min", fun=min)
bccm_temp_5_mean_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_mean", fun=mean)
bccm_temp_5_max_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_max", fun=max)
bccm_temp_5_min_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_min", fun=min)

########
#DO
########
#subset the temp data by these start and end year
do_sub <- BCCM_h_do_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
do_summary_year <- do_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE), min_5_year = min(min_5m, na.rm = TRUE))
#Summarise the do data into the desired climatologies for the entire time period
do_summary <- do_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(do_5_mean = mean(mean_5_year, na.rm = TRUE), do_5_min = mean(min_5_year, na.rm = TRUE))
#Combine the above dataframe with the sshPAR grid and filter out erroneous entries
do_summary_grid <- dplyr::full_join(do_summary, BCCM_h_do_grid, by=c("xi_rho","eta_rho")) %>% na.omit() 
#Convert to spatial points and transform to grid crs
do_points <- sf::st_as_sf(do_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of interest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_do_5_mean_rast <- terra::rasterize(do_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(do_points), res=3000),field = "do_5_mean", fun=mean)
bccm_do_5_min_rast <- terra::rasterize(do_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(do_points), res=3000),field = "do_5_min", fun=mean)


#read in SalishSea cast data
#Annual MinMax is highest/lowest monthly mean within a given year; and Monthly MinMax is highest/lowest daily mean within given month

# dfo NH4, NO3, DO, PAR
#read in the nc file
SSC_dfo_yearmax_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmax_grid_bgc_T_1986-2005.nc")
SSC_dfo_yearmean_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmean_grid_bgc_T_1986-2005.nc")
SSC_dfo_yearmin_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmin_grid_bgc_T_1986-2005.nc")

# get time units
time_units_attr <- ncmeta::nc_atts("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmax_grid_bgc_T_1986-2005.nc", "time_counter") 
time_units <- time_units_attr$value[time_units_attr$name == "units"]
date_origin <- ymd(time_units$units)

#Load in the data and grid
SSC_dfo_yearmax_data <- tidync::hyper_tibble(SSC_dfo_yearmax_nc)%>%
  rename(PARmax = PAR) %>%
  select(-c(dissolved_inorganic_carbon, total_alkalinity, dissolved_oxygen, ammonium, nitrate))

SSC_dfo_yearmean_data <- tidync::hyper_tibble(SSC_dfo_yearmean_nc)%>%
  rename(PARmean = PAR,
         NH4mean = ammonium,
         DOmean = dissolved_oxygen,
         NO3mean = nitrate) %>%
  select(-c(dissolved_inorganic_carbon, total_alkalinity))

SSC_dfo_yearmin_data <- tidync::hyper_tibble(SSC_dfo_yearmin_nc)%>%
  rename(PARmin = PAR,
         DOmin = dissolved_oxygen) %>%
  select(-c(dissolved_inorganic_carbon, total_alkalinity, ammonium, nitrate))

SSC_dfo_year_data <- SSC_dfo_yearmax_data %>%
  full_join(SSC_dfo_yearmean_data, by = join_by(x, y, deptht, time_counter)) %>%
  full_join(SSC_dfo_yearmin_data, by = join_by(x, y, deptht, time_counter)) %>%
  mutate(year = year(as.POSIXct(time_counter, origin = date_origin, tz = "UTC"))) %>%
  filter(deptht > 0 & deptht < 6)

#grid 
SSC_dfo_year_grid <- SSC_dfo_yearmax_nc %>% tidync::activate("D2,D3") %>% tidync::hyper_tibble()

# ubc NH4, NO3, DO, PAR
SSC_ubc_NH4_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Ammonium_Surf_AnnualMean.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_NH4_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Ammonium_4pt5_AnnualMean.csv") %>% mutate(deptht = 4.5)
SSC_ubc_NH4_csv <-bind_rows(SSC_ubc_NH4_sur_csv, SSC_ubc_NH4_5_csv)
SSC_ubc_NO3_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Nitrate_Surf_AnnualMean.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_NO3_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Nitrate_4pt5_AnnualMean.csv") %>% mutate(deptht = 4.5)
SSC_ubc_NO3_csv <-bind_rows(SSC_ubc_NO3_sur_csv, SSC_ubc_NO3_5_csv)
SSC_ubc_DO_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/DO_Surf_AnnualMeanMinMax.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_DO_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/DO_4pt5_AnnualMeanMinMax.csv") %>% mutate(deptht = 4.5)
SSC_ubc_DO_csv <-bind_rows(SSC_ubc_DO_sur_csv, SSC_ubc_DO_5_csv)
SSC_ubc_PAR_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/PAR_Surf_AnnualMeanMinMax.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_PAR_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/PAR_4pt5_AnnualMeanMinMax.csv") %>% mutate(deptht = 4.5)
SSC_ubc_PAR_csv <-bind_rows(SSC_ubc_PAR_sur_csv, SSC_ubc_PAR_5_csv)

dfo_grid <- SSC_dfo_year_grid %>%  rename(longitude = nav_lon, latitude = nav_lat)
SSC_ubc_year_data <- SSC_ubc_NH4_csv %>% 
  left_join(dfo_grid, by = join_by(longitude, latitude)) %>% rename (NH4mean = Ammonium_mean)
#START FROM HERE

  full_join(SSC_ubc_NH4_5_csv, by = join_by(longitude, latitude, time)) %>%
  



# DFO temp and salt
#read in the nc file
SSC_dfo_monthmax_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_monmax_grid_T_1986-2005.nc")
SSC_dfo_monthmean_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_monmean_grid_T_1986-2005.nc")
SSC_dfo_monthmin_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_monmin_grid_T_1986-2005.nc")

# get time units
time_units_attr <- ncmeta::nc_atts("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_monmax_grid_T_1986-2005.nc", "time_counter") 
time_units <- time_units_attr$value[time_units_attr$name == "units"]
date_origin <- ymd(time_units$units)

#Load in the data and grid
SSC_dfo_monthmax_data <- tidync::hyper_tibble(SSC_dfo_monthmax_nc, force = TRUE)%>%
  rename(tempmax = votemper,
         saltmax = vosaline) 

SSC_dfo_monthmean_data <- tidync::hyper_tibble(SSC_dfo_monthmean_nc, force = TRUE)%>%
  rename(tempmean = votemper,
         saltmean = vosaline) 

SSC_dfo_monthmin_data <- tidync::hyper_tibble(SSC_dfo_monthmin_nc, force = TRUE)%>%
  rename(tempmin = votemper,
         saltmin = vosaline) 

SSC_dfo_month_data <- SSC_dfo_monthmax_data %>%
  full_join(SSC_dfo_monthmean_data, by = join_by(x, y, deptht, time_counter)) %>%
  full_join(SSC_dfo_monthmin_data, by = join_by(x, y, deptht, time_counter)) %>%
  mutate(year = year(as.POSIXct(time_counter, origin = date_origin, tz = "UTC"))) %>%
  filter(deptht > 0 & deptht < 6)

#grid 
SSC_dfo_month_grid <- SSC_dfo_monthmax_nc %>% tidync::activate("D2,D3") %>% tidync::hyper_tibble()

#ubc temp and salinity
SSC_ubc_salt_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Salinity_Surf_MonthlyMean.csv")
SSC_ubc_salt_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Salinity_4pt5_MonthlyMean.csv")
SSC_ubc_temp_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Temperature_All_MonthlyMeanMinMax.csv")



## Create climatologies of hindcast variables for prediction to occur onto for present day
#Set the start and end year for the climatology to be calculated based off of
# for seagrass climatologies were created in decadal slices as wanted to match data based on climatologies, not annual data as eelgrass is a perrenial species
#that spreads by colonal growth and beds are usually fairly static even if they expand and contract
# 1993-2002, 2003-2012, 2013-2021. Also did 1993-2021 for prediction grid
start <- 1993
end <- 2002

########
#Annual data
########
#subset the year data by these start and end year, and just surface values
Year_surf_sub <- SSC_dfo_year_data %>% 
  dplyr::filter(year >= start & year <= end,
                deptht < 2)
#Summarise the predictor into the desired climatologies for the entire time period
year_surf_summary <- Year_surf_sub %>% dplyr::group_by(x, y) %>% 
  dplyr::summarise(PAR_mean_sur = mean(PARmean, na.rm = TRUE), PAR_max_sur = mean(PARmax, na.rm = TRUE), PAR_min_sur = mean(PARmin, na.rm = TRUE),
                   NH4_mean_sur = mean(NH4mean, na.rm = TRUE), NO3_mean_sur = mean(NO3mean, na.rm = TRUE),
                   DO_mean_sur = mean(DOmean, na.rm = TRUE), DO_min_sur = mean(DOmin, na.rm = TRUE))

#Combine the above dataframe with the NH4 grid and filter out erroneous entries
year_surf_summary_grid <- dplyr::full_join(year_surf_summary, SSC_dfo_year_grid, by=c("x","y")) %>% na.omit() %>% dplyr::filter(DO_mean_sur > 0)
#Convert to spatial points and transform to grid crs
year_surf_points <- sf::st_as_sf(year_surf_summary_grid, coords=(c("nav_lon", "nav_lat")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 500m equal area grid (size of salishseacast pixels) with mean = akin to project, merge with surface values because loose very shallow areas
ssc_PAR_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "PAR_mean_sur", fun=mean)
ssc_PAR_sur_max_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "PAR_max_sur", fun=mean)
ssc_PAR_sur_min_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "PAR_min_sur", fun=mean)
ssc_NH4_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "NH4_mean_sur", fun=mean)
ssc_NO3_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "NO3_mean_sur", fun=mean)
ssc_DO_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "DO_mean_sur", fun=mean)
ssc_DO_sur_min_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "DO_min_sur", fun=mean)


#subset the year data by these start and end year, and just 5m values
Year_5_sub <- SSC_dfo_year_data %>% 
  dplyr::filter(year >= start & year <= end,
                deptht > 2)
#Summarise the predictor into the desired climatologies for the entire time period
year_5_summary <- Year_5_sub %>% dplyr::group_by(x, y) %>% 
  dplyr::summarise(PAR_mean_5 = mean(PARmean, na.rm = TRUE), PAR_max_5 = mean(PARmax, na.rm = TRUE), PAR_min_5 = mean(PARmin, na.rm = TRUE),
                   NH4_mean_5 = mean(NH4mean, na.rm = TRUE), NO3_mean_5 = mean(NO3mean, na.rm = TRUE), 
                   DO_mean_5 = mean(DOmean, na.rm = TRUE), DO_min_5 = mean(DOmin, na.rm = TRUE))

#Combine the above dataframe with the grid and filter out erroneous entries
year_5_summary_grid <- dplyr::full_join(year_5_summary, SSC_dfo_year_grid, by=c("x","y")) %>% na.omit() %>% dplyr::filter(DO_mean_5 > 0)
#Convert to spatial points and transform to grid crs
year_5_points <- sf::st_as_sf(year_5_summary_grid, coords=(c("nav_lon", "nav_lat")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 500m equal area grid (size of salishseacast pixels) with mean = akin to project, merge with surface values because loose very shallow areas
ssc_PAR_5_mean_rast <- terra::rasterize(year_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_5_points), res=500),field = "PAR_mean_5", fun=mean) %>%
  terra::merge(ssc_PAR_sur_mean_rast, na.rm=TRUE)
ssc_PAR_5_max_rast <- terra::rasterize(year_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_5_points), res=500),field = "PAR_max_5", fun=mean)%>%
  terra::merge(ssc_PAR_sur_max_rast, na.rm=TRUE)
ssc_PAR_5_min_rast <- terra::rasterize(year_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_5_points), res=500),field = "PAR_min_5", fun=mean)%>%
  terra::merge(ssc_PAR_sur_min_rast, na.rm=TRUE)
ssc_NH4_5_mean_rast <- terra::rasterize(year_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_5_points), res=500),field = "NH4_mean_5", fun=mean)%>%
  terra::merge(ssc_NH4_sur_mean_rast, na.rm=TRUE)
ssc_NO3_5_mean_rast <- terra::rasterize(year_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_5_points), res=500),field = "NO3_mean_5", fun=mean)%>%
  terra::merge(ssc_NO3_sur_mean_rast, na.rm=TRUE)
ssc_DO_5_mean_rast <- terra::rasterize(year_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_5_points), res=500),field = "DO_mean_5", fun=mean)%>%
  terra::merge(ssc_DO_sur_mean_rast, na.rm=TRUE)
ssc_DO_5_min_rast <- terra::rasterize(year_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_5_points), res=500),field = "DO_min_5", fun=mean)%>%
  terra::merge(ssc_DO_sur_min_rast, na.rm=TRUE)


########
#monthly data
########
#subset the year data by these start and end year, and just surface values
Month_surf_sub <- SSC_dfo_month_data %>% 
  dplyr::filter(year >= start & year <= end,
                deptht < 2)
#Summarise the predictor by year
month_surf_summary_year <- Month_surf_sub %>% dplyr::group_by(x, y, year) %>% 
  dplyr::summarise(salt_mean_sur_year = mean(saltmean, na.rm = TRUE), salt_min_sur_year = min(saltmin, na.rm = TRUE),
                   temp_mean_sur_year = mean(tempmean, na.rm = TRUE), temp_max_sur_year = max(tempmax, na.rm = TRUE), temp_min_sur_year = min(tempmin, na.rm = TRUE))

#Summarise the predictor into the desired climatologies for the entire time period
month_surf_summary <- month_surf_summary_year %>% dplyr::group_by(x, y) %>% 
  dplyr::summarise(salt_mean_sur = mean(salt_mean_sur_year, na.rm = TRUE), salt_min_sur = mean(salt_min_sur_year, na.rm = TRUE),
                   temp_mean_sur = mean(temp_mean_sur_year, na.rm = TRUE), temp_max_sur = mean(temp_max_sur_year, na.rm = TRUE), temp_min_sur = mean(temp_min_sur_year, na.rm = TRUE))

#Combine the above dataframe with the grid and filter out erroneous entries
month_surf_summary_grid <- dplyr::full_join(month_surf_summary, SSC_dfo_month_grid, by=c("x","y")) %>% na.omit() %>% dplyr::filter(salt_mean_sur > 0)
#Convert to spatial points and transform to grid crs
month_surf_points <- sf::st_as_sf(month_surf_summary_grid, coords=(c("nav_lon", "nav_lat")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 500m equal area grid (size of salishseacast pixels) with mean = akin to project, merge with surface values because loose very shallow areas
ssc_temp_sur_mean_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_mean_sur", fun=mean)
ssc_temp_sur_max_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_max_sur", fun=mean)
ssc_temp_sur_min_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_min_sur", fun=mean)
ssc_salt_sur_mean_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "salt_mean_sur", fun=mean)
ssc_salt_sur_min_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "salt_min_sur", fun=mean)


#subset the year data by these start and end year, and just 5m values
Month_5_sub <- SSC_dfo_month_data %>% 
  dplyr::filter(year >= start & year <= end,
                deptht > 2)

#Summarise the predictor by year
month_5_summary_year <- Month_5_sub %>% dplyr::group_by(x, y, year) %>% 
  dplyr::summarise(salt_mean_5_year = mean(saltmean, na.rm = TRUE), salt_min_5_year = min(saltmin, na.rm = TRUE),
                   temp_mean_5_year = mean(tempmean, na.rm = TRUE), temp_max_5_year = max(tempmax, na.rm = TRUE), temp_min_5_year = min(tempmin, na.rm = TRUE))

#Summarise the predictor into the desired climatologies for the entire time period
month_5_summary <- month_5_summary_year %>% dplyr::group_by(x, y) %>% 
  dplyr::summarise(salt_mean_5 = mean(salt_mean_5_year, na.rm = TRUE), salt_min_5 = mean(salt_min_5_year, na.rm = TRUE),
                   temp_mean_5 = mean(temp_mean_5_year, na.rm = TRUE), temp_max_5 = mean(temp_max_5_year, na.rm = TRUE), temp_min_5 = mean(temp_min_5_year, na.rm = TRUE))

#Combine the above dataframe with the NH4 grid and filter out erroneous entries
month_5_summary_grid <- dplyr::full_join(month_5_summary, SSC_dfo_month_grid, by=c("x","y")) %>% na.omit() %>% dplyr::filter(salt_mean_5 > 0)
#Convert to spatial points and transform to grid crs
month_5_points <- sf::st_as_sf(month_5_summary_grid, coords=(c("nav_lon", "nav_lat")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 500m equal area grid (size of salishseacast pixels) with mean = akin to project, merge with surface values because loose very shallow areas
ssc_temp_5_mean_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "temp_mean_5", fun=mean) %>%
  terra::merge(ssc_temp_sur_mean_rast, na.rm=TRUE)
ssc_temp_5_max_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "temp_max_5", fun=mean)%>%
  terra::merge(ssc_temp_sur_max_rast, na.rm=TRUE)
ssc_temp_5_min_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "temp_min_5", fun=mean)%>%
  terra::merge(ssc_temp_sur_min_rast, na.rm=TRUE)
ssc_salt_5_mean_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "salt_mean_5", fun=mean)%>%
  terra::merge(ssc_salt_sur_mean_rast, na.rm=TRUE)
ssc_salt_5_min_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "salt_min_5", fun=mean)%>%
  terra::merge(ssc_salt_sur_min_rast, na.rm=TRUE)

#Stack all of the layers together
BCCM_layers <- c(bccm_NH4_5_mean_rast, bccm_NO3_5_mean_rast, bccm_salt_5_mean_rast, bccm_salt_5_min_rast, bccm_PAR_5_mean_rast, bccm_PAR_5_min_rast, bccm_PAR_5_max_rast, bccm_temp_s_mean_rast, bccm_temp_s_max_rast, bccm_temp_s_min_rast, bccm_temp_5_mean_rast, bccm_temp_5_max_rast, bccm_temp_5_min_rast, bccm_do_5_mean_rast, bccm_do_5_min_rast)
SSC_layers <- c(ssc_NH4_5_mean_rast, ssc_NO3_5_mean_rast, ssc_salt_5_mean_rast, ssc_salt_5_min_rast, ssc_PAR_5_mean_rast, ssc_PAR_5_min_rast, ssc_PAR_5_max_rast, ssc_temp_s_mean_rast, ssc_temp_s_max_rast, ssc_temp_s_min_rast, ssc_temp_5_mean_rast, ssc_temp_5_max_rast, ssc_temp_5_min_rast, ssc_do_5_mean_rast, ssc_do_5_min_rast)
 c(ssc_salt_5_mean_rast) 
  
BCCM_Resampled <- terra::sapp(BCCM_layers, function(x){
  terra::focal(x,w=3, fun = "mean", na.policy = "only", na.rm = TRUE) %>% terra::project(bathy20m, method="cubicspline") 
})

SSC_Resampled <- terra::sapp(SSC_layers, function(x){
  terra::focal(x,w=7, fun = "mean", na.policy = "only", na.rm = TRUE) %>% terra::project(bathy20m, method="cubicspline") 
})

salt_5_mean_rast<- terra::merge(SSC_Resampled, BCCM_Resampled, na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)

terra::writeRaster(salt_5_mean_rast, "code/output_data/processed_ocean_variables/test2.tif", overwrite = TRUE)










BCCM_Resampled <- terra::sapp(BCCM_layers, function(x){
  terra::focal(x,w=3, fun = "mean", na.policy = "only", na.rm = TRUE) %>% terra::project(bathy20m, method="cubicspline") %>% 
    terra::mask(bathy20m) %>% terra::crop(bathy20m)
})

#Combine BCCM climatologies with the static predictor layers created previously to create  full hindcast predictor layer stack
Predictor_Hindcast_Climatologies <- BCCM_Resampled

names(Predictor_Hindcast_Climatologies) <- c("NH4_5m_mean", "NH4_5m_max", "NH4_5m_min", "NO3_5m_mean", "NO3_5m_max", "NO3_5m_min","salt_5m_mean", "salt_5m_min", "PAR_5m_mean", "PAR_5m_min", "PAR_5m_max", "temp_s_mean", "temp_s_max", "temp_s_min", "temp_5m_mean", "temp_5m_max", "temp_5m_min", "do_5m_mean", "do_5m_min")

#Save the hindcast predictor raster
terra::writeRaster(Predictor_Hindcast_Climatologies, "code/output_data/processed_ocean_variables/BCCM/Predictor_Hindcast_Climatologies_1993-2021.tif", overwrite = TRUE)

