###############################################################################
#
# Authors:      Ashley Park 
# Affiliation:  Fisheries and Oceans Canada (DFO) and University of British Columbia
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca 
# Project:      BC Seagrass SDM
#
#
# Objective:
# ---------
# prepare oceanography data from BCCM, NEP36 and Salish Sea Cast at 20 m resolution 
#
# STILL TO DO
# FIGURE OUT HOW TO APPLY SSH. dOES THIS NEED TO BE APPLIED TO HINDCAST OR JUST PROJECTIONS?? Or do you just subtract the hindcast from the projection and that is how much you change the depth by?
# ADD NEP36 data
###############################################################################

#load packages####
library(sf)
library(tidyverse)
library(terra)
library(tidync)
library(ncmeta)

#read in 20m bathymetry layer####
bathy_hg <- rast("raw_data/envlayers-20m-hg//bathymetry.tif")
bathy_ncc <- rast("raw_data/envlayers-20m-ncc/bathymetry.tif")
bathy_qcs <- rast("raw_data/envlayers-20m-qcs/bathymetry.tif")
bathy_wcvi <- rast("raw_data/envlayers-20m-wcvi/bathymetry.tif")
bathy_ss <- rast("raw_data/envlayers-20m-shelfsalishsea/bathymetry.tif")

bathy20m <- mosaic(bathy_hg, bathy_ncc, bathy_qcs, bathy_wcvi, bathy_ss)


#read in BCCM data 
# the min and max for monthly values were calculated from the 3-day average values

#### NH4 ####
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

#### NO3 ####
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

#### Salinity #### 
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


#### PAR ####
#read in the nc file
BCCM_h_urad_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_1993to2019_Urad.nc")
#Load in the data and grid
BCCM_h_urad_data <- tidync::hyper_tibble(BCCM_h_urad_nc)
BCCM_h_urad_data$years <- as.numeric(BCCM_h_urad_data$years)
BCCM_h_urad_data <- BCCM_h_urad_data %>%
  mutate(date = as.Date('1992-01-01') + years(years),
         year = year(date))
BCCM_h_urad_grid <- BCCM_h_urad_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

#### DO ####
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

#### Temperature ####
#read in the nc file
BCCM_h_temp_nc <- tidync::tidync("raw_data/BCCM/hindcast/bcc42_run4_mon1993to2021_temp.nc")
#Load in the data and grid
BCCM_h_temp_data <- tidync::hyper_tibble(BCCM_h_temp_nc)
BCCM_h_temp_data$months <- as.numeric(BCCM_h_temp_data$months)
BCCM_h_temp_data <- BCCM_h_temp_data %>%
  mutate(date = as.Date('1992-12-01') + months(months),
         year = year(date),
         month = month(date), 
         diff_sur = max_sur - min_sur,
         diff_5m = max_5m - min_5m)
BCCM_h_temp_grid <- BCCM_h_temp_nc %>% tidync::activate("D0,D1") %>% tidync::hyper_tibble()

####read in SalishSea cast data####
#Annual MinMax is highest/lowest monthly mean within a given year; and Monthly MinMax is highest/lowest daily mean within given month

#### Yearly data for NH4, NO3, DO, PAR####
#### dfo NH4, NO3, DO, PAR ####
#read in the nc file
SSC_dfo_yearmax_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmax_grid_bgc_T_1986-2005.nc")
SSC_dfo_yearmean_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmean_grid_bgc_T_1986-2005.nc")
SSC_dfo_yearmin_nc <- tidync::tidync("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmin_grid_bgc_T_1986-2005.nc")

# get time units
time_units_attr <- ncmeta::nc_atts("raw_data/SalishSeaCast/Historical1986-2005/SalishSeaCast-VNR023_yearmax_grid_bgc_T_1986-2005.nc", "time_counter") 
time_units <- time_units_attr$value[time_units_attr$name == "units"]
date_origin <- ymd(time_units$units)

#Load in the data and grid
# Consulting with Amber Holdsworth for the oxygen, a delta needs to be applied to lower the high oxygen values in the model.  The delta was computed from the median bias between the obs and model profiles = subtracted  35.88 mmol/m3 in all simulations.
SSC_dfo_yearmax_data <- tidync::hyper_tibble(SSC_dfo_yearmax_nc)%>%
  rename(PARmax = PAR) %>%
  select(-c(dissolved_inorganic_carbon, total_alkalinity, dissolved_oxygen, ammonium, nitrate))

SSC_dfo_yearmean_data <- tidync::hyper_tibble(SSC_dfo_yearmean_nc)%>%
  mutate(DOmean = dissolved_oxygen - 35.88) %>%
  rename(PARmean = PAR, NH4mean = ammonium, NO3mean = nitrate) %>%
  select(-c(dissolved_inorganic_carbon, total_alkalinity, dissolved_oxygen))

SSC_dfo_yearmin_data <- tidync::hyper_tibble(SSC_dfo_yearmin_nc)%>%
  mutate(DOmin = dissolved_oxygen - 35.88) %>%
  rename(PARmin = PAR) %>%
  select(-c(dissolved_inorganic_carbon, total_alkalinity, ammonium, nitrate, dissolved_oxygen))

SSC_dfo_year_data <- SSC_dfo_yearmax_data %>%
  full_join(SSC_dfo_yearmean_data, by = join_by(x, y, deptht, time_counter)) %>%
  full_join(SSC_dfo_yearmin_data, by = join_by(x, y, deptht, time_counter)) %>%
  mutate(year = year(as.POSIXct(time_counter, origin = date_origin, tz = "UTC"))) %>%
  filter(deptht > 0 & deptht < 6)

#grid 
SSC_dfo_year_grid <- SSC_dfo_yearmax_nc %>% tidync::activate("D2,D3") %>% tidync::hyper_tibble()

#### ubc NH4, NO3, DO, PAR ####
# read in csvs
SSC_ubc_NH4_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Ammonium_Surf_AnnualMean.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_NH4_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Ammonium_4pt5_AnnualMean.csv") %>% mutate(deptht = 4.5)
SSC_ubc_NH4_csv <-bind_rows(SSC_ubc_NH4_sur_csv, SSC_ubc_NH4_5_csv) %>% rename (NH4mean = Ammonium_mean)
SSC_ubc_NO3_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Nitrate_Surf_AnnualMean.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_NO3_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Nitrate_4pt5_AnnualMean.csv") %>% mutate(deptht = 4.5)
SSC_ubc_NO3_csv <-bind_rows(SSC_ubc_NO3_sur_csv, SSC_ubc_NO3_5_csv) %>% rename (NO3mean = Nitrate_mean)
SSC_ubc_DO_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/DO_Surf_AnnualMeanMinMax.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_DO_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/DO_4pt5_AnnualMeanMinMax.csv") %>% mutate(deptht = 4.5)
SSC_ubc_DO_csv <-bind_rows(SSC_ubc_DO_sur_csv, SSC_ubc_DO_5_csv)  %>% rename (DOmean = DO_mean, DOmin = DO_min) %>% select(-c(DO_max))
SSC_ubc_PAR_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/PAR_Surf_AnnualMeanMinMax.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_PAR_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/PAR_4pt5_AnnualMeanMinMax.csv") %>% mutate(deptht = 4.5)
SSC_ubc_PAR_csv <-bind_rows(SSC_ubc_PAR_sur_csv, SSC_ubc_PAR_5_csv) %>% rename (PARmean = PAR_mean, PARmin = PAR_min, PARmax = PAR_max)

# get x and y from dfo grid to be able to join ubc and dfo data
dfo_grid <- SSC_dfo_year_grid %>%  rename(longitude = nav_lon, latitude = nav_lat)
SSC_ubc_year_data <- SSC_ubc_NH4_csv %>% 
  left_join(dfo_grid, by = join_by(longitude, latitude)) %>%
  left_join(SSC_ubc_NO3_csv, by = join_by(longitude, latitude, time, deptht)) %>% 
  left_join(SSC_ubc_DO_csv, by = join_by(longitude, latitude, time, deptht)) %>% 
  left_join(SSC_ubc_PAR_csv, by = join_by(longitude, latitude, time, deptht)) %>%
  mutate(year = year(time), time_counter = NA) %>%
  select(-c(time, latitude, longitude))

#bind ubc and dfo data
SSC_year_data <-bind_rows(SSC_ubc_year_data, SSC_dfo_year_data) 

#### SalishSeaCast Monthly data for temperature and salinity####
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
  rename(tempmax = votemper, saltmax = vosaline) 

SSC_dfo_monthmean_data <- tidync::hyper_tibble(SSC_dfo_monthmean_nc, force = TRUE)%>%
  rename(tempmean = votemper, saltmean = vosaline) 

SSC_dfo_monthmin_data <- tidync::hyper_tibble(SSC_dfo_monthmin_nc, force = TRUE)%>%
  rename(tempmin = votemper, saltmin = vosaline) 

SSC_dfo_month_data <- SSC_dfo_monthmax_data %>%
  full_join(SSC_dfo_monthmean_data, by = join_by(x, y, deptht, time_counter)) %>%
  full_join(SSC_dfo_monthmin_data, by = join_by(x, y, deptht, time_counter)) %>%
  mutate(year = year(as.POSIXct(time_counter, origin = date_origin, tz = "UTC"))) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(saltmax, saltmin))

#grid 
SSC_dfo_month_grid <- SSC_dfo_monthmean_nc %>% tidync::activate("D2,D3") %>% tidync::hyper_tibble()

#ubc temp and salinity
SSC_ubc_salt_sur_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Salinity_Surf_MonthlyMean.csv") %>%  mutate(deptht = 0.5)
SSC_ubc_salt_5_csv <- read.csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Salinity_4pt5_MonthlyMean.csv") %>%  mutate(deptht = 4.5)
SSC_ubc_salt_csv <-bind_rows(SSC_ubc_salt_sur_csv, SSC_ubc_salt_5_csv)  %>% rename (saltmean = salinity) %>% mutate (time =as.Date(time)) 
SSC_ubc_salt_csv$latitude <- round(SSC_ubc_salt_csv$latitude, digits = 5) 
SSC_ubc_salt_csv$longitude<- round(SSC_ubc_salt_csv$longitude, digits = 4) 
SSC_ubc_temp_csv <- readr::read_csv("raw_data/SalishSeaCast/2007-2023 Current (V21-11)/Temperature_All_MonthlyMeanMinMax.csv") %>%  rename(deptht = depth, tempmean=temperature_mean, tempmax= temperature_max, tempmin = temperature_min) %>% filter(deptht < 6)
SSC_ubc_temp_csv$deptht[SSC_ubc_temp_csv$deptht > 1] <- 4.5
SSC_ubc_temp_csv$deptht[SSC_ubc_temp_csv$deptht < 1] <- 0.5
SSC_ubc_temp_csv$latitude <- round(SSC_ubc_temp_csv$latitude, digits = 5) 
SSC_ubc_temp_csv$longitude<- round(SSC_ubc_temp_csv$longitude, digits = 4) 

# get x and y from dfo grid to be able to join ubc and dfo data
dfo_mgrid <- SSC_dfo_month_grid %>%  rename(longitude = nav_lon, latitude = nav_lat)
dfo_mgrid$latitude <- round(dfo_mgrid$latitude, digits = 5) 
dfo_mgrid$longitude<- round(dfo_mgrid$longitude, digits = 4) 

SSC_ubc_month_data <- SSC_ubc_salt_csv %>% 
  left_join(dfo_mgrid, by = join_by(longitude, latitude)) %>% 
  left_join(SSC_ubc_temp_csv, by = join_by(longitude, latitude, time, deptht)) %>% 
  mutate(year = year(time), time_counter = NA) %>%
  select(-c(time, latitude, longitude))

#bind ubc and dfo data
SSC_month_data <-bind_rows(SSC_ubc_month_data, SSC_dfo_month_data) %>%
  mutate(tempdiff = tempmax - tempmin)



## Create climatologies of hindcast variables for prediction to occur onto for present day
# Set the start and end year for the climatology to be calculated based off of
# For seagrass climatologies were created in decadal slices as wanted to match data based on climatologies, not annual data as eelgrass is a perennial species that spreads by colonal growth;
# Beds are usually stable over time even if they expand and contract 
# 1993-2002, 2003-2012, 2013-2023. Also did 1993-2023 for prediction grid
start <- 2003
end <- 2012


#Summarise and make rasters
#### ROMS data####

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

#subset the salt data by these start and end year
salt_sub <- BCCM_h_salt_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year, want freshest month on average, and also monthly sd for each year 
# cv is calculated as the sd for a year divided by the mean of the year *100. Then the mean across the decade is calculated from each year
salt_summary_year <- salt_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE), min_5_year = min(mean_5m, na.rm = TRUE), sd_5_year = sd(mean_5m, na.rm = TRUE)) %>%
  mutate(cv_5_year = sd_5_year/mean_5_year*100)
#Summarise the salt data into the desired climatologies for the entire time period
salt_summary <- salt_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(salt_5_mean = mean(mean_5_year, na.rm = TRUE),  salt_5_min = mean(min_5_year, na.rm = TRUE), salt_5_cv = mean(cv_5_year, na.rm = TRUE))
#Combine the above dataframe with the salt grid and filter out erroneous entries
salt_summary_grid <- dplyr::full_join(salt_summary, BCCM_h_salt_grid, by=c("xi_rho","eta_rho")) %>% na.omit()
#Convert to spatial points and transform to grid crs
salt_points <- sf::st_as_sf(salt_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_salt_5_mean_rast <- terra::rasterize(salt_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(salt_points), res=3000),field = "salt_5_mean", fun=mean)
bccm_salt_5_min_rast <- terra::rasterize(salt_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(salt_points), res=3000),field = "salt_5_min", fun=mean)
bccm_salt_5_cv_rast <- terra::rasterize(salt_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(salt_points), res=3000),field = "salt_5_cv", fun=mean)

#subset the PAR data by these start and end year
urad_sub <- BCCM_h_urad_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year (actually already summarized by year so this doesn't do anything but relabel)
PAR_summary_year <- urad_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE), max_5_year = mean(max_5m, na.rm = TRUE), min_5_year = mean(min_5m, na.rm = TRUE))
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

#subset the do data by these start and end year
do_sub <- BCCM_h_do_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
do_summary_year <- do_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_5_year = mean(mean_5m, na.rm = TRUE), min_5_year = min(mean_5m, na.rm = TRUE))
#Summarise the do data into the desired climatologies for the entire time period
do_summary <- do_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(do_5_mean = mean(mean_5_year, na.rm = TRUE), do_5_min = mean(min_5_year, na.rm = TRUE))
#Combine the above dataframe with the sshPAR grid and filter out erroneous entries
do_summary_grid <- dplyr::full_join(do_summary, BCCM_h_do_grid, by=c("xi_rho","eta_rho")) %>% na.omit() 
#Convert to spatial points and transform to grid crs
do_points <- sf::st_as_sf(do_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of interest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_do_5_mean_rast <- terra::rasterize(do_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(do_points), res=3000),field = "do_5_mean", fun=mean)
bccm_do_5_min_rast <- terra::rasterize(do_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(do_points), res=3000),field = "do_5_min", fun=mean)

#subset the temp data by these start and end year
temp_sub <- BCCM_h_temp_data %>% dplyr::filter(year >= start & year <= end)
#Summarise the predictor by year
temp_summary_year <- temp_sub %>% dplyr::group_by(xi_rho, eta_rho, year) %>% dplyr::summarise(mean_sur_year = mean(mean_sur, na.rm = TRUE), max_sur_year = max(mean_sur, na.rm = TRUE), min_sur_year = min(mean_sur, na.rm = TRUE), mean_5_year = mean(mean_5m, na.rm = TRUE), max_5_year = max(mean_5m, na.rm = TRUE), min_5_year = min(mean_5m, na.rm = TRUE), maxdiff_5_year = max(diff_5m, na.rm = TRUE), maxdiff_surf_year = max(diff_sur, na.rm = TRUE), sd_5_year = sd(mean_5m, na.rm = TRUE), sd_sur_year = sd(mean_sur, na.rm = TRUE)) %>%
  mutate(cv_5_year = sd_5_year/mean_5_year*100, cv_sur_year = sd_sur_year/mean_sur_year*100)
#Summarise the temp data into the desired climatologies for the entire time period
temp_summary <- temp_summary_year %>% dplyr::group_by(xi_rho, eta_rho) %>% dplyr::summarise(temp_s_mean = mean(mean_sur_year, na.rm = TRUE), temp_s_max = mean(max_sur_year, na.rm = TRUE), temp_s_min = mean(min_sur_year, na.rm = TRUE), temp_5_mean = mean(mean_5_year, na.rm = TRUE), temp_5_max = mean(max_5_year, na.rm = TRUE), temp_5_min = mean(min_5_year, na.rm = TRUE), temp_5_diff = mean(maxdiff_5_year, na.rm = TRUE), temp_sur_diff = mean(maxdiff_surf_year, na.rm = TRUE), temp_5_cv = mean(cv_5_year, na.rm = TRUE), temp_sur_cv = mean(cv_sur_year, na.rm = TRUE))
#Combine the above dataframe with the temp grid and filter out erroneous entries
temp_summary_grid <- dplyr::full_join(temp_summary, BCCM_h_temp_grid, by=c("xi_rho","eta_rho")) %>% na.omit() # %>% dplyr::filter(temp_s_min >= 0, temp_10_min >= 0)
#Convert to spatial points and transform to grid crs
temp_points <- sf::st_as_sf(temp_summary_grid, coords=(c("lon_rho", "lat_rho")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of interest onto the 3km equal area grid (size of BCCM pixels) with mean = akin to project
bccm_temp_s_mean_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_s_mean", fun=mean)
bccm_temp_s_max_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_s_max", fun=mean)
bccm_temp_s_min_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_s_min", fun=mean)
bccm_temp_s_diff_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_sur_diff", fun=mean)
bccm_temp_s_cv_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_sur_cv", fun=mean)
bccm_temp_5_mean_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_mean", fun=mean)
bccm_temp_5_max_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_max", fun=mean)
bccm_temp_5_min_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_min", fun=mean)
bccm_temp_5_diff_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_diff", fun=mean)
bccm_temp_5_cv_rast <- terra::rasterize(temp_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(temp_points), res=3000),field = "temp_5_cv", fun=mean)


#### SalishSeaCast Yearly data for NH4, NO3, PAR, DO####

#subset the year data by these start and end year, and just surface values
Year_surf_sub <- SSC_year_data %>% 
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
#Rasterize all  predictors of interest onto the 500m equal area grid (size of salishseacast pixels) with mean = akin to project, merge with surface values because loose very shallow areas
ssc_PAR_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "PAR_mean_sur", fun=mean)
ssc_PAR_sur_max_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "PAR_max_sur", fun=mean)
ssc_PAR_sur_min_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "PAR_min_sur", fun=mean)
ssc_NH4_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "NH4_mean_sur", fun=mean)
ssc_NO3_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "NO3_mean_sur", fun=mean)
ssc_DO_sur_mean_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "DO_mean_sur", fun=mean)
ssc_DO_sur_min_rast <- terra::rasterize(year_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(year_surf_points), res=500),field = "DO_min_sur", fun=mean)

#subset the year data by these start and end year, and just 5m values
Year_5_sub <- SSC_year_data %>% 
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

#### SalishSeaCast Monthly data for temperature and salinity####

#subset the year data by these start and end year, and just surface values
Month_surf_sub <- SSC_month_data %>% 
  dplyr::filter(year >= start & year <= end,
                deptht < 2)
#Summarise the predictor by year
month_surf_summary_year <- Month_surf_sub %>% dplyr::group_by(x, y, year) %>% 
  dplyr::summarise(salt_mean_sur_year = mean(saltmean, na.rm = TRUE), salt_min_sur_year = min(saltmean, na.rm = TRUE), salt_sd_sur_year = sd(saltmean, na.rm = TRUE),
                   temp_mean_sur_year = mean(tempmean, na.rm = TRUE), temp_max_sur_year = max(tempmean, na.rm = TRUE), temp_min_sur_year = min(tempmean, na.rm = TRUE), temp_sd_sur_year = sd(tempmean, na.rm = TRUE),
                   temp_diff_sur_year = max(tempdiff, na.rm = TRUE)) %>%
  mutate(salt_cv_sur_year = salt_sd_sur_year/salt_mean_sur_year*100,
         temp_cv_sur_year = temp_sd_sur_year/temp_mean_sur_year*100)

#Summarise the predictor into the desired climatologies for the entire time period
month_surf_summary <- month_surf_summary_year %>% dplyr::group_by(x, y) %>% 
  dplyr::summarise(salt_mean_sur = mean(salt_mean_sur_year, na.rm = TRUE), salt_min_sur = mean(salt_min_sur_year, na.rm = TRUE), salt_cv_sur = mean(salt_cv_sur_year, na.rm = TRUE),
                   temp_mean_sur = mean(temp_mean_sur_year, na.rm = TRUE), temp_max_sur = mean(temp_max_sur_year, na.rm = TRUE), temp_min_sur = mean(temp_min_sur_year, na.rm = TRUE),
                   temp_diff_sur = mean(temp_diff_sur_year, na.rm = TRUE), temp_cv_sur = mean(temp_cv_sur_year, na.rm = TRUE))

#Combine the above dataframe with the grid and filter out erroneous entries
month_surf_summary_grid <- dplyr::full_join(month_surf_summary, SSC_dfo_month_grid, by=c("x","y")) %>% na.omit() %>% dplyr::filter(salt_mean_sur > 0)
#Convert to spatial points and transform to grid crs
month_surf_points <- sf::st_as_sf(month_surf_summary_grid, coords=(c("nav_lon", "nav_lat")), crs = sf::st_crs(4326)) %>% sf::st_transform(crs = sf::st_crs(3573))
#Rasterize all  predictors of intrerest onto the 500m equal area grid (size of salishseacast pixels) with mean = akin to project, merge with surface values because loose very shallow areas
ssc_temp_sur_mean_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_mean_sur", fun=mean)
ssc_temp_sur_max_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_max_sur", fun=mean)
ssc_temp_sur_min_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_min_sur", fun=mean)
ssc_temp_sur_cv_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_cv_sur", fun=mean)
ssc_temp_sur_diff_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "temp_diff_sur", fun=mean)
ssc_salt_sur_mean_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "salt_mean_sur", fun=mean)
ssc_salt_sur_min_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "salt_min_sur", fun=mean)
ssc_salt_sur_cv_rast <- terra::rasterize(month_surf_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_surf_points), res=500),field = "salt_cv_sur", fun=mean)

#subset the year data by these start and end year, and just 5m values
Month_5_sub <- SSC_month_data %>% 
  dplyr::filter(year >= start & year <= end,
                deptht > 2)

#Summarise the predictor by year
month_5_summary_year <- Month_5_sub %>% dplyr::group_by(x, y, year) %>% 
  dplyr::summarise(salt_mean_5_year = mean(saltmean, na.rm = TRUE), salt_min_5_year = min(saltmean, na.rm = TRUE), salt_sd_5_year = sd(saltmean, na.rm = TRUE),
                   temp_mean_5_year = mean(tempmean, na.rm = TRUE), temp_max_5_year = max(tempmean, na.rm = TRUE), temp_min_5_year = min(tempmean, na.rm = TRUE), temp_sd_5_year = sd(tempmean, na.rm = TRUE),
                   temp_diff_5_year = max(tempdiff, na.rm = TRUE)) %>%
  mutate(salt_cv_5_year = salt_sd_5_year/salt_mean_5_year*100,
         temp_cv_5_year = temp_sd_5_year/temp_mean_5_year*100)
#Summarise the predictor into the desired climatologies for the entire time period
month_5_summary <- month_5_summary_year %>% dplyr::group_by(x, y) %>% 
  dplyr::summarise(salt_mean_5 = mean(salt_mean_5_year, na.rm = TRUE), salt_min_5 = mean(salt_min_5_year, na.rm = TRUE), salt_cv_5 = mean(salt_cv_5_year, na.rm = TRUE),
                   temp_mean_5 = mean(temp_mean_5_year, na.rm = TRUE), temp_max_5 = mean(temp_max_5_year, na.rm = TRUE), temp_min_5 = mean(temp_min_5_year, na.rm = TRUE),
                   temp_diff_5 = mean(temp_diff_5_year, na.rm = TRUE), temp_cv_5 = mean(temp_cv_5_year, na.rm = TRUE))

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
ssc_temp_5_diff_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "temp_diff_5", fun=mean)%>%
  terra::merge(ssc_temp_sur_diff_rast, na.rm=TRUE)
ssc_temp_5_cv_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "temp_cv_5", fun=mean)%>%
  terra::merge(ssc_temp_sur_cv_rast, na.rm=TRUE)
ssc_salt_5_mean_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "salt_mean_5", fun=mean)%>%
  terra::merge(ssc_salt_sur_mean_rast, na.rm=TRUE)
ssc_salt_5_min_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "salt_min_5", fun=mean)%>%
  terra::merge(ssc_salt_sur_min_rast, na.rm=TRUE)
ssc_salt_5_cv_rast <- terra::rasterize(month_5_points, terra::rast(crs ="EPSG:3573", extent = terra::ext(month_5_points), res=500),field = "salt_cv_5", fun=mean)%>%
  terra::merge(ssc_salt_sur_cv_rast, na.rm=TRUE)

#Stack all of the layers together
BCCM_layers <- c(bccm_NH4_5_mean_rast, bccm_NO3_5_mean_rast, bccm_salt_5_mean_rast, bccm_salt_5_min_rast, bccm_PAR_5_mean_rast, bccm_PAR_5_min_rast, 
                 bccm_PAR_5_max_rast, bccm_temp_s_mean_rast, bccm_temp_s_max_rast, bccm_temp_s_min_rast, bccm_temp_5_mean_rast, bccm_temp_5_max_rast, 
                 bccm_temp_5_min_rast, bccm_do_5_mean_rast, bccm_do_5_min_rast, bccm_salt_5_cv_rast, bccm_temp_5_cv_rast, bccm_temp_5_diff_rast,
                 bccm_temp_s_cv_rast, bccm_temp_s_diff_rast)
SSC_layers <- c(ssc_NH4_5_mean_rast, ssc_NO3_5_mean_rast, ssc_salt_5_mean_rast, ssc_salt_5_min_rast, ssc_PAR_5_mean_rast, ssc_PAR_5_min_rast, 
                ssc_PAR_5_max_rast, ssc_temp_sur_mean_rast, ssc_temp_sur_max_rast, ssc_temp_sur_min_rast, ssc_temp_5_mean_rast, ssc_temp_5_max_rast, 
                ssc_temp_5_min_rast, ssc_DO_5_mean_rast, ssc_DO_5_min_rast, ssc_salt_5_cv_rast, ssc_temp_5_cv_rast, ssc_temp_5_diff_rast,
                ssc_temp_sur_cv_rast, ssc_temp_sur_diff_rast)

# resample (and reproject) onto DFO bathy grid
BCCM_Resampled <- terra::sapp(BCCM_layers, function(x){
  terra::focal(x,w=3, fun = "mean", na.policy = "only", na.rm = TRUE) %>% terra::project(bathy20m, method="cubicspline") 
})
#applied larger w to ssc layers otherwise was having bccm layers selected over ssc in nearshore shallow areas
#extrapolation didn't seem too unreasonable 
SSC_Resampled <- terra::sapp(SSC_layers, function(x){
  terra::focal(x,w=7, fun = "mean", na.policy = "only", na.rm = TRUE) %>% terra::project(bathy20m, method="cubicspline") 
})

#merge ssc and bccm layers, chooses ssc layer over bccm if present in same cells. Mask it by bathy layer
NH4_5_mean_rast <- terra::merge(SSC_Resampled[[1]], BCCM_Resampled[[1]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
NO3_5_mean_rast <- terra::merge(SSC_Resampled[[2]], BCCM_Resampled[[2]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
salt_5_mean_rast <- terra::merge(SSC_Resampled[[3]], BCCM_Resampled[[3]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
salt_5_min_rast <- terra::merge(SSC_Resampled[[4]], BCCM_Resampled[[4]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
PAR_5_mean_rast <- terra::merge(SSC_Resampled[[5]], BCCM_Resampled[[5]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
PAR_5_min_rast <- terra::merge(SSC_Resampled[[6]], BCCM_Resampled[[6]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
PAR_5_max_rast <- terra::merge(SSC_Resampled[[7]], BCCM_Resampled[[7]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_s_mean_rast <- terra::merge(SSC_Resampled[[8]], BCCM_Resampled[[8]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_s_max_rast <- terra::merge(SSC_Resampled[[9]], BCCM_Resampled[[9]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_s_min_rast <- terra::merge(SSC_Resampled[[10]], BCCM_Resampled[[10]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_5_mean_rast <- terra::merge(SSC_Resampled[[11]], BCCM_Resampled[[11]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_5_max_rast <- terra::merge(SSC_Resampled[[12]], BCCM_Resampled[[12]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m) 
temp_5_min_rast <- terra::merge(SSC_Resampled[[13]], BCCM_Resampled[[13]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
do_5_mean_rast <- terra::merge(SSC_Resampled[[14]], BCCM_Resampled[[14]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
do_5_min_rast <- terra::merge(SSC_Resampled[[15]], BCCM_Resampled[[15]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
salt_5_cv_rast <- terra::merge(SSC_Resampled[[16]], BCCM_Resampled[[16]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_5_cv_rast<- terra::merge(SSC_Resampled[[17]], BCCM_Resampled[[17]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_5_diff_rast<- terra::merge(SSC_Resampled[[18]], BCCM_Resampled[[18]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_s_cv_rast<- terra::merge(SSC_Resampled[[19]], BCCM_Resampled[[19]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)
temp_s_diff_rast<- terra::merge(SSC_Resampled[[20]], BCCM_Resampled[[20]], na.rm=TRUE) %>% terra::mask(bathy20m) %>% terra::crop(bathy20m)

Predictor_Hindcast_Climatologies <- c(NH4_5_mean_rast, NO3_5_mean_rast, salt_5_mean_rast, salt_5_min_rast, PAR_5_mean_rast, PAR_5_min_rast, 
                 PAR_5_max_rast, temp_s_mean_rast, temp_s_max_rast, temp_s_min_rast, temp_5_mean_rast, temp_5_max_rast, 
                 temp_5_min_rast, do_5_mean_rast, do_5_min_rast, salt_5_cv_rast, temp_5_cv_rast, temp_5_diff_rast, temp_s_cv_rast, temp_s_diff_rast)

names(Predictor_Hindcast_Climatologies) <- c("NH4_5m_mean", "NO3_5m_mean", "salt_5m_mean", "salt_5m_min", "PAR_5m_mean", "PAR_5m_min", 
                                             "PAR_5m_max", "temp_s_mean", "temp_s_max", "temp_s_min", "temp_5m_mean", "temp_5m_max", 
                                             "temp_5m_min", "do_5m_mean", "do_5m_min", "salt_5m_cv", "temp_5m_cv", "temp_5m_diff", "temp_s_cv", "temp_s_diff")

#Save the hindcast predictor raster
terra::writeRaster(Predictor_Hindcast_Climatologies, paste0("code/output_data/processed_ocean_variables/Predictor_Hindcast_Climatologies_", paste(start),"-",paste(end),".tif"), overwrite = TRUE)
#checking individual rasters
#terra::writeRaster(Predictor_Hindcast_Climatologies[[15]], paste0("code/output_data/processed_ocean_variables/Predictor_Hindcast_Climatologies_", names(Predictor_Hindcast_Climatologies[[15]]),".tif"), overwrite = TRUE)

# comparing 2013-2023 data
# PARs all look weird across Salish Sea and BCCM boundaries.
# temp surface min values go below 0, so not going to use those. Matt used values below 0 filtered out to remove some inlets
# Also planning to filter out DO <100 as that removes some inlets that are not predicted well by BCCM
# NH4 is high in the fraser, some inlets (holberg), and in spots around HG; 
#NO3 is high in the Broughton and Holberg Inlet
# salt mean looks good. Salt min only changes around Nass, Skeena, Fraser, and a few major inlets, otherwise looks similar to salt mean

#### Culmulative effects layer from Murray, change from polygon to 20 m raster
cul_eff <- vect(st_read(dsn = "raw_data/anthropogenic/Cumulative_Impacts_Pacfic_Canada.gdb", layer = "Cumulative_Impacts_Pacific_Canada")) # impact score representing impacts from all activities and all habitats in each PU grid cell.
cul_eff_rast20 <- rasterize(cul_eff, bathy20m, field = "Cumul_Impact_ALL", fun = max) %>% 
  terra::focal(w=21, fun = "mean", na.policy = "only", na.rm = TRUE) %>%  # only change cells that are NA to fill in holes and up to coastline
  terra::mask(bathy20m) %>% terra::crop(bathy20m) 
writeRaster(cul_eff_rast20, file.path("code/output_data/culmulative_effects_all_20m.tif"), overwrite=TRUE)
