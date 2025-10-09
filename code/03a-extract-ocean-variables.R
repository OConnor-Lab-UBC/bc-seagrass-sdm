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
# extract oceanography data from BCCM, NEP36, Salish Sea Cast, CHELSA, WorldClim
#
# STILL TO DO
# ADD NEP36 data and atmospheric data
###############################################################################

#load packages####
library(tidyverse)
library(tidync)
library(ncmeta)
library(lubridate)

##########################
###### BCCM data ########
#########################

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

save(BCCM_h_NH4_data, BCCM_h_NH4_grid, BCCM_h_NO3_data, BCCM_h_NO3_grid, BCCM_h_salt_data, BCCM_h_salt_grid, BCCM_h_urad_data, BCCM_h_urad_grid,
     BCCM_h_do_data, BCCM_h_do_grid, BCCM_h_temp_data, BCCM_h_temp_grid, file = "code/output_data/intermediate_ocean_variables/BCCM_data.RData")


###################################
###### Salish Sea Cast data #######
##################################
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

save(SSC_year_data, dfo_grid, file = "code/output_data/intermediate_ocean_variables/SSC_year_data.RData")



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

save(SSC_month_data, dfo_mgrid, file = "code/output_data/intermediate_ocean_variables/SSC_month_data.RData")
 
##########################
###### NEP36 data #######
#########################


#Annual MinMax is highest/lowest monthly mean within a given year; and Monthly MinMax is highest/lowest daily mean within given month

#### yearly data for NH4, just need mean
#read in the nc file
NEP36_NH4_yearmean_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmean_1m_ptrc_T_1996-2024.nc")

time_units_info <- ncmeta::nc_atts("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmean_1m_ptrc_T_1996-2024.nc", "time")
time_units <- time_units_info$value[time_units_info$name == "units"]
origin_string <- sub(".*since ", "", time_units)

NEP36_NH4_yearmean_data <- tidync::hyper_tibble(NEP36_NH4_yearmean_nc) %>%
  mutate(year = year(as.POSIXct(time * 86400,
                                origin = origin_string,
                                tz = "UTC"))) %>%
  rename(NH4mean = NH4) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(time))

###yearly data for alkalinity, DIC, NO3, O2, don't need max for any of these variables
NEP36_yearmean_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmean_1d_mask_T_1996-2024.nc")
NEP36_yearmin_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmin_1d_mask_T_1996-2024.nc")

time_units_info <- ncmeta::nc_atts("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmean_1d_mask_T_1996-2024.nc", "time_counter")
time_units <- time_units_info$value[time_units_info$name == "units"]
origin_string <- sub(".*since ", "", time_units)

NEP36_yearmean_data <- tidync::hyper_tibble(NEP36_yearmean_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time_counter,
                                origin = origin_string,
                                tz = "UTC")))%>%
  rename(DOmean = O2_mask, NO3mean = NO3_mask) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(Alkalini_mask, DIC_mask, salt_mask, temp_mask, time_counter))

NEP36_yearmin_data <- tidync::hyper_tibble(NEP36_yearmin_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time_counter,
                                origin = origin_string,
                                tz = "UTC")))%>%
  rename(DOmin = O2_mask) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(Alkalini_mask, DIC_mask, salt_mask, NO3_mask, temp_mask, time_counter))


### yearly data for PH and PAR
NEP36_PH_PAR_yearmax_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmax_1m_diad_T_1996-2024.nc")
NEP36_PH_PAR_yearmean_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmean_1m_diad_T_1996-2024.nc")
NEP36_PH_PAR_yearmin_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmin_1m_diad_T_1996-2024.nc")

time_units_info <- ncmeta::nc_atts("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-yearmax_1m_diad_T_1996-2024.nc", "time")
time_units <- time_units_info$value[time_units_info$name == "units"]
origin_string <- sub(".*since ", "", time_units)

NEP36_PAR_yearmax_data <- tidync::hyper_tibble(NEP36_PH_PAR_yearmax_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time * 86400,
                                origin = origin_string,
                                tz = "UTC"))) %>%
  rename(PARmax = PAR) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(PH, time))

NEP36_PAR_yearmean_data <- tidync::hyper_tibble(NEP36_PH_PAR_yearmean_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time * 86400,
                                origin = origin_string,
                                tz = "UTC"))) %>%
  rename(PARmean = PAR) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(PH, time))

NEP36_PAR_yearmin_data <- tidync::hyper_tibble(NEP36_PH_PAR_yearmin_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time * 86400,
                                origin = origin_string,
                                tz = "UTC"))) %>%
  rename(PARmin = PAR) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(PH, time))

#all NEP36 yearly data for PAR, O2, No3, NH4
NEP36_year_data <- NEP36_NH4_yearmean_data %>%
  full_join(NEP36_yearmean_data, by = join_by(x, y, deptht, year)) %>%
  full_join(NEP36_yearmin_data, by = join_by(x, y, deptht, year))%>%
  full_join(NEP36_PAR_yearmax_data, by = join_by(x, y, deptht, year))%>%
  full_join(NEP36_PAR_yearmean_data, by = join_by(x, y, deptht, year))%>%
  full_join(NEP36_PAR_yearmin_data, by = join_by(x, y, deptht, year))

#grid 
NEP36_year_grid <- NEP36_yearmean_nc %>% tidync::activate("D2,D3") %>% tidync::hyper_tibble()

save(NEP36_year_data, NEP36_year_grid, file = "code/output_data/intermediate_ocean_variables/NEP36_year_data.RData")

###monthly data for salinity, temperature
NEP36_monmax_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-monmax_1d_mask_T_1996-2024.nc")
NEP36_monmean_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-monmean_1d_mask_T_1996-2024.nc")
NEP36_monmin_nc <- tidync::tidync("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-monmin_1d_mask_T_1996-2024.nc")

time_units_info <- ncmeta::nc_atts("raw_data/NEP36/hindcast/NEP36-CanOE-TKE-monmax_1d_mask_T_1996-2024.nc", "time_counter")
time_units <- time_units_info$value[time_units_info$name == "units"]
origin_string <- sub(".*since ", "", time_units)

NEP36_monmax_data <- tidync::hyper_tibble(NEP36_monmax_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time_counter,
                                origin = origin_string,
                                tz = "UTC")))%>%
  rename(tempmax = temp_mask) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(Alkalini_mask, DIC_mask, NO3_mask, O2_mask, salt_mask))

#save(NEP36_monmax_data, file = "code/output_data/intermediate_ocean_variables/NEP36_monmax_data.RData")

NEP36_monmean_data <- tidync::hyper_tibble(NEP36_monmean_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time_counter,
                                origin = origin_string,
                                tz = "UTC")))%>%
  rename(tempmean = temp_mask, saltmean = salt_mask) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(Alkalini_mask, DIC_mask, NO3_mask, O2_mask))
#save(NEP36_monmean_data, file = "code/output_data/intermediate_ocean_variables/NEP36_monmean_data.RData")

NEP36_monmin_data <- tidync::hyper_tibble(NEP36_monmin_nc, force = TRUE) %>%
  mutate(year = year(as.POSIXct(time_counter,
                                origin = origin_string,
                                tz = "UTC")))%>%
  rename(tempmin = temp_mask) %>%
  filter(deptht > 0 & deptht < 6) %>%
  select(-c(Alkalini_mask, DIC_mask, NO3_mask, O2_mask, salt_mask))
#save(NEP36_monmin_data, file = "code/output_data/intermediate_ocean_variables/NEP36_monmin_data.RData")

NEP36_month_data <- NEP36_monmin_data %>%
  full_join(NEP36_monmean_data, by = join_by(x, y, deptht, time_counter, year)) %>%
  full_join(NEP36_monmax_data, by = join_by(x, y, deptht, time_counter, year)) %>%
  mutate(tempdiff = tempmax - tempmin)

#grid 
NEP36_month_grid <- NEP36_monmin_nc %>% tidync::activate("D2,D3") %>% tidync::hyper_tibble()
save(NEP36_month_data, NEP36_month_grid, file = "code/output_data/intermediate_ocean_variables/NEP36_month_data.RData")


##################
####CHELSA######
#################
# 1km monthly from 1981 to 2019
# Has temp, temp max, temp min, precipitation, surface downwelling shorewave flux, wind speeds
# tas is Daily mean air temperatures at 2 metres from hourly ERA5 data in units K/10. no scale or offset
# pr is "Amount" means mass per unit area. "Precipitation" in the earth's atmosphere means precipitation of water in all phases. units are kg m-2 month-1/100. no scale or offset
# rsds is Surface downwelling shortwave flux in air. Attenuating effects of clouds are accounted for units MJ m-2 d-1. Scale is 0.001 and offset is 0

library(curl)
library(terra)
library(tidyverse)

# would crash so just downloaded seperate variables, sometimes would have to restart a certain year/month
url_file <- "raw_data/CHELSA/envidatS3paths.txt"
urls <- trimws(readLines(url_file))
urls <- urls[nchar(urls) > 0] # remove blank lines

out_csv <- "raw_data/CHELSA/CHELSA_BC_data_PR2.csv"

# If it already exists, remove it (to start clean)
if (file.exists(out_csv)) file.remove(out_csv)

#out_dir  <- "raw_data/CHELSA/BC"
#dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# BC bounding box (xmin, xmax, ymin, ymax)
coastal_bc_bbox  <- ext(-134, -122, 48, 56)

process_url <- function(u) {
  fname <- basename(u)
  
  # detect variable type
  var <- case_when(
    grepl("tas", fname)  ~ "tas",
    grepl("pr", fname)   ~ "pr",
    grepl("rsds", fname) ~ "rsds",
    TRUE                 ~ "unknown"
  )
  
  # parse month and year based on variable
  if (var %in% c("pr", "tas")) {
    # format: CHELSA_<var>_<MM>_<YYYY>_V*.tif
    mo <- as.integer(sub(".*_(\\d{2})_\\d{4}_V.*", "\\1", fname))
    yr <- as.integer(sub(".*_\\d{2}_(\\d{4})_V.*", "\\1", fname))
  } else if (var == "rsds") {
    # format: CHELSA_rsds_<YYYY>_<MM>_V*.tif
    mo <- as.integer(sub(".*_(\\d{2})_V.*", "\\1", fname))
    yr <- as.integer(sub(".*_(\\d{4})_\\d{2}_V.*", "\\1", fname))
  } else {
    mo <- NA
    yr <- NA
  }
  
  message("→ Processing ", var, " ", yr, "-", sprintf("%02d", mo))
  
  tmp <- tempfile(fileext = ".tif")
  tryCatch({
    curl_download(u, tmp)
    r <- rast(tmp)
    
    # crop to coastal BC
    r_bc <- crop(r, coastal_bc_bbox)
    
    # convert to dataframe
    vals <- as.data.frame(r_bc, xy = TRUE, cells = TRUE, na.rm = TRUE)
    names(vals)[4] <- "value"
    
    vals <- vals %>%
      mutate(year = yr,
             month = mo,
             variable = var)
    
    # write directly to CSV
    write.table(vals, out_csv, sep = ",", col.names = !file.exists(out_csv),
                row.names = FALSE, append = TRUE)
    
    unlink(tmp)
  }, error = function(e) {
    message("Failed: ", fname, " (", e$message, ")")
  })
}

for (u in urls) {
  process_url(u)
}

message("Extraction complete. Data saved to: ", out_csv)


#then bring csv together with 
Chelsa_month_pr1 <- read_csv("raw_data/CHELSA/CHELSA_BC_data_PR.csv") 
Chelsa_month_pr2 <- read_csv("raw_data/CHELSA/CHELSA_BC_data_PR2.csv")
unique_year_month1 <- Chelsa_month_pr1 %>%
  distinct(year, month)
unique_year_month2 <- Chelsa_month_pr2 %>%
  distinct(year, month)
Chelsa_month_pr1<- Chelsa_month_pr1 %>% 
  filter(month < 5)
# missing may, download may seperately. 

Chelsa_month_pr3 <- read_csv("raw_data/CHELSA/CHELSA_BC_data_PR3.csv")
Chelsa_month_pr <- rbind(Chelsa_month_pr1, Chelsa_month_pr2, Chelsa_month_pr3)
rm(Chelsa_month_pr1)
rm(Chelsa_month_pr2)
rm(Chelsa_month_pr3)

Chelsa_month_pr <- Chelsa_month_pr %>% 
  rename(precip = value) %>% select(-variable)
unique_year_month <- Chelsa_month_pr %>%
  distinct(year, month)


Chelsa_month_tas1 <- read_csv("raw_data/CHELSA/CHELSA_BC_data_tas.csv")
Chelsa_month_tas2 <- read_csv("raw_data/CHELSA/CHELSA_BC_data_tas2.csv")
unique_year_month1 <- Chelsa_month_tas1 %>%
  distinct(year, month)
unique_year_month2 <- Chelsa_month_tas2 %>%
  distinct(year, month)
# not missing any 

Chelsa_month_tas <- rbind(Chelsa_month_tas1, Chelsa_month_tas2)%>%
  rename(tempmean_air = value) %>% select(-variable)
rm(Chelsa_month_tas1)
rm(Chelsa_month_tas2)
unique_year_month <- Chelsa_month_tas %>%
  distinct(year, month)


Chelsa_month_rsds1 <- read_csv("raw_data/CHELSA/CHELSA_BC_data_rsds.csv")
Chelsa_month_rsds2 <- read_csv("raw_data/CHELSA/CHELSA_BC_data_rsds2.csv")
unique_year_month1 <- Chelsa_month_rsds1 %>%
  distinct(year, month)
unique_year_month2 <- Chelsa_month_rsds2 %>%
  distinct(year, month)
# not missing any


Chelsa_month_rsds <- rbind(Chelsa_month_rsds1, Chelsa_month_rsds2)%>%
  rename(rsds = value) %>% select(-variable)
rm(Chelsa_month_rsds1)
rm(Chelsa_month_rsds2)


CHELSA_month_data <- Chelsa_month_pr %>%
  full_join(Chelsa_month_tas, by = join_by(x, y, month, year)) %>%
  full_join(Chelsa_month_rsds, by = join_by(x, y, month, year)) 

save(CHELSA_month_data, file = "code/output_data/intermediate_ocean_variables/CHELSA_month_data.RData")


