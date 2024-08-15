###############################################################################
#
# Authors:      Ashley Park (code modified from Jessica Nephin and Sarah Friesen)
# Affiliation:  Fisheries and Oceans Canada (DFO)
# Group:        Marine Spatial Ecology and Analysis
# Location:     Institute of Ocean Sciences
# Contact:      e-mail: ashley.park@dfo-mpo.gc.ca
# Project:      UVic turbidity layers
#
# Overview:

#
###############################################################################


# required packages
require(parallel)
require(ncdf4)
require(dplyr)
require(reshape2)
require(rgdal)
require(raster)
require(terra)
require(ggplot2)


# working directory
# setwd('.')

#get input bathy raster files for grid template
# Load bathy raster
bathy <- vrt("raw_data/bathymetry-20m/strait_of_georgia.tif")

# #set the CRS
# geoCRS <- crs(eHG)

#####Create grid
#create empty raster
# ras<-merge(eHG,eNCC,eQCS,eSASIMA,eWCVI)
# ras[]<-NA


#------------------------------------------------------------------------------#


# get input salish sea cast files

salishsea <- tidync("raw_data/salishseacast/SalishSeaCast-VNR023_1d_grid_T_19900901-19900901.nc")

#Load data and grid
salishseaData <- hyper_tibble(salishsea)
salishseaData <- salishseaData %>% filter(deptht<=0.5000003)
salishseaGrid <- salishsea %>% activate("D3,D2") %>% hyper_tibble()
salishseaGrid <- salishseaGrid %>% filter(nav_lat!=0)

#Combine with grid and filter lad value 
salishseav_grid <- full_join(salishseaData,salishseaGrid,by=c("x","y")) %>%
  filter(vosaline>0)

#Convert to spatial points and transform to grid crs
salishseaPoints <- st_as_sf(salishseav_grid, coords=(c("nav_lon","nav_lat")), crs = st_crs(4326))
salishseaPointseq <- st_transform(salishseaPoints,crs = st_crs(3573))

#Rasterize onto 500m equal area 
salishseaSalRast <- rasterize(salishseaPointseq,rast(crs="EPSG:3573",extent=ext(salishseaPointseq),res=500),field="vosaline", fun=mean)
#Allow to extrapolate by 1 pixel
salishseaSalRast <- focal(salishseaSalRast,w=3,fun="mean",na.policy="only")
writeRaster(salishseaSalRast,"raster/salishseaSalinity.tif", overwrite=TRUE)



wrf <- tidync("raw_data/WRF/slp_WRF_CanESM2_d03_1991-03-16.nc")

#Load data and grid
wrfData <- hyper_tibble(wrf)
wrfGrid <- wrf %>% activate("D0,D1") %>% hyper_tibble()
wrfv_grid <- full_join(wrfData,wrfGrid,by=c("x","y"))

#Convert to spatial points and transform to grid crs
wrfPoints <- st_as_sf(wrfv_grid, coords=(c("lon","lat")), crs = st_crs(4326))
wrfPointseq <- st_transform(wrfPoints,crs = st_crs(3573))

#Rasterize onto 3km equal area 
wrfRast <- rasterize(wrfPointseq,rast(crs="EPSG:3573",extent=ext(wrfPointseq),res=3000),field="slp", fun=mean)
#Allow to extrapolate by 1 pixel
wrfRast <- focal(wrfRast,w=3,fun="mean",na.policy="only")
writeRaster(wrfRast,"raster/wrf_slp.tif", overwrite=TRUE)


salishsea_tempsal <- "raw_data/salishseacast/SalishSeaCast-VNR023_1d_grid_T_19900901-19900901.nc"

nc <- nc_open(salishsea_tempsal)

lat <- ncvar_get(nc, "vosaline")


salinity <- rast(salishsea_tempsal, crs = 	"EPSG:4326")
, paste0(c("nav_lat", "nav_lon", "vosaline")))

salinity <- project(salinity, bathy) 

#TSM HG
salinity_ss <- terra::resample(salinity, bathy, method ="bilinear")
writeRaster(salinity_ss, file.path("./raster/salinity.tif"), overwrite=TRUE)

ggplot(salinity_ss)+
  geom_tile(aes(x = X_m, y = Y_m, colour=vosaline, width=20,height=20))+
  scale_colour_viridis_c()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  coord_sf(expand = FALSE)+
  ylab("")+
  xlab("") 




#
#Take mean value across years
BCCMAv <- BCCMData %>%
  group_by(xi_rho,eta_rho)%>%
  summarise(u = mean(u_eastward), v = mean(v_northward), Sal = mean(Salinity), Temp = mean(Temperature))

#Combine with grid and filter lad value 
BCCMAv_grid <- full_join(BCCMAv,BCCMGrid,by=c("xi_rho","eta_rho")) %>%
  filter(Sal>0.52) #Remove land and 1 cell error in model output equal to 0.5 salinity

#Convert to spatial points and transform to grid crs
BCCMPoints <- st_as_sf(BCCMAv_grid, coords=(c("lon_rho","lat_rho")), crs = st_crs(4326))
BCCMPointseq <- st_transform(BCCMPoints,crs = st_crs(3573))


#Rasterize onto 3km equal area (size of BCCM pixels) with mean = akin to project
BCCMSalRast <- rasterize(BCCMPointseq,rast(crs="EPSG:3573",extent=ext(BCCMPointseq),res=3000),field="Sal", fun=mean)
#Allow to extrapolate by 1 pixel
BCCMSalRast <- focal(BCCMSalRast,w=3,fun="mean",na.policy="only")
writeRaster(BCCMSalRast,"raster/BCCMSalinity.tif", overwrite=TRUE)


# 
# #------------------------------------------------------------------------------#
# 
# # list all nc files and combine with season names
# ncfiles <- list.files(path = "netcdf/", pattern = ".nc", full.names = TRUE)
# 
# # get x,y data
# ncxy <- nc_open(filename=ncfiles[1])
# x <- ncvar_get(ncxy, varid="lon") 
# y <- ncvar_get(ncxy, varid="lat") 
# nc_close(ncxy)
# 
# 
# # get mean info TSM
# nctsm_spring <- nc_open(filename=TSMfile_spring)
# nctsm_summer <- nc_open(filename=TSMfile_summer)
# nctsm_fall <- nc_open(filename=TSMfile_fall)
# #extract
# tsm_mean_spring_masked <- ncvar_get(nctsm_spring, varid="tsm_mean_masked") 
# tsm_mean_summer_masked <- ncvar_get(nctsm_summer, varid="tsm_mean_masked")
# tsm_mean_fall_masked <- ncvar_get(nctsm_fall, varid="tsm_mean_masked")
# 
# nc_close(nctsm_spring)
# nc_close(nctsm_summer)
# nc_close(nctsm_fall)
# 
# # get mean info chla
# ncchla_spring <- nc_open(filename=chlafile_spring)
# ncchla_summer <- nc_open(filename=chlafile_summer)
# ncchla_fall <- nc_open(filename=chlafile_fall)
# #extract
# chla_mean_spring_masked <- ncvar_get(ncchla_spring, varid="chla_mean_masked") 
# chla_mean_summer_masked <- ncvar_get(ncchla_summer, varid="chla_mean_masked")
# chla_mean_fall_masked <- ncvar_get(ncchla_fall, varid="chla_mean_masked")
# 
# nc_close(ncchla_spring)
# nc_close(ncchla_summer)
# nc_close(ncchla_fall)
# 
# # get mean info CDOM
# ncCDOM_spring <- nc_open(filename=CDOMfile_spring)
# ncCDOM_summer <- nc_open(filename=CDOMfile_summer)
# ncCDOM_fall <- nc_open(filename=CDOMfile_fall)
# #extract
# CDOM_mean_spring_masked <- ncvar_get(ncCDOM_spring, varid="adg_mean_masked") 
# CDOM_mean_summer_masked <- ncvar_get(ncCDOM_summer, varid="adg_mean_masked")
# CDOM_mean_fall_masked <- ncvar_get(ncCDOM_fall, varid="adg_mean_masked")
# 
# nc_close(ncCDOM_spring)
# nc_close(ncCDOM_summer)
# nc_close(ncCDOM_fall)
# 
# #create spatial dataset from longitude and latitude
# df <- data.frame(longitude = c(x), latitude = c(y), TSM_spring = c(tsm_mean_spring_masked), TSM_summer = c(tsm_mean_spring_masked), TSM_fall = c(tsm_mean_fall_masked), chla_spring = c(chla_mean_spring_masked), chla_summer = c(chla_mean_summer_masked), chla_fall = c(chla_mean_fall_masked), CDOM_spring = c(CDOM_mean_spring_masked), CDOM_summer = c(CDOM_mean_summer_masked), CDOM_fall = c(CDOM_mean_fall_masked)) 
# 
# #make spatial points dataframe
# # projections
# latlon <- "+proj=longlat +datum=WGS84"
# #bcalbers <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
# 
# spdat <- df
# coordinates(spdat) <- ~longitude+latitude
# proj4string(spdat) <- CRS(latlon)
# 
# spdat <- spTransform(spdat, geoCRS)
# 
# # Clean up workspace to free up space
# toRem <- c("eHG", "eNCC", "eQCS", "eSASIMA", "eWCVI", "ncCDOM_fall", "ncCDOM_spring", "ncCDOM_summer", 
#            "ncchla_fall", "ncchla_summer", "ncchla_spring", 
#            "nctsm_fall", "nctsm_summer", "nctsm_spring", "CDOM_mean_fall_masked", "CDOM_mean_spring_masked", "CDOM_mean_summer_masked",
#            "tsm_mean_fall_masked", "tsm_mean_spring_masked", "tsm_mean_summer_masked", "chla_mean_fall_masked", "chla_mean_spring_masked", 
#            "chla_mean_summer_masked")
# toRem <- toRem[toRem %in% ls()]
# rm( list = toRem )
# # Garbage Collection
# gc(reset=TRUE)
# 
# 
# 
# #interpolate turbidity onto sdm grid using thin plate spline model####
# xy <- data.frame(coordinates(spdat))
# attr(xy, "dimnames") <- NULL
# 
# #interpolate chala values
# chla_spring_val <- df[["chla_spring"]]
# tps_chla_spring <- fields::fastTps(xy, chla_spring_val, theta=20000)
# s_chla_spring <- raster::interpolate(ras, tps_chla_spring)
# r_chla_spring <- mask(s_chla_spring, ras)
# names(r_chla_spring) <- "chla_mean_spring"
# out_raster <- paste0("data/chla/Rasters/chla_mean_spring", ".tif")
# writeRaster(rU, out_raster, overwrite=TRUE)
# 
# 
# 
# 
# 
# 
# # remove first and last row and column to match with u and v at rho points
# if ( variable == "curr"){
#   deplev <- deplev[-1,-1,]
#   deplev <- deplev[-nrow(deplev),-ncol(deplev),]
# }
# 
# # list all nc files and combine with season names
# ncfiles <- list.files("ave", pattern = "o_aver", full.names = TRUE)
# 
# # get x,y data
# ncxy <- nc_open(filename=ncfiles[1])
# x <- ncvar_get(ncxy, varid="lon_rho") 
# y <- ncvar_get(ncxy, varid="lat_rho") 
# nc_close(ncxy)
# # remove first and last row and column to match with u and v at rho points
# if ( variable == "curr"){
#   x <- x[-1,-1] 
#   x <- x[-nrow(x),-ncol(x)]
#   y <- y[-1,-1] 
#   y <- y[-nrow(y),-ncol(y)]
# }
# 
# 
# #------------------------------------------------------------------------------#
# # functions
# 
# # open nc, get variable, unstagger for each time step
# u2rho <- function(ncname,dt){
#   require(ncdf4)
#   require(abind)
#   # open nc
#   nc <- nc_open(filename=ncname)
#   # get nc variables, count = x,y,z,t
#   u <- ncvar_get(nc, varid="u", start=c(1,1,1,dt), count=c(235,410,30,1))
#   # close nc
#   nc_close(nc)
#   # Linear interpolation to unstagger u to match to center (rho) lat, long coordinates.
#   ua <- abind(u[-1,,],u[-nrow(u),,], along=4)
#   u_rho <- apply(ua , 1:3, mean, na.rm=T )
#   u_rho <- u_rho[,-1,] # remove first column
#   u_rho <- u_rho[,-ncol(u_rho),] # remove last column
#   # return
#   u_rho
# }
# 
# # open nc, get variables, unstagger for each year file
# v2rho <- function(ncname, dt){
#   require(ncdf4)
#   require(abind)
#   # open nc
#   nc <- nc_open(filename=ncname)
#   # get nc variables, count = x,y,z,t
#   v <- ncvar_get(nc, varid="v", start=c(1,1,1,dt), count=c(236,409,30,1))
#   # close nc
#   nc_close(nc)
#   # Linear interpolation to unstagger v to match to center (rho) lat, long coordinates.
#   va <- abind(v[,-1,],v[,-ncol(v),], along=4)
#   v_rho <- apply(va , 1:3, mean, na.rm=T )
#   v_rho <- v_rho[-1,,] # remove first row
#   v_rho <- v_rho[-nrow(v_rho),,] # remove last row
#   # return
#   v_rho
# }
# 
# # Get variable for each 15 day period
# getVar <- function(ncname, dt){
#   if ( variable == "curr") {
#     # unstagger functions
#     uvar <- u2rho( dt=dt, ncname=ncname)
#     vvar <- v2rho( dt=dt, ncname=ncname)
#     # Calculate speed from u and v
#     var <- sqrt( uvar^2 + vvar^2 )
#   } else {
#     require(ncdf4)
#     # open nc
#     nc <- nc_open(filename=ncname)
#     # get nc variables, count = x,y,z,t
#     var <- ncvar_get(nc, varid=variable, start=c(1,1,1,dt), count=c(236,410,30,1))
#     # close nc
#     nc_close(nc)
#   }
#   # return
#   return(var)
# }
# 
# # Calculate mean for each z level and seperate into depth bins
# getDepth <- function(ncfile, dtlevel){
#   # require
#   require(reshape2)
#   # calculate variables for each z level (1 to 30)
#   varz <- getVar( ncname=ncfile, dt=dtlevel )
#   # group levels into depth bins
#   bs <- deplev; bs[,,30] <- 1; bs[,,1:29] <- 0
#   b10 <- deplev > -10; b10[b10 == FALSE] <- 0; b10[b10 == TRUE] <- 1
#   b25 <- deplev < -10 & deplev > -25
#   b25[b25 == FALSE] <- 0; b25[b25 == TRUE] <- 1
#   b50 <- deplev < -25 & deplev > -50
#   b50[b50 == FALSE] <- 0; b50[b50 == TRUE] <- 1
#   b75 <- deplev < -50 & deplev > -75
#   b75[b75 == FALSE] <- 0; b75[b75 == TRUE] <- 1
#   b100 <- deplev < -75 & deplev > -100
#   b100[b100 == FALSE] <- 0; b100[b100 == TRUE] <- 1
#   b150 <- deplev < -100 & deplev > -150
#   b150[b150 == FALSE] <- 0; b150[b150 == TRUE] <- 1
#   b200 <- deplev < -150 & deplev > -200
#   b200[b200 == FALSE] <- 0; b200[b200 == TRUE] <- 1
#   b250 <- deplev < -200 & deplev > -250
#   b250[b250 == FALSE] <- 0; b250[b250 == TRUE] <- 1
#   b250tb <- deplev < -250; b250tb[b250tb == FALSE] <- 0
#   b250tb[b250tb == TRUE] <- 1
#   bb <- deplev; bb[,,1] <- 1; bb[,,2:30] <- 0
#   # depth bins
#   bins <- c("s","10","25","50","75","100","150","200","250","250tb","b")
#   # get values for depth bins
#   bybin <- function(bin, mat){
#     # Multiply variable matrix by depth bin matrix
#     d <- mat * get(paste0("b",bin))
#     # convert zeros to NA before averaging within depthbins
#     d[d == 0] <- NA
#     var <- apply(d, c(1,2), mean, na.rm=T)
#     # merge with lat, lon data 
#     varz <- data.frame(Latitude=c(y),Longitude=c(x),var=c(var))
#   }
#   varbins <- lapply(bins, FUN=bybin, mat=varz)
#   names(varbins) <- bins
#   # melt
#   mdat <- melt(varbins, id.vars = c("Latitude","Longitude","var"))
#   colnames(mdat)[4] <- "Depthbin" 
#   # return
#   return(mdat)
# }
# 
# # Run bydepth function for each timestep in parallel
# getdt <- function( ncfile ) {
#   # require
#   require(reshape2)
#   require(parallel)
#   ## create cluster object
#   cl <- makeCluster(6)
#   ## make variables and packages available to cluster
#   clusterExport(cl, varlist=c( "getDepth","getVar","v2rho","u2rho",
#                                "deplev","type","variable", "x","y"))
#   # run through all 24 timesteps in each yearly nc file
#   vardt <- parLapply( cl, 1:24, fun=getDepth, ncfile=ncfile )
#   ## stop cluster
#   stopCluster(cl)
#   # Melt
#   mdt <- melt(vardt, id.vars = c("Latitude","Longitude","var", "Depthbin"))
#   colnames(mdt)[5] <- "Timestep"
#   # remove na
#   mdt <- mdt[complete.cases(mdt),]
#   # return 
#   return( mdt )
# }
# 
# 
# #------------------------------------------------------------------------------#
# # run functions
# 
# # start time counter
# start <- Sys.time()
# # Run through all nc files
# ## create cluster object
# clyr <- makeCluster(2)
# ## make variables and packages available to cluster
# clusterExport(clyr, varlist=c( "getdt", "getDepth","getVar","v2rho","u2rho",
#                                "deplev","type","variable", "x","y"))
# # run getdt
# varall <- parLapply( clyr, ncfiles, fun=getdt )
# ## stop cluster
# stopCluster(clyr)
# # Melt
# names(varall) <- 1998:2007
# mall <- melt(varall, id.vars = c("Latitude","Longitude","var", "Depthbin", "Timestep"))
# colnames(mall)[6] <- "Year"
# 
# 
# #------------------------------------------------------------------------------#
# # Calculate min, max, and seasonal means    
# 
# # add season attribute
# mall$Season <- "win"
# mall$Season[mall$Timestep %in% 7:18] <- "sum"
# 
# # get mean values for each season
# sdat <- mall %>% group_by( Latitude, Longitude, Depthbin, Season ) %>%
#   summarise( mean = mean(var) ) %>%
#   as.data.frame()
# 
# mmdat <- mall %>% group_by( Latitude, Longitude, Depthbin ) %>%
#   summarise( min = min(var),
#              max = max(var) ) %>%
#   as.data.frame()
# 
# # combine
# dat <- data.frame( mmdat, 
#                    meanSummer=sdat$mean[sdat$Season == "sum"], 
#                    meanWinter=sdat$mean[sdat$Season == "win"] )
# 
# # Remove zeros
# # Zeros represent areas where the depthbin is deeper than the botton depth
# dat <- dat[!dat$max == 0, ]
# 
# #---------------------------------------------------------------------------------#
# # Export data
# 
# # save as csv
# write.csv(dat, 
#           file=paste0("Data/",variable,"/",type,"_",variable,".csv"), 
#           row.names = FALSE)
# 
# # remove error in salt dataset
# # really low value of one record should be removed
# if( variable == "salt" ) {
#   dat <- dat[!dat$min < 1,]
# }
# 
# # projections
# latlon <- "+proj=longlat +datum=WGS84"
# bcalbers <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
# 
# # export shapefile for each depth range
# # depth bins
# bins <- c("s","10","25","50","75","100","150","200","250","250tb","b")
# for(b in bins){
#   
#   # subset data
#   spdat <- dat[dat$Depthbin == b,]
#   
#   # convert to sp
#   coordinates(spdat) <- ~Longitude+Latitude
#   proj4string(spdat) <- CRS(latlon)
#   spdat <- spTransform(spdat, bcalbers)
#   
#   # save
#   writeOGR(spdat,dsn=paste0("Data/",variable,"/Shapefiles"),
#            layer=paste0(type,"_",variable,"_",b),
#            driver = "ESRI Shapefile", overwrite_layer = TRUE)
# }
# 
# 
# #end time counter
# end <- Sys.time()
# end - start   
# 
# #---------------------------------------------------------------------------------#
# # clean up
# rm(list=ls())
# gc(reset=T)
