###############################################################################
#
# Authors:      Ashley Park, modified code from Jessica Nephin
# Affiliation:  Fisheries and Oceans Canada (DFO)
# Group:        Marine Spatial Ecology and Analysis
# Location:     Institute of Ocean Sciences
# Contact:      e-mail: Ashley.Park@dfo-mpo.gc.ca
#
#
# Objective:
# ---------
# Spatialise Dive Survey Observations for eelgrass and phyllospadix observations
#
#
# Overview:
# --------
# Create spatial lines from extended start and end points, using the coastline
# as a guide to extend lines nearshore. Generates points along extended transect 
# lines at a selected distance. Extracts depth from bathymetry raster at points 
# along lines. Merges points with quadrats using closest depth match. Merges 
# the resultant 'spatialised points' with species data. 
#
#
# Requirements (input data):
# -------------------------
# 1) Must have fields: Survey, Year, Month, Day, HKey, Quadrat, CorDepthM
#    - 'HKey' field is a unique identifier for every transect
#    - 'Quadrat' field is a unique identifier for each quadrat within a transect
#    - 'CorDepthM' field is the chart datum corrected quadrat depth in meters
#
# 2) Lat and long must be decimal degrees represented with fields:  
#     LonShallow, LatShallow, LonDeep, LatDeep
#
# 3) When presence/absence are represented by 0/1 values, 
#    the field containing 0/1 values must be named 'SpNum'


#when spatializing need to remember to not aggregate by years
#what to do with substrate and slope when aggregating
#there are a lot of transects from earlier surveys with no lat/long shallow, need to figure out how to make line without having datapoint
###############################################################################


# Load packages
library(sf)
library(tidyverse)
library(terra)
library(matrixStats)
# library(raster)
library(geosphere)
library(Hmisc)
# library(reshape2)

#load_data####
load("code/output_data/seagrass_data.RData")

coastline_full <- st_read("raw_data/CHS_HWL2015_Coastline.gdb", layer = "Line_CHS_Pacific_HWL_2015_5028437")
boundary <- coastline_full %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_crop(st_bbox(coastline_full) + c(0, 100000, -10000, -80000)) %>%
  st_transform(crs = "EPSG:3005") 

# Load bathy raster
bathy <- vrt("raw_data/Bathymetry/Coastwide_20m_mosaic.tif")

#----------------------------------------------------------------------------#
# Remove quadrats with depth difference > depth_dist
# between corrected quadrat depth and depth from bathymetry raster
depth_dist <- 5

# How are the species data represented?
# If presence/absence is coded using species codes, SpNum == FALSE
# If presence/absence is coded using a numeric column (0/1), SpNum == TRUE
# Name of column in dataset with 0/1 values must be named: 'SpNum'
SpNum <- FALSE

# Name of column with species names
spNames <- "Species"

# Species of interest
# if all species are needed = 'all'
spp <- 'all'


#----------------------------------------------------------------------------#
# Clean input data
# Find missing start and end latitude and longitude values

#Find number of quadrats per HKey
qCnt <- dat %>% 
  group_by(HKey) %>%
  dplyr::summarize( qCnt = length(Quadrat))

# Add number of quadrats to df
dat <- left_join(dat, qCnt, by="HKey")

# Find transects which don't have x,y values for either the start or end and remove
missing_both <- which( apply(is.na(dat[c("LonShallow","LatShallow","LonDeep","LatDeep")]), 1, all) )
no_xy <- dat[missing_both,]
dat <- dat[-c(missing_both),]

#find transects that have same lat and lon and identify which ones may be cliff sites
same_startend <- which( paste0(dat$LonDeep, dat$LatDeep) == paste0(dat$LonShallow, dat$LatShallow) )

# Find transects which only have a single x,y value and they are short transects, so likely cliffs
missing_startorend_cliffs <- which( !complete.cases(dat[, c("LonShallow","LatShallow","LonDeep","LatDeep")]) & dat$Transect_length < 40)
cliffs <- dat[c(same_startend, missing_startorend_cliffs),]
dat <- dat[-c(same_startend, missing_startorend_cliffs),] #remove all potential cliffs from database
cliffs <- cliffs%>% filter(Transect_length < 41) #this is the cliffs data set to be spatialized separately

#find remaining transects that only have one xy coodinate per transect
missing_startorend <- which( !complete.cases(dat[, c("LonShallow","LatShallow","LonDeep","LatDeep")]) )
single_xy <- dat[c(missing_startorend),]
unique(single_xy$Source) #RSU, GSU, Cuke and BHM surveys
dat <- dat[-c(missing_startorend),] #remove those from dat

#likely can only spatialize if have deep end coordinate and transect length is <200m (otherwise other points of land may confuse spatialization)
single_xy <- single_xy %>% 
  filter(Transect_length < 201,
         !is.na(LatDeep),
         !is.na(LonDeep))
### for a later date to figure out if possible to spatialize these!



# Create transect dataset with filtered dat dataset
trans <- dat[!duplicated(dat$HKey), 
             c("Survey","Source", "Year","Month","Day","HKey",
               "LonShallow","LatShallow","LonDeep","LatDeep")]


## add this in the future but not just for point
single_xy_trans <- single_xy[!duplicated(single_xy$HKey),
                             c("Survey","Year","Month","Day","HKey",
                               "LonShallow","LatShallow","LonDeep","LatDeep")]


#----------------------------------------------------------------------------#
# Create spatial lines from shallow to deep transect points
# extended by 'dist_extend' from the deep point
# extended by distance to coastline or 'threshold' from shallow point, threshold decided to be 50m
# (shallow point tends not to reach shoreline)

# For transects with a unique shallow and deep position,
# Calculate distance from shallow point to deep point
shallow.sf <- trans %>% st_as_sf(coords = c("LonShallow", "LatShallow"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")
deep.sf <- trans %>% st_as_sf(coords = c("LonDeep", "LatDeep"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

s2e <- st_distance(shallow.sf, deep.sf, by_element = TRUE)
trans$length <- as.integer(s2e)

# Remove transects with start and end points greater than this distance apart.
# dist_transect apart because either the deep or shallow position is likely wrong
trans <- trans %>%
  mutate(Keep = case_when(length<800 & Source=="RSU_bio" ~ 'Yes',
                          length<800 & Source=="GSU_bio" ~ 'Yes',
                          length<400 & Source=="Cuke_bio" ~ 'Yes',
                          length<2000 & Source=="GDK_bio" ~ 'Yes',
                          length<350 & Source=="BHM" ~ 'Yes',
                          length<350 & Source=="MSEA"~ 'Yes',
                          length<150 & Source=="Multispecies_bio" ~ 'Yes',
                          TRUE ~ 'No')) %>%
  filter(Keep == "Yes") %>%
  select(-Keep)


#Convert shallow points to spatial class with dataset now filtered to remove transects with likely incorrect xy
shallow.sf<- trans %>% st_as_sf(coords = c("LonShallow", "LatShallow"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

nearest <- st_distance(shallow.sf, boundary)
nearest_min <- matrixStats::rowMins(nearest)
trans$distance <- as.integer(nearest_min)

# If distance to coastline is greater than threshold, use threshold, which is 50m
trans$distance[ trans$distance > 50 ] <- 50

# Create featureID and re-assign row names
trans$featureID <- 1:nrow(trans)
row.names(trans) <- trans$featureID

# Distance to extend deep transect point toward deeper depths
dist_extend <- 50

# convert to lines extended past points along transect bearing
lines.list <- list()
for (i in 1:nrow(trans)){
  # get deep, shallow and distance from shallow point to coastline
  deep <- trans[i, c("LonDeep","LatDeep")]
  shallow <- trans[i, c("LonShallow","LatShallow")]
  d <- trans[i, "distance"]
  # if deep and shallow points differ create lines
  if ( ! all(deep == shallow) ){
    # calculate bearing from deep to shallow
    b <- finalBearing(deep, shallow)
    # extend point from shallow (nearshore)
    e.shallow <- destPoint(p=shallow, b=b, d=d) # extend by the nearest distance to shore
    # extend point from deep
    e.deep <- destPoint(p=deep, b=b-180, d=dist_extend) # extend 50m past deep point
    # create line from extended points
    l <- rbind( as.matrix( e.deep ),
                as.matrix( e.shallow ))
    l <- l %>% as.data.frame %>% 
      st_as_sf(coords = c("lon","lat"), crs = "EPSG:4326") %>%
      st_transform(crs = "EPSG:3005") %>%
      group_by() %>% 
      summarize() %>%
      st_cast("LINESTRING")
    # hkey
    hkey <- as.character( trans[i, "HKey"] )
    # add to list
    lines.list[[hkey]] <- st_as_sfc(l)
  }
}

# create spatial lines and spatial lines dataframe
sl<- lines.list %>%
  bind_rows() %>%
  pivot_longer(cols = everything(), names_to = as.character("HKey"), values_to = "Lines") %>%
  st_as_sf(crs = "EPSG:3005")

sldf<- sl %>%
  merge(trans)


#----------------------------------------------------------------------------#
# Generate points along the transect at a certain distance set by 'dist_points'

# Generate points
spdf.list <- list()
for( i in 1:nrow(sl)){
  # Get Hkey
  hkey <- trans$HKey[trans$featureID == i]
  # Get length of line
  linelength <- st_length(sl[i,])
  # Number of points for breaks at 'dist_points' distance
  npts <- ceiling( linelength / 20 ) # 20 becasue 20 m bathy layer
  # Set seed for consistent results
  set.seed(42)
  # Regularly sample points along line
  pts <- st_line_sample(sl[i,], n = npts, type = "regular") # 20 becasue 20 m bathy layer
  ckey <- as.character(hkey)
  spdf.list[[ckey]] <-
    data.frame(st_coordinates(st_as_sf(pts, crs =	"EPSG:3005")),
               featureID = rep(i, npts),
               HKey = rep(hkey, npts)) 
}

# bind all points together
spdf <- do.call("rbind", spdf.list)
spdf<- spdf %>%
  st_as_sf(coords = c("X", "Y"), crs = 	"EPSG:3005")


#----------------------------------------------------------------------------#
# Add points with only a single start or end x,y
# 
# if( nrow(single_xy_trans) > 0 ){
#   # Create featureID and re-assgin row names
#   single_xy_trans$featureID <- nrow(trans)+1:nrow(single_xy_trans)
#   row.names(single_xy_trans) <- single_xy_trans$featureID
#   
#   # Use start x,y unless start is NA, then use end
#   single_xy_trans$Lat <- ifelse( is.na(single_xy_trans$LatDeep), 
#                                  single_xy_trans$LatShallow,
#                                  single_xy_trans$LatDeep )
#   single_xy_trans$Long <- ifelse( is.na(single_xy_trans$LonDeep), 
#                                   single_xy_trans$LonShallow,
#                                   single_xy_trans$LonDeep )
#   
#   # convert to SpatialPointsDataFrame
#   single_spdf <- SpatialPointsDataFrame( coords = coordinates(single_xy_trans[c("Long","Lat")]),
#                                          data = single_xy_trans[c("featureID", "HKey")] ) 
#   # project to albers
#   proj4string(single_spdf) <- CRS("+proj=longlat")
#   single_spdf <- spTransform(single_spdf, geoCRS)
#   
#   # bind single x,y values with generated points along transects
#   spdf <- bind(spdf,single_spdf)
# }
# 

#----------------------------------------------------------------------------#
# Extract depth from bathymetry raster at spdf points

# extract depth from bathy raster
extractdepth <- terra::extract(x = bathy, y = spdf)
names(extractdepth) <-c("ID", "bathy")

# combine with spdf
ptsdat<- cbind(st_drop_geometry(spdf), st_coordinates(st_as_sf(spdf, coords = c("x", "y"), crs = "EPSG:3005")), bathy=extractdepth$bathy)



#----------------------------------------------------------------------------#
# Merge points with quadrats using closest depth match
# If loop throws error increase tol in find.matches (doesn't work when set to Inf)

# Add transect data with only a single start or end x,y back to dat
#dat <- rbind(dat, single_xy)

# Remove HKeys not in ptsdat
qdat <- dat[dat$HKey %in% unique(ptsdat$HKey),]


###NEED TO START EDITING FROM THIS POINT

# Match pts to quadrats by depth
matchSpatial <- function( x ){  
  # require
  require(Hmisc)
  # subset by hkey
  quad <- qdat[qdat$HKey == x,]
  pts <- ptsdat[ptsdat$HKey == x,]
  # match
  matchdepth <- find.matches(quad$CorDepthM, pts$bathy, tol=1000, maxmatch=1)
  # merge quad and pts based of matchdepth
  mdat <- data.frame( quad, pts[matchdepth$matches, c("x", "y", "bathy")],
                      ID = paste(quad$HKey, matchdepth$matches, sep="_"))
  # Calulcate difference between quadrat depth and bathy
  mdat$depthdiff <- abs( mdat$CorDepthM - mdat$bathy )
  # return
  return( mdat )
}
# Run getDist in parallel
cl <- parallel::makeCluster( parallel::detectCores() - 1 )
## make variables available to cluster
parallel::clusterExport( cl, varlist=c("qdat","ptsdat") )
# run on cluster
spatialised.list <- parallel::parLapply( cl, unique(ptsdat$HKey), matchSpatial )
## stop cluster
parallel::stopCluster( cl )

# bind data.frames together
spatialised <- do.call("rbind", spatialised.list)



#----------------------------------------------------------------------------#
# Clean up - Remove points that are likely incorrect

# check depth diff
cat( "\n")
cat( "#----------------------------------------------------------------------#\n")
cat( "Summary of depth difference between dive and bathy depths:","\n")
summary(spatialised$depthdiff)

# Combine transect and quadrat to get unique quadrat identifier
spatialised$QID <- paste(spatialised$HKey, spatialised$Quadrat, sep="_")

# Remove quadrats where bathy is NA
removed_quadrats <- NULL
fsp <- spatialised[ which(complete.cases(spatialised$bathy)), ]
rem_outside <- spatialised[ which(is.na(spatialised$bathy)), ]
if( nrow(rem_outside) > 0 ) rem_outside$reason <- "Outside bathy extent"
removed_quadrats <- rbind(removed_quadrats, rem_outside)
# Number of quadrats removed
cat( "\n")
cat( length(unique(removed_quadrats$QID)), 
     "quadrats removed because outside of bathymetry extent","\n\n")


# Remove quadrats with depth difference > depth_dist
# between corrected quadrat depth and depth from bathymetry raster
dsp <- fsp[which(fsp$depthdiff <= depth_dist),]
rem_depthdiff <- fsp[which(fsp$depthdiff > depth_dist),]
if( nrow(rem_depthdiff) > 0 ) {
  rem_depthdiff$reason <- paste("Depth difference greater than",depth_dist, "m")
}
removed_quadrats <- rbind(removed_quadrats, rem_depthdiff)
# Number of quadrats removed
cat( paste0( length(unique(rem_depthdiff$QID)), 
             " quadrats removed because depth mismatch > ",depth_dist, "m (", 
             round(( length(unique(rem_depthdiff$QID)) / 
                       length(unique(fsp$QID)) )*100), "%)"),"\n\n" )


# Remove NAs in x and y's (usually due to NAs present in CorDepthM)
rem_naxy <- dsp[!complete.cases(dsp[, c("x","y")]),]
if( nrow(rem_naxy) > 0 ) rem_naxy$reason <- "NA values in x,y points along transects"
removed_quadrats <- rbind(removed_quadrats, rem_naxy)
csp <- dsp[complete.cases(dsp[, c("x","y")]),]
# Number of transects removed
cat( length(unique(dsp$QID))-length(unique(csp$QID)), 
     "quadrats removed because NA values in quadrat x and y","\n\n")


# Quadrats removed
write.csv( removed_quadrats, file = file.path("Data",projdir,"Output/Quadrats_Removed.csv") )


# Mean number of quadrats aggregated into a single spatial point
nquads <- aggregate( Quadrat ~ ID, data= csp, function(x) length(unique(x)))
cat( round( mean(nquads$Quadrat), 1 ), "quadrats per spatial point on average", 
     "(range", min(nquads$Quadrat), "to", max(nquads$Quadrat), "quadrats)","\n\n")
names(nquads)[2] <- "NumQuadrats"

# Quadrats retained
write.csv( csp, 
           file = file.path("Data",projdir, 
                            paste0("Output/SpatializedQuadrats_MeltedDataframe_",
                                   filename,".csv")) )



#----------------------------------------------------------------------------#
# Aggregate presence / absence by spatial points
# Create species by site matrix to fill in absences

# Remove empty records with no species names
names(csp)[names(csp) == spNames] <- "Species"
csp <- csp[csp$Species != "",]

# Set value variables for dcast()
if( SpNum ) {
  value.var <- "SpNum"
} else {
  csp$presence <- 1
  value.var <- "presence"
}

# Create to species by site matrix where site is a spatial point
spcast <- dcast(ID ~ Species, value.var = value.var, fun = mean, data=csp)
spcast[is.na(spcast)] <- 0

# Subset just species of interest
if ( spp == "all" ){
  spint <- spcast
} else {
  spint <- spcast[, c("ID", spp)]
}

# Species prevalence
if ( length(spp) > 1 | spp == "all" ){
  prev <- round( apply(spint[,-1], 2, function(x) 100*(length(x[x > 0])/length(x)) ),2)
} else {
  prev  <- round( 100*(length(spint[,-1][spint[,-1] > 0])/length(spint[,-1])) ,2)
}
cat( "\n\n")
cat( "#----------------------------------------------------------------------#\n")
cat( "Species prevalence:\n")
prev

# Merge with nquads
spint <- merge( nquads, spint,  by= "ID" )

# Transect attributes
att <- csp[!duplicated(csp$ID),]
att <- att[c("Survey","Year","Month","Day","HKey","ID" ,"x","y",
             "LonDeep","LatDeep","LonShallow","LatShallow")]
# Quadrat attributes - Mean depth from quadrats aggregated to spatialised points
mean_att <- aggregate( . ~ ID, mean, data = csp[c("ID", "CorDepthM", "bathy", "depthdiff")])
names(mean_att) <- c("ID", "mean_CorDepthM", "mean_bathy", "mean_depthdiff")

# Add back attributes
att <- merge(att, mean_att, by="ID")
spat <- merge(att, spint, by="ID")
spat <- spat[order(spat$HKey),]

# Ensure presence/absence
if( SpNum ) {
  spat[spp][spat[spp] > 0] <- 1
} 

# check
cat( "\n\n")
cat( "#----------------------------------------------------------------------#\n")
cat( "First 5 rows and of spatialised site by species matrix:\n")
head(spat, 5)
cat( "\n\n")



#----------------------------------------------------------------------------#
# Convert to spdf and export

# export as csv
write.csv( spat,
           file = file.path("Data",projdir,
                            paste0("Output/SpatializedQuadrats_SitesvSpeciesMatrix_",
                                   filename,".csv")) )

# as spdf
coordinates(spat) <- ~x+y
proj4string(spat) <- geoCRS

# export as shapefile
# likely to have issues with attribute field names shortening
writeOGR(spat, dsn=file.path("Data",projdir,"Output"), 
         layer=paste0("SpatializedQuadrats_SitesvSpeciesMatrix_",filename),
         driver="ESRI Shapefile", overwrite_layer = T)





#----------------------------------------------------------------------------#
# Print end of file message and elapsed time
cat( "\n\n")
cat( "#----------------------------------------------------------------------#\n")
cat( "Finished: ", sep="" ); print( Sys.time( ) - sTimeBRT )

# Stop sinking output
sink( type = "message" )
sink( )
closeAllConnections()
