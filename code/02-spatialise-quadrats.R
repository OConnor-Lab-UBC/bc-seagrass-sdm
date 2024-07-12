###############################################################################
#
# Authors:      Ashley Park
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



#when spatializing need to decide what to do with aggregating. This code provides out of quadrats aggregated to 20 m and notaggregated
#First test sdmTMB with data not aggregating

###############################################################################


# Load packages
library(sf)
library(tidyverse)
library(terra)
library(matrixStats)
library(geosphere)
library(Hmisc)

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

#find transects that only have one xy coodinate per transect
missing_startorend <- which( !complete.cases(dat[, c("LonShallow","LatShallow","LonDeep","LatDeep")]) )
single_xy <- dat[c(missing_startorend),]
unique(single_xy$Source) #RSU, GSU, Cuke, GDK, and BHM surveys

#find transects that have same lat and lon 
same_startend <- which( paste0(dat$LonDeep, dat$LatDeep) == paste0(dat$LonShallow, dat$LatShallow))
same_xy <- dat[c(same_startend),]

#make cliffs dataset
cliffs <- dat[c(missing_startorend, same_startend),]
cliffs <- cliffs%>% filter(Transect_length < 41) #this is the cliffs data set to be spatialized separately and retain all quadrats
unique(cliffs$HKey) #1497 transects are cliffs

#clean up dat
dat_for_trans <- dat[-c(missing_startorend,same_startend, missing_both),] #remove those from dat, took out  becasue we spatialize some of these. Need to confirm this right to do

single_xy <- single_xy%>% filter(Transect_length > 40) #remove cliffs from single_xy
unique(single_xy$HKey) #2006 transects have single xy
same_xy <- same_xy%>% filter(Transect_length > 40) #remove cliffs from same_xy 
unique(same_xy$HKey) #58 transects have same xy

#likely can only spatialize if have deep end coordinate and transect length is <200m (otherwise other points of land may confuse spatialization)
single_xy_tospatialize <- single_xy %>% 
  filter(Transect_length < 201,
         !is.na(LatDeep),
         !is.na(LonDeep))
unique(single_xy_tospatialize$HKey) #1666 transects have single xy that can be spatialized

#get seperate dataframe of where we can extract some quadrats at the one point
single_xy_topoint <- single_xy %>% 
  filter(Transect_length > 200|Transect_length < 201 & is.na(LatDeep)) %>%
  rbind(same_xy) # add in same_xy transects
 unique(single_xy_topoint$HKey) #396 transects have single xy we can get quadrats from one point

# Create transect dataset with filtered dat dataset
trans <- dat_for_trans[!duplicated(dat_for_trans$HKey), 
             c("Survey","Source", "Year","Month","Day","HKey",
               "LonShallow","LatShallow","LonDeep","LatDeep")]

# Create transect dataset with filtered single_xy to spatialize dataset
single_xy_trans <- single_xy_tospatialize[!duplicated(single_xy_tospatialize$HKey),
                             c("Survey","Source","Year","Month","Day","HKey",
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
      dplyr::summarize() %>%
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
  pts <- st_line_sample(sl[i,], n = npts, type = "regular") 
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


##############----------------------------------------------------------------------------#
# Add points with only a single x,y at end to spatialize
# Create featureID and re-assign row names

singlexy.sf <- single_xy_trans %>% st_as_sf(coords = c("LonDeep", "LatDeep"), crs = "EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

# look at https://stackoverflow.com/questions/51292952/snap-a-point-to-the-closest-point-on-a-line-segment-using-sf

# create spatial lines and remove transects likely not correct due to incorrect xy by distance
sl2<- singlexy.sf %>% 
  mutate(
    Lines = st_nearest_points(geometry, boundary[st_nearest_feature(geometry, boundary, pairwise = TRUE),], pairwise = TRUE),
    closest_point = st_cast(Lines, 'POINT')[seq(2, nrow(.)*2, 2)],
    length = as.integer(st_length(Lines)))

#filter transect database by distance
single_xy_trans$length <- sl2$length

single_xy_trans <- single_xy_trans %>%
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

#filter spatial lines dataset
sl2<- sl2 %>%  
  mutate(Keep = case_when(length<800 & Source=="RSU_bio" ~ 'Yes',
                          length<800 & Source=="GSU_bio" ~ 'Yes',
                          length<400 & Source=="Cuke_bio" ~ 'Yes',
                          length<2000 & Source=="GDK_bio" ~ 'Yes',
                          length<350 & Source=="BHM" ~ 'Yes',
                          length<350 & Source=="MSEA"~ 'Yes',
                          length<150 & Source=="Multispecies_bio" ~ 'Yes',
                          TRUE ~ 'No')) %>%
  filter(Keep == "Yes") %>%
  select(HKey, Lines) %>%
  st_drop_geometry() %>%
  st_as_sf(crs = "EPSG:3005")

single_xy_trans$featureID <- 1:nrow(single_xy_trans)
row.names(single_xy_trans) <- single_xy_trans$featureID
  
# Generate points along the transect at a certain distance set by 'dist_points'
# Generate points
spdf.list2 <- list()
for( i in 1:nrow(sl2)){
  # Get Hkey
  hkey <- single_xy_trans$HKey[single_xy_trans$featureID == i]
  # Get length of line
  linelength <- st_length(sl2[i,])
  # Number of points for breaks at 'dist_points' distance
  npts <- ceiling( linelength / 20 ) # 20 becasue 20 m bathy layer
  # Set seed for consistent results
  set.seed(42)
  # Regularly sample points along line
  pts <- st_line_sample(sl2[i,], n = npts, type = "regular") 
  ckey <- as.character(hkey)
  spdf.list2[[ckey]] <-
    data.frame(st_coordinates(st_as_sf(pts, crs =	"EPSG:3005")),
               featureID = rep(i, npts),
               HKey = rep(hkey, npts)) 
}

# bind all points together
spdf2 <- do.call("rbind", spdf.list2)
spdf2<- spdf2 %>%
  st_as_sf(coords = c("X", "Y"), crs = 	"EPSG:3005")

# extract depth from bathy raster
extractdepth3 <- terra::extract(x = bathy, y = spdf2)
names(extractdepth3) <-c("ID", "bathy")

# combine with spdf
ptsdat2<- cbind(st_drop_geometry(spdf2), st_coordinates(st_as_sf(spdf2, coords = c("x", "y"), crs = "EPSG:3005")), bathy=extractdepth3$bathy)



#----------------------------------------------------------------------------#
# Extract depth from bathymetry raster at spdf points

# extract depth from bathy raster
extractdepth <- terra::extract(x = bathy, y = spdf)
names(extractdepth) <-c("ID", "bathy")

# combine with spdf
ptsdat<- cbind(st_drop_geometry(spdf), st_coordinates(st_as_sf(spdf, coords = c("x", "y"), crs = "EPSG:3005")), bathy=extractdepth$bathy)

ptsdat <- ptsdat %>% rbind(ptsdat2)

#----------------------------------------------------------------------------#
# Merge points with quadrats using closest depth match
# If loop throws error increase tol in find.matches (doesn't work when set to Inf)

# Add transect data with only a single start or end x,y back to dat
#dat <- rbind(dat, single_xy)

# Remove HKeys not in ptsdat
qdat <- dat[dat$HKey %in% unique(ptsdat$HKey),]


# Match pts to quadrats by depth
matchSpatial <- function( x ){  
  # require
  require(Hmisc)
  # subset by hkey
  quad <- qdat[qdat$HKey == x,]
  pts <- ptsdat[ptsdat$HKey == x,]
  # match
  matchdepth <- find.matches(quad$CorDepthM, pts$bathy, tol=1000, maxmatch=1)
  # merge quad and pts based off matchdepth
  mdat <- data.frame( quad, pts[matchdepth$matches, c("X", "Y", "bathy")],
                      ID = paste(quad$HKey, matchdepth$matches, sep="_"))
  # Calculate difference between quadrat depth and bathy
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

# Remove quadrats where bathy is NA becasue outside extent of bathymetry
spatialised <- spatialised[ which(complete.cases(spatialised$bathy)), ]


#### cliffs dataset
# get lat lon of one point
cliffs<- cliffs %>%
  mutate(LatDeep = case_when(!is.na(LatDeep) ~ LatDeep,
                             is.na(LatDeep) ~ LatShallow),
         LonDeep = case_when(!is.na(LonDeep) ~ LonDeep,
                             is.na(LonDeep) ~ LonShallow))%>%
  filter(!is.na(LatDeep))

cliffs.sf<- cliffs %>% 
  st_as_sf(coords = c("LonDeep", "LatDeep"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

cliffs.spdf<- cbind(cliffs, st_coordinates(st_as_sf(cliffs.sf, coords = c("x", "y"), crs = "EPSG:3005")))
cliffs.spdf$bathy <- NA
cliffs.spdf$ID <- paste(cliffs.spdf$HKey, "1", sep="_")
cliffs.spdf$depthdiff <- NA

spatialised<- spatialised%>%
  rbind(cliffs.spdf) # add in cliffs to spatialized

#### sites single xy that are just one point
# get lat lon of one point
single_xy_topoint<- single_xy_topoint %>%
  mutate(LatDeep = case_when(!is.na(LatDeep) ~ LatDeep,
                             is.na(LatDeep) ~ LatShallow),
         LonDeep = case_when(!is.na(LonDeep) ~ LonDeep,
                             is.na(LonDeep) ~ LonShallow))%>%
  filter(!is.na(LatDeep))

singleptxy.sf <- single_xy_topoint %>% st_as_sf(coords = c("LonDeep", "LatDeep"), crs = 	"EPSG:4326") %>%
  st_transform(crs = "EPSG:3005")

# extract depth from bathy raster
extractdepth2 <- terra::extract(x = bathy, y = singleptxy.sf)
names(extractdepth2) <-c("ID", "bathy")

# combine with singleptxy.sf
singleptxy.spdf<- cbind(single_xy_topoint, st_coordinates(st_as_sf(singleptxy.sf, coords = c("x", "y"), crs = "EPSG:3005")), bathy=extractdepth2$bathy)

singleptxy.spdf$ID <- paste(singleptxy.spdf$HKey, "1", sep="_")
singleptxy.spdf$depthdiff <- abs( singleptxy.spdf$CorDepthM - singleptxy.spdf$bathy)

#filter out records with depth diff >5m
singleptxy.spdf <- singleptxy.spdf %>% 
  filter(depthdiff<= 5)

spatialised<- spatialised%>%
  rbind(singleptxy.spdf) # add in single xy point to spatialized

#----------------------------------------------------------------------------#
# Clean up - Remove points that are likely incorrect

# check depth diff
summary(spatialised$depthdiff)

# Combine transect and quadrat to get unique quadrat identifier
spatialised$QID <- paste(spatialised$HKey, spatialised$Quadrat, sep="_")

#not including this one for now as diver depths are more accurate than bathy raster
# Remove quadrats with depth difference > 5
# between corrected quadrat depth and depth from bathymetry raster
#spatialised <- spatialised[which(spatialised$depthdiff <= 5),]

# Remove NAs in x and y's (usually due to NAs present in CorDepthM)
spatialised <- spatialised[complete.cases(spatialised[, c("X","Y")]),]


#Mean number of quadrats aggregated into a single spatial point
nquads <- aggregate( Quadrat ~ ID, data= spatialised, function(x) length(unique(x)))
cat( round( mean(nquads$Quadrat), 1 ), "quadrats per spatial point on average",
     "(range", min(nquads$Quadrat), "to", max(nquads$Quadrat), "quadrats)","\n\n")
names(nquads)[2] <- "NumQuadrats"

# Quadrats retained
#select 
spatialised<-spatialised %>% select("Survey","Year","Month","Day","HKey","ID" , "X", "Y",
             "LonDeep","LatDeep","LonShallow","LatShallow", "CorDepthM", "Slope", "Substrate", "PH", "ZO")

write.csv( spatialised, file="code/output_data/SpatializedQuadrats_notaggregated.csv") 


# Convert to spdf and export
spatialised.spdf <- spatialised %>%
  st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 

# export as shapefile
# likely to have issues with attribute field names shortening
st_write(spatialised.spdf, "code/output_data/SpatializedQuadrats_notaggregated.shp", append=FALSE)



#----------------------------------------------------------------------------#

# Aggregate presence / absence by spatial points, want to test both if need to aggregate, look at spatial autocorrelation
# 
# # Transect attributes
# att <- csp[!duplicated(csp$ID),]
# att <- att[c("Survey","Year","Month","Day","HKey","ID" ,"X","Y",
#              "LonDeep","LatDeep","LonShallow","LatShallow")]
# # Quadrat attributes - Mean depth from quadrats aggregated to spatialised points
# mean_att <- aggregate( . ~ ID, mean, data = csp[c("ID", "CorDepthM", "bathy", "depthdiff", "Slope", "PH", "ZO")])
# names(mean_att) <- c("ID", "mean_CorDepthM", "mean_bathy", "mean_depthdiff", "mean_slope", "PH", "ZO")
# 
# #Quadrat attributes - most common substrate
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }
# 
# mode_att <- aggregate( . ~ ID, Mode, data = csp[c("ID", "Substrate")])
# 
# # Add back attributes
# att <- merge(att, mean_att, by="ID")
# spat <- merge(att, mode_att, by="ID")
# spat <- merge(spat, nquads, by="ID")
# spat <- spat[order(spat$HKey),]
# 
# # Ensure presence/absence
# spat$PH[spat$PH > 0] <- 1
# spat$ZO[spat$ZO > 0] <- 1

# # check
# cat( "\n\n")
# cat( "#----------------------------------------------------------------------#\n")
# cat( "First 5 rows and of spatialised site by species matrix:\n")
# head(spat, 5)
# cat( "\n\n")



#----------------------------------------------------------------------------#
# # export as csv
# write.csv( spat, file="code/output_data/SpatializedQuadrats_SitesvSpeciesMatrix_aggregated.csv") 
# 
# # Convert to spdf and export
# spat <- spat %>%
#   st_as_sf(coords = c("X", "Y"), crs = "EPSG:3005") 
# 
# # export as shapefile
# # likely to have issues with attribute field names shortening
# st_write(spat, "code/output_data/SpatializedQuadrats_aggregated.shp", append=FALSE)
# 
# 
