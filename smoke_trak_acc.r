#!/usr/bin/Rscript
# Smoke Spottr

library(tidyverse) # General data wrangling
library(sf) # Spatial
library(tmap) # Map plotting
library(tmaptools) # Provides the kernel-density smoothing function
library(geosphere) # Calculate distances/bearings on a sphere
library(raster) # Raster files
library(RPostgres)
library(DBI)
setwd("/mnt/R")
source("../config.r")

md = raster("mag/D_Grid_mf_2020.grd")



endTime = as.POSIXct(format(Sys.time(),tz="GMT"),tz="GMT")
startTime = endTime-24*60*60*7

#endTime = as.POSIXct("2018-03-18 02:07:46",tz="GMT")
db_dat = dbGetQuery(con_prod,"select * from smoke_reports")


db_raw = dplyr::select(db_dat,id,lat,lng,created_at,location_type,smoke_intensity,smell_intensity,smell_smoke,see_smoke,picture,alpha)
db_raw$smell_smoke=as.character(db_raw$smell_smoke)
db_raw$see_smoke=as.character(db_raw$see_smoke)
db_raw$picture=as.character(db_raw$picture)
db_raw = st_as_sf(db_raw,coords=c("lng","lat"),crs=4326)
write_sf(db_raw,"tiles/smoke_raw_points.geojson",delete_dsn = TRUE)


db_dat$created_at = as.POSIXct(format(db_dat$created_at,"%Y-%m-%d %H:%M:%S"),tz="GMT")
dat = db_dat %>% filter(created_at > startTime & created_at < endTime & location_type=="dist" & !is.na(alpha))
dat_local = db_dat %>% filter(created_at > startTime & created_at < endTime & location_type=="atl" & see_smoke == TRUE)



# Define array of start- and end- times for moving window analysis
start_array = seq(startTime,endTime,15*60)
end_array = start_array + 30*60

# PolygonList
polygon_list = list()

# PointLIst 
point_list = list()

# VectorList
vector_list = list()

for(i in seq_along(start_array)){
  data = filter(dat,created_at > start_array[i] & created_at < end_array[i]) %>% dplyr::select(created_at, id,lat,lng,alpha,alpha_accuracy)
  if(nrow(data)<3){next}
  data$alpha_accuracy[is.na(data$alpha_accuracy)] = 15
  data_left = data
  data_right = data
  data_left$side = "left"
  data_right$side = "right"
  data_left$alpha = data_left$alpha+data_left$alpha_accuracy*.75
  data_right$alpha = data_right$alpha-data_right$alpha_accuracy*.75
  data$side = "mid"
  data = rbind(data_left,data_right,data)
  
  # Read test data file
  #data = read_csv("smoke_trak/test_data.csv")
  
  # Map of Victoria and roads for plotting
  #vic = read_sf("E:/geodata/Aus_Coastline/australia/cstauscd_r.shp")
  #roads = read_sf("E:/geodata/globalmap2001/Vector/transportation/roadl.shp")
  
  # Convert the input data file to a spatial object with WGS84 lat/longs
  from_data = st_as_sf(data,coords=c("lng","lat"),crs=4326)
  
  # Plot these points on a map of Victoria
  #tm_shape(from_data) + tm_dots(size=0.5,col="green") + tm_shape(vic) + tm_borders() + tm_shape(roads) + tm_lines() + tm_layout(title = "Point Locations")
  
  # Make a copy of the input dataset - we will calculate the end-points of the vectors
  to_points = data
  
  ex_mat = st_coordinates(from_data)
  mDelta = extract(md,ex_mat)
  data$alpha = data$alpha + mDelta
  
  # This function takes lat/longs, a bearing, and a distance in m, and returns the end locations of each vector
  sp = destPoint(cbind(data$lng,data$lat),data$alpha,500)
  data$lng = sp[,1]
  data$lat = sp[,2]
  
  from_data = st_as_sf(data,coords=c("lng","lat"),crs=4326)
  
  dp = destPoint(cbind(data$lng,data$lat),data$alpha,25000)
  
  # Replace the lat/longs in this copy of input dataset with these coordinates
  to_points$lng = dp[,1]
  to_points$lat = dp[,2]
  
  # Convert these end points to a spatial object
  to_data = st_as_sf(to_points,coords=c("lng","lat"),crs=4326)
  
  # We need some unique IDs for the start and end point data frames
  # so they can be grouped together into lines
  from_data$idx = seq_along(from_data$alpha)
  to_data$idx = seq_along(to_data$alpha)
  
  # Join the two data frames into one
  comb = rbind(from_data,to_data)
  combx = filter(comb,side=="mid")
  
  # Group by ID, cast the resulting point groups as spatial lines
  comb = group_by(comb,idx) %>% summarise() %>% st_cast("LINESTRING")
  combx = group_by(combx,created_at,id) %>% summarise() %>% st_cast("LINESTRING")
  
  vector_list[[i]] = combx
  
  #sitelocs = tibble(name=c("Mt Wellington","Bridge","Droughty Pt"),
  #                  lat = c(-42.896574,-42.865191,-42.929398),
  #                  lon = c(147.237183,147.344055,147.417940))
  
  #sitelocs = st_as_sf(sitelocs,coords=c("lon","lat"),crs=4326)
  
  # Plot a map of points, vectors on Melbourne
  #tm_shape(comb) + tm_lines(col="blue") +tm_shape(from_data) + tm_dots(size=0.1,col="green") + tm_shape(vic) + tm_borders() + tm_shape(roads) + tm_lines() + tm_layout(title = "Point Locations") + tm_shape(sitelocs) + tm_symbols(shape=17)
  
  # Calculate intersection of all vectors
  line_is = st_intersection(comb,comb)
  
  # This gave us ALL intersections - both intersection points, and the lines they sit on
  # So keep only the geometries that are points
  line_is = line_is[st_geometry_type(line_is)=="POINT",]
  line_is=st_difference(line_is)
  line_is$idx = NULL
  line_is$idx.1 = NULL
  
    
  data_local = filter(dat_local,created_at > start_array[i] & created_at < end_array[i]) %>% dplyr::select(created_at, lat,lng)
  data_local = st_as_sf(data_local,coords=c("lng","lat"),crs=4326)
  data_local$type = rep("Local",nrow(data_local))
 
  line_is$created_at = rep(end_array[i],nrow(line_is))
  line_is$type = rep("Intersection",nrow(line_is))
  
  line_is = rbind(line_is,data_local)
  
  if(nrow(line_is)==0){next }
  
  point_list[[i]] = line_is
  
  if(nrow(line_is)<2){next }
  
  # Plot a map of these point intersections
  #tm_shape(comb) + tm_lines(col="blue") +tm_shape(from_data) + tm_dots(size=0.5,col="green") + tm_shape(line_is) + tm_squares(size=0.2,col="red") + tm_shape(vic) + tm_borders() + tm_shape(roads) + tm_lines() + tm_layout(title = "Line Intersections")
  
  # To calculate a kernel density, we need to reproject the data into a projected coordinate system
  # I'm using Albers because it covers the whole continent, and accuracy isn't vital
  line_is = st_transform(line_is,crs="+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")
  
  # Calculate a kernel density raster with a 3km bandwidth
  sm = smooth_map(line_is,bandwidth=2,to.Raster=TRUE,nrow=500,ncol=500)
  rast = sm[[1]]
  
  # Plot this raster on the map of Melbourne
  #tm_shape(comb,is.master = TRUE) + tm_lines(alpha=0) + tm_shape(rast) + tm_raster(palette="YlOrRd",title="KDE") + tm_shape(vic) + tm_borders()  + tm_shape(roads,col="black",alpha=0.5) + tm_lines() + tm_layout(title = "Intersection Kernel Density")
  
  # Find the 95th percentile of the raster value
  # And classify a new raster with just those values present
  # This is the bit that we will need real data to figure out if this is the best critera
  # What if we have LOTS of plumes all over the place? Will a 95th percentile make some go missing?
  # What if there is no smoke - will 95th percentile start highlighting smoke that isn't there?
  # May need to identify a useful fixed density instead
  
  # Calculate 95th percentile value
  q=quantile(rast,.95)
  
  # Copy the raster
  trast = rast
  
  # Set values of this raster where original raster is more than quantile value
  values(trast)[values(rast)>q]=1
  values(trast)[values(rast)<=q]=NA
  
  # Plot classified raster
  #tm_shape(comb,is.master = TRUE) + tm_lines(alpha=0) + tm_shape(trast) + tm_raster(palette="YlOrRd",title="KDE") + tm_shape(vic) + tm_borders()  + tm_shape(roads,col="black",alpha=0.5) + tm_lines() + tm_layout(title = "Threshold Density")
  orast = rasterToPolygons(trast,dissolve=TRUE)
  orast = st_as_sf(orast, crs="+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")
  orast = st_transform(orast,crs=4326)
  orast=st_buffer(orast,.000000001)
  orast$created_at = end_array[i]
  polygon_list[[i]]=orast
  
}

polygon_list = compact(polygon_list)
all_poly = do.call(rbind,polygon_list)

vector_list = compact(vector_list)
all_vector = do.call(rbind,vector_list)

point_list = compact(point_list)
all_points = do.call(rbind,point_list)

write_sf(all_poly,"tiles/smoke_poly.geojson",delete_dsn=TRUE)
write_sf(all_vector,"tiles/smoke_lines.geojson",delete_dsn=TRUE)
write_sf(all_points,"tiles/smoke_points.geojson",delete_dsn=TRUE)
