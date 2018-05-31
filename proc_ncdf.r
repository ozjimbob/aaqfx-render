#!/usr/bin/Rscript

library(ncdf4)
library(raster)
library(rgdal)
require(zoo)

#############################
### FUNCTIONS 
#############################
setwd("/mnt/R")

# Get list of species in the file
get_spec_list=function(input){
  spec=ncvar_get(nc=input,varid="name_fspec")
  nvar = nchar(spec)/4
  stseq=seq(from=1,to=nchar(spec),by=4)
  snpfn = function(st){substr(spec,st,st+3)}
  unlist(lapply(stseq,FUN=snpfn))
}

# Drill down a time series for given species and location
drill_spec=function(input,spec,lat,lon){
  sl = get_spec_list(input)
  vno = which(sl %in% spec)
  if(length(vno)==0){
    stop(paste0("Chemical species ",spec," not found."))
  }
  odat=ncvar_get(nc=input,varid="cavg")
  lons=ncvar_get(input,"lon")
  lats=ncvar_get(input,"lat")
  if(lat < min(lat) | lat > max(lat) | lon < min(lon) | lon > max(lon)){
    stop("Coordinates provided are out of range of grid.")
  }
  tlon=which.min(abs(lons-lon))
  tlat=which.min(abs(lats-lat))
  odat=odat[tlon,tlat,vno,]
  odat
}

# Return a map raster for a given hour and species
rast_spec=function(input,spec,hour){
  sl = get_spec_list(input)
  vno = which(sl %in% spec)
  if(length(vno)==0){
    stop(paste0("Chemical species ",spec," not found."))
  }
  if(hour > 24){
    stop(paste0("Hour must be between 1 and 24"))
  }
  odat=ncvar_get(nc=input,varid="cavg")
  lons=ncvar_get(input,"lon")
  lats=ncvar_get(input,"lat")
  lodiff=mean(diff(lons))
  ladiff=mean(diff(lats))
  
  odat=odat[,,vno,hour]
  rotm <- function(x) t(apply(x, 2, rev))
  odat=rotm(rotm(rotm(odat)))
  ra=raster(odat,xmn=min(lons)-(lodiff/2),xmx=max(lons)+(lodiff/2),ymn=min(lats)-(ladiff/2),ymx=max(lats)+(ladiff/2))
  proj4string(ra)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
  ra
}


# Return a 3D matrix for a given hour and species
cube_spec=function(input,spec){
  sl = get_spec_list(input)
  vno = which(sl %in% spec)
  if(length(vno)==0){
    stop(paste0("Chemical species ",spec," not found."))
  }
  
  odat=ncvar_get(nc=input,varid="cavg",c(1,1,vno,1),c(-1,-1,1,-1))
 
    odat
}


# Return a 3D matrix for a given hour and species
cube_met=function(input,spec){
  
  odat=ncvar_get(nc=input,varid=spec,c(1,1,1),c(-1,-1,-1))
  
  odat
}

cube_layer=function(input,spec,layer){
  
  odat=ncvar_get(nc=input,varid=spec,c(1,1,layer,1),c(-1,-1,-1,-1))
  
  odat
}

difun=function(x){rollmean(x,2)}



#########################
####  Examples ##########
#########################


data <- nc_open("aqfx/vtas.nc")

unlink("aqfx_done/vtas_01.nc")
file.rename("aqfx_done/vtas_02.nc","aqfx_done/vtas_01.nc")
file.rename("aqfx_done/vtas_03.nc","aqfx_done/vtas_02.nc")
file.rename("aqfx_done/vtas_04.nc","aqfx_done/vtas_03.nc")
file.rename("aqfx_done/vtas_05.nc","aqfx_done/vtas_04.nc")
file.rename("aqfx_done/vtas_06.nc","aqfx_done/vtas_05.nc")
file.rename("aqfx_done/vtas_07.nc","aqfx_done/vtas_06.nc")


#PM25
datPM25 = cube_met(data,"PM25")[,,1:24]
datOzone = cube_met(data,"Ozone")[,,1:24]
datLevo = cube_met(data,"Levo")[,,1:24]
datAOD = cube_met(data,"AOD")[,,1:24]






times = 0:23
tdate = Sys.Date()
year = as.numeric(format(tdate,"%Y"))
month= as.numeric(format(tdate,"%m"))
day = as.numeric(format(tdate,"%d"))
adate = as.numeric(as.POSIXct(sprintf("%d-%d-%d %d:00:00",year,month,day,times)))
  
lats = ncvar_get(data,"lat")
lons = ncvar_get(data,"lon")

dimTime <- ncdim_def("time","seconds since 1970-1-1 00:00:00",adate,longname="time",unlim=TRUE)
dimLat <- ncdim_def("lat","degrees_north",lats,longname="latitude")
dimLon <- ncdim_def("lon","degrees_east",lons,longname="longitude")
  
varPM25<- ncvar_def("PM25","ugm3",list(dimLon,dimLat,dimTime),longname="PM25 Concentration",prec="float")
varOzone<- ncvar_def("Ozone","ppb",list(dimLon,dimLat,dimTime),longname="Ozone Concentration",prec="float")
varLevo<- ncvar_def("Levo","ugm3",list(dimLon,dimLat,dimTime),longname="Levoglucosan Concentration",prec="float")
varAOD<- ncvar_def("AOD","ExC",list(dimLon,dimLat,dimTime),longname="Aerosol Optical Depth",prec="float")

  
  
new_file = nc_create("aqfx_done/vtas_07.nc",list(varPM25,varOzone,varLevo, varAOD))
  
ncatt_put(new_file,0,attname="CDI",attval="Climate Data Interface version 1.7.2 (http://mpimet.mpg.de/cdi")
ncatt_put(new_file,0,attname="CDO",attval="Climate Data Operators version 1.7.2 (http://mpimet.mpg.de/cdo)")
ncatt_put(new_file,0,attname="Conventions",attval="CF-1.4")
  
  
ncvar_put(new_file,varPM25,datPM25)
ncvar_put(new_file,varOzone,datOzone)
ncvar_put(new_file,varLevo,datLevo)
ncvar_put(new_file,varAOD,datAOD)

nc_close(new_file)
 

cmd="cdo -O mergetime aqfx_done/vtas_*.nc aqfx_done/merge_vtas.nc"
system(cmd)

source("../config.r")
system(geopw)


