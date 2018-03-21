#!/usr/bin/Rscript

library(ncdf4)
library(raster)
library(rgdal)
require(zoo)

#############################
### FUNCTIONS 
#############################


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

# List of dates to access, should be today's date minus a few/plus a few?
# Will see what martin comes up with

date_list=c("05","06","07")
day_list=c("01","02","03")

# For each day
for(id in seq_along(date_list)){
  
  # To do: will need to download the file from Raijin, process, and delete and end of each loop
  # as files are ~ 12gb in size.

  # Define the file location
  tfile<-sprintf("vtas_201704%s_cb05_aer2_0.030dg_s4.nc",date_list[id])
  
  # SCP copy...

  # Load the net cdf
  data<-nc_open(tfile)
  
  #PM25
  datPM25 = cube_spec(data,"EC25")
  prcp = cube_met(data,"prcp")[,,2:25]
  temp_a = cube_layer(data,"temp_a",1)[,,2:25]
  pres = cube_layer(data,"pres",1)[,,2:25]
  
  # wind vectors
  
  u=ncvar_get(nc=input,varid="u_comp",c(1,1,1,1),c(-1,-1,-1,-1))
  u = apply(u,c(2,3),difun)  
  v=ncvar_get(nc=input,varid="v_comp",c(1,1,1,1),c(-1,-1,-1,-1))
  v = apply(v,c(1,3),difun)  
  v=aperm(v, c(2,1,3))
  ws = sqrt(u**2+v**2)
  wd=90-(atan2(v,u)*180/pi)
  wd[wd<0]=wd[wd<0]+360
  
  ws=ws[,,2:25]
  wd=wd[,,2:25]
  
  times = 0:23
  day = ncvar_get(data,"day")
  month = ncvar_get(data,"month")
  year = ncvar_get(data,"year")
  adate = as.numeric(as.POSIXct(sprintf("%d-%d-%d %d:00:00",year,month,day,times)))
  
  lats = ncvar_get(data,"lat")
  lons = ncvar_get(data,"lon")
  
  
  dimTime <- ncdim_def("time","seconds since 1970-1-1 00:00:00",adate,longname="time",unlim=TRUE)
  dimLat <- ncdim_def("lat","degrees_north",lats,longname="latitude")
  dimLon <- ncdim_def("lon","degrees_east",lons,longname="longitude")
  
  varEC25 <- ncvar_def("EC25","ugm3",list(dimLon,dimLat,dimTime),longname="PM25 Concentration",prec="float")
  varprcp <- ncvar_def("prcp","mm",list(dimLon,dimLat,dimTime),longname="Precipitation",prec="float")
  vartemp_a<- ncvar_def("temp_a","degK",list(dimLon,dimLat,dimTime),longname="Ambient Temp",prec="float")
  varpres<- ncvar_def("pres","hPa",list(dimLon,dimLat,dimTime),longname="Atmospheric Pressure",prec="float")
  var_ws<- ncvar_def("ws","m/s",list(dimLon,dimLat,dimTime),longname="Wind Speed",prec="float")
  var_wd<- ncvar_def("wd","degrees",list(dimLon,dimLat,dimTime),longname="Wind Direction",prec="float")

  
  
  # Create one file per day of modelling
  # With file name indicating day of sequence?
  new_file = nc_create(sprintf("test_%s.nc",day_list[id]),list(varEC25,varprcp,vartemp_a, varpres,var_ws,var_wd))
  
  ncatt_put(new_file,0,attname="CDI",attval="Climate Data Interface version 1.7.2 (http://mpimet.mpg.de/cdi")
  ncatt_put(new_file,0,attname="CDO",attval="Climate Data Operators version 1.7.2 (http://mpimet.mpg.de/cdo)")
  ncatt_put(new_file,0,attname="Conventions",attval="CF-1.4")
  
  
  ncvar_put(new_file,varEC25,datPM25)
  ncvar_put(new_file,varprcp,prcp)
  ncvar_put(new_file,vartemp_a,temp_a)
  ncvar_put(new_file,varpres,pres)
  ncvar_put(new_file,var_ws,ws)
  ncvar_put(new_file,var_wd,wd)
  
  nc_close(new_file)

  # Delete original ncdf

}

system("cdo mergetime test_*.nc /mnt/geoserver/data/aqfx/merge.nc")i
unlink(paste0("test_",day_list,".nc"))

