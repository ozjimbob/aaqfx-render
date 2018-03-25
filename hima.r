#!/usr/bin/Rscript
library(ncdf4)
library(raster)
library(rgdal)
require(zoo)
library(RCurl)

setwd("/mnt/R")
source("../config.r")

#date="20170416212000"

ttime = Sys.time()
ttime = as.POSIXct(format(ttime,tz="UTC"),tz="UTC")-15*60 
last10 = paste0(substr(format(ttime,"%M"),1,1),"0")
date = paste0(format(ttime,"%Y%m%d%H"),last10)

working = "hima/"
# download sub-files
url_root = sprintf("ftp://ftp-reg.cloud.bom.gov.au/gms/IDE00218.%s.nc",date)
fn = getBinaryURL(url_root,userpwd=bompwd)
writeBin(fn, paste0(working,"bandb.nc"))

url_root = sprintf("ftp://ftp-reg.cloud.bom.gov.au/gms/IDE00219.%s.nc",date)
fn = getBinaryURL(url_root,userpwd=bompwd)
writeBin(fn, paste0(working,"bandg.nc"))

url_root = sprintf("ftp://ftp-reg.cloud.bom.gov.au/gms/IDE00220.%s.nc",date)
fn = getBinaryURL(url_root,userpwd=bompwd)
writeBin(fn, paste0(working,"bandr.nc"))

url_root = sprintf("ftp://ftp-reg.cloud.bom.gov.au/gms/IDE00221.%s.nc",date)
fn = getBinaryURL(url_root,userpwd=bompwd)
writeBin(fn, paste0(working,"bandir.nc"))

run1 = sprintf('%sgdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %sbandb.nc %stemp.tif',exeroot,working,working)
run2 = sprintf('%sgdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sbandb.tif',exeroot,working,working)
system(run1)
system(run2)



run1 = sprintf('%sgdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %sbandg.nc %stemp.tif',exeroot,working,working)
run2 = sprintf('%sgdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sbandg.tif',exeroot,working,working)
system(run1)
system(run2)

run1 = sprintf('%sdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %sbandr.nc %stemp.tif',exeroot,working,working)
run2 = sprintf('%sgdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sbandr.tif',exeroot,working,working)
system(run1)
system(run2)

run1 = sprintf('%sgdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %sbandir.nc %stemp.tif',exeroot,working,working)
run2 = sprintf('%sgdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sbandir.tif',exeroot,working,working)
system(run1)
system(run2)

blue=raster(paste0(working,"bandb.tif"))
gr1=raster(paste0(working,"bandg.tif"))
red=raster(paste0(working,"bandr.tif"))
gr2=raster(paste0(working,"bandir.tif"))

blue = resample(blue,red)
gr1 = resample(gr1,red)
gr2 = resample(gr2,red)


red = red*255
blue= blue*255
green = gr1 * 0.86 + gr2 * (1-0.86)
green = green * 255


# Scale Red
x = c(0, 24, 65, 125, 199, 255)
y = c(0, 65, 116, 179, 232, 255)
k=approxfun(x,y)
values(red)=k(values(red))


# Scale green
x = c(0, 5, 33, 75, 139, 208, 255)
y = c(0, 0, 64, 116, 178, 235, 255)
k=approxfun(x,y)
values(green)=k(values(green))

# Scale blue
x = c(0, 10, 40, 85, 149, 208, 255)
y = c(0, 3, 63, 116, 185, 238, 255)
k=approxfun(x,y)
values(blue)=k(values(blue))


red=floor(red)
green=floor(green)
blue=floor(blue)

t1 = brick(red,green,blue)
plotRGB(t1)


writeRaster(t1,"hima/stack_vis.tif",overwrite=TRUE,datatype="INT1U")




url_root = sprintf("ftp://ftp-reg.cloud.bom.gov.au/gms/IDE00213.%s.nc",date)
fn = getBinaryURL(url_root,userpwd=bompwd)
writeBin(fn, paste0(working,"bandct.nc"))

run1 = sprintf('%sgdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %sbandct.nc %stemp.tif',exeroot,working,working)
run2 = sprintf('%sgdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sbandct.tif',exeroot,working,working)
system(run1)
system(run2)


# Load IR Raster
IR=raster(paste0(working,"bandct.tif"))



ir2= 255 -  pmax( 2*(values(IR)-180), pmin(0, 255))
values(IR) = ir2

IR = disaggregate(IR,2,"bilinear")
# Define colour ramp
colv = colorRampPalette(c("black","grey60","white"),1)(256)

# Set colortable
colortable(IR)=colv
values(IR)[values(IR)<3]=3
# Set colortable

writeRaster(IR,"hima/stack_IR.tif",overwrite=TRUE,datatype="INT1U")


unlink(paste0(working,c("bandb.nc","bandct.nc","bandg.nc","bandir.nc","bandr.nc","temp.tif","bandb.tif","bandct.tif","bandg.tif","bandir.tif",  "bandr.tif")))

cmd=paste0(exeroot,"gdal_edit.py")
args=c("-unsetnodata","hima/stack_vis.tif")
system2(cmd,args)

if(sum(is.na(values(t1)))<100000){
	cmd=paste0(exeroot,"gdal2tiles.py")
	args=c("hima/stack_vis.tif","-a 255,255,255,255","-z 5-11",sprintf("tiles/hima_%s_vis",date))
	system2(cmd,args)
}

# Expand indexed tif to a rgba vrt
cmd=paste0(exeroot,"gdal_translate")
args=c("-of vrt","-expand rgba","hima/stack_IR.tif","hima/stack_IR.vrt")
system2(cmd,args)

cmd=paste0(exeroot,"gdal_edit.py")
args=c("-unsetnodata","hima/stack_IR.vrt")
system2(cmd,args)

cmd=paste0(exeroot,"gdal2tiles.py")
args=c("hima/stack_IR.vrt","-a 1","-z 5-11",sprintf("tiles/hima_%s_IR",date))
system2(cmd,args)

unlink(paste0(working,c("stack_IR.tif","stack_vis.tif","stack_IR.vrt")))

ttime = Sys.time()
ttime = as.POSIXct(format(ttime,tz="UTC"),tz="UTC")-24*60*60


flist = file.info(list.files("tiles",pattern="hima_",full.names=TRUE))
todel = rownames(flist)[flist$ctime < ttime]
unlink(todel,recursive=TRUE)

source("cat.r")


