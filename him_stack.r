#!/usr/bin/Rscript
library(ncdf4)
library(raster)
library(rgdal)
require(zoo)


# Get current time as string suitable for Himawari filename substitution
#ttime = Sys.time()
#ttime = as.POSIXct(format(ttime,tz="UTC"),tz="UTC")
#last10 = paste0(substr(format(ttime,"%M"),1,1),"0")
#date = paste0(format(ttime,"%Y%m%d%H"),last10,"00")

date="20170416212000"

# Set working directory
working = "/mnt/hima/"

# NOTDONE
# Download most recent files
#



# Format filenames for the colour bands
tfile_1=sprintf("%s-P1S-ABOM_BRF_B01-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc",date)
tfile_2=sprintf("%s-P1S-ABOM_BRF_B02-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc",date)
tfile_3=sprintf("%s-P1S-ABOM_BRF_B03-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc",date)
tfile_4=sprintf("%s-P1S-ABOM_BRF_B04-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc",date)


# Translate Band 1 to GeoTiff
run1 = sprintf('gdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %s%s %stemp.tif',working,tfile_1,working)
run2 = sprintf('gdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sband1.tif',working,working)
system(run1)
system(run2)

# Translate Band 2 to Geotiff
run1 = sprintf('gdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %s%s %stemp.tif',working,tfile_2,working)
run2 = sprintf('gdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sband2.tif',working,working)
system(run1)
system(run2)

# Translate Band 3 to Geotiff
run1 = sprintf('gdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %s%s %stemp.tif',working,tfile_3,working)
run2 = sprintf('gdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sband3.tif',working,working)
system(run1)
system(run2)

# Translate Band 4 to Geotiff
run1 = sprintf('gdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %s%s %stemp.tif',working,tfile_4,working)
run2 = sprintf('gdalwarp -overwrite -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sband4.tif',working,working)
system(run1)
system(run2)

# Load bands as rasters
blue=raster(paste0(working,"band1.tif"))
gr1=raster(paste0(working,"band2.tif"))
red=raster(paste0(working,"band3.tif"))
gr2=raster(paste0(working,"band4.tif"))

# Convert to 0-255 color range and adjust green
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

# Convert to integers
red=floor(red)
green=floor(green)
blue=floor(blue)

t1 = brick(red,green,blue)

# Write TIF
writeRaster(t1,sprintf("%sstack_%s.tif",working,date),overwrite=TRUE,datatype="INT1U")

# Render Tiles
cmd="/usr/bin/gdal2tiles.py"
args=c(sprintf("%sstack_%s.tif",working,date),"-z 4-10","-r bilinear",sprintf("%s/h8_vis_%s",working,date))
system2(cmd,args)


# Read IR 2km band
tfile_IR=sprintf("%s-P1S-ABOM_OBS_B13-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc",date)

# Transform IR band to geotiff
run1 = sprintf('gdal_translate -a_srs "+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0=140.7 +no_defs" -a_ullr -5500000 5500000 5500000 -5500000 %s%s %stemp.tif',working,tfile_IR,working)
run2 = sprintf('gdalwarp -overwrite -r cubic -t_srs "+proj=latlong +ellps=WGS84" -wo SOURCE_EXTRA=100 %stemp.tif -te 139 -44 153 -33 %sstack_IR.tif',working,working)
system(run1)
system(run2)

# Read ID band
IR=raster(paste0(working,"stack_IR.tif"))


# Normalize values
values(IR) <- 255-pmax( 0, pmin( values(IR-180), 255))

# Define colour ramp
colv = colorRampPalette(c("black","grey20","grey","yellow","purple","blue"),.3)(256)

# Set colortable
colortable(IR)=colv
plot(IR)

# Write IR band
writeRaster(IR,sprintf("%sIR_%s.tif",working,date),overwrite=TRUE,datatype="INT1U")

# Convert to RGB tif
cmd="/usr/bin/gdal_translate"
args=c("-of vrt","-expand rgb",sprintf("%sIR_%s.tif",working,date),sprintf("%sIR_%s.vrt",working,date))
system2(cmd,args)

# Render tiles
cmd="/usr/bin/gdal2tiles.py"
args=c(sprintf("%sIR_%s.vrt",working,date),"-z 4-10","-r bilinear",sprintf("%s/h8_IR_%s",working,date))
system2(cmd,args)


unlink(paste0(working,c("band1.tif","band2.tif","band3.tif","band4.tif","temp.tif","stack_IR.tif",sprintf("stack_%s.tif",date),sprintf("IR_%s.tif",date),sprintf("IR_%s.vrt",date))))

