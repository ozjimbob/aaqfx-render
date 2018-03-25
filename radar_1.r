#!/usr/bin/Rscript
# Radar interpretation
library(tidyverse)
library(stringr)
library(gtools)
library(raster)
library(geosphere)
library(RCurl)
setwd("/mnt/R")
source("../config.r")

#02 = Melbourne
#49 = Yarrawonda

rad_tab=data.frame(rad=c("02","49"),time=c(6,10))
this_rad="02"

ss_rad = rad_tab[rad_tab$rad==this_rad,]
interval = ss_rad$time

date = as.POSIXct(format(Sys.time(),tz="UTC"),tz="UTC") - (interval) * 60
mins = formatC((as.numeric(format(date,"%M")) %/% interval) * interval,width=2,flag="0")
fdate = paste0(format(date,"%Y%m%d%H"),mins)
#2018 03 04 15 18



url_root = sprintf("ftp://ftp-reg.cloud.bom.gov.au/radar/IDR%sPOL.%s.txt",this_rad,fdate)

dat=getURLContent(url_root,userpwd=bompwd,binary=FALSE)
dat=strsplit(dat,"\n")[[1]]

delta_table = c("!","[","a","b","c","]","@",
                "/","d","e","f","g","h","\\",
                "i","j","k","<","l","m","n",
                "o","p","-",".","+","q","r",
                "s","t","u",">","v","w","x",
                "(","y","S","T","U","V",")",
                "$","{","W","X","Y","}","&")
delta_table = asc(delta_table)
delta_matrix = matrix(delta_table,nrow=7,byrow=TRUE)

delta_lookup = function(x){
  pos = which(delta_matrix==x)
  col=ceiling(pos / 7)
  vcol=col-4
  
  row = which(x == delta_matrix[,col])
  row=row-4
  return(c(vcol,row))
}

line_lat = dat[str_detect(dat,"LATITUDE")] %>% str_split(" ",simplify=TRUE)
LAT = -as.numeric(line_lat[2])

line_stn = dat[str_detect(dat,"STNID")] %>% str_split(" ",simplify=TRUE)
STNID = as.numeric(line_stn[2])

line_lat = dat[str_detect(dat,"LATITUDE")] %>% str_split(" ",simplify=TRUE)
LAT = -as.numeric(line_lat[2])

line_lon = dat[str_detect(dat,"LONGITUDE")] %>% str_split(" ",simplify=TRUE)
LON = as.numeric(line_lon[2])

line_RNGRES =  dat[str_detect(dat,"RNGRES")] %>% str_split(" ",simplify=TRUE)
RNGRES = as.numeric(line_RNGRES[2])

line_ANGRES =  dat[str_detect(dat,"ANGRES")] %>% str_split(" ",simplify=TRUE)
ANGRES = as.numeric(line_ANGRES[2])

line_STARTRNG =  dat[str_detect(dat,"STARTRNG")] %>% str_split(" ",simplify=TRUE)
STARTRNG = as.numeric(line_STARTRNG[2])

line_ENDRNG =  dat[str_detect(dat,"ENDRNG")] %>% str_split(" ",simplify=TRUE)
ENDRNG = as.numeric(line_ENDRNG[2])

line_VIDRES =  dat[str_detect(dat,"VIDRES")] %>% str_split(" ",simplify=TRUE)
VIDRES = as.numeric(line_VIDRES[2])

start_dat = which(str_detect(dat,"COPYRIGHT"))+1
end_dat = which(str_detect(dat,"END RADAR IMAGE"))-1

dat = dat[start_dat:end_dat]

nbins = (ENDRNG - STARTRNG) / RNGRES


out_mat = matrix(0,nrow=360,ncol=nbins)

for(i in dat){
  ang = as.numeric(substr(i,2,4))
  #print(ang)
  # Grab encoding vector
  avec = substr(i,5,nchar(i))
  
  # Convert to ascii
  avec = asc(avec)
  
  # Generate flagmode vector
  mvec = character(length=length(avec))
  
  # Mark ABS absolute encoding values
  mvec[avec >= strtoi(41,16) & avec <= strtoi(50,16)] = "ABS"
  
  # Mark RLE run-length encoding values
  mvec[avec >= strtoi(30,16) & avec <= strtoi(39,16)] = "RLE"
  
  buffer=rep(NA,nbins)
  output_position = 1
  input_position = 1
  RLE_cache = 0
  while(input_position < length(mvec)){
    # If we have an absolute value, just store it
    if(mvec[input_position]=="ABS"){
      buffer[output_position] = abs(avec[input_position])-65
      RLE_cache = abs(avec[input_position])-65
      input_position = input_position+1
      output_position = output_position+1
      next
    }
    # If we have an RLE - count digits, then insert RLE_cache previous value
    if(mvec[input_position]=="RLE"){
      if(mvec[input_position+1]!="RLE"){
        RLE_count = as.numeric(chr(avec[input_position]))
        skip=0
      }
      
      if(mvec[input_position+1]=="RLE" & mvec[input_position+2]!="RLE"){
        RLE_count = as.numeric(paste0(chr(avec[input_position]),chr(avec[input_position+1])))
        skip=1
      }
      
      if(mvec[input_position+1]=="RLE" & mvec[input_position+2]=="RLE"){
        RLE_count = as.numeric(paste0(chr(avec[input_position]),chr(avec[input_position+1]),chr(avec[input_position+2])))
        skip=2
      }
      buffer[output_position:(output_position+(RLE_count-1))]=RLE_cache
      output_position = output_position + RLE_count
      input_position = input_position + skip + 1
      next
    }
    
    if(mvec[input_position] == ""){
      deltas = delta_lookup(avec[input_position])
      RLE_cache = RLE_cache+ deltas[1]
      buffer[output_position]=RLE_cache 
      RLE_cache = RLE_cache+ deltas[2]
      buffer[output_position+1]=RLE_cache 
      output_position = output_position + 2
      input_position = input_position + 1
    }
  }
  
  
  
  out_mat[ang,] = buffer
  
  
}
out_mat[is.na(out_mat)]=0


## Make translation matrix (radar-dependent)

#d2r=function(x){x*pi/180}
#r2d=function(x){x*(180/pi)}


## Define grid

grid_extent=extent(LON-2,LON+2,LAT-2,LAT+2)
grid_layer=raster(grid_extent,nrows=800,ncols=800)

#dist_template=grid_layer
#ang_template=grid_layer

#for(this_cell in 1:(800*800)){
  
#  cp=xyFromCell(dist_template,this_cell)
  
#  cell_lon=cp[1]
#  cell_lat=cp[2]
  
#  d=distCosine(p1=c(cell_lon,cell_lat),p2=c(LON,LAT),r=6371000)
  
#  ag=bearing(p1=c(LON,LAT),p2=c(cell_lon,cell_lat))
  
#  dist_template[this_cell]=d
#  ang_template[this_cell]=ag
#}

#ang_template2 =floor(ang_template %% 361)+1


#dist_template2 = (dist_template-STARTRNG)/RNGRES
#dist_template2 = dist_template2+1
#dist_template2[dist_template2<1]=1
#dist_template2 = floor(dist_template2)

#pos_vec = values(dist_template2-1) * 360 + values(ang_template2)

pos_vec = readRDS(paste0("radar/posvec_",this_rad,".RDS"))
values(grid_layer) = out_mat[pos_vec]
proj4string(grid_layer)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

values(grid_layer)[values(grid_layer)==0]=NA

colortable(grid_layer)=c("#FFFFFF00", "#CCCCFFDD", "#9999FFEE", "#6565FFff", "#3232FFff", "#0000FFff", "#3300CBff" ,"#660098ff" ,"#990065ff" ,"#CC0032ff","#FF0000ff", "#CC0000ff" ,"#980000ff" ,"#650000ff", "#320000ff", "#000000ff",rep("#FFFFFF00",256-16))
writeRaster(grid_layer,sprintf("radar/rad%s_%s.tif",fdate,this_rad),overwrite=TRUE,options=c("TIFFTAG_GDAL_NODATA=0"))


# Expand indexed tif to a rgba vrt
cmd=paste0(exeroot,"gdal_translate")
args=c("-of vrt","-expand rgba",sprintf("radar/rad%s_%s.tif",fdate,this_rad),sprintf("radar/rad%s_%s.vrt",fdate,this_rad))
system2(cmd,args)

cmd=paste0(exeroot,"gdal2tiles.py")
args=c(sprintf("radar/rad%s_%s.vrt",fdate,this_rad),"-a 0,0,0,0","-z 5-11",sprintf("tiles/rad_%s_%s",fdate,this_rad))
system2(cmd,args)

# Remove tifs

unlink(sprintf("radar/rad%s_%s.vrt",fdate,this_rad))
unlink(sprintf("radar/rad%s_%s.tif",fdate,this_rad))

# Remove old tiles

date = as.POSIXct(format(Sys.time(),tz="UTC"),tz="UTC") - (interval-1) * 60 - 120*60
flist = file.info(list.files("tiles",pattern="rad_",full.names=TRUE))
todel = rownames(flist)[flist$ctime < date]
unlink(todel,recursive=TRUE)

source("cat.r")
