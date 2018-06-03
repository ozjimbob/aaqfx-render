#!/usr/bin/Rscript


tdate = format(Sys.Date(),"%Y%m%d")
setwd("/mnt/R")
a=1
cnt = 1
while(a==1){
  cmd = sprintf("scp gxw581@raijin.nci.org.au:/short/en0/share/aqfx/Latest/AQFx_NCoutput_%s.tar aqfx/aqfx.tar",tdate)
  a=system(cmd)
  cnt = cnt + 1
  if(cnt = 6 * 12){
    quit()
  }
  if(a==0){
    break
  }
  Sys.sleep(10*60)
  
  
}


cmd="tar -xf aqfx/aqfx.tar -C aqfx"
system(cmd)
unlink("aqfx/aqfx.tar")

unlink(sprintf("aqfx/aus_%sCavg.nc",tdate))
unlink(sprintf("aqfx/se_%sCavg.nc",tdate))
unlink(sprintf("aqfx/vtas_%sCavg.nc",tdate))
#full_vtas = sprintf("aqfx/vtas_%sCavg.nc",tdate)
#full_se = sprintf("aqfx/se_%sCavg.nc",tdate)
#full_aus = sprintf("aqfx/aus_%sCavg.nc",tdate)

#file.rename(full_vtas,"aqfx/full_vtas.nc")
#file.rename(full_se,"aqfx/full_se.nc")
#file.rename(full_aus,"aqfx/full_aus.nc")

vtas = sprintf("aqfx/%s_PM25plus_vtas.nc",tdate)
se = sprintf("aqfx/%s_PM25plus_se.nc",tdate)
aus = sprintf("aqfx/%s_PM25plus_aus.nc",tdate)

file.rename(vtas,"aqfx/vtas.nc")
file.rename(se,"aqfx/se.nc")
file.rename(aus,"aqfx/aus.nc")

source("proc_ncdf.r")


