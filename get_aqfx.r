#!/usr/bin/Rscript
tdate = format(Sys.Date(),"%Y%m%d")
setwd("/mnt/R")
cmd = sprintf("scp gxw581@raijin.nci.org.au:/short/en0/share/aqfx/Latest/AQFx_NCoutput_%s.tar aqfx/aqfx.tar",tdate)
system(cmd)
cmd="tar -xf aqfx/aqfx.tar -C aqfx"
system(cmd)
unlink("aqfx/aqfx.tar")

unlink(sprintf("aqfx/aus_%sCavg.nc",tdate))
unlink(sprintf("aqfx/se_%sCavg.nc",tdate))
unlink(sprintf("aqfx/vtas_%sCavg.nc",tdate))



vtas = sprintf("aqfx/%s_PM25plus_vtas.nc",tdate)
se = sprintf("aqfx/%s_PM25plus_se.nc",tdate)
aus = sprintf("aqfx/%s_PM25plus_aus.nc",tdate)

file.rename(vtas,"aqfx/vtas.nc")
file.rename(se,"aqfx/se.nc")
file.rename(aus,"aqfx/aus.nc")


