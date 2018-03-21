#!/usr/bin/Rscript
library(tidyverse)
library(rjson)

fl = list.files("tiles")

of = list()

rad_tab=data.frame(rad=c("02","49"),time=c(6,10))


rdr = fl[substr(fl,1,3)=="rad"]

idx = 1
for(i in rdr){
	time = substr(i,5,16)
	site = substr(i,18,19)
	interval = rad_tab[rad_tab$rad==as.character(site),]$time
	filename = i
	ftime = as.POSIXct(time,format="%Y%m%d%H%M")
	stime = ftime - interval*60 + 1
	thisdf = data.frame(prod="radar",site=site,end_time=as.character(ftime),start_time=as.character(stime),layer=filename,stringsAsFactors=FALSE)
	of[[idx]]=thisdf
	idx=idx+1
}

of = bind_rows(of)

ol2 = list()
idx=1
unq_prods = unique(of$prod)
for(tprod in unq_prods){
	print(tprod)
	tof = filter(of,prod==tprod)
	ol2[[tprod]]=list()
	unq_sites = unique(tof$site)
	idx=1
	for(tsite in unq_sites){
		print(tsite)
		sub_site = filter(tof,site==tsite)
		ol2[[tprod]][[idx]]=list()
		ol2[[tprod]][[idx]]$id = tsite
		ol2[[tprod]][[idx]]$start_time = sub_site$start_time
		ol2[[tprod]][[idx]]$end_time = sub_site$end_time
		ol2[[tprod]][[idx]]$layer = sub_site$layer
		idx=idx+1
	}
}


# Himawari
of=list()
rdr = fl[substr(fl,1,4)=="hima"]

idx = 1
for(i in rdr){
	time = substr(i,6,17)
	site = substr(i,19,nchar(i))
	print(site)
	print(time)
	interval = 10
	filename = i
	ftime = as.POSIXct(time,format="%Y%m%d%H%M")
	stime = ftime - interval*60 + 1
	print(filename)
	thisdf = data.frame(prod="satellite",site=site,end_time=as.character(ftime),start_time=as.character(stime),layer=filename,stringsAsFactors=FALSE)
	of[[idx]]=thisdf
	idx=idx+1
}

of = bind_rows(of)

ol3 = list()
idx=1
unq_prods = unique(of$prod)
for(tprod in unq_prods){
	print(tprod)
	tof = filter(of,prod==tprod)
	ol3[[tprod]]=list()
	unq_sites = unique(tof$site)
	idx=1
	for(tsite in unq_sites){
		print(tsite)
		sub_site = filter(tof,site==tsite)
		ol3[[tprod]][[idx]]=list()
		ol3[[tprod]][[idx]]$id = tsite
		ol3[[tprod]][[idx]]$start_time = sub_site$start_time
		ol3[[tprod]][[idx]]$end_time = sub_site$end_time
		ol3[[tprod]][[idx]]$layer = sub_site$layer
		idx=idx+1
	}
}

ol2 = c(ol2,ol3)
out=toJSON(as.list(ol2))
cat(out,file="tiles/layers.json")

