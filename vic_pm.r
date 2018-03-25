#!/usr/bin/Rscript
library(tidyverse)
library(rvest)
library(RPostgres)
library(DBI)

source("../config.r")

setwd("/mnt/R")

vic_page=read_html("http://www.epa.vic.gov.au/Air/Bulletins/aqbhour.asp")
tbls <- html_nodes(vic_page, "table") %>%
  .[2] %>%
  html_table(fill = TRUE)

tbls=tbls[[1]]

tbls=tbls[4:nrow(tbls),]
names(tbls)=c("Region","Name","CO","O3","NO2","SO2","PM25","PM10","Vis","AQI","Summ",'AQI2')
tbls$Name=gsub("\\*","",tbls$Name)
tbls$PM25=gsub("\\*","",tbls$PM25)

pm_df = tbls

stations=read_csv("pm/vic_stations.csv")
names(stations) = c("Code","Name","Ele","Lat","Lon")
pm_df=left_join(pm_df,stations)


pm_df$CO = as.numeric(pm_df$CO)
pm_df$O3 = as.numeric(pm_df$O3)
pm_df$NO2 = as.numeric(pm_df$NO2)
pm_df$SO2 = as.numeric(pm_df$SO2)
pm_df$PM25 = as.numeric(pm_df$PM25)
pm_df$PM10 = as.numeric(pm_df$PM10)
pm_df$Vis = as.numeric(pm_df$Vis)
pm_df$AQI = as.numeric(pm_df$AQI)

pm_df = filter(pm_df,!is.na(Ele))
