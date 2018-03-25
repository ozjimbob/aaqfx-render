library(tidyverse)
library(RCurl)
library(xml2)

url_root = "ftp://ftp-reg.cloud.bom.gov.au/fwo/IDV60920.xml"
dat=getURLContent(url_root,userpwd="bom893:bmaT94jN",binary=FALSE)
dat=read_xml(dat)
dl = xml_find_all(dat,".//observations")
dl = xml_find_all(dl,".//station")

st_id = xml_attr(dl,"wmo-id")
st_name = xml_attr(dl,"stn-name")
st_lat = xml_attr(dl,"lat")
st_lon = xml_attr(dl,"lon")

get_var = function(dl,varname){
  x=map(dl,xml_find_all,xpath=sprintf("period/level/element[@type='%s']",varname))
  x=map(x,xml_text)
  x[sapply(x, is_empty)] <- NA
  as.numeric(unlist(x))
}


df = tibble(id = st_id,
            name = st_name,
            lat=as.numeric(st_lat),
            lon=as.numeric(st_lon),
            apparent_temp = get_var(dl,"apparent_temp"),
            air_temperature = get_var(dl,"air_temperature"),
            dew_point = get_var(dl,"dew_point"),
            pres = get_var(dl,"pres"),
            rel_humidity = get_var(dl,"rel-humidity"),
            wind_dir_deg = get_var(dl,"wind_dir_deg"),
            wind_spd_kmh = get_var(dl,"wind_spd_kmh"),
            rainfall_24hr = get_var(dl,"rainfall_24hr")
            )

library(tmap)
library(sf)
dat = st_as_sf(df,coords=c("lon","lat"),crs=4326)

m = tm_shape(dat) + tm_text("air_temperature")
