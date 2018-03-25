library(tidyverse)
library(httr)
library(sf)
kml = GET("https://sentinel.ga.gov.au/geoserver/public/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=public:hotspot_current&outputFormat=application%2Fvnd.google-earth.kml%2Bxml",write_disk("hs.kml"))

dat = read_sf("hs.kml")
library(mapview)
