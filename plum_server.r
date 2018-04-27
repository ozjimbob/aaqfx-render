#!/usr/bin/Rscript

library(plumber)
r<- plumb("smoke_trak_server.r")  # Where 'myfile.R' is the location of the file shown above
r$run(port=1235)

