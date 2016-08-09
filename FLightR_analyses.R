# Geolocation analysis with Open Source Tools
# 2016 North American Ornithological Congress, Washington D.C.

## Sponsored by: 

## Migrate Technology LLC.-- www.migratetech.co.uk

## The Cooper Ornithological Society

## The National Science Foundation

#--------------------------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#--------------------------------------------------------------

#FLightR analyses.
#Note that we have already performed some of the preparatory steps

library(SGAT)
library(BAStag)
library(GeoLight)

# install latest FLightR
library(devtools)
install_github('eldarrak/FlightR')
library(FLightR)

#read in pre-processed data and calibration data
Proc.data<-get.tags.data("data/A2_FLightR_twl.csv") #opens and formats data straight from TAGS formatted csv file
load("data/FLightR_calibration") #loads object called Calibration

# Establish a spatial grid and rules for possible migration paths
# Default grid resolution is is 50 X 50 km.
# Terms "left," "right," "bottom," and "top" define your bounding box. 
# 
# distance.from.land.allowed.to.use should be vector with length of two, 
#   first number is negative distance allowed to use while over land 
#   (restricts birds to flying only over coastlines and water) 
#   and second is distance from land allowed to use while over water 
#   (restricts birds to flying only over coastlines and land). 
#   
# distance.from.land.allowed.to.stay should be vector with length of two, 
#   first number is negative distance where bird is allowed to be stationary, 
#   (restricts birds to landing only on coastlines and land)
#   and second is distance from land allowed to fly over during twilight while 
#   over water. (restricts birds to landing only on coastlines and water)

Grid<-make.grid(left=-14, bottom=30, right=13, top=57,
                distance.from.land.allowed.to.use=c(-Inf, Inf),  #Use infinity to withold any restrictions on migration paths
                distance.from.land.allowed.to.stay=c(-Inf, Inf))

# CREATE A PROPOSAL
# Here we create an array of settings and data that incorporates all the objects
# created at earlier steps: 
#    - the light data with the detected twilight events (Proc.data), 
#    - the spatial parameters (Grid), 
#    - geographic coordinates of the starting location (start)
#    - and the calibration parameters (Calibration).
# This can take a while.

all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration)

#Save this if you want
#save(all.in, file = "data/FLightR_alldata.Rdata", compress = T)
#load(file = "data/FLightR_alldata.Rdata")

# RUN THE PARTICLE FILTER

# Here is where the results are are calculated (coordinates, behavior, stationarity).
# Within the function run.particle.filter, the following parameters can be preset:
#   -number of particles (1e4 is recommended for test and 1e6 for the analysis) 
#   -known.last = TRUE if you know the track ends where it began  (FALSE is the default) 
#   -check.outliers = TRUE, for the "on a fly" discard of outliers (only recommended to make pretty maps).

nParticles=1e4     #just a quick trial
a= Sys.time()       #This lets you measure the analysis time

Result<-run.particle.filter(all.in, threads=-1,
                            nParticles=nParticles, known.last=TRUE,
                            precision.sd=25, check.outliers=F)
b= Sys.time()
b-a                 #how long did it take?

#Now save your results are as an RData object.
#save(Result, file="data/A2_FLightR_results.RData")
load("data/A2_FLightR_results.RData")

#PLOTTING
#Plot a simple map
map.FLightR.ggmap(Result)

#PLOTTING
#Plot and save a simple map
map.FLightR.ggmap(Result, save.options = list(filename = "data/FLightR.map.pdf"))

#Plot lon lat graph
plot.lon.lat(Result)
