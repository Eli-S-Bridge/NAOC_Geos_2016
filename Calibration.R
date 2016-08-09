#Geolocation analysis with Open Source Tools
#2016 North American Ornithological Congress, Washington D.C.

#Sponsored by: 

#Migrate Technology LLC.-- www.migratetech.co.uk

#The Cooper Ornithological Society

#The National Science Foundation

#--------------------------------------------------------------

#Calibrating for geolocator analyses

library(SGAT)
library(BAStag)
library(FLightR)

#start by reading in a data set consisting of twilight times
twl <- read.csv("data/749_twl.csv")

#format the date columns
twl$tFirst <- as.POSIXct(twl$tFirst, "UTC")   
twl$tSecond <- as.POSIXct(twl$tSecond, "UTC")

#Lets say we know the bird stayed put for about 3 weeks after tag deployment.
#We can establish a subset of the twilight data for just those three weeks,
#from the start date (2011-06-15) to 2011-07-06
calib1 <- subset(twl, tFirst < as.POSIXct("2011-07-06", "UTC"))

#Define the coordinates of the tag for the calibration period (lon and lat)
start <- c(-80.46,	42.62)

#You can use getElevation() in GeoLight to get the median sun angle for
#the twilights in your calibration period
library(GeoLight)
cal1 <- getElevation(twl = calib1, known.coord = start, plot = T)

cal1

#You can now use GeoLight functions to get a first look at location data on a map
#Calculate coordinates for a track and run tripmap
track <- coord(twl = twl, degElevation = cal1, tol = 0,
               method = "NOAA", note = TRUE)
tripMap(track, equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)





#You can also choose a calibration period from a light image

library("BAStag")                                 #open the BAStag package
d.lux <- readLig("data/749_000.lig", skip = 0)         #read the data into a dataframe called d.lux
d.lux <- subset(d.lux,select=c("Date","Light"))   #reduce the dataframe to just Date and Light

offset = 19
### Plot your light data with predicted twilight times assuming bird was stationary at deployment site for entire year
lightImage(d.lux, offset = offset)

### Leave plot from above open. The locator function can be used click on the outer boundaries of the calibration 
### period. Once you run code below, you get a target cursor. First click defines the start 
### of the calibration period, second click defines the end. User can define calibration period however they choose. 
### Choice will affect resulting twilight error distribution.
tm.calib <- as.POSIXct(locator(n=2)$x, origin = "1970-01-01")

### See results of your selection
tm.calib
calib1 <- subset(twl, tFirst > tm.calib[1] & tFirst < tm.calib[2])

#You can now repeat the steps above to get a sun angle and a quick track

cal1 <- getElevation(twl = calib1, known.coord = start, plot = T)
cal1
track <- coord(twl = twl, degElevation = cal1, tol = 0,
               method = "NOAA", note = TRUE)
tripMap(track, equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)





#What if you are not sure what data are good for calibrating?
#Here's a process that uses FlightR to delineate good calibration periods

Proc.data<-get.tags.data("data/A2_FLightR_twl.csv") #opens and formats data straight from TAGS formatted csv file

start=c(5.43, 52.93)  #tracking orgin, start location longitude and latitude

plot.slopes.by.location(Proc.data=Proc.data, location=start)

#Use abline to visualize potential calibration periods
abline(v=as.POSIXct("2013-08-22")) # end of first calibration period
abline(v=as.POSIXct("2014-04-05")) # start of the second calibration period

# Next create a data.frame with a separate line is designated for each calibration period. 
# The columns are: 
#     (1) start of the calibration period, 
#     (2) end of the calibration period, 
#     (3) longitude of the calibration location and, 
#     (4) latitude of the calibration location.

Calibration.periods<-data.frame(   
  calibration.start=as.POSIXct(c(NA, "2014-05-05")),   #This will create two lines
  calibration.stop=as.POSIXct(c("2013-08-23", NA)),
  lon=5.43, lat=52.93) #use c() also for the geographic coordinates, if you have more than one calibration location (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))

#View results
Calibration.periods

#create a calibration object 
Calibration<-make.calibration(Proc.data, Calibration.periods)

#save it for later use
save(Calibration, file = "data/FLightR_calibration")


