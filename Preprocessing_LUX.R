#Geolocation analysis with Open Source Tools
#2016 North American Ornithological Congress, Washington D.C.

#Sponsored by: 

#Migrate Technology LLC.-- www.migratetech.co.uk

#The Cooper Ornithological Society

#The National Science Foundation

####################################################################
#### Read in data from a .lux file from a Migrate Technology tag ###
####################################################################

# If necessary tell R where to find files.
# setwd("~/NAOC_Geos_2016")

#First reduce the data down to just datestamps and light levels
#use the readMTlux function in TwGeos to read in the data

library("TwGeos")
d.lux<-readMTlux("data/A2_raw_data.lux")     #read the data into a dataframe called d.lux
head(d.lux)         #view the first few lines
max(d.lux$Light)    #note the maximum value
min(d.lux$Light)    #note the minimum value

#Lets view the light data.
#You can use the plot function to look at small pieces of the dataset.
plot(d.lux$Date[7000:8000], d.lux$Light[7000:8000], type = "l")

#In lux files light values go very high, so we can log transform data before selecting twilights.
d.lux$Light<-log(d.lux$Light)

#Now try the simple plot again
plot(d.lux$Date[7000:8000], d.lux$Light[7000:8000], type = "l")

#For a more complete view use the lightimage() function in the TwGeos package.
#In this graph each vertical line is a day (24 hours) of data.
#Low light levels are shown in dark shades of gray and high levels are light shades.
#The offset value of 12 adjusts the y-axis to put night (dark shades) in the middle.
lightImage(d.lux, offset = 12, zlim = c(0, max(d.lux$Light)), dt = 120) 
#Note the dark areas near the beginning and end of the dataset. 
#These are probably associated with nesting and the pre/post deployment period.

#------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#------------------------------------------

#Options for editing twilights.
#   preprocessLight() from the TwGeos package
#   findTwilights() from the TwGeos package
#   twilightCalc() from the GeoLight package
#   TAGS - a web-based editor

#Establish a threshold for determining twilights (what light value separates day and night?)
#The best choice is the lowest value that is consistently above any noise in the nighttime light levels
#For this log transformed LUX tag a good choice appears to be 1.5 which is about the same as log(4.5)
threshold = log(4.5)

#------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#------------------------------------------

## preprocessLight() is an interactive function for editing light data and deriving twilights
## Unfortunately, it does not work with an RStudio web interface (i.e. a virtual machine)
## Note: if you are working on a Mac set gr.Device to X11 and install Quartz first (https://www.xquartz.org)
## See help file for details on the interactive process.

## for pc
twl <- preprocessLight(d.lux, threshold = threshold, offset = 18, lmax = 12, gr.Device = "default")

## for mac
twl <- preprocessLight(d.lux, threshold = threshold, offset = 18, lmax = 12, gr.Device = "x11")


#------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#------------------------------------------

## The findTwilights function in TwGeos package just finds the twilight times 
## (without individual insepction and without editing)
## A so called 'seed' is required, This is just a known date and time when you know it is night
## (see help file: "?findTwilights()").
## You can establish the seed by graphing the data and clicking on a nightime period.

plot(d.lux$Date[3000:5000], d.lux$Light[3000:5000], type = "o", pch = 16, cex = 0.5)
seed <- as.POSIXct(locator(n=1)$x, origin  = "1970-01-01", tz = "GMT")
twl  <- findTwilights(d.lux, threshold, include = seed)

## See if it worked
lightImage(d.lux, offset = 12, zlim = c(0, 12), dt = 120)
tsimagePoints(twl$Twilight, offset = 12, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

## The function twilightEdit may help to find outliers and either remove or edit them according to one rule:
## rule: if a twilight time is 'outlier.mins' (e.g. 45) minutes different to its sourrounding twilight times, 
## defined by 'window' - the number of sourrounding twilight times (e.g. 4), and the sourrounding twiligth times 
## are within 'stationary.mins' (e.g. 15) minutes, the outlayer will be moved in between the two neighboring twilights
## or deleted if the sourrounding twilights are > 'stationary.mins'.

## This allows fast and easily reproducable definition of twilight times.
twl <- twilightEdit(twl, window = 4, outlier.mins = 45, stationary.mins = 25, plot = T)

## Plot: grey points are either deleted (crossed) or edited (moved) twiligth times

## Plot the edited data on ligthImage
lightImage(d.lux, offset = 12, zlim = c(0, 12), dt = 120)
tsimagePoints(twl$Twilight[!twl$Deleted], offset = 12, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise[!twl$Deleted], "dodgerblue", "firebrick"))


## adjust twilight data since this tag records maximum values over time (in this case every 5 minutes)
twl <- twilightAdjust(twl, 5*60)

#------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#------------------------------------------

#You can use the twilightCalc() function in GeoLight to go through the data a bit at a time.
#This can take a few minutes.
#twilightCalc tries to choose the most likely twilights, but it can be wrong.
#If you want to trust the function and not look at the data (not recommended) the set "ask" to "FALSE."
library(GeoLight)    #activate the GeoLight package
twl <- twilightCalc(datetime = d.lux$Date, light = d.lux$Light,
                    LightThreshold = threshold, ask = TRUE, preSelection = TRUE, 
                    nsee = 3000, allTwilights = FALSE)

head(twl)          #view the twilight data
class(twl$tFirst)  #make sure the dates are in POSIXct format (they are)
#These data are already formatted for GeoLight, so lets save it.
write.csv(twl, file = "data/A2_twl_GeoL.csv", quote = FALSE, row.names = FALSE)

#To format for FlightR we need the data in TAGS format and not log transformed
#First read in the raw, untransformed light data and make some formatting changes
d.lux<-readMTlux("data/A2_raw_data.lux")    
d.lux <- data.frame(datetime = format(d.lux$Date, format = "%Y-%m-%dT%H:%M:%S.000Z"),
                    light = d.lux$Light, twilight = 0, interp = FALSE, excluded = FALSE)

#Get the twilight data into the same format
twl <- data.frame(datetime = format(twl$tFirst, format = "%Y-%m-%dT%H:%M:%S.000Z"), 
                  light = 4.5, twilight = twl$type, interp = TRUE, excluded = FALSE)

#Bind the rows of the two data frames together and sort them 
twl <- rbind(d.lux,twl)                        
twl <- twl[order(as.character(twl$datetime)),]

#Lets remove the useless data from the beginning and end of the data set
datetime <- strptime(twl$datetime, "%Y-%m-%dT%H:%M:%OSZ", "GMT")  #vector of datestamps in a format R can work with
twl <- twl[datetime < as.POSIXct("2014-05-16", "GMT") & datetime > as.POSIXct("2013-06-11", "GMT"),]    #Remove data logged after Apr 20, 2012

#Lets save it
write.csv(twl$allTwilights, file = "data/A2_twl_FlightR.csv", quote = FALSE, row.names = FALSE)

#------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#------------------------------------------

# To pre- process the dataset with TAGS it may be necessary to compress the data (i.e. remove repeated light levels)
source("Geolocator_compression.R")
d.lux<-readMTlux("data/A2_raw_data.lux")    #restore d.lux to orginal state (with large numbers))
d.lux <- geolocator_compress(datetime = d.lux$Date, light = d.lux$Light, max = 15)
# Unfortunately, TAGS rounds light data, which is not good for FLightR analyses.
# A work around is to multiply the light values by 10 for processing.
# We'll undo this by dividing by ten later on
d.lux$light <- d.lux$light * 10

# Save the data to a text file that TAGS can read
write.table(x = d.lux, file = "data/A2_short.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = c("datetime", "light"))

#Process the data with TAGS and read it back into R.

#The file "A2_TAGS_twl.txt" has been processed in TAGS with a threshold of 49.9 (roughly 4.5 times 10)
twlx <- read.csv("data/A2_TAGS_twl.txt")       #read in the twilight data

#FLightR makes use of all that continuous light level from a Migration Technology tag
#So we should restore that data to the output from TAGS
twl <- twl[twl$twilight != 0 & twl$excluded == FALSE,]  #get just the twilights from the TAGS output (twilight not equal to 0)
twl$light <- twl$light/10       #divide by 10 to get back to the original light values (almost)

d.lux <- readMTlux("data/A2_raw_data.lux")      #restore d.lux to orginal state (with large numbers))
d.lux <- data.frame(datetime = format(d.lux$Date, format = "%Y-%m-%dT%H:%M:%S.000Z"),
                    light = d.lux$Light, twilight = 0, interp = FALSE, excluded = FALSE)
twl <- rbind(d.lux,twl)
twl <- twl[order(as.character(twl$datetime)),]

#As we did above, lets remove the useless data from the beginning and end of the data set
datetime <- strptime(twl$datetime, "%Y-%m-%dT%H:%M:%OSZ", "GMT")  #Get datestamps into a format R can work with
twl <- twl[datetime < as.POSIXct("2014-05-16", "GMT") & datetime > as.POSIXct("2013-06-11", "GMT"),]           #Remove data logged after Apr 20, 2012

#Save a csv file for later use
write.csv(twl, file = "data/A2_FlightR_twl.csv", quote = FALSE, row.names = FALSE)

