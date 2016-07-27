#Geolocation analysis with Open Source Tools
#2016 North American Ornithological Congress, Washington D.C.

#Sponsored by: 

#Migrate Technology LLC.-- www.migratetech.co.uk

#The Cooper Ornithological Society

#The National Science Foundation

#--------------------------------------------------------------

#Preprocessing a .lig file from a Migrate Technology tag

#First reduce the data down to just datestamps and light levels
#use the readLig function in BAStag to read in the data

library("BAStag")                                 #open the BAStag package
d.lig <- readLig("data/749_000.lig", skip = 0)         #read the data into a dataframe called d.lig
d.lig <- subset(d.lig,select=c("Date","Light"))   #reduce the dataframe to just Date and Light

#Lets view the data.
#You can use the plot function to look at small pieces of the dataset.
plot(d.lig$Date[3000:5000], d.lig$Light[3000:5000], type = "l")

#For a more complete view use the lightimage() function in the BAStag package.
#In this graph each vertical line is a day (24 hours) of data.
#Low light levels are shown in dark shades of gray and high levels are light shades.
#The offset value of 17 adjusts the y-axis to put night (dark shades) in the middle.
lightImage(d.lig, offset = 17, zlim = c(0, 64), dt = 120) 
#Note the dark spots near the end of the dataset. These are probably associated with nesting (in a dark box).

#Options for editing twilights.
#   preprocessLight() from the BAStag package
#   twilightCalc() from the GeoLight package
#   TAGS - a web-based editor

#Establish a threshold for determining twilights (what light value separates day and night?)
#The best choice is the lowest value that is consistently above any noise in the nighttime light levels
#For many BAS data sets, 2.5 is a good threshold
threshold = 2.5

#preprocessLight() is an interactive function for editing light data and deriving twilights
#But it will not work with Rstudio sever.
#You can try it on your own machine, but it is tricky with a Mac
twl <- preprocessLight(d.lig, threshold = threshold)

#You can use the twilightCalc() function in GeoLight to go through the data a bit at a time.
#This can take a few minutes.
#twilightCalc tries to choose the most likely twilights, but it can be wrong.
#If you want to trust the function and not look at the data (not recommended) the set "ask" to "FALSE."
library(GeoLight)    #activate the GeoLight package
twl <- twilightCalc(datetime = d.lig$Date, light = d.lig$Light,
                    LightThreshold = threshold, ask = T, preSelection = TRUE, 
                    nsee = 5000, allTwilights = F)
#
#The allTwilight=TRUE option causes the function to save two dataframes in a list
#The first dataframe  ("allTwilights") contains all the light data with twilight instances added in (Useful for FlightR and SGAT)
#The second ("consecTwilights") is just the twilight events and is useful for GeoLight

#Let's truncate the data to cut off the messy data at the end
#and save it to a csv file for later use.
datetime <- as.POSIXct(twl$tFirst, "UTC")  #Get datestamps into a format R can work with
twl <- twl[datetime < as.POSIXct("2012-05-20", "UTC"),] 
write.csv(twl, file = "data/749_twl.csv", quote = FALSE, row.names = FALSE)


# You can also use the online data editing tools in TAGS
# at tags.animalmigration.org
# To pre- process the dataset with TAGS it may be necessary 
# to compress the data (i.e. remove repeated light levels)
# We've provided a function to let you do this compression
source("Geolocator_compression.R")
d.lig <- geolocator_compress(datetime = d.lig$Date, light = d.lig$Light)
#save the data to a text file that TAGS can read
write.table(x = d.lig, file = "data/749_000_short.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = c("datetime", "light"))

#Process the data with TAGS and read it back into R.

twl <- read.csv("data/749_TAGS_twl.txt")
# This data set was only edited up until April 20, 2012. 
# (By then the birds were back on the breeding grounds and going in and out of next boxes)
# So lets remove the unedited data (everything later than April 20).
datetime <- strptime(twl$datetime, "%Y-%m-%dT%H:%M:%OSZ", "GMT")  #Get datestamps into a format R can work with
twl <- twl[datetime < as.POSIXct("2012-05-20", "GMT"),]                  #Remove data logged after Apr 20, 2012
# These data are in TAGS format, which FlightR can work with. So lets save them.
write.csv(twl, file = "data/749_twl_FlightR.csv", quote = FALSE, row.names = FALSE)






