#Geolocation analysis with Open Source Tools
#2016 North American Ornithological Congress, Washington D.C.

#Sponsored by: 

#Migrate Technology LLC.-- www.migratetech.co.uk

#The Cooper Ornithological Society

#The National Science Foundation

#--------------------------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#--------------------------------------------------------------

#Geolight analyses.

library(GeoLight)

#For a fast and simple map of locations:
#Call up the twlight data
twl <- read.csv("data/749_twl.csv")

#make sure the datetime columns are in the correct format
twl$tFirst <- as.POSIXct(twl$tFirst, "UTC")   
twl$tSecond <- as.POSIXct(twl$tSecond, "UTC")

#provide a sun angle for the light threshold value 
cal1 = -5.52

#generate tracking data
track <- coord(twl = twl, degElevation = cal1, tol = 0,
               method = "NOAA", note = TRUE)

#plot a map
tripMap(track, equinox = TRUE, xlim = c(-90,-70), ylim = c(10,60), legend = TRUE)

#You can remove some obvious outliers with the distance filter. 
#SIMEON, I DON'T THINK THE UNITS OPTION FOR "DAY" IS WORKING.
d.filt <- distanceFilter(twl=twl, degElevation = cal1, distance = 35, units = "hour")
tripMap(track[d.filt,], equinox = TRUE, xlim = c(-90,-70), ylim = c(10,60), legend = TRUE)

#You can also remove outliers by applying a loess smoother to filter outlying twilights.
l.filt = loessFilter(twl, k = 3, plot = TRUE)
tripMap(track[l.filt,], equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)

#To use both filters, just multiply them (logical AND)
b.filt <- l.filt & d.filt
tripMap(track[b.filt,], equinox = TRUE, xlim = c(-90,-70), ylim = c(-30,50), legend = TRUE)

#TO REMOVE EQUINOX DATA
#The 2011 fall equinox was Sept 23
fall_eq <- as.POSIXct("2011-09-23")
#The 2012 spring equinox was March 20
spring_eq <- as.POSIXct("2012-03-20")
#how big a window
window <- 20
fall_eq_start <- fall_eq - 86400*window/2 # 86400 seconds in one day
fall_eq_end <- fall_eq + 86400*window/2 # 86400 seconds in one day
#The 2012 spring equinox was March 20
spring_eq_start <- spring_eq - 86400*window/2
spring_eq_end <- spring_eq + 86400*window/2

fe.filt <- (twl$tFirst < fall_eq_start | twl$tSecond > fall_eq_end)
se.filt <- (twl$tFirst < spring_eq_start | twl$tSecond > spring_eq_end)

#multiply (logical AND) all the filters together
all.filt <- l.filt & d.filt & fe.filt & se.filt
tripMap(track[all.filt,], equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)

#To discern movement and stationary periods you can employ a changepoint analysis.
#Let's apply it to the filtered data
cl <- changeLight(twl[all.filt,], quantile = 0.95, rise.prob = NA, set.prob = NA, days = 5, plot = TRUE, summary = TRUE)
siteMap(track[all.filt,], cl$site, type = "cross", hull = T, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)

#The mergeSites() function lets you group data according to your changeLight results
#SIMEON THESE RESULTS SEEM OFF.
ms <- mergeSites(twl[all.filt,], site = cl$site, degElevation = cal1, distThreshold = 350)
siteMap(track[all.filt,], ms$site, type = "cross", hull = F, xlim = c(-90,-70), ylim = c(10,60), legend = TRUE)
lines(x=ms$summary$Lon, y=ms$summary$Lat)

