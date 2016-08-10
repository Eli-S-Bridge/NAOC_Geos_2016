#Geolocation analysis with Open Source Tools
#2016 North American Ornithological Congress, Washington D.C.

#Sponsored by: 

#Migrate Technology LLC.-- www.migratetech.co.uk

#The Cooper Ornithological Society

#The National Science Foundation

#--------------------------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ Yeahh I like this snake!
#--------------------------------------------------------------

## Geolight analyses.
## Based on GeoLight_2.01

# library(devtools)
# install_github("SLisovski/GeoLight", ref = "Update_2.01", force = T)

library(GeoLight)
library(TwGeos)

## For a fast and simple map of locations:
## Call up the twlight data
twl <- read.csv("data/749_twl.csv")

## make sure the datetime columns are in the correct format
twl$Twilight <- as.POSIXct(twl$Twilight, "GMT")   
  
## transformt into GeoLight format
twl <- export2GeoLight(twl)


## provide a sun angle for the light threshold value 
cal1 = -5.52

## generate tracking data
track <- coord(twl = twl, degElevation = cal1, tol = 0,
               method = "NOAA", note = TRUE)

## plot a map
tripMap(track, equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)

## You can filter outliers by applying a loess smoother to filter outlying twilights.
l.filt = loessFilter(twl, k = 3, plot = TRUE)
tripMap(track[l.filt,], equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)


 ## If you want to show a trip Map like this is is probably a good idea to remove the highly
 ## unrealistic positions around the equinox

    ## The 2011 fall equinox was Sept 23
    fall_eq <- as.POSIXct("2011-09-23")
    ## The 2012 spring equinox was March 20
    spring_eq <- as.POSIXct("2012-03-20")

    ## how big a window (days)
    window <- 25

    fall_eq_start <- fall_eq - (window/2)*24*60*60 # POSXct dates are in seconds
    fall_eq_end <- fall_eq + (window/2)*24*60*60
    
    spring_eq_start <- spring_eq - (window/2)*24*60*60
    spring_eq_end <- spring_eq + (window/2)*24*60*60

    fe.filt <- (twl$tFirst < fall_eq_start | twl$tSecond > fall_eq_end)
    se.filt <- (twl$tFirst < spring_eq_start | twl$tSecond > spring_eq_end)

    ## multiply (logical AND) all the filters together
    all.filt <- l.filt & fe.filt & se.filt
    tripMap(track[all.filt,], equinox = TRUE, xlim = c(-85,-72.5), ylim = c(15, 50), legend = TRUE)

    
    ## color code the locations over time (green to purple)
    cols <- rainbow(nrow(track), start = 0.2, end = 0.8)
    tripMap(track[all.filt,], equinox = FALSE, legend = FALSE, xlim = c(-85,-72.5), ylim = c(15, 50),
            pch = "", col = "grey50")
    points(track[all.filt,], pch = 16, col = cols[all.filt])    


    ## Latitude and Longitude over time
    opar <- par(mfrow = c(2,1), mar = c(3,4.5,1,1))
    plot(twl[,1], track[,1], type = "o", pch = 16, ylim = c(-85,-72.5), col = ifelse(all.filt, cols, "grey80"),
         xlab = "", ylab = "Longitude")
    plot(twl[,1], track[,2], type = "o", pch = 16, ylim = c(15, 50), col = ifelse(all.filt, cols, "grey80"),
         xlab = "", ylab = "Latitude")
    par(opar)
    


### Movement analysis #####    
    
## GeoLight offers tools to discern movement and stationary periods using a changepoint analysis on the twilight times.
## Forthermore, GeoLight allows you to simplify the movement path to stationary sites only.
    
## 1) changeLight
cl <- changeLight(twl[all.filt, ], quantile = 0.9, days = 5, plot = TRUE, summary = TRUE)
    siteMap(track[all.filt, ], cl$site, type = "cross", hull = T, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)

## 2) mergeSites
ms <- mergeSites(twl[all.filt,], site = cl$site, degElevation = cal1, distThreshold = 75)
  siteMap(track[all.filt,], ms$site, type = "cross", hull = F, xlim = c(-84,-75), ylim = c(20,47), legend = TRUE,
          pch = 16, cex  = 2, col = rainbow(6, start = 0.2, end = 0.8))

  arrows(x0 = c(-80.36, -78.26, -77.09, -80.85), y0 = c(42.4, 38.57, 35.77, 25.74), 
         x1 = c(-78.66, -77.13, -80.42, -80.66), y1 = c(39.23, 36.70, 24.92, 38.689))
  

## Schedule extracts the time of 'arrival' and 'departure' at the sites
## NOTE: the dates have to be treated very carefully and should only be used as a rough estimate.
schedule(twl[all.filt,1], twl[all.filt,2], ms$site)
