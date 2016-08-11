## Geolocation analysis with Open Source Tools
## 2016 North American Ornithological Congress, Washington D.C.

## Sponsored by: 

## Migrate Technology LLC.-- www.migratetech.co.uk

## The Cooper Ornithological Society

## The National Science Foundation

#--------------------------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ 
#--------------------------------------------------------------

## SGAT analysis WITHOUT a land mask.

library(SGAT)
library(MASS)  #needed for fitting distributions
library(TwGeos)
library(maptools)
  data(wrld_simpl)
# library(GeoLight)

## Start by getting a fast and simple map of locations:
## Call up the twlight data

d.lig <- readLig("data/749_000.lig", skip = 0)    # read the data into a dataframe called d.lig
d.lig <- subset(d.lig,select=c("Date","Light"))   # reduce the dataframe to just Date and Light
d.lig$Date  <- as.POSIXct(d.lig$Date, "GMT")

## We've already defined twilights... 
twl <- read.csv("data/749_twl.csv", col.names = c("Twilight", "Rise"))
twl$Twilight <- as.POSIXct(twl$Twilight, "GMT")

lightImage(d.lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

## Let's truncate the data using point and click
# trnc <- as.POSIXct(locator(n=2)$x, origin = "1970-01-01", tz = "GMT")
trnc <- as.POSIXct(c("2011-06-26", "2012-05-11"), tz = "GMT")
abline(v  = trnc, col  ="orange")

## Subset the twlight data.
twl <- twl[twl$Twilight > trnc[1] & twl$Twilight < trnc[2],]


### Select Calibration Periods ###

## To see which periods correspond to calibration sites, you can add sunrise/sunset times
## of the known position.

# Calibration Coordinates
lon.calib <- -80.46
lat.calib <- 42.62

tm <- seq(twl$Twilight[1], twl$Twilight[nrow(twl)], by = "day")
rise <- rep(c(TRUE, FALSE), length(tm))

c.dat <- data.frame(Twilight = twilight(rep(tm, each = 2), lon = lon.calib, lat = lat.calib, 
                                        rise = rise, zenith = 96), Rise = rise)

lightImage(d.lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))
tsimagePoints(c.dat$Twilight, offset = 19, pch = 16, cex = 0.25,
              col = "orange")


## Select one or more periods
# tm.calib1 <- as.POSIXct(locator(n=2)$x, origin = "1970-01-01", tz = "GMT")
tm.calib1 <- as.POSIXct(c("2011-06-26", "2011-07-12"), tz = "GMT") # results of previous line
# tm.calib2 <- as.POSIXct(locator(n=2)$x, origin = "1970-01-01", tz = "GMT")
tm.calib2 <- as.POSIXct(c("2012-04-08", "2012-05-08"), tz = "GMT")

abline(v = c(tm.calib1, tm.calib2), col  ="orange")

## Calibration data subset
d.calib <- subset(twl, (Twilight>=tm.calib1[1] & Twilight<=tm.calib1[2]) |
                       (Twilight>=tm.calib2[1] & Twilight<=tm.calib2[2]))


### Calibration (The twiligth model) ####

sun  <- solar(d.calib[,1]) #calculate some solar parameters
z    <- refracted(zenith(sun, lon.calib, lat.calib)) #adjust solar zenigh angles for atmispheric refraction

## To use SGAT, we need to provide some paramters that describe
## the expected error in the twilight times. We can use the calibration
## data to estimate these parameters.
## First calculate actual twilights based on date and location
twl_t <- twilight(d.calib[,1], lon.calib, lat.calib, rise = d.calib[,2], zenith = max(z)+0.15) 
## Then get the distribution of differences (errors) between the calculated and estimated twilights 
twl_dev <- ifelse(d.calib$Rise, as.numeric(difftime(d.calib[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, d.calib[,1], units = "mins")))

## View the distribution of errors
hist(twl_dev, freq = F, xlim = c(0, 50), breaks = 15, main = "", xlab = "Twilight error (mins)")

## Fit a log-Normal distribution to these errors.
fitml <- fitdistr(twl_dev, "log-Normal") 

## plot the fitted distribution
lines(1:50, dlnorm(1:50, fitml$estimate[1], fitml$estimate[2]), col = "firebrick", lwd = 3, lty = 2)

## If it looks good save the mean and standard deviation for the model
alpha <- c(fitml$estimate[1], fitml$estimate[2]) ## Twilight model parameters

## derive a zenth angle for sunrise (angle in degrees measured from directly overhead)
zenith0 <- median(z) ## used only for simple threshold map
zenith0
zenith  <- quantile(z, prob = 1)  #use the ~90-100% of the zenith angles.
zenith 

path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith0, tol = 0.15)
## have a quick look to see that it makes sense.
tripMap(path$x, equinox = TRUE, xlim = c(-85,-72.5), ylim = c(15,45), legend = TRUE)
abline(h = lat.calib, lty = 2)

## We can use this initial path to define a prior set of 
## coordinates for the SGAT model to start from.
## SGAT estimates locations for sunrise and sunset (x0) 
## as well as intermediate locations (z0).
x0 <- path$x
z0 <- trackMidpts(x0)

## Now you have to specify a movement model or a distribution of movement speeds
## This is not an ideal model for e.g. passerine migration. You just want to specify 
## something realistic for your species. A gamma model is good for specifying a 
## high frequency of short movements and rare long-distance movements (the distribution
## needs to reflect the speeds for the entire track!).

beta = c(0.45, 0.05)  # mean and standard deviation for your model
## plot a gamma distribution
plot(0:80, dgamma(0:80, beta[1], beta[2]), type = "l", col="red",
     xlab = "Speed (km/h)", ylab = "Probability")


## Establish fixed positions for start and end (if both are known)
fixedx <- rep(F, nrow(x0))                     # Just a list of "FALSE" the same length as x0
fixedx[c(1:3, (nrow(x0)-2):nrow(x0))] <- T     # Make the first three and last three "TRUE."
x0[fixedx, 1] <- lon.calib                     # fix the longitude for the known locations
x0[fixedx, 2] <- lat.calib                     # fix the latitude for the known locations

z0 <- trackMidpts(x0) # update z0 positions



#### Finally! THE ESTELLE MODEL!! ####

## The threshold.model function requires the following

  # the twilight times and whether they are sunrise or sunset
  # a modified Log Normal model for the twilight errors (with relaxed assumptions to tune the proposals)
  # the parameters of the distribution of twilight errors (alpha)
  # the parameters of the distribution of travel speed (beta)
  # optional: the spatial mask
  # the initial x and z, and
  # the zenith angle that defines twilight.

model <- thresholdModel(twl$Twilight, twl$Rise,
                        twilight.model = "ModifiedLogNormal",
                        alpha = alpha, beta = beta,
                        x0 = x0,z0 = z0,zenith = zenith,fixedx = fixedx)

proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))  #specify a multivariate normal distribution for the sampler
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

## The function estelle.metropolis draws samples from the posterior distribution of the
## locations defined by the model. The sampler is run with the modified model
## until all the constraints of the full model are met.

## Do a short run with a "forgiving" model to get intial x and z data that make sense
fit <- estelleMetropolis(model,proposal.x,proposal.z,iters=500,thin=20,chains=1)

## The result is a multidimentional list that contains model parameters, x estimates, and z estimates.

## The next step is to tune the proposal distributions. 
## The model and proposals are redefined using the last set of locations 
## from the previous run to initialize

x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

## Note that the model uses a not-so-forgiving LogNormal distribution now.
model <- thresholdModel(twl$Twilight, twl$Rise, twilight.model = "LogNormal",
                        alpha = alpha, beta = beta, x0 = x0, z0 = z0, zenith = zenith, 
                        fixedx = fixedx)

proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))


fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         iters=300,thin=20,chains=1)

## Perform tuning three more times
for(k in 1:3) {
  proposal.x <- mvnorm(chainCov(fit$x),s=0.2)
  proposal.z <- mvnorm(chainCov(fit$z),s=0.2)
  fit <- estelleMetropolis(model,proposal.x,proposal.z,
                           x0=chainLast(fit$x),
                           z0=chainLast(fit$z),
                           iters=300,thin=20,chains=1)
}

## Final Run
## Now draw a larger sample
## We can now define the samplers based on covariance estimates 
proposal.x <- mvnorm(chainCov(fit$x),s=0.25)
proposal.z <- mvnorm(chainCov(fit$z),s=0.25)

fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         x0=chainLast(fit$x),
                         z0=chainLast(fit$z),
                         iters=2000,thin=20,chains=1)



## summarize the location data in fit
s <- locationSummary(fit$z,time=model$time,collapse=T)
# save(s, file = "data/SGAT_results1.Rdata")
# load(file = "data/SGAT_results1.Rdata")

fixedz <- fixedx[-length(fixedx)] > 0 & fixedx[-length(fixedx)]==fixedx[-1]
dt <- ifelse(fixedz,0,model$dt)
im <- locationImage(fit$z,xlim=xlim,ylim=ylim,nx=4*diff(xlim),ny=4*diff(ylim),
                    weight=dt,collapse=TRUE)

opar <- par(mar=c(2,2,2,2)+0.1)
plot(NA, xlim = range(im$x), ylim = range(im$y), xlab = "Longitude", ylab = "Latitude")
plot(wrld_simpl,col= "grey90",border="grey10", add = T)
image(im$x,im$y,im$W,xlab="",ylab="",cex.axis=0.7, add = T, col = c("transparent", rev(topo.colors(200))))

plot(wrld_simpl, add = T)

lines(s$Lon.mean, s$Lat.mean,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.5))
points(s$Lon.mean, s$Lat.mean,pch=16,cex=0.5,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.2))
points(lon.calib, lat.calib, pch = 22, bg = "white")
par(opar)


## Let's view the data by month as 2D histograms.
## First create an empty background raster
r <- raster(nrows = 80, ncols = 100, xmn = xlim[1], xmx = xlim[2], ymn = ylim[1], ymx = ylim[2])
## The slices function in SGAT can cut your data up into days, weeks, months, etc.
s <- slices(type = "intermediate", breaks = "month", mcmc = fit, grid = r)

opar <- par(mfrow = c(4,3), mar = c(0,0,0,0), oma = c(0.4,0.4,0.4,0.4))
for (k in sliceIndices(s)) {   #loop through the slices to plot the data
  tm <- sliceInterval(s, k)    
  sk <- slice(s, k)
  plot(NA, xlim = range(im$x), ylim = range(im$y), xaxt = "n", yaxt = "n", xlab = "", ylab  = "")
  plot(wrld_simpl, col = "grey95", add = T)
  plot(sk, add = T, col = rev(heat.colors(64, alpha = 0.7)), useRaster = F, legend = F)

  mtext(sprintf("%s - %s", as.Date(tm[1]), as.Date(tm[2])), 1, line = -2, cex = 0.8)
}
par(opar)

## Just plot the mean locations (medians might be better when a mask is employed)
plot(s$Lon.mean, s$Lat.mean, type = "n", xlab = "Longitude", ylab = "Latitude")
plot(wrld_simpl, col = "grey95", add = T)
lines(s$Lon.mean, s$Lat.mean,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.5))
points(s$Lon.mean, s$Lat.mean,pch=16,cex=0.5,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.2))