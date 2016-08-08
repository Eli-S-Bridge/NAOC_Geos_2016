#Geolocation analysis with Open Source Tools
#2016 North American Ornithological Congress, Washington D.C.

#Sponsored by: 

#Migrate Technology LLC.-- www.migratetech.co.uk

#The Cooper Ornithological Society

#The National Science Foundation

#--------------------------------------------------------------

#SGAT analyses.

library(SGAT)
library(MASS)  #needed for fitting distributions
library(BAStag)
library(GeoLight)

#Start by getting a fast and simple map of locations:
#Call up the twlight data
                             #open the BAStag package
d.lig <- readLig("data/749_000.lig", skip = 0)         #read the data into a dataframe called d.lig
d.lig <- subset(d.lig,select=c("Date","Light"))   #reduce the dataframe to just Date and Light
d.lig$Date  <- as.POSIXct(d.lig$Date, "GMT")

#We've already defined twilights... 
twl <- read.csv("data/749_twl.csv")
#...but they need to be reformatted
twl <- data.frame(Twilight = as.POSIXct(twl$tFirst, "GMT"), Rise = twl$type %% 2, stringsAsFactors = F) #use modulus operator (%%) to transform 2 to 0 and 1 to 1

lightImage(d.lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

#Let's truncate the data using point and click
trnc <- as.POSIXct(locator(n=2)$x, origin = "1970-01-01")

#view the time points selected to make sure they make sense
trnc

#Subset the twlight data.
twl <- twl[twl$Twilight > trnc[1] & twl$Twilight < trnc[2],]

#View the data again and select a calibration period
lightImage(d.lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))
tm.calib <- as.POSIXct(locator(n=2)$x, origin = "1970-01-01")

### See results of your selection
tm.calib

#Calibration data subset
d.calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])

# Calibration Coordinates
lon.calib <- -80.46
lat.calib <- 42.62

sun  <- solar(d.calib[,1]) #calculate some solar parameters
z    <- refracted(zenith(sun, lon.calib, lat.calib)) #adjust solar zenigh angles for atmispheric refraction

#To use SGAT, we need to provide some paramters that describe
#the expected error in the twilight times. We can use the calibration
#data to estimate these parameters.
#First calculate actual twilights based on date and location
twl_t <- twilight(d.calib[,1], lon.calib, lat.calib, rise = d.calib[,2], zenith = max(z)+0.25) 
#Then get the distribution of differences (errors) between the calculated and estimated twilights 
twl_dev <- ifelse(d.calib$Rise, as.numeric(difftime(d.calib[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, d.calib[,1], units = "mins")))

#View the distribution of errors
hist(twl_dev, freq = F, ylim = c(0, 0.1), xlim = c(0, 50), main = "", xlab = "Twilight error (mins)")

#Fit a log-Normal distribution to these errors.
fitml <- fitdistr(twl_dev, "log-Normal") 

#plot the fitted distribution
lines(1:100, dlnorm(1:100, fitml$estimate[1], fitml$estimate[2]), col = "firebrick", lwd = 3, lty = 2)

#If it looks good save the mean and standard deviation for the model
alpha <- c(fitml$estimate[1], fitml$estimate[2]) ## Twilight model parameters


#SIMEON
#Need some help here. The zenith angle needs to be about 97 to get reasonable results
zenith0 <- median(z) #used only for simple threshold map
zenith0
zenith  <- quantile(z, prob = .05)  #use the 95% of the zenith angles. Why??? This is too low a value
zenith <- zenith0  #still too low a value
zenith <- 97  #Maybe a bit too high. 
zenith <- median(z)

path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol = 0.15)
#have a quick look to see that it makes sense.
tripMap(path$x, equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)

#We can use this initial path to define a prior set of 
#coordinates for the SGAT model to start from.
#SGAT estimates locations for noon and midnight (x0) 
#as well as intermediate locations (z0).
x0 <- path$x
z0 <- trackMidpts(x0)

#Now you have to specify a movement model or a distribution of movement speeds
#This is not an ideal model. You just want to specify something realistic for 
#your species. A gamma model is good for specifying a high frequency of short
#movements and rare long-distance movements.
beta = c(0.7, 0.05)  #mean and standard deviation for your model
#plot a gamma distribution
plot(0:80, dgamma(0:80, beta[1], beta[2]), type = "l", col="red",
     xlab = "Speed (km/h)", ylab = "Probability")

#Establish fixed positions for start and end (if both are known)
fixedx <- rep(F, nrow(x0))            #Just a list of "FALSE" the same length as x0
fixedx[1:3] <- T                      #Make the first three "TRUE."
fixedx[(nrow(x0)-2):nrow(x0)] <- T    #Make the last three "TRUE."
x0[fixedx, 1] <- lon.calib            #fix the longitude for the known locations
x0[fixedx, 2] <- lat.calib            #fix the latitude for the known locations

z0 <- trackMidpts(x0) # update z0 positions


#Finally! THE MODEL!!

#The threshold.model function requires the following

  # the twilight times and whether they are sunrise or sunset
  # a modified Log Normal model for the twilight errors
  # the parameters of the distribution of twilight errors (alpha)
  # the parameters of the distribution of travel speed (beta)
  # the initial x and z, and
  # the zenith angle that defines twilight.

model <- thresholdModel(twl$Twilight, twl$Rise,
                        twilight.model="ModifiedLogNormal",
                        alpha=alpha, beta=beta,
                        x0=x0,z0=z0,zenith=zenith,fixedx=fixedx)


proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))  #specify a multivariate normal distribution for the sampler
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

# The function estelle.metropolis draws samples from the posterior distribution of the
# locations defined by the model. The sampler is run with the modified model
# until all the constraints of the full model are met.

#Do a short run with a "forgiving" model to get intial x and z data that make sense
fit <- estelleMetropolis(model,proposal.x,proposal.z,iters=200,thin=20,chains=1)

#The result is a multidimentional list that contains model parameters, x estimates, and z estimates.

#The next step is to tune the proposal distributions. 
#The model and proposals are redefined using the last set of locations 
#from the previous run to initialize

x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

#Note that the model uses a not-so-forgiving LogNormal distribution now.
model <- thresholdModel(twl$Twilight, twl$Rise, twilight.model = "LogNormal", 
                        alpha = alpha, beta = beta, x0 = x0, z0 = z0, zenith = zenith, 
                        fixedx = fixedx)

proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))


fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         iters=300,thin=20,chains=2)

#Perform tuning three more times
for(k in 1:3) {
  proposal.x <- mvnorm(chainCov(fit$x),s=0.2)
  proposal.z <- mvnorm(chainCov(fit$z),s=0.2)
  fit <- estelleMetropolis(model,proposal.x,proposal.z,
                           x0=chainLast(fit$x),
                           z0=chainLast(fit$z),
                           iters=300,thin=20,chains=2)
}

#Final Run
#Now draw a larger sample
#We can now define the samplers based on covariance estimates 
proposal.x <- mvnorm(chainCov(fit$x),s=0.25)
proposal.z <- mvnorm(chainCov(fit$z),s=0.25)

fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         x0=chainLast(fit$x),
                         z0=chainLast(fit$z),
                         iters=2000,thin=20,chains=2)

#summarize the location data in fit
s <- locationSummary(fit$x,time=model$time,collapse=T)
#save(s, file = "data/SGAT_results1.Rdata")
load(file = "data/SGAT_results1.Rdata")

#Have a look at the results
library(maptools)
data(wrld_simpl)

#limits for map drawing
xlim = c(-90,-70)
ylim = c(10,50)

fixedz <- fixedx[-length(fixedx)] > 0 & fixedx[-length(fixedx)]==fixedx[-1]
dt <- ifelse(fixedz,0,model$dt)
im <- locationImage(fit$z,xlim=xlim,ylim=ylim,nx=4*diff(xlim),ny=4*diff(ylim),
                    weight=dt,collapse=TRUE)

opar <- par(mar=c(2,2,2,2)+0.1)
plot(wrld_simpl,col= "grey90",border="grey10" , xlim = range(im$x), ylim = range(im$y))
image(im$x,im$y,im$W,xlab="",ylab="",cex.axis=0.7, add = T, col = c("transparent", rev(topo.colors(200))))

plot(wrld_simpl, add = T)

lines(s$Lon.mean, s$Lat.mean,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.5))
points(s$Lon.mean, s$Lat.mean,pch=16,cex=0.5,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.2))
box()


#Let's do it again but with a land mask. Contrain the birds to land surfaces.

#Create a function for the land mask

land.mask <- function(xlim, ylim, n = 4, land = TRUE) {
  #Create a raster grid of "NA" that fills the x and y limits. Use a lat/lon projections
  #resolution is 1/n degrees
  r <- raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1], 
              xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  #Replace the NA values with 1 where there is land
  #This code is kind of complicated because it allows for grids that wrap around the date line
  r <- cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE), 
             rasterize(wrld_simpl, r, 1, silent = TRUE), 
             rasterize(elide(wrld_simpl, shift = c(360, 0)), r, 1, silent = TRUE))
 
   #make the raster a matrix with column order reversed and NA set to TRUE
  r <- as.matrix(is.na(r))[nrow(r):1, ]

  if (land) #reverse the TRUE/FALSE if land is set to TRUE
    r <- !r
  
  #Define all the x and y bins for the matrix so you can look up a particular value
  xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
  ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)
  
  #This function just fits p into the appropriat bins in the grid and returns TRUE or FALSE. 
  function(p) {
    r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
  }
}

#Now use the function above to create a particular mask with defined extent and land set to T or F.
is.land <- land.mask(xlim = xlim, ylim = ylim, n = 4, land = T)

#Try some coordinates in the carribean sea
is.land(cbind(-80, 20))

#Try some coordinates in the eastern US
is.land(cbind(-80, 40))

#Try some coordinates outside the x and y limits
is.land(cbind(-60, 40))

#The location estimates derived by Estelle can effectively excluded 
#from the land by imposing a prior on the x (and z) locations so that 
#locations on the land have a vanishingly small probability of occurrence. 
#The prior is defined on the log scale. Here, we don’t want to exlude 
#them but give location estimates on land a higher prior.

#Make another function that define the log prior for x and z coordinates
log.prior <- function(p) {
  f <- is.land(p)                       # f will be TRUE, FALSE, or NA.
  ifelse(f | is.na(f), log(2), log(1))  # If f is TRUE or NA, then the prior is log(2) (note that off-grid is treated like land),
}                                       # otherwise it;s log(1) or 0. You can change these values as needed


#INITIALIZATION
x0 <- path$x
fixedx <- rep(F, nrow(x0))
fixedx[1:2] <- T
fixedx[nrow(x0) - 1] <- T
fixedx[nrow(x0)] <- T
x0[fixedx, 1] <- lon.calib
x0[fixedx, 2] <- lat.calib
z0 <- trackMidpts(x0) # update z0 positions

model <- thresholdModel(twl$Twilight,twl$Rise, 
                        twilight.model="ModifiedLogNormal",
                        alpha=alpha,beta=beta,
                        logp.x = log.prior, logp.z = log.prior, #mask is applied here
                        x0=x0,z0=z0,zenith=zenith,fixedx=fixedx)

proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

fit <- estelleMetropolis(model,proposal.x,proposal.z,iters=200,thin=20,chains=1)

#TUNING

x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- thresholdModel(twl$Twilight, twl$Rise, twilight.model = "LogNormal", 
                        alpha = alpha, beta = beta, x0 = x0, z0 = z0, zenith = zenith, 
                        logp.x = log.prior, logp.z = log.prior,
                        fixedx = fixedx)

proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         iters=300,thin=20,chains=2)

for(k in 1:3) {
  proposal.x <- mvnorm(chainCov(fit$x),s=0.2)
  proposal.z <- mvnorm(chainCov(fit$z),s=0.2)
  fit <- estelleMetropolis(model,proposal.x,proposal.z,
                           x0=chainLast(fit$x),
                           z0=chainLast(fit$z),
                           iters=300,thin=20,chains=2)
}

#FINAL RUN

proposal.x <- mvnorm(chainCov(fit$x),s=0.25)
proposal.z <- mvnorm(chainCov(fit$z),s=0.25)

fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         x0=chainLast(fit$x),
                         z0=chainLast(fit$z),
                         iters=2000,thin=20,chains=2)

#Let's view the data month by month as 2D histograms.
#First create an empty background raster
r <- raster(nrows = 80, ncols = 100, xmn = xlim[1], xmx = xlim[2], ymn = ylim[1], ymx = ylim[2])
#The slices function in SGAT can cut your data up into days, weeks, months, etc.
s <- slices(type = "primary", breaks = "month", mcmc = fit, grid = r)
for (k in sliceIndices(s)) {   #loop through the slices to plot the data
  tm <- sliceInterval(s, k)    
  sk <- slice(s, k)
  plot(wrld_simpl, col = "grey95", xlim = xlim, ylim = ylim)
  plot(sk, add = T, col = heat.colors(64, alpha = 0.7), useRaster = F)
  axis(1)
  axis(2)
  box()
  title(sprintf("%s - %s", as.Date(tm[1]), as.Date(tm[2])), cex.main = 0.8)
}

#Just plot the mean locations (medians might be better when a mask is employed)
s <- locationSummary(fit$x,time=model$time,collapse=T)
plot(wrld_simpl, col = "grey95", xlim = xlim, ylim = ylim)
lines(s$Lon.mean, s$Lat.mean,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.5))
points(s$Lon.mean, s$Lat.mean,pch=16,cex=0.5,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.2))


# What if we know the birds are stationary from 1 December to 1 February?
# It is possible to do a grouped analysis and specify that there is no movement 
# during this period.

stat <- which(twl$Twilight > as.POSIXct("2011-12-01 00:00:01", "GMT") & twl$Twilight < as.POSIXct("2012-02-01 00:00:01", "GMT"))
min(stat) ; max(stat)
grp <- c(1:(min(stat)-1), rep(min(stat),length(stat)), 
         seq(from = min(stat)+1, by = 1, length.out = nrow(twl)-(max(stat))))

#INITIALIZATION
x0 <- path$x
fixedx <- rep(F, nrow(x0))
fixedx[1:2] <- T
fixedx[nrow(x0) - 1] <- T
fixedx[nrow(x0)] <- T
x0[fixedx, 1] <- lon.calib
x0[fixedx, 2] <- lat.calib
z0 <- trackMidpts(x0) # update z0 positions



#night only grouping
grp2 <- rep(1:nrow(twl), each = 2, length.out = nrow(twl) + 1)[-1]
#simplest grouping (no grouping at all)
grp3 <- 1:658

model <- groupedThresholdModel(twl$Twilight,twl$Rise, group = grp,
                        twilight.model="ModifiedLogNormal",
                        alpha=alpha,beta=beta,
                        logp.x = log.prior, logp.z = log.prior, #mask is applied here
                        x0=x0,z0=z0,zenith=zenith,fixedx=fixedx)

proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# THIS FAILS!!!
fit <- estelleMetropolis(model,proposal.x,proposal.z,iters=200,thin=20,chains=1)







#TUNING

x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- groupedThresholdModel(twl$Twilight,twl$Rise, group = grp,
                               twilight.model="LogNormal",
                               alpha=alpha,beta=beta,
                               logp.x = log.prior, logp.z = log.prior, #mask is applied here
                               x0=x0,z0=z0,zenith=zenith,fixedx=fixedx)

proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         iters=300,thin=20,chains=2)

for(k in 1:3) {
  proposal.x <- mvnorm(chainCov(fit$x),s=0.2)
  proposal.z <- mvnorm(chainCov(fit$z),s=0.2)
  fit <- estelleMetropolis(model,proposal.x,proposal.z,
                           x0=chainLast(fit$x),
                           z0=chainLast(fit$z),
                           iters=300,thin=20,chains=2)
}

#FINAL RUN

proposal.x <- mvnorm(chainCov(fit$x),s=0.25)
proposal.z <- mvnorm(chainCov(fit$z),s=0.25)

fit <- estelleMetropolis(model,proposal.x,proposal.z,
                         x0=chainLast(fit$x),
                         z0=chainLast(fit$z),
                         iters=2000,thin=20,chains=2)

#Other possible model components and constraints?

# Missing twilight data can be specified
#    missing	 <- integer vector indicating which twilights were unobserved and why.
# Available travel time can be specified
#    dt  <-  time intervals for speed calculation in hours.
# Twlights can be grouped to specify stationary periods
#    group	<- integer vector that defines the twilight groups. If codegroup[k]==j then the k-th twilight occurs at location j.

