
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

