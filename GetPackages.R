# Check to make sure the required packages are installed on your machine
# If not, they will be installed

reqPackages <- c("devtools","GeoLight","rages","geosphere","raster","fields","forecast",
                  "circular","truncnom","parallel","bit","rgdal","CircStats","Rcpp","RcppArmadillo",
				  "ggmap","ggsn","sp","maptools","rgeos","MASS")


get.packages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]

if(length(get.packages)>0) install.packages(get.packages,repos = "https://cloud.r-project.org/")

# Install necessary packages from Github using the devtools library #

library(devtools)
install_github("SWotherspoon/SGAT")
install_github("SLisovski/TwGeos")
install_github("SWotherspoon/BAStag") 
install_github("eldarrak/FLightR")


