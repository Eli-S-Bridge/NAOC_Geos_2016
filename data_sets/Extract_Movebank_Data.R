# Geolocation analysis with Open Source Tools
# 2016 North American Ornithological Congress, Washington D.C.

## Sponsored by: 

## Migrate Technology LLC.-- www.migratetech.co.uk

## The Cooper Ornithological Society

## The National Science Foundation

#--------------------------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#--------------------------------------------------------------

library(TwGeos)

#extract data from a Movebank download
#!!! First set your working directory !!!#
#!!! Do Session -> Set Working Directory -> To Source File Location !!!#

d.mb <- read.csv("Movebank_BLPW_A_lon-72.82_lat44.53.csv")

#View data
head(d.mb)

#Reduce to Date and Light columns with correct date-time format
d.mb <- data.frame(Date = d.mb$timestamp, Light = d.mb$gls.light.level)
d.mb$Date <- strptime(d.mb$Date, "%Y-%m-%d %H:%M:%S.000", "GMT")  #vector of datestamps in a format R can work with

lightImage(d.mb, offset = 13, zlim = c(0, 64), dt = 120) 

#save the data
write.csv(d.mb, file = "BLPW_A.csv", quote = FALSE, row.names = FALSE)


