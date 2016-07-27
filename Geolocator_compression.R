#Geolocation analysis with Open Source Tools
#2016 North American Ornithological Congress, Washington D.C.

#Sponsored by: 

#Migrate Technology LLC.-- www.migratetech.co.uk

#The Cooper Ornithological Society

#The National Science Foundation

#--------------------------------------------------------------

#General data compression function for Geolocator files

geolocator_compress <- function(datetime, light, min = 0, max = 9999999999) {
  light[light > max] <- max  #trim high values to designated maximum
  light[light < min] <- min  #raise low values to designated minimum
  l1 <- light[2:(length(light)-1)] #all but first and last light levels
  l2 <- light[3:length(light)]     #light levels subsequent to each in l1
  l3 <- light[1:(length(light)-2)] #light levels previous to each in l1
  keep <- c(TRUE, l1!=l2 | l1!=l3, TRUE) #True/False vector corresponding with lt. True means there is a change in light levels (keep these rows)
  df <- data.frame(datetime = datetime[keep], light = light[keep]) #new dataframe with just the changing light levels
  return(df)
}