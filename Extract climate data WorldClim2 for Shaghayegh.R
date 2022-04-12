#download the dataset you need from WorldClim 2 (https://worldclim.org/version2)
#they are zipped .tif data (basically big maps of the world with each pixel being a data point)
#other datasets have .asc files, but the approach is the same
#you first make the data accessible with "raster" and then extract the data points you need with "extract" - is that easy once you figure it out :)


install.packages("raster")
install.packages("sp")
install.packages("rgdal")

library(raster)
library(sp)
library(rgdal)


#load you population data file
ALL_pops <- read.csv("Path/to/Population_info.csv", header=TRUE)

#extract only longitude and latitude data (obviously the factors need to be called "Longitude" and "Latitude" in your original dataset)
ALL_latlon <- ALL_pops[,c('Longitude','Latitude')]

#this converts you coordinate to spatial points - that doesn't really seem to change anything, so I skipped it, but you might want to double check it
#ALL_latlon_converted <- SpatialPoints(ALL_latlon, proj4string = r_Tmean_April@crs)

#there are 12 files for each variables, one for each month. So if you wan to extract average T for April and May...
r_Tmean_April <- raster("wc2.0_30s_tavg/wc2.0_30s_tavg_04.tif")
r_Tmean_May <- raster("wc2.0_30s_tavg/wc2.0_30s_tavg_05.tif")

#precipitation
r_prec_April <- raster("wc2.0_30s_prec/wc2.0_30s_prec_04.tif")
r_prec_May <- raster("wc2.0_30s_prec/wc2.0_30s_prec_05.tif")

#radiation
r_Radiation_April <- raster("wc2.0_30s_srad/wc2.0_30s_srad_04.tif")
r_Radiation_May <- raster("wc2.0_30s_srad/wc2.0_30s_srad_05.tif")

#and so on...
#fianlly extract the values for your coordinates
avTemp_April <- extract(r_Tmean_April,ALL_latlon)
avTemp_May <- extract(r_Tmean_May,ALL_latlon)

prec_April <- extract(r_prec_April,ALL_latlon)
prec_May <- extract(r_prec_May,ALL_latlon)

radiation_April <- extract(r_Radiation_April,ALL_latlon)
radiation_May <- extract(r_Radiation_May,ALL_latlon)


#and bind them all together (PopID being teh column in your initial dataset that contains the population identifier)
ALLpops_climate <- cbind.data.frame(ALL_pops$PopID,A,avTemp_April, avTemp_May, prec_April, prec_May, radiation_April, radiation_May)

write.csv(ALLpops_climate, "ALLpops_climate.csv")
