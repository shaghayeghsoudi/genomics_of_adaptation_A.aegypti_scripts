library(raster)

### read in your climate rasters (assuming you want to extract data from multiple layers)

climateFiles<-list.files(path="pathToClimateRasters",full.names=TRUE,pattern=".asc") 
# NOTE: I think they are .asc. If not you may have to modifty the pattern to the correct file extension

climate<-stack(climateFiles)

# if a single raster layer you can you use
# climate<-raster("rasterfilename")

### read in your locality data
## This assumes the locality data are in decimal degress with columns x, y 


sites<-read.csv("pathToFile/localities.csv",header=TRUE,stringsAsFactor=FALSE)

# if columns are lat long, change y+x to long+lat
coordinates(sites)<-~y+x
proj4string(sites)<-CRS("+proj=longlat +ellps=GRS80 +datum=WGS84")


### extract climate data 

values<-extract(climate,sites,cellnumbers=TRUE)


### make dataframe

sites_climate<-cbind(as.data.frame(sites),values)


### save to file

write.csv(sites_climate,"myClimateValues.csv",row.names=FALSE)