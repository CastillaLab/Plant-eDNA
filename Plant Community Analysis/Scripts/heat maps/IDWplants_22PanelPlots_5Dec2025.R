#IDW script for plant dataset - JDR November 2025
#Plot surfaces for all lakes
#All plants, terrestrial, aquatic, wetland plants
setwd("~/Desktop/PLOTS-22panel_5Dec2025")

#Load necessary libraries
library(sf)
library(raster)
library(prevR)
library(colorRamps)
library(prettymapr)  #!#Added 5/3/21 for scale bar and north arrow
library(scales)

#Read in the Lake polygons
lakes <- st_read("~/Downloads/Lake_Polygons", layer = "Lake_Polygons")

#!#Adding rivers to the plots - JDR 5/3/21
#!# Downloaded from https://gis-michigan.opendata.arcgis.com/datasets/hydrography-lines-v17a?geometry=-85.547%2C45.183%2C-84.549%2C45.352
rivs <- st_read("~/Downloads/Hydrography_Lines_(v17a)")
rivs <- st_transform(rivs, crs=crs(lakes))

#Load in the per-sample data for species richness - total, aquatic, wetland, terrestrial
srdat <- read.csv("SRcounts_byHab.csv", header = TRUE)

#List of lakes... make this a loop
mylake <- c("Fourth",
			"Ocqueoc",
			"Haithco",
			"Dumont",
			"Holloway",
			"Pentwater",
			"Walloon",
			"Five.Channels",
			"George",
			"Mullett",
			"Long",
			"Cass",
			"Austin",
			"Thompson",
			"Wycamp",
			"Torch",
			"Manistique",
			"Brev",
			"Houghton",
			"Higgins",
			"Kimball",
			"Pickerel")


#Loop over the variables (habitats)
names(srdat)[7:10] #These are the columns we want to plot...

for(var in 7:10) {
pdf(paste0("IDW_SR-",names(srdat)[var],".pdf"), width=15, height=15)
par(mfrow = c(5,5), mar=c(3,3,3,3))

  #Loop over the lakes
  for(lk in mylake) {
	  print(lk)
	  #ID the polygon in the shapefile associated with the lake of interest
	  thisone <- c()
	  tmptest <- vector("list",length(lakes))
	  tmpdat <- srdat[grep(lk, srdat$lake_name),]
	  sampCoords <- st_as_sf(data.frame(lon=tmpdat$longitude, lat=tmpdat$latitude), 
	                         coords=c("lon","lat"), crs=4326)
	  sampCoords <- st_transform(sampCoords, crs(lakes))
	  for(x in 1:nrow(lakes)) {
	    tmptest[[x]] <- st_intersects(lakes[x,], sampCoords, sparse=FALSE)
		  if(sum(tmptest[[x]]) > 0) {
			  thisone <- c(thisone,x)
		  }
	  }

	  rm(tmptest, tmpdat)
	  tmplake <- lakes[thisone,]

		keepfeat <- st_intersects(rivs, st_buffer(tmplake, dist=10), sparse = FALSE)
	  if(sum(keepfeat) > 0) {
	    tmprivs <- rivs[which(keepfeat == TRUE),]
	    tmprivs <- tmprivs[which(tmprivs$FCC != "H21"),]
	  } else {
	    tmprivs <- NULL
  	}
	
	  dat <- srdat[grep(lk, srdat$lake_name),]
	
	  #Reduce name length for Five Channels
	  if(dat$lake_name[1] == "Five.Channels.Impoundment") {
		  dat$lake_name <- "Five.Channels"
	  }

	  dat <- dat[!is.na(dat$longitude),]

	  lakedat <- data.frame(long = dat$longitude, lat = dat$latitude, sr = dat[,var])
	  lakedat <- st_as_sf(lakedat, coords=c("long","lat"), crs=4326)
	  lakedat <- st_transform(lakedat, crs(tmplake))
	  grid <- st_sample(tmplake, size=100000, type="regular")
		
	  #Using an IDW power of 2 for all surfaces
	  finalPOW <- 2

	  #Fit the IDW model
	  idwmodel <- idw(lakedat$sr~1, locations = lakedat, 
	                  newdata = grid, idp = finalPOW, 
	                  maxdist = Inf)

	  #Plot the lake and surface
	  plot(st_geometry(tmplake), main = gsub("\\.", " ", dat$lake_name[1]))
	  my.bbox <- st_bbox(st_transform(tmplake, crs=4326))
	  xyats <- st_as_sf(
	    data.frame(lat=seq(my.bbox["xmin"], my.bbox["xmax"], length.out = 3),
	               lon=seq(my.bbox["ymin"], my.bbox["ymax"], length.out = 3)), 
	    coords = c("lat","lon"), crs=4326)
	  CRSxyats <- st_transform(xyats, crs(tmplake))
	  axis(1, at = st_coordinates(CRSxyats)[,1], labels = round(st_coordinates(xyats)[,1], 3)) 
	  axis(2, at = st_coordinates(CRSxyats)[,2], labels = round(st_coordinates(xyats)[,2], 3))
	  if(length(tmprivs) > 0) {
	    plot(st_geometry(tmprivs), add = TRUE, col = "blue", lwd = 2)
	  }
	  idwDF <- as.data.frame(st_coordinates(idwmodel["var1.pred"]))
	  idwDF$Z <- idwmodel["var1.pred"]$var1.pred
	  rasIDW <- rasterFromXYZ(idwDF)
	  col <- matlab.like(50)
	  plot(rasIDW, col = col, add = TRUE, alpha =0.75, legend.shrink = 0.4, legend.mar = 4.5, legend.cex = 0.5)
	  if(length(which(tmplake$LAKE_TYPE == "Island")) > 0) {
	    plot(st_geometry(tmplake[which(tmplake$LAKE_TYPE == "Island"),]), 
	                    add=TRUE, col = "white")
	  }
	  points(lakedat, pch=19, cex=0.65, col=alpha("black",0.65))
	
	  #!#Adding North arrow & scale bar - JDR 5/3/21
	  addnortharrow("topright", scale = 0.4)
	  addscalebar(htin=0.075, label.cex = 0.75)
	
	  rm(CRSxyats, col, idwmodel, idwDF, rasIDW, lakedat, finalPOW, grid)
  }
	dev.off()
}

