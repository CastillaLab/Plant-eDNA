#IDW script for plant dataset - JDR 3/24/2025
#Plot surfaces for all lakes
#All plants, terrestrial, aquatic, wetland plants
setwd("~/Downloads/PLOTS-22panel")

#Load necessary libraries
library(rgdal)
library(sp)
library(raster)
library(gstat)
library(prevR)
library(colorRamps)
library(prettymapr)  #!#Added 5/3/21 for scale bar and north arrow
library(rgeos)  #!# Added 5/3/21 for river additions

#Read in the Lake polygons
#Uses readOGR from the 'rgdal' package
#lakes <- readOGR(dsn="~/Google Drive/My Drive/MSU/GLRI_eDNA/GLRI_2017/ACTUAL/FIELD/Sampling/Lake_Polygons", layer = "Lake_Polygons")
lakes <- readOGR(dsn="~/Downloads/Lake_Polygons", layer="Lake_Polygons")
#!#Adding rivers to the plots - JDR 5/3/21
#!# Downloaded from https://gis-michigan.opendata.arcgis.com/datasets/hydrography-lines-v17a?geometry=-85.547%2C45.183%2C-84.549%2C45.352
#rivs <- readOGR(dsn="~/Google Drive/My Drive/MSU/GLRI_eDNA/GLRI_2017/MANUSCRIPT/fish/heatmaps/REVISED_May2021/Hydrography_Lines_(v17a)")
rivs <- readOGR(dsn="~/Downloads/Hydrography_Lines_(v17a)")

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

#Loop over the lakes
#for(lk in mylake) {
#	pdf(paste0("IDW_alienSR-",lk,".pdf"), width = 14, height = 3.5)
#	par(mfrow = c(1,4), mar = c(3,3,3,3))


#Loop over the variables (habitats)
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
	for(x in 1:length(lakes)) {
		tmptest[[x]] <- point.in.SpatialPolygons(tmpdat$longitude, tmpdat$latitude, lakes[x,])
		if(sum(tmptest[[x]]) > 0) {
			thisone <- c(thisone,x)
		}
	}

	rm(tmptest, tmpdat)
	tmplake <- lakes[thisone,]

	#####
	#!#NEW STUFF ADDED TO IDENTIFY RIVERS HERE!!!
	#!#JDR - 5/3/2021
	#####
	tmplake2 <- spTransform(tmplake, crs(rivs))
	keepfeat <- intersect(rivs, buffer(tmplake2, width = 10))
	if(length(keepfeat) > 0) {
		keepfeat <- keepfeat$OBJECTID[keepfeat$FCC != "H21"]
		if(length(keepfeat) > 0) {
			tmprivs <- rivs[rivs$OBJECTID %in% keepfeat,]
			tmprivs <- spTransform(tmprivs, crs(lakes))
		} else {
			tmprivs <- NULL
		}
	}
	#####
	#####

	dat <- srdat[grep(lk, srdat$lake_name),]
	
	#Reduce name length for Five Channels
	if(dat$lake_name[1] == "Five.Channels.Impoundment") {
		dat$lake_name <- "Five.Channels"
	}

	dat <- dat[!is.na(dat$longitude),]

	#Loop over the variables to plot - SR for total, aquatic, wetland, and terrestrial (columns 7-10)
	#for(var in c(7:10)) {
		lakedat <- data.frame(long = dat$longitude, lat = dat$latitude, sr = dat[,var])
		coordinates(lakedat) <- ~long+lat
		grid <- spsample(tmplake, n = 100000, type = "regular")
		proj4string(lakedat) <- proj4string(grid)

		#Using an IDW power of 2 for all surfaces
		finalPOW <- 2

		#Fit the IDW model
		idwmodel <- idw(lakedat$sr~1, locations = lakedat, newdata = grid, idp = finalPOW, maxdist = Inf)
		idwDF <- as.data.frame(idwmodel)[,1:3]
		rasIDW <- rasterFromXYZ(idwDF)

		#Then make the plot
		plot(tmplake, main = gsub("\\.", " ", dat$lake_name[1]))

		xats <- round(seq(bbox(tmplake)[1,1], bbox(tmplake)[1,2], length.out = 4),3)
		yats <- round(seq(bbox(tmplake)[2,1], bbox(tmplake)[2,2], length.out = 4),3)
		axis(1, at = xats)
		axis(2, at = yats) 
		if(length(tmprivs) > 0) {
			plot(tmprivs, add = TRUE, col = "blue", lwd = 2)
		}
		points(lakedat$long, lakedat$lat, pch = 19, cex = 0.65)
		col <- matlab.like(50)
		plot(rasIDW, col = col, add = TRUE, alpha =0.75, legend.shrink = 0.4, legend.mar = 4.5, legend.cex = 0.5)
				
		#!#Adding North arrow & scale bar - JDR 5/3/21
		addnortharrow("topright", scale = 0.4)
		addscalebar(htin=0.075, label.cex = 0.75)
		rm(xats, yats, col, idwmodel, idwDF, rasIDW, lakedat, finalPOW, grid)
	}
	dev.off()
}

