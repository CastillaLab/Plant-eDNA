library(raster)
elevationmapdata <- getData("GADM", country="US", level=0)

### each ma.dem line pulls a piece of the map

ma.dem1 <- getData("SRTM", lon=-103.81, lat=29.52)
ma.dem2 <- getData("SRTM", lon=-103.9, lat=29.55)
ma.dem3 <- getData("SRTM", lon=-103.65, lat=30)
ma.dem4 <- getData("SRTM", lon=-95, lat=30)
ma.dem5 <- getData("SRTM", lon=-115, lat=30)
ma.dem6 <- getData("SRTM", lon=-107, lat=30)
ma.dem7 <- getData("SRTM", lon=-97, lat=30)
ma.dem8 <- getData("SRTM", lon=-93, lat=30)
ma.dem9 <- getData("SRTM", lon=-93, lat=34)
ma.dem10 <- getData("SRTM", lon=-115, lat=34)
ma.dem11 <- getData("SRTM", lon=-107, lat=34)
ma.dem12 <- getData("SRTM", lon=-101, lat=32)

##putting the pieces together

ma.dem <- mosaic(ma.dem1, ma.dem2, ma.dem3, ma.dem4, ma.dem5, ma.dem6, ma.dem7, ma.dem8, ma.dem9, ma.dem10, ma.dem11, ma.dem12, fun=mean)
plot(ma.dem)

###cropping the map ma.dem.c is the close map, ma.dem.f is the zoomed out map

ma.dem.c <- crop(ma.dem, extent(-103.815, -103.77, 29.540, 29.557))
plot(ma.dem.c)
plot(elevationmapdata, add=TRUE)


ma.dem.f <- crop(ma.dem, extent(-108, -100, 28, 31))
plot(ma.dem.f)
plot(elevationmapdata, add=TRUE)

p.lon <- c(-103.80400,
           -103.80500,
           -103.80600,
           -103.80600,
           -103.80600,
           -103.80600,
           -103.80700,
           -103.80700,
           -103.80800,
           -103.80900,
           -103.79600,
           -103.79800,
           -103.80000,
           -103.80200,
           -103.78700,
           -103.78700,
           -103.78700,
           -103.78700,
           -103.78700,
           -103.78700,
           -103.78700,
           -103.78600,
           -103.78500,
           -103.78400,
           -103.78300,
           -103.78300,
           -103.78300,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.79400,
           -103.78700,
           -103.79400,
           -103.78700,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79300,
           -103.79200,
           -103.79200,
           -103.79200,
           -103.79100,
           -103.79100,
           -103.79100,
           -103.79300,
           -103.79300,
           -103.79400,
           -103.79400,
           -103.78600,
           -103.78600,
           -103.78600,
           -103.78600,
           -103.78600,
           -103.78600,
           -103.77900,
           -103.77900,
           -103.77900,
           -103.78000,
           -103.78000,
           -103.78000,
           -103.78000,
           -103.78000,
           -103.78000,
           -103.79200,
           -103.79200)
p.lat <- c(29.54600,
           29.54741,
           29.54766,
           29.54734,
           29.54739,
           29.54738,
           29.54717,
           29.54719,
           29.54634,
           29.54527,
           29.55052,
           29.54801,
           29.54697,
           29.54626,
           29.54976,
           29.54952,
           29.54944,
           29.54933,
           29.54971,
           29.54977,
           29.54965,
           29.54974,
           29.55069,
           29.55203,
           29.55202,
           29.55235,
           29.55264,
           29.55424,
           29.55424,
           29.55424,
           29.55225,
           29.55216,
           29.55215,
           29.55204,
           29.55226,
           29.55213,
           29.55210,
           29.55202,
           29.55203,
           29.55182,
           29.55181,
           29.55177,
           29.55173,
           29.55171,
           29.55163,
           29.54970,
           29.55159,
           29.54968,
           29.55210,
           29.55208,
           29.55205,
           29.55204,
           29.55200,
           29.55193,
           29.55192,
           29.55185,
           29.55186,
           29.55187,
           29.55189,
           29.55186,
           29.55154,
           29.55164,
           29.55195,
           29.55139,
           29.55192,
           29.55178,
           29.55075,
           29.55001,
           29.55448,
           29.55317,
           29.55278,
           29.55253,
           29.55245,
           29.55252,
           29.55286,
           29.55282,
           29.55263,
           29.55264,
           29.55264,
           29.55250,
           29.55246,
           29.55213,
           29.55208,
           29.55134,
           29.55129)
##adding points based on coordinates in .csv
points(p.lon, p.lat, pch=21, col="black", bg="white", cex=2)


##now I want to make the points colored by species

oaks<-read.csv("DDRS.oaks.csv",header=TRUE)
poi.c <- as.character(oaks$species)


poi.c[poi.c=="Q.grisea"] <- "yellowgreen"
poi.c[poi.c=="Q.pungens"] <- "gold"
points(oaks$longitude, oaks$latitude, pch=21, bg=poi.c, cex=1.3)

legend("bottomright", legend=c("Q.grisea", "Q.pungens"), pch=c(21, 21), pt.bg=c("yellowgreen", "gold"), pt.cex=1.3, bg="white")

