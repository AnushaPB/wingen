
m <- matrix(nrow = 100, ncol = 100)
lyr <- raster(m)
extent(lyr) <- extent(0,100,0,100)

image <- jpeg::readJPEG("sims/data/arda_resample.jpg")
lyr[] <- image
plot(lyr)

write.table(image, "sims/data/arda_100.csv", sep=",",  col.names = FALSE, row.names = FALSE)


lyr50 <- aggregate(lyr,2)
lyr50 <- as.matrix(lyr50)
write.table(lyr50, "sims/data/arda_50.csv", sep=",",  col.names = FALSE, row.names = FALSE)

# test

lyr50 <- read.csv("arda_50.csv")
lyr50 <- raster(as.matrix(lyr50))
plot(lyr50)
