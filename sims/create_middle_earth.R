library("raster")

# function to rescale raster from 0 to 1
rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = cellStats(x, "min", na.rm = TRUE)
  if(is.null(x.max)) x.max = cellStats(x, "max", na.rm = TRUE)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

#DATA SOURCE: https://scholarworks.wm.edu/asoer/3/

if(!file.exists("sims/data/DEM_middle_earth.tif")){
  # load data
  e1 <- raster("sims/data/DEM_50m_Quad1.tif")
  e2 <- raster("sims/data/DEM_50m_Quad2.tif")
  e3 <- raster("sims/data/DEM_50m_Quad3.tif")
  e4 <- raster("sims/data/DEM_50m_Quad4.tif")

  # stich together
  stitched <- merge(e1, e2, e3, e4)
  plot(stitched)
  writeRaster(stitched, "sims/data/DEM_middle_earth.tif", overwrite = TRUE)
}

# Aggregate stitched layer
agg1 <- aggregate(stiched, 50)

# Crop stiched layer to remove weird warped part
crp <- crop(agg1, c(2527364, 4200000, 1500000, 3150000))
# Change extent so it is easier to figure out crop math
extent(crp) <- c(0, 4200000 - 2527364, 0, 3150000 - 1500000)
# Crop to a square size
crp2 <- crop(crp, c(0, 1650000, 0, 1650000))

# Aggregate to get a resolution close to 100 (aim high)
agg2 <- aggregate(crp2,6)
# Change extent to make math easier to crop to 100 x 100
extent(agg2) <- c(0,110,0,110)
# Crop to 100 x 100
crp3 <- crop(agg2, c(0,100,0,100))

# Rescale from 0 to 1
rms <- rescale(crp3)
# Turn into matrix
m <- as.matrix(rms)
# change NA values to 0 (gnx doesn't like NA values)
final_mat[is.na(final_mat)] <- 0

# write out matrix for gnx
write.table(final_mat, "sims/data/middle_earth.csv", sep = ",", row.names = FALSE, col.names = FALSE)

# read in matrix and plot to confirm it looks right
rmr <- read.csv("sims/data/middle_earth.csv", header = FALSE)
rmr <- as.matrix(rmr)
plot(raster(rmr))
