lyr <- read.csv(here("sims", "data", "middle_earth.csv"), header = FALSE)
lyr <- raster(as.matrix(lyr))
extent(lyr) <- extent(0,100,-100,0)

bottleneck_lyr <- lyr/2


write.table(as.matrix(bottleneck_lyr), here("sims", "data", "bottleneck_lyr.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
lyr_bc <- read.csv( here("sims", "data", "bottleneck_lyr.csv"), header = FALSE)
lyr_bc <- raster(as.matrix(lyr_bc))
extent(lyr_bc) <- extent(0,100,-100,0)

plot(lyr, zlim = c(0,1))
plot(lyr_bc, zlim = c(0,1))
