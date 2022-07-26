lyr <- read.csv(here("sims", "data", "middle_earth.csv"), header = FALSE)

# increase initial pop size
lyr_init <- lyr*2
lyr_init[lyr_init > 1] <- 1

# on diagonal decrease pop size
lyr_bc <- lyr_init
tri <- lower.tri(lyr_bc)
tri <- tri[, ncol(tri):1]
lyr_bc[tri] <- lyr_init[tri]/4

write.table(as.matrix(lyr_bc), here("sims", "data", "bottleneck_lyr.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
lyr_bc <- read.csv( here("sims", "data", "bottleneck_lyr.csv"), header = FALSE)
lyr_bc <- raster(as.matrix(lyr_bc))
extent(lyr_bc) <- extent(0,100,-100,0)

write.table(as.matrix(lyr_init), here("sims", "data", "init_lyr.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
lyr_init <- read.csv( here("sims", "data", "init_lyr.csv"), header = FALSE)
lyr_init <- raster(as.matrix(lyr_init))
extent(lyr_init) <- extent(0,100,-100,0)

plot(lyr_init, zlim = c(0,1), col = magma(100))
plot(lyr_bc, zlim = c(0,1), col = magma(100))

cellStats(lyr_init, sum)
cellStats(lyr_bc, sum)
cellStats(lyr_init, sum)/cellStats(lyr_bc, sum)
