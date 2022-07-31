library(here)
library(raster)

wdir <- here("paperex", "simex")

lyr <- read.csv(here(wdir, "data", "middle_earth.csv"), header = FALSE)
lyr <- raster(as.matrix(lyr))
extent(lyr) <- extent(0,100,-100,0)
bkg <- lyr
bkg[bkg < 0.01] <- NA


for (n in c(1000, 5000, 10000)){

  geo <- read.csv(here(wdir, "data", paste0("mod-sim_params_it-0_t-", n, "_spp-spp_0.csv")))
  coords <- geo[,c("idx","x","y")]
  coords$y <- -coords$y
  kde <- raster(MASS::kde2d(coords$x, coords$y, h = c(10,10), n = 100, lims = c(0,100,-100,0)))

  plot_gd(kde, breaks = 10)
  writeRaster(kde, here(wdir, "outputs", paste0(n,"_kde.tif")), overwrite = TRUE)
}







