## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(wingen)
library(terra)
library(raster)
library(viridis)
library(sf)
library(ggplot2)

## ----load example data, fig.width = 5, fig.height = 5-------------------------
load_middle_earth_ex()

# Genetic data
lotr_vcf

# Coordinates
head(lotr_coords)

# Raster data
lotr_lyr

# Map of data
plot(lotr_lyr, col = magma(100), axes = FALSE, box = FALSE)
points(lotr_coords, col = mako(1, begin = 0.8), pch = 3, cex = 0.5)

## ----fig.width = 5, fig.height = 5--------------------------------------------
ex_raster1 <- coords_to_raster(lotr_coords, buffer = 1, plot = TRUE)
ex_raster2 <- coords_to_raster(lotr_coords, buffer = 1, agg = 2, plot = TRUE)
ex_raster3 <- coords_to_raster(lotr_coords, buffer = 1, disagg = 4, plot = TRUE)
ex_raster4 <- coords_to_raster(lotr_coords, buffer = 1, res = 10, plot = TRUE)

## -----------------------------------------------------------------------------
# First, we create example latitude-longitude coordinates and rasters

## Example raster:
lyr_longlat <- rast(
  ncols = 40, nrows = 40, xmin = -110, xmax = -90, ymin = 40, ymax = 60,
  crs = "+proj=longlat +datum=WGS84"
)

## Example coordinates:
coords_df <- data.frame(x = c(-110, -90), y = c(40, 60))
coords_longlat <- st_as_sf(coords_df, coords = c("x", "y"), crs = "+proj=longlat")

# Next, the coordinates and raster can be projected to an equal area projection, in this case the Goode Homolosine projection (https://proj.org/operations/projections/goode.html):
coords_eq <- st_transform(coords_longlat, crs = "+proj=goode")
lyr_eq <- project(lyr_longlat, "+proj=goode")

# Coordinates can be switched back to latitude-longitude the same way by replacing "goode" with "longlat"

## ----fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE-------------
preview_gd(lotr_lyr,
  lotr_coords,
  method = "window",
  wdim = 7,
  fact = 3,
  sample_count = TRUE,
  min_n = 2
)

## ----moving window, fig.width = 5, fig.height = 5, cache = TRUE, warning = TRUE, message = FALSE----
wgd <- window_gd(lotr_vcf,
  lotr_coords,
  lotr_lyr,
  stat = "pi",
  wdim = 7,
  fact = 3,
  rarify = TRUE,
  rarify_n = 2,
  rarify_nit = 5,
  L = 100
)

par(pty = "s")
# The plot_gd function plots the genetic diversity layer
plot_gd(wgd, bkg = lotr_range, main = "Moving window pi")
# Add coordinates to map
points(lotr_coords, pch = 16)

# The plot_count function plots the sample count layer
plot_count(wgd, main = "Moving window sample counts")

## ----fig.width = 5, fig.height = 5--------------------------------------------
ggplot_gd(wgd, bkg = lotr_range) +
  coord_sf() +
  geom_point(data = lotr_coords, aes(x = x, y = y)) +
  ggtitle("Moving window pi")

ggplot_count(wgd) +
  coord_sf() +
  ggtitle("Moving window sample counts")

## ----fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE-------------
# First, create a new aggregated raster
lotr_lyr3 <- aggregate(lotr_lyr, 3)

# Then calculate the distance matrix
distmat <- get_geodist(coords = lotr_coords, lyr = lotr_lyr3)

preview_gd(lotr_lyr3,
  lotr_coords,
  method = "circle",
  maxdist = 10,
  distmat = distmat,
  sample_count = TRUE,
  min_n = 2
)

## ----circle, fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE, message = FALSE----
cgd <- circle_gd(lotr_vcf,
  lotr_coords,
  lotr_lyr3,
  stat = "pi",
  maxdist = 10,
  distmat = distmat,
  rarify = FALSE,
  L = 100
)

par(pty = "s")
# The plot_gd function plots the genetic diversity layer
plot_gd(cgd, bkg = lotr_range, main = "Circle moving window pi")
# The plot_count function plots the sample count layer
plot_count(cgd, main = "Circle moving window sample counts")

## ----variable circle, fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE, message = FALSE----
vcgd <- circle_gd(lotr_vcf,
  lotr_coords,
  lotr_lyr3,
  stat = "pi",
  maxdist = lotr_lyr3 * 100,
  distmat = distmat,
  rarify = FALSE,
  L = 100
)

par(pty = "s")
# The plot_gd function plots the genetic diversity layer
plot_gd(vcgd, bkg = lotr_range, main = "Circle moving window pi")
# The plot_count function plots the sample count layer
plot_count(vcgd, main = "Circle moving window sample counts")

## ----cache = TRUE, warning = FALSE, message = FALSE---------------------------
lotr_distmat <- get_resdist(coords = lotr_coords, lyr = aggregate(lotr_lyr, 3))

## ----fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE, message = FALSE----
preview_gd(lotr_lyr,
  lotr_coords,
  method = "resist",
  maxdist = 60,
  distmat = lotr_distmat,
  fact = 3,
  sample_count = TRUE,
  min_n = 2
)

## ----resist, fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE, message = FALSE----
rgd <- resist_gd(lotr_vcf,
  lotr_coords,
  lotr_lyr,
  distmat = lotr_distmat,
  stat = "pi",
  maxdist = 60,
  fact = 3,
  rarify = FALSE,
  L = 100
)

par(pty = "s")
# The plot_gd function plots the genetic diversity layer
plot_gd(rgd, bkg = lotr_range, main = "Resistance moving window pi")
# The plot_count function plots the sample count layer
plot_count(rgd, main = "Resistance moving window sample counts")

## ----fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE-------------
# create layer
resist_lyr <- rast(aggregate(lotr_lyr, 5))

# create an impassable area down the middle coded as NA
resist_lyr[] <- 1
resist_lyr[ext(40, 60, -100, 0)] <- NA

# get resistance distances
distmat <- get_resdist(coords = lotr_coords, lyr = resist_lyr)

# plot preview
preview_gd(resist_lyr, lotr_coords, method = "resist", distmat = distmat, maxdist = 10)

# create moving window map
resist_gd <- resist_gd(
  lotr_vcf,
  lotr_coords,
  resist_lyr,
  maxdist = 10,
  distmat = distmat
)

plot_gd(resist_gd)

## ----krige results, fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE----
# Note: this step can take a little while
# index = 1 kriges the first layer in wgd (the genetic diversity layer)
kgd <- krig_gd(wgd, index = 1, grd = lotr_lyr, disagg_grd = 2)

par(pty = "s")
plot_gd(kgd, main = "Kriged pi")

## ----universal krige results, fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE----
kgd_ordinary <-
  krig_gd(
    wgd,
    index = 1,
    grd = lotr_lyr,
    disagg_grd = 2,
    krig_method = "ordinary"
  )

par(pty = "s")
plot_gd(kgd_ordinary, main = "Kriged pi")

## ----fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE-------------
kgd_autoKrige <-
  krig_gd(
    wgd,
    index = 1,
    grd = lotr_lyr,
    disagg_grd = 2,
    autoKrige_output = TRUE
  )

summary(kgd_autoKrige)

## ----mask results 1, fig.width = 5, fig.height = 5, warning = FALSE-----------
# Disaggregate lotr_lyr to make it the same resolution as kgd before masking
## Note: lotr_lyr is a RasterLayer which we convert to a SpatRaster with rast()
mask_lyr <- disagg(rast(lotr_lyr), 2)
mgd <- mask_gd(kgd, mask_lyr, minval = 0.1)

par(pty = "s")
plot_gd(mgd, main = "Kriged & carrying capacity masked pi")

## ----mask results 2, fig.width = 5, fig.height = 5, warning = FALSE-----------
mgd <- mask_gd(kgd, lotr_range)

par(pty = "s")
plot_gd(mgd, main = "Kriged & range masked pi")

## ----eval = FALSE-------------------------------------------------------------
#  # First, install and load the SpatialKDE package
#  if (!require("SpatialKDE", quietly = TRUE)) install.packages("SpatialKDE")
#  library(SpatialKDE)

## ----mask results 3, fig.width = 5, fig.height = 5, eval = FALSE--------------
#  # Note: this code is not evaluated when building the vignette as it requires the SpatialKDE package
#  
#  # Spatial KDE requires an sf data.frame containing only POINTS that is in a projected coordinate system
#  # The simulated coordinates are not projected, but for the purpose of this example, we'll pretend they are projected under the Goode Homolosine projection
#  # This kind of arbitrary setting of crs should not be done for real data (see above example for how to properly project coordinates)
#  lotr_sf <- st_as_sf(lotr_coords, coords = c("x", "y")) %>% st_set_crs("+proj=goode")
#  
#  # The grid layer must also be a RasterLayer for SpatialKDE
#  # Here, we use the kriged raster as our grid to get a KDE layer of the same spatial extent and resolution
#  grid <- raster(kgd)
#  
#  # See the SpatialKDE package for more options and details about using KDE
#  kde_lyr <- kde(
#    lotr_sf,
#    kernel = "quartic",
#    band_width = 15,
#    decay = 1,
#    grid = grid,
#  )
#  
#  par(pty = "s")
#  # Visualize KDE layer
#  plot_count(kde_lyr, main = "KDE")
#  
#  # Mask with mask_gd
#  mgd <- mask_gd(kgd, kde_lyr, minval = 1)
#  
#  plot_gd(mgd, main = "Kriged & sample density masked pi")

## ----range map background, fig.width = 5, fig.height = 5, warning = FALSE-----
par(pty = "s")
plot_gd(mgd, bkg = lotr_range, main = "Kriged & masked pi")

## ----parallelization, fig.width = 5, fig.height = 5, eval = FALSE, message = FALSE----
#  # Note: this code is not evaluated when building the vignette as it spawns multiple processes
#  
#  # setup parallel session
#  future::plan("multisession", workers = 2)
#  
#  wgd <- window_gd(lotr_vcf,
#    lotr_coords,
#    lotr_lyr,
#    stat = "pi",
#    wdim = 7,
#    fact = 3,
#    rarify_n = 2,
#    rarify_nit = 5,
#    rarify = TRUE
#  )
#  
#  # end parallel session
#  future::plan("sequential")

## ----cache = TRUE, warning = FALSE, warning = FALSE, message = FALSE----------
multistat_wgd <- window_gd(lotr_vcf,
  lotr_coords,
  lotr_lyr,
  stat = c("pi", "Ho", "biallelic_richness"),
  wdim = 7,
  fact = 3,
  rarify = FALSE,
  L = 100
)

## ----fig.width = 5, fig.height = 5, warning = FALSE---------------------------
par(pty = "s")
plot_gd(multistat_wgd, bkg = lotr_range)

## ----warning = FALSE, message = FALSE-----------------------------------------
hwe_wgd <- window_gd(lotr_vcf[1:10, ],
  lotr_coords,
  lotr_lyr,
  stat = "hwe",
  wdim = 3,
  fact = 5,
  rarify = FALSE,
  L = 100,
  sig = 0.10
)

bs_wgd <- window_gd(lotr_vcf[1:10, ],
  lotr_coords,
  lotr_lyr,
  stat = "basic_stats",
  wdim = 3,
  fact = 5,
  rarify = FALSE,
  L = 100
)

## ----fig.width = 5, fig.height = 5, warning = FALSE---------------------------
par(pty = "s")
plot_gd(hwe_wgd, bkg = lotr_range)

## ----fig.width = 5, fig.height = 5, warning = FALSE---------------------------
par(pty = "s")
plot_gd(bs_wgd, bkg = lotr_range)

## ----fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE, message = FALSE----
# First, we extract the raster values at those coordinates
vals <- extract(lotr_lyr, lotr_coords)

# Next, we run the window_general function with the env vector and set the `stat` to mean
# Note: we can also provide additional arguments to functions, such as na.rm = TRUE
we <- window_general(vals,
  coords = lotr_coords,
  lyr = lotr_lyr,
  stat = mean,
  wdim = 7,
  fact = 3,
  rarify_n = 2,
  rarify_nit = 5,
  rarify = TRUE,
  na.rm = TRUE
)

par(pty = "s")
plot_gd(we, main = "Mean raster value")

## ----fig.width = 5, fig.height = 5, cache = TRUE, warning = FALSE, message = FALSE----
# We use the vcfR package to convert vcf to genind for our example
library(vcfR)

# Convert existing vcf example file into genind object:
genind <- vcfR2genind(lotr_vcf)

# Run window_general with no rarefaction
we_gi <- window_general(genind,
  coords = lotr_coords,
  lyr = lotr_lyr,
  stat = "allelic_richness",
  wdim = 7,
  fact = 3,
  na.rm = TRUE
)

# Plot results
par(pty = "s")
plot_gd(we_gi, main = "Allelic richness")

