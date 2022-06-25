# Code to create example data

ex_vcf <- vcfR::read.vcfR("inst/extdata/ex_vcf.vcf")
usethis::use_data(ex_vcf, overwrite = TRUE)

ex_coords <- read.csv("inst/extdata/ex_coords.csv")
usethis::use_data(ex_coords, overwrite = TRUE)

ex_lyr <- raster::raster("inst/extdata/ex_layer.tif")
usethis::use_data(ex_lyr, overwrite = TRUE)

# MIDDLE EARTH EXAMPLE DATA

# load coords
middle_earth_coords <- read.csv("inst/extdata/mod-sim_params_it-0_t-500_spp-spp_0.csv") %>%
  dplyr::select(idx, x, y) %>%
  dplyr::mutate(y = -y)

# get subsample
samples <- read.csv("inst/extdata/samples_seed42.csv")
middle_earth_coords <- middle_earth_coords[samples$inds, ]

# load genetic data
vcf <- vcfR::read.vcfR("inst/extdata/mod-sim_params_it-0_t-500_spp-spp_0.vcf")

# subsample loci and individuals
# note: first column is FORMAT, hence c(1, samples + 1)
middle_earth_vcf <- vcf[sample(1:nrow(vcf@gt), 1000), c(1, samples$inds + 1)]
# check to make sure order and IDs are the same
stopifnot(colnames(middle_earth_vcf@gt)[-1] == middle_earth_coords$idx)

# save results
middle_earth_coords <- middle_earth_coords %>% dplyr::select(x, y)
usethis::use_data(middle_earth_coords, overwrite = TRUE)
usethis::use_data(middle_earth_vcf, overwrite = TRUE)

# load and save raster layer
lyr <- read.csv("inst/extdata/middle_earth.csv", header = FALSE)
middle_earth_lyr <- raster::raster(as.matrix(lyr))
raster::extent(middle_earth_lyr) <- raster::extent(0,100,-100,0)
usethis::use_data(middle_earth_lyr, overwrite = TRUE)
