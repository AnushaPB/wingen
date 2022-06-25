# Code to create example data

ex_vcf <- vcfR::read.vcfR("inst/extdata/ex_vcf.vcf")
usethis::use_data(ex_vcf, overwrite = TRUE)

ex_coords <- read.csv("inst/extdata/ex_coords.csv")
usethis::use_data(ex_coords, overwrite = TRUE)

ex_lyr <- raster::raster("inst/extdata/ex_layer.tif")
usethis::use_data(ex_lyr, overwrite = TRUE)

# MIDDLE EARTH EXAMPLE DATA
vcf <- vcfR::read.vcfR("inst/extdata/mod-sim_params_it-0_t-500_spp-spp_0.vcf")
middle_earth_vcf <- vcf[sample(1:nrow(vcf@gt), 1000),]
usethis::use_data(middle_earth_vcf, overwrite = TRUE)

middle_earth_coords <- read.csv("inst/extdata/mod-sim_params_it-0_t-500_spp-spp_0.csv") %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(y = -y)
usethis::use_data(middle_earth_coords, overwrite = TRUE)

lyr <- read.csv("inst/extdata/middle_earth.csv", header = FALSE)
middle_earth_lyr <- raster::raster(as.matrix(lyr))
raster::extent(middle_earth_lyr) <- raster::extent(0,100,-100,0)
usethis::use_data(middle_earth_lyr, overwrite = TRUE)
