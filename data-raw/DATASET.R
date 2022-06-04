
ex_vcf <- vcfR::read.vcfR("inst/extdata/ex_vcf.vcf")
usethis::use_data(ex_vcf, overwrite = TRUE)

ex_coords <- read.csv("inst/extdata/ex_coords.csv")
usethis::use_data(ex_coords, overwrite = TRUE)

ex_lyr <- raster::raster("inst/extdata/ex_layer.tif")
usethis::use_data(ex_lyr, overwrite = TRUE)
