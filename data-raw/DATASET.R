# Code to create middle earth example data -------------------------------------------------------------

# load coords
lotr_coords <- read.csv("inst/extdata/mod-sim_params_it-0_t-1000_spp-spp_0.csv") %>%
  dplyr::select(idx, x, y) %>%
  dplyr::mutate(y = -y)

# get subsample
samples <- read.csv("inst/extdata/samples_seed42.csv")
lotr_coords <- lotr_coords[samples$inds, ]

# load genetic data
vcf <- vcfR::read.vcfR("inst/extdata/mod-sim_params_it-0_t-1000_spp-spp_0.vcf")

# subsample loci and individuals
# note: first column is FORMAT, hence c(1, samples + 1)
lotr_vcf <- vcf[sample(1:nrow(vcf@gt), 1000), c(1, samples$inds + 1)]
# check to make sure order and IDs are the same
stopifnot(colnames(lotr_vcf@gt)[-1] == lotr_coords$idx)

# save results
lotr_coords <- lotr_coords %>% dplyr::select(x, y)
usethis::use_data(lotr_coords, overwrite = TRUE)
usethis::use_data(lotr_vcf, overwrite = TRUE)

# load and save raster layer
lyr <- read.csv("inst/extdata/middle_earth.csv", header = FALSE)
lotr_lyr <- raster::raster(as.matrix(lyr))
raster::extent(lotr_lyr) <- raster::extent(0,100,-100,0)
usethis::use_data(lotr_lyr, overwrite = TRUE)

# Code to create tiny example dataset ------------------------------------------------------------------
mini_lyr <- raster::aggregate(lotr_lyr, 10)
mini_vcf <- lotr_vcf[1:10,1:11]
mini_coords <- lotr_coords[1:10,]

usethis::use_data(mini_lyr, overwrite = TRUE)
usethis::use_data(mini_vcf, overwrite = TRUE)
usethis::use_data(mini_coords, overwrite = TRUE)

