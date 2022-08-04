# Code to create middle earth example data -------------------------------------------------------------

# load and save raster layer
lyr <- read.csv("inst/extdata/middle_earth.csv", header = FALSE)
lotr_lyr <- raster::raster(as.matrix(lyr))
raster::extent(lotr_lyr) <- raster::extent(0, 100, -100, 0)
usethis::use_data(lotr_lyr, overwrite = TRUE)

# make a fake range map
lotr_lyr[lotr_lyr < 0.01] <- NA
lotr_lyr <- lotr_lyr * 0
lotr_range <- raster::rasterToPolygons(lotr_lyr, dissolve = TRUE, n = 16)
usethis::use_data(lotr_range, overwrite = TRUE)

# load coords
lotr_coords <- read.csv("inst/extdata/mod-sim_params_it-0_t-1000_spp-spp_0.csv") %>%
  dplyr::select(idx, x, y) %>%
  dplyr::mutate(y = -y)

# get subsample
# use lotr layer as probability so that sampling is more even across the landscape
p <- extract(lotr_lyr, lotr_coords[, c("x", "y")])
set.seed(42)
samples <- sample(nrow(lotr_coords), 100, prob = 1 / p)
lotr_coords <- lotr_coords[samples, ]

# load genetic data
vcf <- vcfR::read.vcfR("inst/extdata/mod-sim_params_it-0_t-1000_spp-spp_0.vcf")

# subsample loci and individuals
# note: first column is FORMAT, hence c(1, samples + 1)
lotr_vcf <- vcf[sample(1:nrow(vcf@gt), 100), c(1, samples + 1)]
# check to make sure order and IDs are the same
stopifnot(colnames(lotr_vcf@gt)[-1] == lotr_coords$idx)

# save results
lotr_coords <- lotr_coords %>% dplyr::select(x, y)
usethis::use_data(lotr_coords, overwrite = TRUE)
usethis::use_data(lotr_vcf, overwrite = TRUE)


# Code to create tiny example dataset ------------------------------------------------------------------
mini_lyr <- raster::aggregate(lotr_lyr, 10)
mini_vcf <- lotr_vcf[1:10, 1:11]
mini_coords <- lotr_coords[1:10, ]

usethis::use_data(mini_lyr, overwrite = TRUE)
usethis::use_data(mini_vcf, overwrite = TRUE)
usethis::use_data(mini_coords, overwrite = TRUE)
