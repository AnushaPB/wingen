# Code to create middle earth example data -------------------------------------------------------------

# load and save raster layer
lyr <- read.csv("inst/extdata/middle_earth.csv", header = FALSE)
lotr_lyr <- raster::raster(as.matrix(lyr))
raster::extent(lotr_lyr) <- raster::extent(0, 100, -100, 0)
usethis::use_data(lotr_lyr, overwrite = TRUE)

# make a fake range map
lotr_range <- lotr_lyr
lotr_range[lotr_range < 0.01] <- NA
lotr_range <- lotr_range * 0
lotr_range <- raster::rasterToPolygons(lotr_range, dissolve = TRUE, n = 16)
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
# check if file exists locally and if not download it from figshare
file <- "inst/extdata/mod-sim_params_it-0_t-1000_spp-spp_0.vcf"
if(!file.exists(file)){
  download.file("https://figshare.com/ndownloader/files/36617433?private_link=7f6783de9b3d7a0ed897", file)
}
vcf <- vcfR::read.vcfR(file)

# subsample loci and individuals
# note: first column is FORMAT, hence c(1, samples + 1)
lotr_vcf <- vcf[, c(1, samples + 1)]
# retain only variant sites
lotr_vcf <- lotr_vcf[is.polymorphic(lotr_vcf), ]
# subsample loci
set.seed(42)
lotr_vcf <- lotr_vcf[sample(1:nrow(lotr_vcf@gt), 100), ]
# check to make sure order and IDs are the same
stopifnot(colnames(lotr_vcf@gt)[-1] == lotr_coords$idx)

# save results
lotr_coords <- lotr_coords %>% dplyr::select(x, y)
usethis::use_data(lotr_coords, overwrite = TRUE)
usethis::use_data(lotr_vcf, overwrite = TRUE)


# Code to create tiny example dataset ------------------------------------------------------------------
mini_lyr <- raster::aggregate(lotr_lyr, 10)

mini_vcf <- lotr_vcf[, 1:11]
mini_vcf <- mini_vcf[is.polymorphic(mini_vcf),]
mini_vcf <- mini_vcf[1:10, ]

mini_coords <- lotr_coords[1:10, ]

# add NAs to vcf for testing
mini_vcf_NA <- mini_vcf
mini_vcf_NA@gt[, 10] <- NA
mini_vcf_NA@gt[9, ] <- NA
mini_vcf_NA@gt[8, c(1, 2, 3)] <- NA
mini_vcf_NA@gt[c(2, 3, 4), 8] <- NA
# complete format column
mini_vcf_NA@gt[, 1] <- "GT"

usethis::use_data(mini_lyr, overwrite = TRUE)
usethis::use_data(mini_vcf, overwrite = TRUE)
usethis::use_data(mini_coords, overwrite = TRUE)
usethis::use_data(mini_vcf_NA, overwrite = TRUE)
