test_that("circle_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(cpi <- circle_gd(mini_vcf_NA, mini_coords, mini_lyr, maxdist = 50))
  expect_s4_class(cpi, "SpatRaster")
  expect_equal(terra::nlyr(cpi), 2)

  # check against expected values
  vals <- terra::global(cpi, fun = "mean", na.rm = TRUE)
  expect_equal(0.3658471, vals["pi", ], tolerance = 0.000001)
  expect_equal(5.12, vals["sample_count", ])

  # test general functions
  distmat <- get_geodist(mini_coords, mini_lyr)
  expect_warning(expect_warning(data <- check_data(mini_vcf_NA, mini_coords)))
  mini_vcf_NA <- data$vcf
  mini_coords <- data$coords
  capture_warnings(wp <- circle_general(vcf_to_dosage(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "pi", rarify = FALSE))
  capture_warnings(wh <- circle_general(vcf_to_het(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "Ho", rarify = FALSE))
  capture_warnings(wb <- circle_general(vcf_to_dosage(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = FALSE))
  capture_warnings(wbr <- circle_general(vcf_to_dosage(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = TRUE))
  capture_warnings(wa <- circle_general(vcfR::vcfR2genind(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "allelic_richness", rarify = FALSE))

  # examples with custom functions
  toy <- vcf_to_dosage(mini_vcf_NA)
  # test on vector
  capture_warnings(wm <- circle_general(toy[, 1], mini_coords, mini_lyr, maxdist = 10, distmat = distmat, stat = mean, na.rm = TRUE))
  # test on matrix
  capture_warnings(wm <- circle_general(toy, mini_coords, mini_lyr, maxdist = 10, distmat = distmat, stat = mean, na.rm = TRUE))
  # test custom functions
  foo <- function(x) var(apply(x, 2, var, na.rm = TRUE), na.rm = TRUE)
  capture_warnings(wm <- circle_general(toy, mini_coords, mini_lyr, maxdist = 10, distmat = distmat, stat = foo))
  foo <- function(x, na.rm = TRUE) var(apply(x, 2, var))
  capture_warnings(wm <- circle_general(toy, mini_coords, mini_lyr, maxdist = 10, distmat = distmat, stat = foo, na.rm = TRUE))
  foo <- function(x, na.rm = TRUE, silly = 2) sd(apply(x, 2, var)) * silly
  capture_warnings(wm <- circle_general(toy, mini_coords, mini_lyr, maxdist = 10, distmat = distmat, stat = foo, na.rm = TRUE, silly = 3))
})
