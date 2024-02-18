test_that("circle_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(cpi <- circle_gd(mini_vcf_NA, mini_coords, mini_lyr, maxdist = 50))
  expect_s4_class(cpi, "SpatRaster")
  expect_equal(terra::nlyr(cpi), 2)

  # check against expected values
  vals <- terra::global(cpi, fun = "mean", na.rm = TRUE)
  expect_equal(0.3733736, vals["pi", ], tolerance = 0.000001)
  expect_equal(4.44, vals["sample_count", ])

  # test general functions
  load_mini_ex(quiet = TRUE)
  expect_warning(expect_warning(data <- check_data(mini_vcf_NA, mini_coords)))
  mini_vcf_NA <- data$vcf
  mini_coords <- data$coords
  capture_warnings(distmat <- get_geodist(mini_coords, mini_lyr))
  capture_warnings(wp <- circle_general(vcf_to_dosage(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "pi", rarify = FALSE))
  capture_warnings(wh <- circle_general(vcf_to_het(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "Ho", rarify = FALSE))
  capture_warnings(wb <- circle_general(vcf_to_dosage(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = FALSE))
  capture_warnings(wbr <- circle_general(vcf_to_dosage(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = TRUE))
  capture_warnings(wa <- circle_general(vcfR::vcfR2genind(mini_vcf_NA), maxdist = 10, distmat = distmat, mini_coords, mini_lyr, stat = "allelic_richness", rarify = FALSE))

  # test with custom maxdist
  mini_lyr[] <- 40
  mini_lyr[1:50] <- 0
  capture_warnings(wp <- circle_general(vcf_to_dosage(mini_vcf_NA), maxdist = mini_lyr, distmat = distmat, mini_coords, mini_lyr, stat = "pi", rarify = FALSE))

  # examples with custom functions
  toy <- vcf_to_dosage(mini_vcf_NA) * 0 + 1
  # check if additional custom arguments provided work
  foo <- function(x, silly) sum(x * silly, na.rm = TRUE)
  capture_warnings(cg_1 <-
    circle_general(
      toy[, 3],
      mini_coords,
      mini_lyr,
      stat = foo,
      maxdist = 50,
      silly = 2
    )[[1]])

  capture_warnings(cg_2 <-
    circle_general(
      toy[, 3],
      mini_coords,
      mini_lyr,
      stat = foo,
      maxdist = 50,
      silly = 1
    )[[1]])

  expect_equal(terra::values(cg_1), terra::values(cg_2) * 2)
})



test_that("circle_gd works for different spatial types", {
  load_mini_ex(quiet = TRUE)

  # sf coords
  sf_coords <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  capture_warnings(wpi_sf <- circle_gd(mini_vcf, sf_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # matrix coords
  mat_coords <- as.matrix(mini_coords)
  capture_warnings(wpi_mat <- circle_gd(mini_vcf, mat_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # df coords
  capture_warnings(wpi_df <- circle_gd(mini_vcf, mini_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # spatvec coords
  vect_coords <- terra::vect(sf_coords)
  capture_warnings(wpi_vect <- circle_gd(mini_vcf, vect_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # compare rasters
  # expect_true(terra::all.equal(wpi_df, wpi_mat))
  # expect_true(terra::all.equal(wpi_df, wpi_sf))
  # expect_true(terra::all.equal(wpi_df, wpi_vect))
  expect_equal(terra::values(wpi_df), terra::values(wpi_mat))
  expect_equal(terra::values(wpi_df), terra::values(wpi_sf))
  expect_equal(terra::values(wpi_df), terra::values(wpi_vect))
})
