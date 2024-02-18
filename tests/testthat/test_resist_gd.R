test_that("resist_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(distmat <- get_resdist(mini_coords, mini_lyr))
  capture_warnings(
    rpi <- resist_gd(
      mini_vcf,
      mini_coords,
      mini_lyr,
      rarify = FALSE,
      maxdist = quantile(distmat, 0.05, na.rm = TRUE),
      distmat = distmat
    )
  )

  expect_s4_class(rpi, "SpatRaster")
  expect_equal(terra::nlyr(rpi), 2)

  # check against expected values
  vals <- terra::global(rpi, fun = "mean", na.rm = TRUE)
  expect_equal(0.2182051, vals["pi", ], tolerance = 0.000001)
  expect_equal(0.5, vals["sample_count", ])

  # check without premaking matrix
  capture_warnings(
    rpi2 <- resist_gd(
      mini_vcf,
      mini_coords,
      mini_lyr,
      rarify = FALSE,
      maxdist = 50,
      distmat = NULL
    )
  )
})

test_that("resist_general returns expected output", {
  load_mini_ex(quiet = TRUE)

  # examples with custom functions
  toy <- vcf_to_dosage(mini_vcf_NA) * 0 + 1
  toy[1:2, ] <- NA

  # check if additional custom arguments provided work
  foo <- function(x, silly) sum(x * silly, na.rm = TRUE)

  capture_warnings(rg_1 <-
    resist_general(
      toy[, 3],
      mini_coords,
      mini_lyr,
      stat = foo,
      maxdist = 100,
      silly = 2
    )[[1]])

  capture_warnings(rg_2 <-
    resist_general(
      toy[, 3],
      mini_coords,
      mini_lyr,
      stat = foo,
      maxdist = 100,
      silly = 1
    )[[1]])

  expect_equal(terra::values(rg_1), terra::values(rg_2) * 2)
})


test_that("resist_gd works for different spatial types", {
  load_mini_ex(quiet = TRUE)

  # sf coords
  sf_coords <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  capture_warnings(wpi_sf <- resist_gd(mini_vcf, sf_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # matrix coords
  mat_coords <- as.matrix(mini_coords)
  capture_warnings(wpi_mat <- resist_gd(mini_vcf, mat_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # df coords
  capture_warnings(wpi_df <- resist_gd(mini_vcf, mini_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # spatvec coords
  vect_coords <- terra::vect(sf_coords)
  capture_warnings(wpi_vect <- resist_gd(mini_vcf, vect_coords, mini_lyr, maxdist = 50, rarify = FALSE))

  # compare rasters
  # expect_true(terra::all.equal(wpi_df, wpi_mat))
  # expect_true(terra::all.equal(wpi_df, wpi_sf))
  # expect_true(terra::all.equal(wpi_df, wpi_vect))
  expect_equal(terra::values(wpi_df), terra::values(wpi_mat))
  expect_equal(terra::values(wpi_df), terra::values(wpi_sf))
  expect_equal(terra::values(wpi_df), terra::values(wpi_vect))
})
