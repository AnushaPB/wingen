test_that("resist_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(distmatp <- get_resdist(mini_coords, mini_lyr, parallel = TRUE, ncores = 2))
  capture_warnings(distmat <- get_resdist(mini_coords, mini_lyr, parallel = FALSE))
  expect_equal(distmat, distmatp)

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

  capture_warnings(
    rpi <- resist_general(
      vcf_to_dosage(mini_vcf),
      mini_coords,
      mini_lyr,
      stat = "pi",
      rarify = FALSE,
      maxdist = 50,
      fact = 2,
      distmat = NULL
    )
  )

  expect_s4_class(rpi, "SpatRaster")

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
