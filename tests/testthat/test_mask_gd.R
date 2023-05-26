test_that("mask_gd returns expected output", {
  load_mini_ex(quiet = TRUE)

  capture_warnings(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE))
  capture_warnings(kpi <- krig_gd(wpi, mini_lyr))
  mpi <- mask_gd(kpi, mini_lyr, minval = 2)

  expect_s4_class(wpi, "SpatRaster")
  expect_equal(terra::nlyr(wpi), 2)
  expect_true(all(terra::values(is.na(mpi)) == terra::values(mini_lyr < 2)))

  capture_warnings(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE))
  capture_warnings(kpi <- krig_gd(wpi, mini_lyr))
  mpi <- mask_gd(kpi, mini_lyr, maxval = 2)

  expect_s4_class(wpi, "SpatRaster")
  expect_equal(terra::nlyr(wpi), 2)
  expect_true(all(terra::values(is.na(mpi)) == terra::values(mini_lyr > 2)))
})

test_that("scale mismatch produces error", {
  load_mini_ex(quiet = TRUE)

  x <- mini_lyr
  mask <- terra::aggregate(mini_lyr, 2)

  # check for error if invalid resample arg is supplied
  expect_error(mse <- mask_gd(x, mask))
})

test_that("different objects can be used for masking", {
  data("lotr_lyr")
  data("lotr_range")

  # SpatVector
  vect_range <- terra::vect(lotr_range)
  m1 <- mask_gd(lotr_lyr, lotr_range)

  # sf object
  sf_range <- sf::st_as_sf(lotr_range)
  m2 <- mask_gd(lotr_lyr, sf_range)

  # sp object
  sp_range <- sf::as_Spatial(sf_range)
  m3 <- mask_gd(lotr_lyr, sp_range)

  # expect_true(terra::all.equal(m1, m2))
  # expect_true(terra::all.equal(m1, m3))
  expect_equal(terra::values(m1), terra::values(m2))
  expect_equal(terra::values(m1), terra::values(m3))
})
