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


