test_that("mask_gd returns expected output", {
  load_mini_ex(quiet = TRUE)

  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  expect_warning(kpi <- krig_gd(wpi, mini_lyr))
  mpi <- mask_gd(kpi, mini_lyr, minval = 2)

  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
  expect_true(raster::cellStats(is.na(mpi) == (mini_lyr < 2), all))

  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  expect_warning(kpi <- krig_gd(wpi, mini_lyr))
  mpi <- mask_gd(kpi, mini_lyr, maxval = 2)

  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
  expect_true(raster::cellStats(is.na(mpi) == (mini_lyr > 2), all))
})

test_that("resampling occurs correctly", {
  load_mini_ex(quiet = TRUE)

  x <- mini_lyr
  mask <- raster::aggregate(mini_lyr, 2)

  # resample to mask to match x
  expect_error(msm <- mask_gd(x, mask, resample = "mask"), NA)
  expect_true(raster::compareRaster(msm, x))
  # resample x to match mask
  expect_error(msx <- mask_gd(x, mask, resample = "x"), NA)
  expect_true(raster::compareRaster(msx, mask))

  # check for error if invalid resample arg is supplied
  expect_error(mse <- mask_gd(x, mask, resample = "grid"))
})
