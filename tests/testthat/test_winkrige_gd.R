test_that("winkrige_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  names(mini_lyr) <- "test"

  # Basic kriging
  capture_warnings(kpi <- winkrige_gd(mini_lyr, grd = mini_lyr))
  expect_s4_class(kpi, "SpatRaster")
  expect_equal(terra::nlyr(kpi), 1)

  # Check model_output = TRUE returns full list
  capture_warnings(kpi_list <- winkrige_gd(mini_lyr, grd = mini_lyr, model_output = TRUE))
  expect_true(is.list(kpi_list))
  expect_s4_class(kpi_list$raster, "SpatRaster")
  expect_s3_class(kpi_list$variogram, "gstatVariogram")
  expect_s3_class(kpi_list$model, "variogramModel")

  # Check extents match
  expect_true(terra::ext(mini_lyr) == terra::ext(kpi))
})

test_that("winkrige_gd with weights produces expected results", {
  load_mini_ex(quiet = TRUE)

  # Create dummy weights (e.g., all sample counts = 10)
  weight_r <- mini_lyr
  weight_r[] <- 10

  capture_warnings(kpi <- winkrige_gd(mini_lyr, grd = mini_lyr, weight_r = weight_r))
  expect_s4_class(kpi, "SpatRaster")

  # Confirm error if weights have values > 0
  weight_r[1] <- 0
  expect_error(
    capture_warnings(winkrige_gd(mini_lyr, grd = mini_lyr, weight_r = weight_r)),
    "All values of weight_r must be > 0"
  )

  # Confirm error if weights have NA values
  weight_r[1] <- NA
  expect_error(
    capture_warnings(winkrige_gd(mini_lyr, grd = mini_lyr, weight_r = weight_r)),
    "All values of weight_r must be non-NA"
  )
})

test_that("winkrige_gd allows custom starting values", {
  load_mini_ex(quiet = TRUE)
  
  capture_warnings(kpi <- winkrige_gd(
    mini_lyr,
    grd = mini_lyr,
    psill_start = 0.5,
    nugget_start = 0.1,
    range_start = 1,
    model_output = TRUE
  ))
  expect_true(is.list(kpi))
  expect_s3_class(kpi$model, "variogramModel")
})

test_that("winkrige_gd can handle multiple models", {
  load_mini_ex(quiet = TRUE)

  # Test with multiple models
  capture_warnings(kpi <- winkrige_gd(
    mini_lyr,
    grd = mini_lyr,
    models = c("Sph", "Exp", "Gau"),
    model_output = TRUE
  ))

  # Test with one model
  capture_warnings(kpi_single <- winkrige_gd(
    mini_lyr,
    grd = mini_lyr,
    models = "Sph",
    model_output = TRUE
  ))
  
  expect_true(inherits(kpi, "list"))
})