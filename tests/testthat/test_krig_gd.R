test_that("krig_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  names(mini_lyr) <- "test"

  # check basic kriging
  capture_warnings(kpi <- krig_gd(mini_lyr, mini_lyr))
  expect_s4_class(kpi, "SpatRaster")
  expect_equal(terra::nlyr(kpi), 1)

  # check kriging of stack
  capture_warnings(expect_warning(kpi <- krig_gd(raster::stack(mini_lyr, mini_lyr), index = 1:2, mini_lyr)))
  expect_s4_class(kpi, "SpatRaster")
  expect_equal(raster::nlayers(mini_lyr), 1)

  # check extents match
  expect_true(terra::ext(mini_lyr) == terra::ext(kpi))
})

test_that("krig_gd returns warning when not provided grd", {
  load_mini_ex(quiet = TRUE)
  expect_error(warnings <- capture_warnings(kpi <- krig_gd(mini_lyr, grd = NULL)), NA)
  expect_equal(warnings[2], "No grd provided, defaults to using first raster layer to create grd")
})

test_that("coord kriging works", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(kpi1 <- krig_gd(mini_lyr, grd = mini_lyr, coords = mini_coords))

  # sf coords
  sf_coords <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  capture_warnings(kpi2 <- krig_gd(mini_lyr, grd = mini_lyr, coords = sf_coords))

  # expect_true(terra::all.equal(kpi1, kpi2))
  expect_equal(terra::values(kpi1), terra::values(kpi2))
})

test_that("krige_gd returns error when provided bad grd", {
  load_mini_ex(quiet = TRUE)

  capture_warnings(expect_error(kpi <- krig_gd(mini_lyr, grd = mini_coords)))
})

test_that("raster_transform transformations are correct", {
  data("mini_lyr")

  r <- terra::rast(mini_lyr)
  grd <- r
  r <- terra::aggregate(r, 2)

  # test pass through
  noTransform <- raster_transform(r, grd)
  # expect_true(terra::all.equal(noTransform[[1]], r))
  # expect_true(terra::all.equal(noTransform[[2]], grd))
  expect_equal(terra::values(noTransform[[1]]), terra::values(r))
  expect_equal(terra::values(noTransform[[2]]), terra::values(grd))

  # test resampling of r
  resample_r <- raster_transform(r, grd, resample = "r")
  # remove values by rast() because we don't expect them to be the same
  # expect_true(terra::all.equal(terra::rast(resample_r[[1]]), terra::rast(grd)))
  # expect_true(terra::all.equal(terra::rast(resample_r[[2]]), terra::rast(grd)))
  capture_warnings(expect_equal(terra::values(terra::rast(resample_r[[1]])), terra::values(terra::rast(grd))))
  capture_warnings(expect_equal(terra::values(terra::rast(resample_r[[2]])), terra::values(terra::rast(grd))))

  # test aggregation of r
  agg_r <- raster_transform(r, grd, agg_r = 2)
  # expect_true(terra::all.equal(agg_r[[1]], raster::aggregate(r, 2)))
  # expect_true(terra::all.equal(agg_r[[2]], grd))
  expect_equal(terra::values(agg_r[[1]]), terra::values(raster::aggregate(r, 2)))
  expect_equal(terra::values(agg_r[[2]]), terra::values(grd))

  # test dissaggregation of r
  disagg_r <- raster_transform(r, grd, disagg_r = 2)
  # expect_true(terra::all.equal(disagg_r[[1]], terra::disagg(r, 2)))
  # expect_true(terra::all.equal(disagg_r[[2]], grd))
  expect_equal(terra::values(disagg_r[[1]]), terra::values(terra::disagg(r, 2)))
  expect_equal(terra::values(disagg_r[[2]]), terra::values(grd))

  # test resampling of grd
  resample_grd <- raster_transform(r, grd, resample = "grd")
  # remove values by rast() because we don't expect them to be the same
  # expect_true(terra::all.equal(terra::rast(resample_grd[[1]]), terra::rast(r)))
  # expect_true(terra::all.equal(terra::rast(resample_grd[[2]]), terra::rast(r)))
  capture_warnings(expect_equal(terra::values(terra::rast(resample_grd[[1]])), terra::values(terra::rast(r))))
  capture_warnings(expect_equal(terra::values(terra::rast(resample_grd[[2]])), terra::values(terra::rast(r))))

  # test aggregation of grd
  agg_grd <- raster_transform(r, grd, agg_grd = 2)
  # expect_true(terra::all.equal(agg_grd[[1]], r))
  # expect_true(terra::all.equal(agg_grd[[2]], terra::aggregate(grd, 2)))
  expect_equal(terra::values(agg_grd[[1]]), terra::values(r))
  expect_equal(terra::values(agg_grd[[2]]), terra::values(terra::aggregate(grd, 2)))

  # test dissaggregation of grd
  disagg_grd <- raster_transform(r, grd, disagg_grd = 2)
  # expect_true(terra::all.equal(disagg_grd[[1]], r))
  # expect_true(terra::all.equal(disagg_grd[[2]], terra::disagg(grd, 2)))
  expect_equal(terra::values(disagg_grd[[1]]), terra::values(r))
  expect_equal(terra::values(disagg_grd[[2]]), terra::values(terra::disagg(grd, 2)))

  # check error if both disagg and agg are provided
  expect_error(raster_transform(r, grd, agg_r = 2, disagg_r = 2))
  expect_error(raster_transform(r, grd, agg_grd = 2, disagg_grd = 2))
})

test_that("krig_method argument works", {
  data("mini_lyr")
  capture_warnings(expect_output(kpi <- krig_gd(mini_lyr, mini_lyr, krig_method = "ordinary"), "[using ordinary kriging]", fixed = TRUE))
  capture_warnings(expect_output(kpi <- krig_gd(mini_lyr, mini_lyr, krig_method = "universal"), "[using universal kriging]", fixed = TRUE))
  capture_warnings(expect_error(kpi <- krig_gd(mini_lyr, mini_lyr, krig_method = "invalid"), "invalid krig_method specified"))
})

test_that("autoKrige_output argument works", {
  data("mini_lyr")
  capture_warnings(kpi <- krig_gd(mini_lyr, mini_lyr, autoKrige_output = TRUE))

  expect_true(is.list(kpi))
  expect_true(inherits(kpi[["raster"]], "SpatRaster"))
  expect_true(inherits(kpi[["autoKrige_output"]], "autoKrige"))

  capture_warnings(kpi <- krig_gd(raster::stack(mini_lyr, mini_lyr), index = 1:2, mini_lyr, autoKrige_output = TRUE))
  capture_warnings(kpi <- krig_gd(mini_lyr, mini_lyr, autoKrige_output = FALSE))
})

test_that("raster transform check", {
  data("mini_lyr")

  mini_lyr <- terra::rast(mini_lyr)
  bad_stack <- c(mini_lyr, mini_lyr)
  expect_error(kpi <- raster_transform(bad_stack, grd = mini_lyr), ">1 layer provided for r")
  expect_error(kpi <- raster_transform(mini_lyr, grd = bad_stack), ">1 layer provided for grd")

  # check resample_first arg
  expect_error(kpi <- raster_transform(mini_lyr, mini_lyr, resample_first = FALSE), NA)
  expect_error(kpi <- raster_transform(mini_lyr, mini_lyr, resample_first = TRUE), NA)
})

test_that("bound check", {
  data("mini_lyr")
  load_mini_ex(quiet = TRUE)
  capture_warnings(kpi <- krig_gd(mini_lyr, mini_lyr, upper_bound = TRUE, lower_bound = TRUE))
  capture_warnings(kpi <- krig_gd(mini_lyr, mini_lyr, upper_bound = FALSE, lower_bound = FALSE))
  capture_warnings(kpi <- krig_gd(mini_lyr, mini_lyr, upper_bound = 0.5, lower_bound = 0))
  expect_true(min(terra::values(kpi)) >= 0)
  expect_true(max(terra::values(kpi)) <= 0.5)
})
