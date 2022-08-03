test_that("krig_gd returns expected output", {
  library("raster")
  load_mini_ex(quiet = TRUE)
  expect_warning(kpi <- krig_gd(mini_lyr, mini_lyr))
  expect_s4_class(mini_lyr, "RasterLayer")
  expect_equal(raster::nlayers(kpi), 1)

  expect_warning(expect_warning(kpi <- krig_gd(raster::stack(mini_lyr, mini_lyr), index = 1:2, mini_lyr)))
  expect_s4_class(kpi, "RasterStack")
  expect_equal(raster::nlayers(mini_lyr), 1)
})

test_that("krig_gd returns warning when not provided grd", {
  library("raster")
  load_mini_ex(quiet = TRUE)
  expect_error(warnings <- capture_warnings(kpi <- krig_gd(mini_lyr, grd = NULL)), NA)
  # clean this up later:
  expect_equal(warnings[1], "no grd provided, defaults to using first raster layer to create grd")
})

test_that("coord kriging works", {
  library("raster")
  load_mini_ex(quiet = TRUE)
  expect_error(kpi <- krig_gd(mini_lyr, grd = lotr_lyr, coords = lotr_coords), NA)
})

test_that("grd kriging works", {
  library("raster")
  load_mini_ex(quiet = TRUE)

  grd <- raster_to_grid(mini_lyr)

  expect_error(expect_warning(kpi <- krig_gd(mini_lyr, grd = grd)), NA)
})


test_that("krige_gd returns error when provided bad grd", {
  library("raster")
  load_mini_ex(quiet = TRUE)

  expect_error(kpi <- krig_gd(mini_lyr, grd = mini_coords))
  expect_error(kpi <- krig_gd(mini_lyr, grd = sp::SpatialPoints(mini_coords)), "unable to find an inherited method for type of grd provided")
})


test_that("krig_gd returns warning when provided crs", {
  library("raster")
  load_mini_ex(quiet = TRUE)
  # clean this up later:
  grd <- raster_to_grid(mini_lyr)
  raster::crs(grd) <- "+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62"

  warnings <- capture_warnings(kpi <- krig_gd(mini_lyr, grd))
  expect_equal(warnings[1], "the provided raster and grid have different crs")
  expect_equal(warnings[2], "NaNs produced")
})


test_that("raster_transform transformations are correct", {
  data("mini_lyr")
  grd <- mini_lyr
  r <- raster::aggregate(mini_lyr, 2)

  # test pass through
  noTransform <- raster_transform(r, grd)
  expect_true(raster::compareRaster(noTransform[[1]], r, stopiffalse = FALSE))
  expect_true(raster::compareRaster(noTransform[[2]], grd, stopiffalse = FALSE))

  # test resampling of r
  resample_r <- raster_transform(r, grd, resample = "r")
  expect_true(raster::compareRaster(resample_r[[1]], grd, stopiffalse = FALSE))
  expect_true(raster::compareRaster(resample_r[[2]], grd, stopiffalse = FALSE))

  # test aggregation of r
  agg_r <- raster_transform(r, grd, agg_r = 2)
  expect_true(raster::compareRaster(agg_r[[1]], raster::aggregate(r, 2), stopiffalse = FALSE))
  expect_true(raster::compareRaster(agg_r[[2]], grd, stopiffalse = FALSE))

  # test dissaggregation of r
  disagg_r <- raster_transform(r, grd, disagg_r = 2)
  expect_true(raster::compareRaster(disagg_r[[1]], raster::disaggregate(r, 2), stopiffalse = FALSE))
  expect_true(raster::compareRaster(disagg_r[[2]], grd, stopiffalse = FALSE))

  # test resampling of grd
  resample_grd <- raster_transform(r, grd, resample = "grd")
  expect_true(raster::compareRaster(resample_grd[[1]], r, stopiffalse = FALSE))
  expect_true(raster::compareRaster(resample_grd[[2]], r, stopiffalse = FALSE))

  # test aggregation of grd
  agg_grd <- raster_transform(r, grd, agg_grd = 2)
  expect_true(raster::compareRaster(agg_grd[[1]], r, stopiffalse = FALSE))
  expect_true(raster::compareRaster(agg_grd[[2]], raster::aggregate(grd, 2), stopiffalse = FALSE))

  # test dissaggregation of grd
  disagg_grd <- raster_transform(r, grd, disagg_grd = 2)
  expect_true(raster::compareRaster(disagg_grd[[1]], r, stopiffalse = FALSE))
  expect_true(raster::compareRaster(disagg_grd[[2]], raster::disaggregate(grd, 2), stopiffalse = FALSE))

  # check error if both disagg and agg are provided
  expect_error(raster_transform(r, grd, agg_r = 2, disagg_r = 2))
  expect_error(raster_transform(r, grd, agg_grd = 2, disagg_grd = 2))
})

test_that("xy argument works", {
  library("raster")
  data("mini_lyr")
  expect_warning(kpi <- krig_gd(mini_lyr, mini_lyr, xy = TRUE))

  expect_warning(kpi <- krig_gd(mini_lyr, mini_lyr, xy = FALSE))
})


test_that("raster transform check", {
  library("raster")
  data("mini_lyr")

  bad_stack <- raster::stack(mini_lyr, mini_lyr)
  expect_error(kpi <- raster_transform(bad_stack, grd = mini_lyr), ">1 layer provided for r")
  expect_error(kpi <- raster_transform(mini_lyr, grd = bad_stack), ">1 layer provided for grd")

  # check resample_first arg
  expect_error(kpi <- raster_transform(mini_lyr, mini_lyr, resample_first = FALSE), NA)
  expect_error(kpi <- raster_transform(mini_lyr, mini_lyr, resample_first = TRUE), NA)
})
