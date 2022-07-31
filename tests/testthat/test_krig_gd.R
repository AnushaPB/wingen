test_that("krig_gd returns expected output", {
  library("raster")
  load_mini_ex()
  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  expect_warning(kpi <- krig_gd(wpi, mini_lyr))
  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
})

test_that("krig_gd returns warning when not provided grd", {
  library("raster")
  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)

  warnings <- capture_warnings(kpi <- krig_gd(wpi))
  # clean this up later:
  expect_equal(warnings[1], "no grd provided, defaults to using first raster layer to create grd")
})


test_that("krig_gd returns warning when provided crs", {
  library("raster")
  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  # clean this up later:
  grd <- raster_to_grid(mini_lyr)
  raster::crs(grd) <- "+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62"

  warnings <- capture_warnings(kpi <- krig_gd(wpi, grd))
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

