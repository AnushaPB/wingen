test_that("check dimensions of raster produced", {
  load_mini_ex()

  r <- coords_to_raster(mini_coords)

  buffer <- 0
  xmin <- min(mini_coords$x, na.rm = TRUE) - buffer
  xmax <- max(mini_coords$x, na.rm = TRUE) + buffer
  ymin <- min(mini_coords$y, na.rm = TRUE) - buffer
  ymax <- max(mini_coords$y, na.rm = TRUE) + buffer

  expect_true(all(dim(r) == c(round(ymax - ymin, 0), round(xmax - xmin, 0), 1)))

  r <- coords_to_raster(mini_coords, buffer = 1)

  buffer <- 1
  xmin <- min(mini_coords$x, na.rm = TRUE) - buffer
  xmax <- max(mini_coords$x, na.rm = TRUE) + buffer
  ymin <- min(mini_coords$y, na.rm = TRUE) - buffer
  ymax <- max(mini_coords$y, na.rm = TRUE) + buffer

  expect_true(all(dim(r) == c(round(ymax - ymin, 0), round(xmax - xmin, 0), 1)))
})

test_that("check resolution of raster produced", {
  load_mini_ex()

  r <- coords_to_raster(mini_coords, res = 5)
  # check resolution is about 5
  expect_equal(round(raster::res(r), 0), c(5, 5))

  r <- coords_to_raster(mini_coords, res = c(5, 4))
  # check resolution is about 4 , 5 (also confirm order of x, y)
  expect_equal(round(raster::res(r), 0), c(5, 4))

  expect_error(r <- coords_to_raster(mini_coords, res = c(5, 4, 3)), "invalid res provided")
})

test_that("make sure it works if coords are in different formats", {
  load_mini_ex()

  rdf <- coords_to_raster(mini_coords)
  rmat <- coords_to_raster(as.matrix(mini_coords))
  expect_equal(rdf, rmat)
})

test_that("aggregation and disaggregation produce correct rasters", {
  load_mini_ex()

  r0 <- coords_to_raster(mini_coords)

  ra <- coords_to_raster(mini_coords, agg = 2)
  expect_true(raster::compareRaster(ra, raster::aggregate(r0, 2)))


  rd <- coords_to_raster(mini_coords, disagg = 2)
  expect_true(raster::compareRaster(rd, raster::disaggregate(r0, 2)))

  expect_warning(rad <- coords_to_raster(mini_coords, agg = 2, disagg = 2))
  expect_true(raster::compareRaster(rad, r0))
})

test_that("plots without errors", {
  load_mini_ex()

  expect_error(r <- coords_to_raster(mini_coords, plot = TRUE), NA)
})
