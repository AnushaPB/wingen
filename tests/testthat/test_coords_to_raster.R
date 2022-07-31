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

test_that("make sure it works if coords are in different formats", {
  load_mini_ex()
  rdf <- coords_to_raster(mini_coords)
  rmat <- coords_to_raster(as.matrix(mini_coords))
  expect_equal(rdf, rmat)
})
