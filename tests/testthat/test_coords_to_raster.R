test_that("check dimensions of raster produced", {
  load_mini_ex(quiet = TRUE)

  r <- coords_to_raster(mini_coords)

  r1 <- coords_to_raster(mini_coords, buffer = 1)

  expect_true(all(range(terra::ext(r1)) - range(terra::ext(r)) == 2))
})

test_that("check resolution of raster produced", {
  load_mini_ex(quiet = TRUE)

  r <- coords_to_raster(mini_coords, res = 5)
  # check resolution is 5
  expect_equal(terra::res(r), c(5, 5))

  r <- coords_to_raster(mini_coords, res = c(5, 4))
  # check resolution is about 4 , 5 (also confirm order of x, y)
  expect_equal(round(terra::res(r), 0), c(5, 4))

  expect_error(r <- coords_to_raster(mini_coords, res = c(5, 4, 3)), "invalid res provided")
})

test_that("make sure it works if coords are in different formats", {
  load_mini_ex(quiet = TRUE)

  # sf object
  sf_coords <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  rsf <- coords_to_raster(sf_coords)

  # SpatVector
  vect_coords <- terra::vect(sf_coords)
  rvect <- coords_to_raster(sf_coords)

  # mat
  rmat <- coords_to_raster(as.matrix(mini_coords))

  # df
  rdf <- coords_to_raster(mini_coords)

  # expect_true(terra::all.equal(rvect, rmat))
  # expect_true(terra::all.equal(rvect, rdf))
  # expect_true(terra::all.equal(rvect, rsf))
  expect_equal(terra::values(rvect), terra::values(rmat))
  expect_equal(terra::values(rvect), terra::values(rdf))
  expect_equal(terra::values(rvect), terra::values(rsf))
})

test_that("aggregation and disaggregation produce correct rasters", {
  load_mini_ex(quiet = TRUE)

  r0 <- coords_to_raster(mini_coords)

  ra <- coords_to_raster(mini_coords, agg = 2)
  expect_true(terra::compareGeom(ra, terra::aggregate(r0, 2)))

  rd <- coords_to_raster(mini_coords, disagg = 2)
  expect_true(terra::compareGeom(rd, terra::disagg(r0, 2)))

  expect_warning(rad <- coords_to_raster(mini_coords, agg = 2, disagg = 2))
  r0_rad <- terra::disagg(terra::aggregate(r0, 2), 2)
  expect_true(terra::compareGeom(rad, r0_rad))
})

test_that("plots without errors", {
  load_mini_ex(quiet = TRUE)

  expect_error(r <- coords_to_raster(mini_coords, plot = TRUE), NA)
})
