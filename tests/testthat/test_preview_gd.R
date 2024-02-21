test_that("preview_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  pw <- preview_gd(
    mini_lyr,
    mini_coords,
    method = "window",
    wdim = 3,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  expect_equal(terra::nlyr(pw), 1)

  pw <- preview_gd(
    mini_lyr,
    mini_coords,
    wdim = 3,
    fact = 3,
    sample_count = FALSE,
    min_n = 2
  )

  expect_equal(pw, NULL)
})

test_that("preview_gd works for different coord types", {
  load_mini_ex(quiet = TRUE)

  coords_sf <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  coords_mat <- as.matrix(mini_coords)
  coords_vect <- terra::vect(coords_sf)

  pw1 <- preview_gd(
    mini_lyr,
    coords_sf,
    wdim = 3,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  pw2 <- preview_gd(
    mini_lyr,
    coords_mat,
    wdim = 3,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  pw3 <- preview_gd(
    mini_lyr,
    coords_vect,
    wdim = 3,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  # expect_true(terra::all.equal(pw1, pw2))
  # expect_true(terra::all.equal(pw1, pw3))
  expect_equal(terra::values(pw1), terra::values(pw2))
  expect_equal(terra::values(pw1), terra::values(pw3))
})


test_that("preview_gd works for circle method with all different coordinate types", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(distmat <- get_geodist(mini_coords, mini_lyr, fact = 3))
  pw <- preview_gd(
    mini_lyr,
    mini_coords,
    method = "circle",
    maxdist = 50,
    distmat = distmat,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  coords_sf <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  coords_mat <- as.matrix(mini_coords)
  coords_vect <- terra::vect(coords_sf)

  pw1 <- preview_gd(
    mini_lyr,
    coords_sf,
    method = "circle",
    maxdist = 50,
    distmat = distmat,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  pw2 <- preview_gd(
    mini_lyr,
    coords_mat,
    method = "circle",
    maxdist = 50,
    distmat = distmat,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  pw3 <- preview_gd(
    mini_lyr,
    coords_vect,
    method = "circle",
    maxdist = 50,
    distmat = distmat,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  expect_true(terra::all.equal(pw1, pw2))
  expect_true(terra::all.equal(pw1, pw3))
})


test_that("preview_gd works for resist method", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(distmat <- get_resdist(mini_coords, mini_lyr, fact = 5))
  capture_warnings(pw <- preview_gd(
    mini_lyr,
    mini_coords,
    method = "resist",
    maxdist = 50,
    distmat = distmat,
    fact = 5,
    sample_count = TRUE,
    min_n = 2
  ))

  expect_equal(terra::nlyr(pw), 1)
})
