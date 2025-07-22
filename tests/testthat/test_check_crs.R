test_that("CRS are handled correctly by window_gd", {
  capture_warnings(library(raster))
  load_mini_ex(quiet = TRUE)

  # NO CRS
  nocrs_lyr <- mini_lyr
  nocrs_coords <- mini_coords

  # CRS
  crs_lyr <- terra::rast(mini_lyr)
  terra::crs(crs_lyr) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
  crs_coords <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  sf::st_crs(crs_coords) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"

  # no CRS lyr and coords
  cw <- capture_warnings(window_gd(mini_vcf, nocrs_coords, nocrs_lyr, rarify = FALSE))
  expect_true(any(cw == "No CRS found for the provided coordinates. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)\n"))
  expect_true(any(cw == "No CRS found for the provided raster. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)\n"))

  # no CRS lyr and CRS coords
  cw <- capture_warnings(window_gd(mini_vcf, crs_coords, nocrs_lyr, rarify = FALSE))
  expect_true(any(cw == "No CRS found for the provided raster. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)\n"))

  # CRS lyr and no CRS coords
  cw <- capture_warnings(window_gd(mini_vcf, nocrs_coords, crs_lyr, rarify = FALSE))
  expect_true(any(cw == "No CRS found for the provided coordinates. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)\n"))

  # CRS lyr and CRS coords (no warnings expected)
  wg <- window_gd(mini_vcf, crs_coords, crs_lyr, rarify = FALSE)
  expect_true(sf::st_crs(wg) == sf::st_crs(crs_lyr))

  # mismatched CRS
  crs2_lyr <- mini_lyr
  terra::crs(crs2_lyr) <- "+proj=longlat +ellps=GRS80 +datum=NAD83
  +no_defs +towgs84=0,0,0"
  expect_error(window_gd(mini_vcf, crs_coords, crs2_lyr, rarify = FALSE), "CRS of the provided coordinates and raster do not match")
})


test_that("CRS are handled correctly by winkrige_gd (r vs grd)", {
  capture_warnings(library(raster))
  load_mini_ex(quiet = TRUE)
  mini_lyr <- terra::rast(mini_lyr)

  # NO CRS
  nocrs_lyr <- mini_lyr
  # CRS
  crs_lyr <- mini_lyr
  terra::crs(crs_lyr) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"

  # no CRS lyr and coords
  cw <- capture_warnings(winkrige_gd(r = nocrs_lyr, grd = nocrs_lyr))
  expect_true(any(cw == "No CRS found for the provided raster (r). Make sure that r and grd have the same projection.\n"))
  expect_true(any(cw == "No CRS found for the provided raster (grd). Make sure that r and grd have the same projection.\n"))

  # no CRS lyr and CRS coords
  cw <- capture_warnings(winkrige_gd(r = nocrs_lyr, grd = crs_lyr))
  expect_true(any(cw == "No CRS found for the provided raster (r). Make sure that r and grd have the same projection.\n"))

  # CRS lyr and no CRS coords
  cw <- capture_warnings(winkrige_gd(r = crs_lyr, grd = nocrs_lyr))
  expect_true(any(cw == "No CRS found for the provided raster (grd). Make sure that r and grd have the same projection.\n"))

  # CRS lyr and CRS coords (no warnings expected)
  capture_warnings(kg <- winkrige_gd(r = crs_lyr, grd = crs_lyr))
  expect_true(sf::st_crs(kg) == sf::st_crs(crs_lyr))

  # mismatched CRS
  crs2_lyr <- mini_lyr
  terra::crs(crs2_lyr) <- "+proj=longlat +ellps=GRS80 +datum=NAD83
  +no_defs +towgs84=0,0,0"
  expect_error(winkrige_gd(r = crs2_lyr, grd = crs_lyr))
  expect_error(winkrige_gd(grd = crs2_lyr, r = crs_lyr))
})
