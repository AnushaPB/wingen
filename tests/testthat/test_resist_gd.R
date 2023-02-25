
test_that("resist_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  distmat <- get_resdist(mini_coords, mini_lyr, ncores = 2)
  capture_warnings(rpi <- resist_gd(mini_vcf, mini_coords, mini_lyr,
                                    rarify = FALSE,
                                    maxdist = quantile(distmat, 0.05, na.rm = TRUE),
                                    distmat = distmat))
  expect_s4_class(rpi, "SpatRaster")
  expect_equal(terra::nlyr(rpi), 2)

  # check against expected values
  vals <- terra::global(rpi, fun = "mean", na.rm = TRUE)
  expect_equal(2.737778, vals["pi",], tolerance = 0.000001)
  expect_equal(0.45, vals["sample_count",])
})
