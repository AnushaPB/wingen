
test_that("circle_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(cpi <- circle_gd(mini_vcf, mini_coords, mini_lyr, maxdist = 50))
  expect_s4_class(cpi, "SpatRaster")
  expect_equal(terra::nlyr(cpi), 2)

  # check against expected values
  vals <- terra::global(cpi, fun = "mean", na.rm = TRUE)
  expect_equal(3.658471, vals["pi",], tolerance = 0.000001)
  expect_equal(5.12, vals["sample_count",])
})
