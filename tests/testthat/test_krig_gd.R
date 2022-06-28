test_that("krig_gd returns expected output", {
  library("raster")
  load_mini_ex()
  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  kpi <- krig_gd(wpi, mini_lyr)
  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
})
