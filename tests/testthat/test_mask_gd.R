test_that("mask_gd returns expected output", {
  library("raster")
  load_mini_ex()
  min_n <- 2
  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  kpi <- krig_gd(wpi, mini_lyr)
  mpi <- mask_gd(kpi, min_n)

  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
  expect_true(cellStats(is.na(mpi) == (kpi[[2]] < min_n), all))
})
