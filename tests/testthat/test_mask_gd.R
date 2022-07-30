test_that("mask_gd returns expected output", {
  load_mini_ex()
  min_n <- 2
  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  expect_warning(kpi <- krig_gd(wpi, mini_lyr))
  mpi <- mask_gd(kpi, mini_lyr, minval = 2)

  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
  expect_true(cellStats(is.na(mpi) == (mini_lyr < 2), all))
})
