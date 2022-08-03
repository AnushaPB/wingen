test_that("preview_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  pw <- preview_gd(mini_lyr,
    mini_coords,
    wdim = 3,
    fact = 3,
    sample_count = TRUE,
    min_n = 2
  )

  expect_equal(raster::nlayers(pw), 2)

  pw <- preview_gd(mini_lyr,
    mini_coords,
    wdim = 3,
    fact = 3,
    sample_count = FALSE,
    min_n = 2
  )

  expect_equal(raster::nlayers(pw), 1)
  expect_equal(raster::unique(pw), c(0, 1, 2))
})
