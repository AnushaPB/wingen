test_that("check that background plotting works", {
  data("mini_lyr")

  # check that there is NO error (NA)
  expect_error(ggplot_gd(mini_lyr, bkg = mini_lyr), NA)
})

test_that("check that index works", {
  data("mini_lyr")

  stk <- raster::stack(mini_lyr, mini_lyr)

  # check that there is NO error (NA)
  expect_error(ggplot_gd(stk, index = 1, bkg = mini_lyr), NA)
  expect_error(ggplot_gd(stk, index = 2, bkg = mini_lyr), NA)
})

test_that("check that plot count function works", {
  data("mini_lyr")

  stk <- raster::stack(mini_lyr, mini_lyr)

  # check that there is NO error (NA)
  expect_error(ggplot_count(mini_lyr), NA)
  expect_error(ggplot_count(stk), NA)
  expect_error(ggplot_count(stk, index = 1), NA)
})
