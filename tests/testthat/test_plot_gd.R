test_that("check that background plotting works", {
  data("mini_lyr")

  # check that there is NO error (NA)
  expect_error(plot_gd(mini_lyr, bkg = mini_lyr), NA)
})

test_that("check that index works", {
  data("mini_lyr")

  stk <- raster::stack(mini_lyr, mini_lyr)

  # check that there is NO error (NA)
  expect_error(plot_gd(stk, index = 1, bkg = mini_lyr), NA)
  expect_error(plot_gd(stk, index = 2, bkg = mini_lyr), NA)
  expect_error(plot_gd(stk, index = 1:3, bkg = mini_lyr))
})

test_that("check that background plot function works", {
  data("mini_lyr")

  stk <- raster::stack(mini_lyr, mini_lyr)

  # check that there is NO error (NA)
  expect_error(plot_bkg(index = 1, stk, bkg = mini_lyr), NA)
  expect_error(plot_bkg(index = 2, stk, bkg = mini_lyr), NA)
  expect_error(plot_bkg(index = 1:2, stk, bkg = mini_lyr), NA)
})

test_that("check that plot count function works", {
  data("mini_lyr")

  stk <- raster::stack(mini_lyr, mini_lyr)

  # check that there is NO error (NA)
  expect_error(plot_count(mini_lyr), NA)
  expect_error(plot_count(stk), NA)
  expect_error(plot_count(stk, index = 1), NA)
})
