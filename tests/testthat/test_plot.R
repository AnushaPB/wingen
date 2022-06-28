test_that("incorrect inputs to plot produces correct errors", {
  data("mini_lyr")
  bad_stack <- raster::stack(mini_lyr, mini_lyr, mini_lyr)

  expect_error(plot_gd(bad_stack))
  expect_error(plot_count(bad_stack))
})
