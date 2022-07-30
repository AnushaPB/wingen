
test_that("window_gd returns expected output", {
  load_mini_ex()
  wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE)
  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
})

test_that("wdim_check fixes wdim", {
  expect_error(wdim_check(1))
  expect_error(wdim_check(c(2, 3)))

  expect_equal(wdim_check(c(3, 3)), c(3, 3))
  expect_warning(expect_equal(wdim_check(c(3, 4)), c(3, 5)))
  expect_warning(expect_equal(wdim_check(c(4, 3)), c(5, 3)))
  expect_warning(expect_warning(expect_equal(wdim_check(c(4, 4)), c(5, 5))))

  expect_warning(wdim_check(c(3, 4)))
  expect_warning(wdim_check(4))
})

test_that("returns matrix with only one zero", {
  n <- wdim_to_mat(c(3,5))

  # only ones and zeroes
  expect_equal(unique(as.vector(n)), c(1,0))

  # only one zero
  expect_equal(sum(as.vector(n) == 0), 1)

  # zero is at center
  center <- n[nrow(n)/2 + 0.5,ncol(n)/2 + 0.5]
  expect_equal(center, 0)

  n <- wdim_to_mat(3)

  # only ones and zeroes
  expect_equal(unique(as.vector(n)), c(1,0))

  # only one zero
  expect_equal(sum(as.vector(n) == 0), 1)

  # zero is at center
  center <- n[nrow(n)/2 + 0.5,ncol(n)/2 + 0.5]
  expect_equal(center, 0)

})

test_that("biallelic richness is calculated correctly", {
  expected <- c(1, 2, 2, 2, 2, 2, 2, 2, 1)
  all_possible_combos <- t(expand.grid(0:2, 0:2))
  ar_vals <- apply(all_possible_combos, 2, helper_calc_biar)
  expect_equal(ar_vals, expected)
  expect_equal(calc_mean_biar(all_possible_combos), mean(expected))

  expect_error(calc_mean_biar(matrix(c(0:4), nrow = 1)), "to calculate biallelic richness, all values in genetic matrix must be NA, 0, 1 or 2")
})

