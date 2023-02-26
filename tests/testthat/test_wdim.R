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

test_that("wdim_to_mat returns correct errors", {
  expect_error(wdim_to_mat(c(3, 3, 3)), "wdim must be a single integer or a vector two integers")
})

test_that("wdim returns matrix with only one zero", {
  n <- wdim_to_mat(c(3, 5))

  # only ones and zeroes
  expect_equal(unique(as.vector(n)), c(1, 0))

  # only one zero
  expect_equal(sum(as.vector(n) == 0), 1)

  # zero is at center
  center <- n[nrow(n) / 2 + 0.5, ncol(n) / 2 + 0.5]
  expect_equal(center, 0)

  n <- wdim_to_mat(3)

  # only ones and zeroes
  expect_equal(unique(as.vector(n)), c(1, 0))

  # only one zero
  expect_equal(sum(as.vector(n) == 0), 1)

  # zero is at center
  center <- n[nrow(n) / 2 + 0.5, ncol(n) / 2 + 0.5]
  expect_equal(center, 0)
})
