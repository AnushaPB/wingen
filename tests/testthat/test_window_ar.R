

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
