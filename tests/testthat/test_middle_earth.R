test_that("middle earth is loaded succesfully", {
  expect_message(load_middle_earth_ex())
  expect_message(load_middle_earth_ex(quiet = TRUE), NA)
  expect_true(class(lotr_vcf)[1] == "vcfR")
  expect_true(class(lotr_coords) == "data.frame")
  expect_true(class(lotr_lyr)[1] == "RasterLayer")
})

test_that("mini middle earth ex is loaded succesfully", {
  expect_message(load_mini_ex())
  expect_message(load_mini_ex(quiet = TRUE), NA)
  expect_true(class(mini_vcf)[1] == "vcfR")
  expect_true(class(mini_coords) == "data.frame")
  expect_true(class(mini_lyr)[1] == "RasterLayer")
})
