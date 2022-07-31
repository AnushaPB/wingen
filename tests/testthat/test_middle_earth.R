
test_that("middle earth is loaded succesfully", {
  expect_message(load_middle_earth_ex())
  expect_true(class(lotr_vcf)[1] == "vcfR")
  expect_true(class(lotr_coords) == "data.frame")
  expect_true(class(lotr_lyr)[1] == "RasterLayer")
})
