
test_that("check vcf check throws correct errors and return correct object",{
  data("mini_vcf")

  expect_error(vcf_check("BAD/PATH"))
  expect_error(vcf_check(1))
  expect_s4_class(vcf_check(mini_vcf), "vcfR")
})

test_that("conversion functions return correct object types and errors",{
  data("mini_vcf")

  expect_error(vcf_to_dosage("BAD/PATH"))
  expect_error(vcf_to_dosage(1))
  expect_type(vcf_to_dosage(mini_vcf), "integer")

  expect_error(vcf_to_genind("BAD/PATH"))
  expect_error(vcf_to_genind(1))
  expect_warning(expect_warning(expect_s4_class(vcf_to_genind(mini_vcf), "genind")))

})

test_that("pop warnings for genind are correct",{
  data("mini_vcf")

  # good cases
  expect_s4_class(vcf_to_genind(mini_vcf, pops = FALSE), "genind")
  expect_warning(expect_warning(expect_s4_class(vcf_to_genind(mini_vcf, pops = NULL), "genind")))
  expect_s4_class(vcf_to_genind(mini_vcf, pops = 1:10), "genind")
  expect_s4_class(vcf_to_genind(mini_vcf, pops = rep(c("a","b"),5)), "genind")

  # bad cases
  expect_error(vcf_to_genind(mini_vcf, pops = 1:2), "length of pops does not match number of individuals in genind")
  expect_error(vcf_to_genind(mini_vcf, pops = c("a","b","c")), "length of pops does not match number of individuals in genind")

})

test_that("check vcf to dosage matrix conversion is correct", {
  data("mini_vcf")
  dos <- vcf_to_dosage(mini_vcf)
  # get genotype matrix, remove FORMAT col (first col), and transpose so rows are individuals and cols are loci
  gt <- t(mini_vcf@gt[,-1])

  gt0 <- which(gt == "0|0", arr.ind = TRUE)
  gt1 <- rbind(which(gt == "0|1", arr.ind = TRUE), which(gt == "1|0", arr.ind = TRUE))
  gt2 <- which(gt == "1|1", arr.ind = TRUE)

  dos0 <- which(dos == 0, arr.ind = TRUE)
  dos1 <- which(dos == 1, arr.ind = TRUE)
  dos2 <- which(dos == 2, arr.ind = TRUE)

  # Reorder heterozygotes to match (since rowbind was used to create gt1)
  dos1 <- dos1[order(rownames(dos1), dos1[,1], dos1[,2]),]
  gt1 <- gt1[order(rownames(gt1), gt1[,1], gt1[,2]),]

  expect_equal(gt0, dos0)
  expect_equal(gt1, dos1)
  expect_equal(gt2, dos2)
})

test_that("check that vcf to genind conversion is correct", {
  data("mini_vcf")
  expect_warning(genind <- vcf_to_genind(mini_vcf))
  # get table and retain only columns with a .0 in it (for genind objects each allele is a column,
  # so by removing every other you get basically a dosage matrix)
  alleles <- colnames(genind@tab)
  allele0 <- stringr::str_split_fixed(alleles, "[.]", n = 2)
  keep <- alleles[allele0[,2] == "0"]

  tab <- genind@tab[, keep]
  gt <- t(mini_vcf@gt[,-1])

  gt0 <- which(gt == "0|0", arr.ind = TRUE)
  gt1 <- rbind(which(gt == "0|1", arr.ind = TRUE), which(gt == "1|0", arr.ind = TRUE))
  gt2 <- which(gt == "1|1", arr.ind = TRUE)

  # NOTE: genind is coded opposite to dosage (e.g. 0 is 2 and 2 is 0)
  tab0 <- which(tab == 2, arr.ind = TRUE)
  tab1 <- which(tab == 1, arr.ind = TRUE)
  tab2 <- which(tab == 0, arr.ind = TRUE)

  # Reorder heterozygotes to match (since rowbind was used to create gt1)
  tab1 <- tab1[order(rownames(tab1), tab1[,1], tab1[,2]),]
  gt1 <- gt1[order(rownames(gt1), gt1[,1], gt1[,2]),]

  expect_equal(gt0, tab0)
  expect_equal(gt1, tab1)
  expect_equal(gt2, tab2)
})

