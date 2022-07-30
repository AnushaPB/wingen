
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
  # get table and retain only unique columns (for genind objects each allele is a column,
  # so by removing every other you get basically a dosage matrix)
  alleles <- colnames(genind@tab)
  allele0 <- stringr::str_split_fixed(alleles, "[.]", n = 2)
  keep <- alleles[!duplicated(allele0[,1])]

  # check number of loci match
  expect_equal(length(keep), nrow(mini_vcf@gt))

  # get tables of genotypes
  tab <- genind@tab[, keep]
  gt <- t(mini_vcf@gt[,-1])

  # since for genind a homozygote can either be 2 or 0 depending on the ref allele selected make
  # two dosage matrices for each possible ref allele
  gt_dos1 <- gt
  gt_dos1[which(gt == "0|0", arr.ind = TRUE)] <- 2
  gt_dos1[rbind(which(gt == "0|1", arr.ind = TRUE), which(gt == "1|0", arr.ind = TRUE))] <- 1
  gt_dos1[which(gt == "1|1", arr.ind = TRUE)] <- 0

  gt_dos2 <- gt
  gt_dos2[which(gt == "0|0", arr.ind = TRUE)] <- 0
  gt_dos2[rbind(which(gt == "0|1", arr.ind = TRUE), which(gt == "1|0", arr.ind = TRUE))] <- 1
  gt_dos2[which(gt == "1|1", arr.ind = TRUE)] <- 2

  # turn into character for comparison
  tab <- apply(tab, 2, as.character)

  # make logical comparison of possible reference alleles (0 = both FALSE, 1 = at least one TRUE, 2 = both true)
  sum_log <- (gt_dos1 == tab) + (gt_dos2 == tab)

  # checks that in at least one comparison the result is true
  expect_true(all(sum_log > 0))

  # checks that all of the cases where true occurs twice are heterozygotes (e.g. 1)
  expect_true(all((sum_log == 2) == (tab == "1")))
})

