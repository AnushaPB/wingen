

test_that("window_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  expect_error(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE), NA)
  expect_s4_class(wpi, "RasterStack")
  expect_equal(raster::nlayers(wpi), 2)
})

test_that("all stats and parallel works", {
  load_mini_ex(quiet = TRUE)
  expect_warning(expect_warning(wp <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat = "pi", rarify = FALSE)))
  expect_warning(expect_warning(wh <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat = "het", rarify = FALSE)))
  expect_warning(expect_warning(wb <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat = "biallelic.richness", rarify = FALSE, rarify_alleles = FALSE)))
  expect_warning(expect_warning(wbr <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat = "biallelic.richness", rarify = FALSE, rarify_alleles = TRUE)))
  expect_warning(expect_warning(wa <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat = "allelic.richness", rarify = FALSE)))

  # check parallel
  expect_warning(expect_warning(wpp <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat = "pi", rarify = FALSE, parallel = TRUE, ncores = 2)))
  names(wpp) <- names(wp) <- NULL
  expect_equal(wpp, wp)
})

test_that("all stats work with just one locus", {
  load_mini_ex(quiet = TRUE)
  expect_warning(wp <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "pi", rarify = FALSE))
  expect_warning(wh <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "het", rarify = FALSE))
  expect_warning(wb <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "biallelic.richness", rarify = FALSE, rarify_alleles = FALSE))
  expect_warning(wb <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "biallelic.richness", rarify = FALSE, rarify_alleles = TRUE))
  expect_warning(wa <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "allelic.richness", rarify = FALSE))
})

test_that("get error with just one individual", {
  load_mini_ex(quiet = TRUE)
  # note: the first col of a vcf is the format col
  expect_error(wp <- window_gd(mini_vcf[, 1:2], mini_coords[1, ], mini_lyr, stat = "pi", rarify = FALSE), "cannot run window_gd with only one individual")
})


test_that("error gets returned for mismatch vcf and coords", {
  load_mini_ex(quiet = TRUE)
  expect_error(wp <- window_gd(mini_vcf, mini_coords[1:2, ], mini_lyr, stat = "pi", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
  expect_error(wa <- window_gd(mini_vcf, mini_coords[1:2, ], mini_lyr, stat = "allelic.richness", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
  expect_error(wp <- window_gd(mini_vcf[, 1:3], mini_coords, mini_lyr, stat = "pi", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
  expect_error(wa <- window_gd(mini_vcf[, 1:3], mini_coords, mini_lyr, stat = "allelic.richness", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
})


test_that("L argument works", {
  load_mini_ex(quiet = TRUE)
  expect_error(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE, L = "nvariants"), NA)
  expect_error(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE, L = 100), NA)
  expect_error(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE, L = NULL), NA)


  expect_error(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, stat = "allelic.richness", rarify = FALSE, L = "nvariants"), NA)
  expect_error(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, stat = "allelic.richness", rarify = FALSE, L = 100), NA)
  expect_error(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, stat = "allelic.richness", rarify = FALSE, L = NULL), NA)
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

test_that("biallelic richness is calculated correctly for all possible combos (including NAs)", {
  all_possible_combos <- t(expand.grid(c(0:2, NA), c(0:2, NA)))
  expected <- c(1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 1, NA)

  ar_vals <- apply(all_possible_combos, 2, helper_calc_biar, rarify_alleles = FALSE)
  expect_equal(ar_vals, expected)
  expect_equal(calc_mean_biar(all_possible_combos, rarify_alleles = FALSE), mean(expected, na.rm = TRUE))

  expect_error(calc_mean_biar(matrix(c(0:4), nrow = 1)), "to calculate biallelic richness, all values in genetic matrix must be NA, 0, 1 or 2")

  # check rarefaction
  min.n <- get_minn(all_possible_combos)
  expect_equal(min.n, 2)

  expected_rar <- c(1, 3 / 2, 5 / 3, 1, 3 / 2, 5 / 3, 1.5, 2, 5 / 3, 3 / 2, 1, 1, 1, 2, 1, NA)
  ar_vals_rar <- apply(all_possible_combos, 2, helper_calc_biar, rarify_alleles = TRUE, min.n = min.n)
  expect_equal(ar_vals_rar, expected_rar)
})


test_that("allelic richness is calculated correctly for dataset with NAs (and rarefaction)", {
  load_mini_ex()

  expect_warning(expect_warning(gen <- vcf_to_genind(mini_vcf_NA)))
  observed_ar <- helper_calc_ar(gen)

  expect_warning(expect_warning(data <- check_data(mini_vcf_NA, mini_coords)))
  dos <- vcf_to_dosage(data$vcf)
  observed_bar <- apply(dos, 2, helper_calc_biar, min.n = get_minn(dos))

  expect_true(all(observed_bar == observed_ar))

  # check rasters
  set.seed(22)
  expect_warning(
    expect_warning(
      trab <- window_gd(mini_vcf_NA,
        mini_coords,
        mini_lyr,
        stat = "biallelic.richness",
        wdim = 3,
        fact = 3,
        rarify_n = 2,
        rarify_nit = 5,
        rarify = TRUE,
        parallel = FALSE,
        rarify_alleles = TRUE
      )
    )
  )


  set.seed(22)
  expect_warning(
    expect_warning(
      tra <- window_gd(mini_vcf_NA,
        mini_coords,
        mini_lyr,
        stat = "allelic.richness",
        wdim = 3,
        fact = 3,
        rarify_n = 2,
        rarify_nit = 5,
        rarify = TRUE,
        parallel = FALSE
      )
    )
  )


  names(tra) <- names(trab) <- NULL
  expect_equal(trab, tra)
})

test_that("allelic richness is calculated correctly for dataset with no NAs", {
  load_mini_ex()

  gen <- vcf_to_genind(mini_vcf)
  observed_ar <- helper_calc_ar(gen)

  data <- check_data(mini_vcf, mini_coords)
  dos <- vcf_to_dosage(data$vcf)
  observed_bar_rar <- apply(dos, 2, helper_calc_biar, min.n = get_minn(dos), rarify_alleles = TRUE)
  observed_bar_norar <- apply(dos, 2, helper_calc_biar, min.n = get_minn(dos), rarify_alleles = FALSE)

  expect_true(all(observed_bar_rar == observed_ar))
  expect_true(all(observed_bar_norar == observed_ar))
  expect_true(all(observed_bar_norar == observed_bar_rar))

  # check rasters
  set.seed(22)
  trab_rar <- window_gd(mini_vcf,
    mini_coords,
    mini_lyr,
    stat = "biallelic.richness",
    wdim = 3,
    fact = 3,
    rarify_n = 2,
    rarify_nit = 5,
    rarify = TRUE,
    parallel = FALSE,
    rarify_alleles = TRUE
  )

  # check rasters
  set.seed(22)
  trab_norar <- window_gd(mini_vcf,
    mini_coords,
    mini_lyr,
    stat = "biallelic.richness",
    wdim = 3,
    fact = 3,
    rarify_n = 2,
    rarify_nit = 5,
    rarify = TRUE,
    parallel = FALSE,
    rarify_alleles = FALSE
  )

  set.seed(22)
  tra <- window_gd(mini_vcf,
    mini_coords,
    mini_lyr,
    stat = "allelic.richness",
    wdim = 3,
    fact = 3,
    rarify_n = 2,
    rarify_nit = 5,
    rarify = TRUE,
    parallel = FALSE
  )

  names(tra) <- names(trab_rar) <- names(trab_norar) <- NULL
  expect_equal(tra, trab_rar)
  expect_equal(tra, trab_norar)
  expect_equal(trab_norar, trab_rar)
})

test_that("vcf path works", {
  load_mini_ex()
  vcfR::write.vcf(mini_vcf, file = "test_temp.vcf")
  vcfpath <- "test_temp.vcf"
  expect_error(wpi <- window_gd(vcfpath, mini_coords, mini_lyr, rarify = FALSE), NA)
  expect_true(file.remove("test_temp.vcf"))
})

test_that("error if bad vcf is given", {
  vcfpath <- "badpath"
  expect_error(wpi <- window_gd(vcfpath, mini_coords, mini_lyr, rarify = FALSE), "Cannot open file: No such file or directory")
  expect_error(wpi <- window_gd(mini_coords, mini_coords, mini_lyr, rarify = FALSE), "Input is expected to be an object of class 'vcfR' or a path to a .vcf file")
})

test_that("return_stat returns correct functions", {
  expect_equal(return_stat("pi"), calc_pi)
  expect_equal(return_stat("het"), calc_mean_het)
  expect_equal(return_stat("biallelic.richness"), calc_mean_biar)
  expect_equal(return_stat("allelic.richness"), calc_mean_ar)
})

test_that("raref works", {
  allele.counts <- c(1, 3)
  expect_equal(raref(allele.counts, min.n = 2), 1.5)
})

test_that("countgen works", {
  all_possible_combos <- t(expand.grid(c(0:2, NA), c(0:2, NA)))
  actual <- apply(all_possible_combos, 2, countgen)
  expected <- c(2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, NA)

  expect_equal(expected, actual)

  actual <- sapply(all_possible_combos[1, ], countgen)
  expected <- c(1, 1, 1, NA, 1, 1, 1, NA, 1, 1, 1, NA, 1, 1, 1, NA)

  expect_equal(expected, actual)
})
