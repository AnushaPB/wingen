test_that("window_gd returns expected output", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE))
  expect_s4_class(wpi, "SpatRaster")
  expect_equal(terra::nlyr(wpi), 2)

  # check against expected values
  vals <- terra::global(wpi, fun = "mean", na.rm = TRUE)
  expect_equal(0.2675758, vals["pi", ], tolerance = 0.000001)
  expect_equal(0.87, vals["sample_count", ])
})

test_that("all stats and parallel works", {
  load_mini_ex(quiet = TRUE)

  stat <- c("pi", "Ho", "biallelic_richness", "allelic_richness", "basic_stats", "hwe")

  capture_warnings(wg <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat, rarify = FALSE))

  # check parallel
  future::plan("multisession", workers = 2)
  capture_warnings(wgp <- window_gd(mini_vcf_NA, mini_coords, mini_lyr, stat, rarify = FALSE))
  future::plan("sequential")

  # expect_true(terra::all.equal(wgp, wg))
  expect_equal(terra::values(wgp), terra::values(wg))
})

test_that("rarifaction works for all options", {
  load_mini_ex(quiet = TRUE)

  capture_warnings(window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE, rarify_nit = 100))
  capture_warnings(window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE, rarify_nit = 10))
  capture_warnings(wpi1 <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE, rarify_nit = 1))
  capture_warnings(wpiall <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE, rarify_nit = "all"))

  # test doesn't matter, just checking if above lines don't error
  expect_true(TRUE)
})

test_that("check that setting the seed produces the same results", {
  load_mini_ex(quiet = TRUE)

  set.seed(42)
  capture_warnings(wg1 <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE))
  set.seed(42)
  capture_warnings(wg2 <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE))
  # expect_true(terra::all.equal(wg1, wg1))
  expect_equal(terra::values(wg1), terra::values(wg2))


  future::plan("multisession", workers = 2)
  set.seed(42)
  capture_warnings(wg1p <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE))
  set.seed(42)
  capture_warnings(wg2p <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE))
  future::plan("sequential")
  # expect_true(terra::all.equal(wg1p, wg1p))
  expect_equal(terra::values(wg1p), terra::values(wg2p))
})

test_that("all stats work with just one site", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(wp <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "pi", rarify = FALSE))
  capture_warnings(wh <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "Ho", rarify = FALSE))
  capture_warnings(wb <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = FALSE))
  capture_warnings(wb <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = TRUE))
  capture_warnings(wa <- window_gd(mini_vcf_NA[8, ], mini_coords, mini_lyr, stat = "allelic_richness", rarify = FALSE))

  expect_true(TRUE)
})

test_that("get error with just one individual", {
  load_mini_ex(quiet = TRUE)
  # note: the first col of a vcf is the format col
  expect_error(wp <- window_gd(mini_vcf[, 1:2], mini_coords[1, ], mini_lyr, stat = "pi", rarify = FALSE), "cannot run window_gd with only one individual")
})

test_that("window_gd returns expected value", {
  load_mini_ex()
  vcf <- mini_vcf[c(1, 4, 8), 1:3]
  coords <- mini_coords[1:2, ]

  capture_warnings(wg <- window_gd(vcf, coords, mini_lyr, stat = "pi", wdim = 5, min_n = 2))
  dos <- vcf_to_dosage(vcf)
  expect_true(calc_pi(dos, L = ncol(dos)) == unique(na.omit(raster::values(wg[[1]]))))

  capture_warnings(wg <- window_gd(vcf, coords, stat = "Ho", mini_lyr, wdim = 5, min_n = 2))
  het <- vcf_to_het(vcf)
  expect_true(calc_mean_het(het) == unique(na.omit(raster::values(wg[[1]]))))

  capture_warnings(wg <- window_gd(vcf, coords, stat = "biallelic_richness", wdim = 5, mini_lyr, min_n = 2, rarify_alleles = TRUE))
  expect_true(calc_mean_biar(dos, rarify_alleles = TRUE) == unique(na.omit(raster::values(wg[[1]]))))

  capture_warnings(wg <- window_gd(vcf, coords, stat = "biallelic_richness", wdim = 5, mini_lyr, min_n = 2, rarify_alleles = FALSE))
  expect_true(calc_mean_biar(dos, rarify_alleles = FALSE) == unique(na.omit(raster::values(wg[[1]]))))

  capture_warnings(wg <- window_gd(vcf, coords, stat = "allelic_richness", wdim = 5, mini_lyr, min_n = 2))
  genind <- vcfR::vcfR2genind(vcf)
  expect_true(calc_mean_ar(genind) == unique(na.omit(raster::values(wg[[1]]))))
})

test_that("error gets returned for mismatch vcf and coords", {
  load_mini_ex(quiet = TRUE)
  expect_error(wp <- window_gd(mini_vcf, mini_coords[1:2, ], mini_lyr, stat = "pi", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
  expect_error(wa <- window_gd(mini_vcf, mini_coords[1:2, ], mini_lyr, stat = "allelic_richness", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
  expect_error(wp <- window_gd(mini_vcf[, 1:3], mini_coords, mini_lyr, stat = "pi", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
  expect_error(wa <- window_gd(mini_vcf[, 1:3], mini_coords, mini_lyr, stat = "allelic_richness", rarify = FALSE), "number of samples in coords data and number of samples in gen data are not equal")
})


test_that("L argument works", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(wpi_Lnv <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE, L = "nvariants"))
  capture_warnings(wpi_L1k <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE, L = 1000))
  capture_warnings(wpi_LNULL <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE, L = NULL))

  mean_Lnv <- mean(terra::values(wpi_Lnv[[1]]), na.rm = TRUE)
  mean_L1k <- mean(terra::values(wpi_L1k[[1]]), na.rm = TRUE)
  mean_LNULL <- mean(terra::values(wpi_LNULL[[1]]), na.rm = TRUE)

  expect_equal(mean_Lnv, mean_L1k * 100)
  expect_equal(mean_Lnv, mean_LNULL / 10)
})


test_that("biallelic richness is calculated correctly for all possible combos (including NAs)", {
  all_possible_combos <- t(expand.grid(c(0:2, NA), c(0:2, NA)))
  expected <- c(1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 1, NA)

  ar_vals <- apply(all_possible_combos, 2, helper_calc_biar, rarify_alleles = FALSE)
  expect_equal(ar_vals, expected)
  expect_true(calc_mean_biar(all_possible_combos, rarify_alleles = FALSE) == mean(expected, na.rm = TRUE))

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

  expect_warning(expect_warning(gen <- vcfR::vcfR2genind(mini_vcf_NA)))
  observed_ar <- helper_calc_ar(gen)

  expect_warning(expect_warning(data <- check_data(mini_vcf_NA, mini_coords)))
  dos <- vcf_to_dosage(data$vcf)
  observed_bar <- apply(dos, 2, helper_calc_biar, min.n = get_minn(dos))

  expect_true(all(observed_bar == observed_ar))

  # check rasters
  set.seed(22)
  capture_warnings(
    trab <- window_gd(mini_vcf_NA,
      mini_coords,
      mini_lyr,
      stat = "biallelic_richness",
      wdim = 3,
      fact = 3,
      rarify_n = 2,
      rarify_nit = 5,
      rarify = TRUE,
      rarify_alleles = TRUE
    )
  )


  set.seed(22)
  capture_warnings(
    tra <- window_gd(mini_vcf_NA,
      mini_coords,
      mini_lyr,
      stat = "allelic_richness",
      wdim = 3,
      fact = 3,
      rarify_n = 2,
      rarify_nit = 5,
      rarify = TRUE
    )
  )

  names(tra) <- names(trab)


  # expect_true(terra::all.equal(trab, tra))
  expect_equal(terra::values(trab), terra::values(tra))
})

test_that("allelic richness is calculated correctly for dataset with no NAs", {
  load_mini_ex()

  gen <- vcfR::vcfR2genind(mini_vcf)
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

  capture_warnings(
    trab_rar <- window_gd(mini_vcf,
      mini_coords,
      mini_lyr,
      stat = "biallelic_richness",
      wdim = 3,
      fact = 3,
      rarify_n = 2,
      rarify_nit = 5,
      rarify = TRUE,
      rarify_alleles = TRUE
    )
  )

  # check rasters
  set.seed(22)
  capture_warnings(
    trab_norar <- window_gd(mini_vcf,
      mini_coords,
      mini_lyr,
      stat = "biallelic_richness",
      wdim = 3,
      fact = 3,
      rarify_n = 2,
      rarify_nit = 5,
      rarify = TRUE,
      rarify_alleles = FALSE
    )
  )

  set.seed(22)
  capture_warnings(
    tra <- window_gd(mini_vcf,
      mini_coords,
      mini_lyr,
      stat = "allelic_richness",
      wdim = 3,
      fact = 3,
      rarify_n = 2,
      rarify_nit = 5,
      rarify = TRUE
    )
  )

  names(tra) <- names(trab_rar) <- names(trab_norar)

  # expect_true(terra::all.equal(tra, trab_rar))
  expect_equal(terra::values(tra), terra::values(trab_rar))

  # expect_true(terra::all.equal(tra, trab_norar))
  expect_equal(terra::values(tra), terra::values(trab_norar))

  # expect_true(terra::all.equal(trab_norar, trab_rar))
  expect_equal(terra::values(trab_norar), terra::values(trab_rar))
})

test_that("vcf path works", {
  load_mini_ex()
  vcfR::write.vcf(mini_vcf, file = "test_temp.vcf")
  vcfpath <- "test_temp.vcf"
  capture_warnings(wpi <- window_gd(vcfpath, mini_coords, mini_lyr, rarify = FALSE))
  expect_true(file.remove("test_temp.vcf"))
})

test_that("error if bad vcf is given", {
  vcfpath <- "badpath"
  expect_error(capture_warnings(wpi <- window_gd(vcfpath, mini_coords, mini_lyr, rarify = FALSE)), "Cannot open file: No such file or directory")
  expect_error(capture_warnings(wpi <- window_gd(mini_coords, mini_coords, mini_lyr, rarify = FALSE)), "Input is expected to be an object of class 'vcfR' or a path to a .vcf file")
})

test_that("return_stat returns correct functions", {
  expect_equal(return_stat("pi"), calc_pi)
  expect_equal(return_stat("Ho"), calc_mean_het)
  expect_equal(return_stat("biallelic_richness"), calc_mean_biar)
  expect_equal(return_stat("allelic_richness"), calc_mean_ar)

  # check custom functions work
  x <- c(1, 2, NA)
  foo <- return_stat(mean, na.rm = TRUE)
  expect_equal(foo(x), 1.5)
  foo <- return_stat(mean, na.rm = FALSE)
  expect_true(is.na(foo(x)))
  foo_null <- function(x) NULL
  foo <- return_stat(foo_null, na.rm = FALSE)
  expect_error(foo(x))
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

test_that("invariant warning is given", {
  data("mini_vcf_NA")
  # check for one site
  ## no NA
  invariant_vcf <- mini_vcf_NA[3, c(1:5)]
  expect_warning(check_data(invariant_vcf), "invariant sites found in vcf")
  # some NA
  invariant_vcf <- mini_vcf_NA[3, c(1:8)]
  expect_warning(expect_warning(check_data(invariant_vcf), "invariant sites found in vcf"))
  # all NA
  invariant_vcf <- mini_vcf_NA[9, ]
  expect_warning(expect_warning(check_data(invariant_vcf), "invariant sites found in vcf"))
})

test_that("custom functions with window general work", {
  load_mini_ex(quiet = TRUE)
  expect_warning(expect_warning(data <- check_data(mini_vcf_NA, mini_coords)))
  mini_vcf_NA <- data$vcf
  mini_coords <- data$coords
  capture_warnings(wp <- window_general(vcf_to_dosage(mini_vcf_NA), mini_coords, mini_lyr, stat = "pi", rarify = FALSE))
  capture_warnings(wh <- window_general(vcf_to_het(mini_vcf_NA), mini_coords, mini_lyr, stat = "Ho", rarify = FALSE))
  capture_warnings(wb <- window_general(vcf_to_dosage(mini_vcf_NA), mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = FALSE))
  capture_warnings(wbr <- window_general(vcf_to_dosage(mini_vcf_NA), mini_coords, mini_lyr, stat = "biallelic_richness", rarify = FALSE, rarify_alleles = TRUE))
  capture_warnings(wa <- window_general(vcfR::vcfR2genind(mini_vcf_NA), mini_coords, mini_lyr, stat = "allelic_richness", rarify = FALSE))

  # examples with custom functions
  toy <- vcf_to_dosage(mini_vcf_NA) * 0 + 1
  toy[1:2, ] <- NA

  # test on vector
  capture_warnings(wm <- window_general(toy[, 1], mini_coords, mini_lyr, stat = mean, na.rm = TRUE))

  # test on matrix
  capture_warnings(wm <- window_general(toy, mini_coords, mini_lyr, stat = mean, na.rm = TRUE))

  # check NA removal
  foo <- function(x) sum(apply(x, 2, sum, na.rm = TRUE), na.rm = TRUE)
  capture_warnings(wm_1 <- window_general(toy, mini_coords, mini_lyr, stat = foo))
  foo <- function(x, na.rm = FALSE) sum(apply(x, 2, sum, na.rm = na.rm), na.rm = na.rm)
  capture_warnings(wm_2 <- window_general(toy, mini_coords, mini_lyr, stat = foo, na.rm = TRUE))
  expect_true(terra::all.equal(wm_1, wm_2))

  # check if additional custom arguments provided work
  foo <- function(x, silly) sum(x * silly, na.rm = TRUE)
  capture_warnings(wm_1 <- window_general(toy[, 3], mini_coords, mini_lyr, stat = foo, silly = 2)[[1]])
  capture_warnings(wm_2 <- window_general(toy[, 3], mini_coords, mini_lyr, stat = foo, silly = 1)[[1]])
  expect_equal(terra::values(wm_1), terra::values(wm_2) * 2)
})

test_that("get_adj works", {
  load_mini_ex(quiet = TRUE)
  data("mini_lyr")
  data("mini_coords")

  mini_lyr <- terra::rast(mini_lyr)
  nmat <- wdim_to_mat(5)
  coord_cells <- raster::extract(mini_lyr, mini_coords, cell = TRUE)[, "cell"]
  adj <- get_adj(i = coord_cells[1], mini_lyr, nmat, coord_cells)

  # fill in window of raster layer
  adjc <- raster::adjacent(mini_lyr, cells = coord_cells[1], directions = nmat, include = TRUE)
  adjci <- purrr::map_dbl(adjc, 1, ~ seq(.x[1], .x[2]))
  mini_lyr[] <- 0
  mini_lyr[adjci] <- 1

  # left in just to help with a manual visual check
  raster::plot(mini_lyr)
  points(mini_coords)
  points(mini_coords[adj, ], col = "red")

  # check that all coords within the window are counted and all outside the window are not
  expect_true(all(terra::extract(mini_lyr, mini_coords[adj, ], ID = FALSE) == 1))
  expect_true(all(terra::extract(mini_lyr, mini_coords[-adj, ], ID = FALSE) != 1))
})


test_that("window_gd works for different spatial types", {
  load_mini_ex(quiet = TRUE)

  # sf coords
  sf_coords <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  capture_warnings(wpi_sf <- window_gd(mini_vcf, sf_coords, mini_lyr, rarify = FALSE))

  # matrix coords
  mat_coords <- as.matrix(mini_coords)
  capture_warnings(wpi_mat <- window_gd(mini_vcf, mat_coords, mini_lyr, rarify = FALSE))

  # df coords
  capture_warnings(wpi_df <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = FALSE))

  # spatvec coords
  vect_coords <- terra::vect(sf_coords)
  capture_warnings(wpi_vect <- window_gd(mini_vcf, vect_coords, mini_lyr, rarify = FALSE))

  # compare rasters
  # expect_true(terra::all.equal(wpi_df, wpi_mat))
  # expect_true(terra::all.equal(wpi_df, wpi_sf))
  # expect_true(terra::all.equal(wpi_df, wpi_vect))
  expect_equal(terra::values(wpi_df), terra::values(wpi_mat))
  expect_equal(terra::values(wpi_df), terra::values(wpi_sf))
  expect_equal(terra::values(wpi_df), terra::values(wpi_vect))
})

test_that("CRS are handled correctly", {
  load_mini_ex(quiet = TRUE)

  # NO CRS
  nocrs_lyr <- mini_lyr
  nocrs_coords <- mini_coords

  # CRS
  crs_lyr <- mini_lyr
  terra::crs(crs_lyr) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
  crs_coords <- sf::st_as_sf(mini_coords, coords = c("x", "y"))
  sf::st_crs(crs_coords) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"

  # no CRS lyr and coords
  cw <- capture_warnings(window_gd(mini_vcf, nocrs_coords, nocrs_lyr, rarify = FALSE))
  expect_equal(length(cw), 2)

  # no CRS lyr and CRS coords
  expect_warning(window_gd(mini_vcf, crs_coords, nocrs_lyr, rarify = FALSE))

  # CRS lyr and no CRS coords
  expect_warning(window_gd(mini_vcf, nocrs_coords, crs_lyr, rarify = FALSE))

  # CRS lyr and CRS coords (no warnings expected)
  wg <- window_gd(mini_vcf, crs_coords, crs_lyr, rarify = FALSE)
  expect_true(sf::st_crs(wg) == sf::st_crs(crs_lyr))

  # mismatched CRS
  crs2_lyr <- terra::rast(mini_lyr)
  terra::crs(crs2_lyr) <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0"
  expect_error(window_gd(mini_vcf, crs_coords, crs2_lyr, rarify = FALSE), "CRS of the provided coordinates and raster do not match")
})


test_that("edge cropping is performed correctly", {
  load_mini_ex(quiet = TRUE)
  wdim <- c(3, 5)
  capture_warnings(wg <- window_gd(mini_vcf, mini_coords, mini_lyr, wdim = wdim, rarify = FALSE, crop_edges = TRUE))

  expect_true((dim(mini_lyr)[1] - dim(wg)[1]) / 2 == (wdim[2] - 1) / 2)
  expect_true((dim(mini_lyr)[2] - dim(wg)[2]) / 2 == (wdim[1] - 1) / 2)
})

test_that("informative error for large min_n", {
  load_mini_ex(quiet = TRUE)
  capture_warnings(expect_error(wg <- window_gd(mini_vcf, mini_coords, mini_lyr, min_n = 100)))
})
