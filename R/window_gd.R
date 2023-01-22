

#' Create a moving window map of genetic diversity
#'
#' Generate a continuous raster map of genetic diversity using moving windows
#'
#' @param vcf object of type vcf or a path to a vcf file (*note:* order matters! The coordinate and genetic data should be in the same order; there are currently no checks for this)
#' @param coords coordinates of samples as sf points, a two-column matrix, or a data.frame representing x and y coordinates. Should be in a Euclidean system (i.e., not longitude latitude) or the window cell height and width will not be equal (see details).
#' @param lyr SpatRaster or RasterLayer to slide the window across. Should be in a Euclidean system (i.e., not longitude latitude) or the window cell height and width will not be equal (see details).
#' @param stat genetic diversity statistic to calculate (can either be `"pi"` for nucleotide diversity (default), `"Ho"` for average observed heterozygosity across all sites, `"allelic_richness"` for average number of alleles across all sites, or `"biallelic_richness"` to get average allelic richness across all sites for a biallelic dataset (this option faster than `"allelic_richness"`))
#' @param wdim dimensions (height x width) of window; if only one value is provided, a square window is created (defaults to 3 x 3 window)
#' @param fact aggregation factor to apply to `lyr` (defaults to 0; *note:* increasing this value reduces computational time)
#' @param rarify if rarify = TRUE, rarefaction is performed (defaults to FALSE)
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction (defaults to 2)
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction (defaults to 5). Can also be set to `"all"` to use all possible combinations of samples of size `rarify_n` within the window.
#' @param min_n min number of samples to use in calculations (any focal cell with a window containing less than this number of samples will be assigned a value of NA; equal to rarify_n if rarify = TRUE, otherwise defaults to 2)
#' @param fun function to use to summarize rarefaction results (defaults to mean, must take `na.rm = TRUE` as an argument)
#' @param L for calculating pi, L argument in \link[hierfstat]{pi.dosage} function. Return the average nucleotide diversity per nucleotide given the length L of the sequence. The wingen defaults is L = "nvariants" which sets L to the number of variants in the VCF. If L = NULL, returns the sum over SNPs of nucleotide diversity (*note:* L = NULL is the \link[hierfstat]{pi.dosage} default which wingen does not use)
#' @param rarify_alleles for calculating biallelic_richness, whether to perform rarefaction of allele counts as in \link[hierfstat]{allelic.richness} (defaults to TRUE)
#' @param parallel whether to parallelize the function (defaults to FALSE)
#' @param ncores if parallel = TRUE, number of cores to use for parallelization (defaults to total available number of cores minus 1)
#' @param crop_edges whether to remove cells on the edge of the raster where the window is incomplete (defaults to FALSE)
#' @details Coordinates and rasters should be in a Euclidean coordinate system (i.e., UTM coordinates) such that raster cell width and height are equal distances.
#' As such, longitude-latitude systems should be transformed before using window_gd. Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#' @export
#'
#' @examples
#'
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE)
#' plot_gd(wpi, main = "Window pi")
#' plot_count(wpi)
#'
window_gd <- function(vcf, coords, lyr, stat = "pi", wdim = 3, fact = 0,
                      rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                      fun = mean, L = "nvariants", rarify_alleles = TRUE,
                      parallel = FALSE, ncores = NULL, crop_edges = FALSE) {
  # check that the input file is a vcf or a path to a vcf object
  vcf <- vcf_check(vcf)

  # check that coords and vcf align and reformat data, if necessary
  # note: list2env adds the new, corrected vcf and coords back to the environment
  list2env(check_data(vcf, coords), envir = environment())

  # convert vcf based on statistic being calculated
  x <- convert_vcf(vcf, stat)

  # run moving window
  result <- window_general(
    x = x,
    coords = coords,
    lyr = lyr,
    stat = stat,
    wdim = wdim,
    fact = fact,
    rarify = rarify,
    rarify_n = rarify_n,
    rarify_nit = rarify_nit,
    min_n = min_n,
    fun = fun,
    L = L,
    rarify_alleles = rarify_alleles,
    parallel = parallel,
    ncores = ncores,
    crop_edges = crop_edges
  )

  return(result)
}

#' General function for making moving window maps
#'
#' Generate a continuous raster map using moving windows. While \link[wingen]{window_gd} is built specifically for making moving window maps of genetic diversity from vcfs,
#' `window_general` can be used to make moving window maps from different data inputs. Unlike `window_gd`, `window_general` will not convert your data into
#' the correct format for calculations of different diversity metrics. To calculate `pi` or `biallelic_richness`, `x` must be a dosage matrix with values of 0, 1, or 2. To calculate
#' `Ho`, `x` must be a heterozygosity matrix where values of 0 = homozygosity and values of 1 = heterozygosity. To calculate `allelic_richness`, `x` must be a `genind` type object.
#' Users can set `x` to a vector and create moving window maps with any function that can be applied to a vector (e.g., `stat = mean`, `var`, `sum`, etc.).
#' `x` can also be a matrix or data frame (where rows are individuals), and then `stat` can be any function that takes a matrix or data frame and outputs a
#' single numeric value (e.g., a function that produces a custom diversity index); however, this should be attempted with caution since this functionality has
#'  not have been tested extensively and may produce errors.
#'
#' @param x data to be summarized by the moving window (*note:* order matters! `coords` should be in the same order, there are currently no checks for this). The class of `x` required depends on the statistic being calculated (see the `stat` argument and the function description for more details)
#' @param stat moving window statistic to calculate (can either be `pi` for nucleotide diversity (`x` must be a dosage matrix), `Ho` for average observed heterozygosity across all loci (`x` must be a heterozygosity matrix) , "allelic_richness" for average allelic richness across all loci (`x` must be a `genind` type object), "biallelic_richness" to get average allelic richness across all loci for a biallelic dataset (`x` must be a dosage matrix). `stat` can also be set to any function that will take `x`as input and return a single numeric value (for example, `x` can be a vector and `stat` can be set equal to a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `stat`, additional arguments to pass to the `stat` function (e.g. if `stat = mean`, users may want to set `na.rm = TRUE`)
#' @inheritParams window_gd
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#'
#' @export
window_general <- function(x, coords, lyr, stat, wdim = 3, fact = 0,
                           rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                           fun = mean, L = "nvariants", rarify_alleles = TRUE,
                           parallel = FALSE, ncores = NULL, crop_edges = FALSE, ...) {
  # check layers and coords (only lyr is modified and returned)
  lyr <- layer_coords_check(lyr, coords)

  # check wdim
  wdim <- wdim_check(wdim)

  # set L if pi is being calculated
  if (!is.null(L) & !is.numeric(L)) if (L == "nvariants") L <- ncol(x)

  # Get function to calculate the desired statistic
  stat_function <- return_stat(stat, ...)

  # check that coords and x align and reformat data, if necessary
  # note: list2env adds the new, corrected x and coords back to the environment
  list2env(check_data(x, coords), envir = environment())

  # make neighbor matrix
  nmat <- wdim_to_mat(wdim)

  # make aggregated raster
  if (fact == 0) lyr <- lyr * 0 else lyr <- terra::aggregate(lyr, fact, fun = mean) * 0

  # get cell index for each coordinate
  coord_cells <- terra::extract(lyr, coords, cell = TRUE)[, "cell"]

  # run sliding window calculations
  if (parallel) {
    # currently, terra uses a C++ pointer which means SpatRasters cannot be directly passed to nodes on a computer cluster
    # instead of saving the raster layer to a file, I am converting it to a RasterLayer temporarily (it will get switched back)
    lyr <- raster::raster(lyr)

    if (is.null(ncores)) ncores <- future::availableCores() - 1

    future::plan(future::multisession, workers = ncores)

    rast_vals <- furrr::future_map_dfr(1:terra::ncell(lyr), window_helper,
      lyr = lyr, x = x, coord_cells = coord_cells, nmat = nmat,
      stat_function = stat_function, rarify = rarify, rarify_n = rarify_n, rarify_nit = rarify_nit,
      min_n = min_n, fun = fun, L = L, rarify_alleles = rarify_alleles,
      .options = furrr::furrr_options(seed = TRUE, packages = c("wingen", "terra", "raster", "adegenet"))
    )

    # convert back to SpatRast
    lyr <- terra::rast(lyr)
  } else {
    rast_vals <- purrr::map_dfr(1:terra::ncell(lyr), window_helper,
      lyr = lyr, x = x, coord_cells = coord_cells, nmat = nmat,
      stat_function = stat_function, rarify = rarify, rarify_n = rarify_n, rarify_nit = rarify_nit,
      min_n = min_n, fun = fun, L = L, rarify_alleles = rarify_alleles
    )
  }

  # format resulting raster values
  result <- vals_to_lyr(lyr, rast_vals, stat)

  # crop resulting raster
  if (crop_edges) result <- edge_crop(result, wdim)

  return(result)
}

#' Helper function for window calculations
#'
#' @param i cell index
#' @param coord_cells cell indices for each coordinate
#' @param nmat neighborhood matrix
#'
#' @inheritParams window_general
#'
#' @return genetic diversity and counts for a single cell
#'
#' @noRd
window_helper <- function(i, lyr, x, coord_cells, nmat, stat_function,
                          rarify, rarify_n, rarify_nit, min_n,
                          fun, L = NULL, rarify_alleles = TRUE) {
  # convert RasterLayer back to SpatRaster (necessary for parallelized task)
  # if (inherits(lyr, "RasterLayer")) lyr <- terra::rast(lyr)

  # if rarify = TRUE, min_n = rarify_n (i.e. minimum defaults to rarify_n)
  if (rarify) min_n <- rarify_n

  # skip if raster value is NA
  if (is.na(lyr[i])) {
    return(data.frame(gd = NA, ns = NA))
  }

  # get sample indices in window
  sub <- get_adj(i, lyr, nmat, coord_cells)

  # if there are too few samples in that window assign the cell value NA
  if (length(sub) < min_n) {
    gd <- NA
  } else if (rarify) {
    gd <- rarify_helper(x, sub, rarify_n, rarify_nit, stat_function, fun, L = L, rarify_alleles = rarify_alleles)
  } else {
    gd <- sample_gd(x, sub, stat_function, L = L, rarify_alleles = rarify_alleles)
  }

  # count the number of samples in the window
  ns <- length(sub)

  return(data.frame(gd = gd, ns = ns))
}



#' Rarefaction helper function
#'
#' @inheritParams window_general
#'
#' @noRd
#'
#' @return genetic diversity statistic for a rarified subsample
#'
#' @noRd
rarify_helper <- function(x, sub, rarify_n, rarify_nit, stat_function,
                          fun = mean, L = NULL, rarify_alleles = TRUE) {
  # if number of samples is less than rarify_n, assign the value NA
  if (length(sub) < rarify_n) {
    gd <- NA
  }

  # if number of samples is greater than rarify_n, rarify
  if (length(sub) > rarify_n) {
    gd <- rarify_gd(x, sub, rarify_nit = rarify_nit, rarify_n = rarify_n, stat_function = stat_function, fun = fun, L = L, rarify_alleles = rarify_alleles)
  }

  # if the number of samples is equal to rarify_n, calculate stat
  if (length(sub) == rarify_n) {
    gd <- sample_gd(x, sub, stat_function, L = L, rarify_alleles = rarify_alleles)
  }

  return(gd)
}


#' Helper function to rarify subsample and calculate genetic diversity
#'
#' @inheritParams window_general
#'
#' @return rarified genetic diversity statistic
#'
#' @noRd
rarify_gd <- function(x, sub, rarify_nit = 5, rarify_n = 4, stat_function,
                      fun, L = NULL, rarify_alleles = TRUE) {
  # check to make sure sub is greater than rarify_n
  if (!(length(sub) > rarify_n)) {
    stop("rarify_n is less than the number of samples provided")
  }

  # define subsample to rarify
  if (rarify_nit == "all") {
    # get all possible combos (transpose so rows are unique combos)
    cmb <- t(utils::combn(sub, rarify_n))
  } else if (choose(length(sub), rarify_n) < rarify_nit) {
    # (note: this combo step is done so when the number of unique combos < rarify_nit, extra calcs aren't performed)
    # get all possible combos (transpose so rows are unique combos)
    cmb <- t(utils::combn(sub, rarify_n))
  } else {
    # randomly sample subsets of size rarify_nit (transpose so rows are unique combos)
    # note: replace is set to FALSE so the same individual cannot be drawn multiple times within the same sample
    # however, individuals can be drawn multiple times across different samples
    cmb <- t(replicate(rarify_nit, sample(sub, rarify_n, replace = FALSE), simplify = TRUE))
  }

  # for each of the possible combos get gendiv stat
  gdrar <- apply(cmb, 1, sample_gd, x = x, stat_function = stat_function, L = L, rarify_alleles = rarify_alleles)

  # summarize rarefaction results
  gd <- fun(gdrar, na.rm = TRUE)

  return(gd)
}


#' Helper function to calculate genetic diversity of a sample
#'
#' @inheritParams window_general
#'
#' @return mean allelic richness of a subsample
#'
#' @noRd
sample_gd <- function(x, sub, stat_function, L = NULL, rarify_alleles = TRUE) {
  if (identical(stat_function, calc_mean_biar)) {
    gd <- stat_function(x[sub, ], rarify_alleles)
  } else if (is.null(L) | !identical(stat_function, calc_pi)) {
    gd <- stat_function(x[sub, ])
  } else {
    gd <- stat_function(x[sub, ], L)
  }
  return(gd)
}


#' Calculate mean allelic richness
#'
#' @param genind genind
#'
#' @return allelic richness averaged across all loci
#'
#' @noRd
calc_mean_ar <- function(genind) {
  ar <- helper_calc_ar(genind)
  gd <- mean(ar, na.rm = TRUE)
  return(gd)
}

#' Helper function to calculate allelic richness
#'
#' @param genind genind object
#'
#' @return allelic richness
#'
#' @noRd
helper_calc_ar <- function(genind) {
  # get number of individuals
  nind <- nrow(genind@tab)

  # assign pops so that the whole sample is treated as one pop
  genind$pop <- rep(factor(1), nind)

  # note: min.n is the The number of alleles down to which the number of alleles should be rarefied.
  # The default is the minimum number of individuals genotyped (times 2 for diploids). However, if there
  # are NA values then it doesn't count those as genotypes. Therefore, to ensure that rarefaciton DOES NOT
  # OCCUR (since we have our own rarefaction step) min.n is set to the total number of individuals
  # (including those with NAs) times two (assuming diploids)

  # note: [,1] references the first column which is AR for each site across all inds (nrow(AR) == L)
  # ar <- hierfstat::allelic.richness(genind, min.n = nind * 2)$Ar[, 1]
  ar <- hierfstat::allelic.richness(genind)$Ar[, 1]
  return(ar)
}


#' Calculate mean heterozygosity
#'
#' @param hetmat matrix of heterozygosity (0/FALSE = homozygote, 1/TRUE = heterozygote)
#'
#' @return heterozygosity averaged across all individuals and then all loci
#'
#' @noRd
calc_mean_het <- function(hetmat) {
  gd <- mean(hetmat, na.rm = TRUE)
  return(gd)
}

#' Calculate nucleotide diversity (pi) from dosage data
#'
#' Wrapper for \link[hierfstat]{pi.dosage} function
#'
#' @param dos a ni X nl dosage matrix containing the number of derived/alternate alleles each individual carries at each SNP
#' @param L length of the sequence (*note:* defaults to number of loci in the provided dosage matrix; TODO: COME BACK AND FIX THIS)
#'
#' @return nucleotide diversity (pi)
#'
#' @noRd
calc_pi <- function(dos, L = NULL) {
  gd <- hierfstat::pi.dosage(dos, L = L)
  return(gd)
}

#' Calculate mean allelic richness for biallelic data
#'
#' @param dos dosage matrix
#'
#' @return allelic richness averaged across all loci
#'
#' @noRd
calc_mean_biar <- function(dos, rarify_alleles = TRUE) {
  if (!all(dos %in% c(0, 1, 2, NA))) {
    stop("to calculate biallelic richness, all values in genetic matrix must be NA, 0, 1 or 2")
  }

  # if rarify_alleles = TRUE, get min.n for rarefaction
  # min.n is set to the number of non-NA genotypes * 2 (for diploid)
  if (rarify_alleles) {
    min.n <- get_minn(dos)
  } else {
    min.n <- NULL
  }

  # if null dimensions (e.g., only one value provided), use helper_calc_biar directly
  # otherwise apply across columns (i.e., loci)
  if (is.null(dim(dos))) {
    ar_by_locus <- sapply(dos, helper_calc_biar, rarify_alleles, min.n)
  } else {
    ar_by_locus <- apply(dos, 2, helper_calc_biar, rarify_alleles, min.n)
  }

  mean_ar <- mean(ar_by_locus, na.rm = TRUE)
  return(mean_ar)
}

#' Helper function to calculate allelic richness for a biallelic locus
#'
#' @param loc genotypes at a biallelic locus (must have values of 0, 1, or 2)
#'
#' @return biallelic richness value
#'
#' @noRd
helper_calc_biar <- function(loc, rarify_alleles = TRUE, min.n = NULL) {
  # check if all NA and if so return NA
  if (all(is.na(loc))) {
    return(NA)
  }

  if (rarify_alleles) {
    # omit NAs before counting alleles
    loc_NArm <- stats::na.omit(loc)

    # make df of counts of reference and alternate
    counts <- c(R = sum(loc_NArm), A = 2 * length(loc_NArm) - sum(loc_NArm, na.rm = TRUE))

    # rarefied counts calculation taken from hierfstat::allelic.richness function
    AR <- raref(counts, min.n = min.n)
  } else {
    # calculate number of unique alleles
    # note: has to be na.omit (na.rm is not an argument for unique)
    uq <- unique(stats::na.omit(loc))
    if (1 %in% uq) {
      AR <- 2
    } else if (0 %in% uq & 2 %in% uq) {
      AR <- 2
    } else {
      AR <- 1
    }
  }

  return(AR)
}

#' Rarify allele counts (based on \link[hierfstat]{allelic.richness}code)
#'
#' @param x allele counts
#' @param min.n the number of alleles down to which the number of alleles should be rarefied.
#'
#' @noRd
raref <- function(x, min.n) {
  nn <- sum(x)
  dum <- exp(lchoose(nn - x, min.n) - lchoose(nn, min.n))
  dum[is.na(dum)] <- 0
  return(sum(1 - dum))
}

#' Helper function to get adjacent cells to a given cell index
#'
#' @param i cell index
#' @param r SpatRast
#' @param n neighborhood matrix
#' @param coord_cells cell numbers of coordinates
#'
#' @return indices of coordinates that are adjacent to the given cell
#'
#' @noRd
get_adj <- function(i, r, n, coord_cells) {
  # get adjacent cells to cell i
  adjc <- terra::adjacent(r, i, directions = n, include = TRUE)
  # get indices of adjacent cells
  adjci <- purrr::map_dbl(adjc, 1, ~ seq(.x[1], .x[2]))
  # get list of indices of coords in that set of cells
  sub <- which(coord_cells %in% adjci)

  return(sub)
}


#' Calculate min.n
#'
#' @param dos dosage matrix
#'
#' @noRd
get_minn <- function(dos) {
  if (is.null(nrow(dos))) {
    # if nrow is NULL then only one sample is included so min.n must be 2
    min.n <- 2
  } else {
    min.n <- 2 * min(apply(dos, 2, countgen), na.rm = TRUE)
  }
  return(min.n)
}

#' Count not NA genotypes and return NA if all are NA
#'
#' @param x dosage matrix
#'
#' @noRd
countgen <- function(x) {
  notNA <- sum(!is.na(x))
  if (notNA == 0) notNA <- NA
  return(notNA)
}

#' Check coordinate and genetic data
#'
#' Check that the number of individuals in each data set align
#'
#' @param x moving window data
#' @param coords coordinates
#'
#' @noRd
#'
check_data <- function(x, coords = NULL) {
  # if x is a vector, convert to a dataframe
  if (is.vector(x)) x <- data.frame(x)

  # check and format coords
  if (!is.null(coords)) {
    if (nrow(coords) == 1) stop("cannot run window_gd with only one individual")
  }

  # check number of samples
  if (inherits(x, "genind")) nind <- nrow(x$tab)

  if (inherits(x, "vcfR")) nind <- (ncol(x@gt) - 1)

  if (inherits(x, "data.frame") | inherits(x, "matrix")) nind <- nrow(x)

  # check coords
  if (!is.null(coords)) {
    if (nind != nrow(coords)) {
      stop("number of samples in coords data and number of samples in gen data are not equal")
    }
  }

  # check for rows or columns with missing data in a vcf and give warning if there are invariant sites
  if (inherits(x, "vcfR")) {
    return(check_vcf_NA(x, coords))
  } else {
    return(list(x = x, coords = coords))
  }
}

#' Check vcf for loci and individuals with all NAs and return corrected vcf and coords
#'
#' @param vcf vcfR
#' @param coords coordinates
#'
#' @noRd
check_vcf_NA <- function(vcf, coords = NULL) {
  # check for mismatch before indexing
  if (!is.null(coords)) {
    if ((ncol(vcf@gt) - 1) != nrow(coords)) {
      stop("number of samples in coords data and number of samples in vcf are not equal")
    }
  }

  if (nrow(vcf@fix) == 1) {
    NA_col <- get_allNA(vcf@gt[-1])
    NA_row <- FALSE
  } else {
    NA_col <- get_allNA(vcf@gt[, -1], MARGIN = 2)
    NA_row <- get_allNA(vcf@gt[, -1], MARGIN = 1)
  }

  if (any(NA_row)) {
    warning("Markers with no scored alleles have been removed")
    vcf <- vcf[!NA_row, ]
  }

  if (any(NA_col)) {
    warning("Individuals with no scored loci have been removed")
    vcf <- vcf[, c(TRUE, !NA_col)]
  }

  # check for invariant sites
  if (any(!vcfR::is.polymorphic(vcf, na.omit = TRUE))) warning("invariant sites found in vcf")

  # make results
  if (is.null(coords)) {
    result <- vcf
  } else {
    coords <- coords[!NA_col, ]
    result <- list(vcf = vcf, coords = coords)
  }

  return(result)
}

#' Helper function to get NA values
#'
#' @inheritParams base::array
#'
#' @noRd
get_allNA <- function(x, MARGIN = NULL) {
  if (is.null(dim(x))) allNA <- is.na(x)
  if (!is.null(dim(x))) {
    allNA <- apply(x, MARGIN, function(x) {
      all(is.na(x))
    })
  }
  return(allNA)
}

#' Convert vcf to correct format based on stat
#'
#' @param vcf vcfR
#' @param stat genetic diversity statistic
#'
#' @noRd
convert_vcf <- function(vcf, stat) {
  if (stat == "allelic_richness") {
    return(vcf_to_genind(vcf))
  }

  if (stat == "Ho") {
    return(vcf_to_het(vcf))
  }

  if (stat == "pi" | stat == "biallelic_richness") {
    return(vcf_to_dosage(vcf))
  }

  stop(paste0(stat, " is an invalid arugment for stat"))
}

#' Rename results from window_gd
#'
#' @param x SpatRaster produced by window_gd
#' @param stat genetic diversity statistic
#'
#' @noRd
name_results <- function(x, stat) {
  names(x[[2]]) <- "sample_count"

  if (is.character(stat)) names(x[[1]]) <- stat else names(x[[1]]) <- "custom"

  return(x)
}

#' Helper function to get genetic diversity functions
#'
#' @param stat moving window statistic to calculate (can either be `pi` for nucleotide diversity, `Ho` for average observed heterozygosity across all loci, "allelic_richness" for average number of alleles across all loci, "biallelic_richness" to get average number of alleles across all loci for a biallelic dataset. `stat` can also be set to functions that will return a single numeric value from the input data (for example a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `x`, additional arguments to pass to the `x` function (e.g. if `x = mean`, users may want to set `na.rm = TRUE`)
#'
#' @return function corresponding with desired statistic
#'
#' @noRd
return_stat <- function(stat, ...) {
  if (inherits(stat, "function")) {
    return(purrr::partial(stat, ...))
  }

  if (stat == "pi") {
    return(calc_pi)
  }

  if (stat == "biallelic_richness") {
    return(calc_mean_biar)
  }

  if (stat == "allelic_richness") {
    return(calc_mean_ar)
  }

  if (stat == "Ho") {
    return(calc_mean_het)
  }

  stop(paste(stat, "is an invalid argument for stat"))
}

#' Helper function to check lyr and coords
#'
#' @param lyr RasterLayer or SpatRaster
#' @param coords sf points, data frame, or matrix representing coordinates
#'
#' @return SpatRaster
#'
#' @noRd
layer_coords_check <- function(lyr, coords) {
  # check coords and lyr
  crs_check_window(lyr, coords)

  # convert to terra
  if (inherits(lyr, "RasterLayer") | inherits(lyr, "RasterStack")) lyr <- terra::rast(lyr)

  # check number of layers
  nlayers <- terra::nlyr(lyr)
  if (nlayers > 1) {
    warning(paste0(nlayers, " provided, but only one is need. Defaults to using the first layer."))
    lyr <- lyr[[1]]
  }

  return(lyr)
}

#' Convert values into new raster layers
#'
#' @param lyr SpatRaster
#' @param rast_vals dataframe of gd and ns
#'
#' @return SpatRaster
#'
#' @noRd
vals_to_lyr <- function(lyr, rast_vals, stat) {
  # make copies of rasters
  alyr <- lyr
  nsagg <- lyr

  # assign values to rasters
  alyr[] <- rast_vals[, "gd"]
  nsagg[] <- rast_vals[, "ns"]

  # stack rasters
  results <- c(alyr, nsagg)

  # set raster layer names based on stat
  results <- name_results(results, stat)

  return(results)
}


#' Crop edge off raster
#'
#' @param x SpatRaster
#' @param wdim window dimensions
#'
#' @return SpatRaster
#'
#' @noRd
edge_crop <- function(x, wdim) {
  if (length(wdim) == 1) wdim <- c(wdim, wdim)

  # get extent
  x_ext <- terra::ext(x)

  # calculate x edge buffer
  x_edge_size <- res(x)[1] * ((wdim[1] - 1) / 2)
  xmin <- x_ext$xmin + x_edge_size
  xmax <- x_ext$xmax - x_edge_size

  # calculate y edge buffer
  y_edge_size <- res(x)[2] * ((wdim[2] - 1) / 2)
  ymin <- x_ext$ymin + y_edge_size
  ymax <- x_ext$ymax - y_edge_size

  # crop raster
  x_crop <- terra::crop(x, terra::ext(xmin, xmax, ymin, ymax))

  return(x_crop)
}
