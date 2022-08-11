

#' Sliding window map of genetic diversity
#'
#' @param vcf object of type vcf ( (*note:* order matters! the coordinate and genetic data should be in the same order, there are currently no checks for this.))
#' @param stat genetic diversity stat to calculate (can either be "pi" for nucleotide diversity, "het" for average heterozygosity across all loci, "allelic.richness" for average allelic richness across all loci, or "biallelic.richness" to get average allelic richness across all loci for a biallelic dataset (this option faster than "allelic.richness"))
#' @param lyr RasterLayer to slide window across
#' @param wdim dimensions (height x width) of window, if only one value is provided a square window is created
#' @param fact aggregation factor to apply to the RasterLayer (*note:* increasing this value reduces computational time)
#' @param rarify if rarify = TRUE, rarefaction is performed
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction
#' @param min_n min number of samples to use in calculations (any focal cell with a window containing less than this number of samples will be assigned a value of NA; equal to rarify_n if rarify = TRUE, otherwise defaults to 2)
#' @param fun function to use to summarize data in window (defaults to base R mean)
#' @param parallel whether to parallelize the function (see vignette for setting up a cluster to do so)
#' @param L for calculating pi, L argument in \link[hierfstat]{pi.dosage} function. Return the average nucleotide diversity per nucleotide given the length L of the sequence. The wingen defaults is L = "nvariants" which sets L to the number of variants in the VCF. If L = NULL, returns the sum over SNPs of nucleotide diversity (note: L = NULL is the \link[hierfstat]{pi.dosage} default which wingen does not to use).
#' @param ncores if parallel = TRUE, number of cores to use for parallelization (defaults to total available number of cores minus 1)
#'
#' @return RasterStack that includes a raster of genetic diversity and a raster of the number of samples within the window for each cell
#' @export
#'
#' @examples
#' library("raster")
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE)
#' plot_gd(wpi, main = "Window pi")
#' plot_count(wpi)
#'
window_gd <- function(vcf, coords, lyr, stat = "pi", wdim = 5, fact = 0, rarify = FALSE, rarify_n = 4, rarify_nit = 5, min_n = 2, fun = mean, parallel = FALSE, L = "nvariants", ncores = NULL) {

  # check that the input file is a vcf or a path to a vcf object
  if (!inherits(vcf, "vcfR") & is.character(vcf)) {
    vcf <- vcfR::read.vcfR(vcf)
  } else if (!inherits(vcf, "vcfR") & !is.character(vcf)) {
    stop("Input is expected to be an object of class 'vcfR' or a path to a .vcf file")
  }

  # check to make sure coords and gen align
  check_data(vcf, coords)

  # calc stats
  if (stat == "allelic.richness") {
    # convert from vcf to genind
    gen <- vcf_to_genind(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = stat, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel, ncores)

    names(results[[1]]) <- "allelic_richness"
  }

  if (stat == "het" | stat == "heterozygosity") {
    # convert from vcf to heterozygosity matrix
    gen <- vcf_to_het(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = stat, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel, ncores)

    names(results[[1]]) <- "heterozygosity"
  }

  if (stat == "pi") {
    # convert from vcf to dosage matrix
    gen <- vcf_to_dosage(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = stat, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel, L, ncores)

    names(results[[1]]) <- "pi"
  }

  if (stat == "biallelic.richness") {
    # convert vcf to dosage matrix
    gen <- vcf_to_dosage(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = stat, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel, ncores)

    names(results[[1]]) <- "biallelic_richness"
  }

  names(results[[2]]) <- "sample_count"

  return(results)
}
#' Helper function for window_gd
#'
#' @param gen genetic data (*note:* order matters! the coordinate and genetic data should be in the same order, there are currently no checks for this.)
#' @inheritParams window_gd
#'
#' @return RasterStack that includes a raster of genetic diversity and a raster of the number of samples within the window for each cell
#' @export
#'
#' @keywords internal
#' @noRd
#'
window_gd_general <- function(gen, coords, lyr, stat = "pi", wdim = 3, fact = 0, rarify = FALSE, rarify_n = 2, rarify_nit = 10, min_n = 2, fun = mean, parallel = FALSE, L = "nvariants", ncores = NULL) {

  # set L if pi is being calculated
  if (stat == "pi" & !is.null(L) & !is.numeric(L)) {
    if (L == "nvariants") {
      L <- ncol(gen)
    }
  }

  # replace stat with function to calculate diversity statistic
  stat <- return_stat(stat)

  # format coords
  coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")

  # confirm that coords and gen align
  check_data(gen, coords)

  nmat <- wdim_to_mat(wdim)

  # make aggregated raster
  if (fact == 0) {
    lyr <- lyr * 0
  } else {
    lyr <- raster::aggregate(lyr, fact, fun = mean) * 0
  }

  # get cell index for each coordinate
  coord_cells <- raster::extract(lyr, coords, cell = TRUE)[, "cells"]

  if (parallel) {
    if (is.null(ncores)) ncores <- future::availableCores() - 1

    future::plan(future::multisession, workers = ncores)

    rast_vals <- furrr::future_map_dfr(1:raster::ncell(lyr),
      window_helper, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, L,
      .options = furrr::furrr_options(seed = TRUE, packages = c("raster", "purrr", "hierfstat", "stats", "adegenet"))
    )
  } else {
    rast_vals <- purrr::map_dfr(
      1:raster::ncell(lyr),
      window_helper, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, L
    )
  }


  # make copies of rasters
  alyr <- lyr
  nsagg <- lyr
  # assign values to rasters
  alyr[] <- rast_vals[, "gd"]
  nsagg[] <- rast_vals[, "ns"]

  results <- raster::stack(alyr, nsagg)

  return(results)
}

#' Helper function for window calculations
#'
#' @param i cell index
#' @param coord_cells cell indices for each coordinate
#' @param nmat neighborhood matrix
#'
#' @inheritParams window_gd_general
#'
#' @keywords internal
#' @noRd
#'
#' @return genetic diversity and counts for a single cell
#' @export
#'
window_helper <- function(i, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, L = NULL) {

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
    gd <- rarify_helper(gen, sub, rarify_n, rarify_nit, stat, fun, L)
  } else {
    gd <- sample_gd(gen, sub, stat, L)
  }

  # count the number of samples in the window
  ns <- length(sub)

  return(data.frame(gd = gd, ns = ns))
}



#' Rarefaction helper function
#'
#' @inheritParams window_gd_general
#'
#' @keywords internal
#' @noRd
#'
#' @return genetic diversity statistic for a rarified subsample
#' @export
#'
rarify_helper <- function(gen, sub, rarify_n, rarify_nit, stat, fun = mean, L = NULL) {
  # if number of samples is less than rarify_n, assign the value NA
  if (length(sub) < rarify_n) {
    gd <- NA
  }

  # if number of samples is greater than rarify_n, rarify
  if (length(sub) > rarify_n) {
    gd <- rarify_gd(gen, sub, rarify_nit = rarify_nit, rarify_n = rarify_n, stat = stat, fun = fun, L = L)
  }

  # if the number of samples is equal to rarify_n, calculate stat
  if (length(sub) == rarify_n) {
    gd <- sample_gd(gen, sub, stat, L)
  }

  return(gd)
}


#' Helper function to rarify subsample and calculate genetic diversity
#'
#' @inheritParams window_gd_general
#'
#' @return rarified genetic diversity statistic
#' @export
#'
#' @keywords internal
#' @noRd
#'
rarify_gd <- function(gen, sub, rarify_nit = 10, rarify_n = 4, stat, fun, L = NULL) {

  # check to make sure sub is greater than rarify_n
  if (!(length(sub) > rarify_n)) {
    stop("rarify_n is less than the number of samples provided")
  }

  # define subsample to rarify
  # (note: this combo step is done so when the number of unique combos < rarify_nit, extra calcs aren't performed)
  if (choose(length(sub), rarify_n) < rarify_nit) {
    # get all possible combos (transpose so rows are unique combos)
    cmb <- t(utils::combn(sub, rarify_n))
  } else {
    # random sample subsets of size rarify_nit (transpose so rows are unique combos)
    cmb <- t(replicate(rarify_nit, sample(sub, rarify_n), simplify = TRUE))
  }

  # for each of the possible combos get gendiv stat
  gdrar <- apply(cmb, 1, sample_gd, gen = gen, stat = stat, L = L)

  # summarize rarefaction results
  gd <- stats::na.omit(fun(gdrar))

  return(gd)
}


#' Helper function to calculate genetic diversity of a sample
#'
#' @inheritParams window_gd_general
#'
#' @return mean allelic richness of a subsample
#' @export
#'
#' @keywords internal
#' @noRd
#'
sample_gd <- function(gen, sub, stat, L = NULL) {
  if (is.null(L) | !identical(stat, calc_pi)) {
    gd <- stat(gen[sub, ])
  } else {
    gd <- stat(gen[sub, ], L)
  }
  return(gd)
}


#' Calculate mean allelic richness
#'
#' @param genind genind
#'
#' @return allelic richness averaged across all loci
#' @export
#'
#' @keywords internal
#' @noRd
#'
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
#' @export
#'
#' @keywords internal
#' @noRd
helper_calc_ar <- function(genind) {
  genind$pop <- rep(factor(1), nrow(genind$tab))
  # note [,1] references the first column which is AR for each locus across all inds (nrow(AR) == L)
  ar <- hierfstat::allelic.richness(genind)$Ar[, 1]
  return(ar)
}

#' Calculate mean heterozygosity
#'
#' @param hetmat matrix of heterozygosity (0/FALSE = homozygote, 1/TRUE = heterozygote)
#'
#' @return heterozygosity averaged across all individuals and then all loci
#' @export
#'
#' @keywords internal
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
#' @export
#'
#' @keywords internal
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
#' @export
#'
#' @keywords internal
#' @noRd
calc_mean_biar <- function(dos) {
  if (!all(dos %in% c(0, 1, 2, NA))) {
    stop("to calculate biallelic richness, all values in genetic matrix must be NA, 0, 1 or 2")
  }

  # if null dimensions (e.g., only one value provided), use helper_calc_biar directly
  # otherwise apply across columns (i.e., loci)
  if (is.null(dim(dos))){
    ar_by_locus <- helper_calc_biar(dos)
  } else {
    ar_by_locus <- apply(dos, 2, helper_calc_biar)
  }

  mean_ar <- mean(ar_by_locus, na.rm = TRUE)
  return(mean_ar)
}

#' Helper function to calculate allelic richness for a biallelic locus
#'
#' @param loc genotypes at a biallelic locus (must have values of 0, 1, or 2)
#'
#' @return biallelic richness value
#' @export
#'
#' @keywords internal
#' @noRd
helper_calc_biar <- function(loc) {
  uq <- unique(loc, na.rm = TRUE)
  if (1 %in% uq) {
    return(2)
  } else if (0 %in% uq & 2 %in% uq) {
    return(2)
  } else {
    return(1)
  }
}

#' Check coordinate and genetic data
#'
#' Check that the number of individuals in each data set align
#'
#' @param gen genetic data
#' @param coords coordinates
#'
#' @keywords internal
#' @noRd
#'
#' @export
#'
check_data <- function(gen, coords) {

  # check number of samples
  if (inherits(gen, "genind")) {
    nind <- nrow(gen$tab)
  }

  if (inherits(gen, "vcfR")) {
    nind <- (ncol(gen@gt) - 1)
  }

  if (inherits(gen, "data.frame") | inherits(gen, "matrix")) {
    nind <- nrow(gen)
  }


  # check to make sure coords and gen align
  if (nind != nrow(coords)) {
    stop("number of samples in coords data and number of samples in gen data are not equal")
  }

}

#' Helper function to get adjacent cells to a given cell index
#'
#' @param i cell index
#' @param r RasterLayer
#' @param n neighborhood matrix
#' @param coord_cells cell numbers of coordinates
#'
#' @return indices of coordinates that are adjacent to the given cell
#' @export
#'
#' @keywords internal
#' @noRd
get_adj <- function(i, r, n, coord_cells) {
  # get adjacent cells to cell i
  adjc <- raster::adjacent(r, i, directions = n, include = TRUE, sorted = TRUE)
  # get indices of adjacent cells
  adjci <- purrr::map_dbl(adjc, 1, function(x) {
    seq(x[1], x[2])
  })
  # get list of indices of coords in that set of cells
  sub <- which(coord_cells %in% adjci)

  return(sub)
}

#' Helper function to get genetic diversity functions
#'
#' @param x genetic diversity statistic
#'
#' @return function corresponding with desired statistic
#' @export
#'
#' @keywords internal
#' @noRd
return_stat <- function(x) {
  if (x == "pi") stat <- calc_pi
  if (x == "biallelic.richness") stat <- calc_mean_biar
  if (x == "allelic.richness") stat <- calc_mean_ar
  if (x == "het") stat <- calc_mean_het
  return(stat)
}
