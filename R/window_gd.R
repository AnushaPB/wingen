#' Create a moving window map of genetic diversity
#'
#' Generate a continuous raster map of genetic diversity using moving windows.
#'
#' @param gen genetic data either as an object of type vcf or a path to a vcf file (*note:* order matters! The coordinate and genetic data should be in the same order; there are currently no checks for this)
#' @param coords coordinates of samples as sf points, a two-column matrix, or a data.frame representing x and y coordinates (see Details for important information about projections)
#' @param lyr SpatRaster or RasterLayer to slide the window across (see Details for important information about projections)
#' @param stat genetic diversity statistic(s) to calculate (see Details, defaults to `"pi"`). Can be a single statistic or a vector of statistics
#' @param wdim dimensions (height x width) of window; if only one value is provided, a square window is created (defaults to 3 x 3 window)
#' @param fact aggregation factor to apply to `lyr` (defaults to 0; *note:* increasing this value reduces computational time)
#' @param rarify if rarify = TRUE, rarefaction is performed (defaults to FALSE)
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction (defaults to min_n)
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction (defaults to 5). Can also be set to `"all"` to use all possible combinations of samples of size `rarify_n` within the window.
#' @param min_n minimum number of samples to use in calculations (any focal cell with a window containing less than this number of samples will be assigned a value of NA; defaults to 2)
#' @param fun function to use to summarize rarefaction results (defaults to mean, must take `na.rm = TRUE` as an argument)
#' @param L for calculating `"pi"`, L argument in \link[hierfstat]{pi.dosage} function. Return the average nucleotide diversity per nucleotide given the length L of the sequence. The wingen default is L = "nvariants" which sets L to the number of variants in the VCF. If L = NULL, returns the sum over SNPs of nucleotide diversity (*note:* L = NULL is the \link[hierfstat]{pi.dosage} default which wingen does not use)
#' @param rarify_alleles for calculating `"biallelic_richness"`, whether to perform rarefaction of allele counts as in \link[hierfstat]{allelic.richness} (defaults to TRUE)
#' @param sig for calculating `"hwe"`, significance threshold (i.e., alpha level) to use for hardy-weinberg equilibrium tests (defaults to 0.05)
#' @param crop_edges whether to remove cells on the edge of the raster where the window is incomplete (defaults to FALSE)
#' @param ... [deprecated] this was intended to be used to pass additional arguments to the `stat` function, however now formal arguments are used instead (see `L`, `rarify_alleles`, and `sig`). Passing additional arguments using `...` is still possible with the `*_general()` functions.
#' @details
#'
#' Coordinates and rasters should be in a projected (planar) coordinate system such that raster cells are of equal sizes.
#' Therefore, spherical systems (including latitute-longitude coordinate systems) should be projected prior to use.
#' Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' Current genetic diversity metrics that can be specified with `stat` include:
#' - `"pi"` for nucleotide diversity (default) calculated using `hierfstat` \link[hierfstat]{pi.dosage}
#' - `"Ho"` for average observed heterozygosity across all sites
#' - `"allelic_richness"` for average number of alleles across all sites
#' - `"biallelic_richness"` for average allelic richness across all sites for a biallelic dataset (this option is faster than `"allelic_richness"`)
#' - `"hwe"` for the proportion of sites that are not in Hardy–Weinberg equilibrium, calculated using `pegas` \link[pegas]{hw.test} at the 0.05 level (other alpha levels  can be specified by adding the sig argument; e.g., `sig = 0.10`).
#' - `"basic_stats"` for a series of statistics produced by `hierfstat` \link[hierfstat]{basic.stats} including
#' mean observed heterozygosity (same as Ho), mean gene diversities within population (Hs),
#' Gene diversities overall (Ht), and Fis following Nei (1987). Population-based statistics (e.g., FST) normally reported by \link[hierfstat]{basic.stats}
#' are not included as they are not meaningful within the individual-based moving windows.
#'
#' @return SpatRaster that includes raster layers of genetic diversity and a raster layer of the number of samples within the window for each cell
#' @export
#'
#' @examples
#'
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE)
#'
window_gd <- function(gen, coords, lyr, stat = "pi", wdim = 3, fact = 0,
                      rarify = FALSE, rarify_n = NULL, rarify_nit = 5, min_n = 2,
                      fun = mean, L = "nvariants", rarify_alleles = TRUE, sig = 0.05,
                      crop_edges = FALSE, ...) {
  # run moving window
  result <-
    purrr::map(
      stat,
      \(stat)
      window_gd_stats(
        gen = gen,
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
        sig = sig,
        crop_edges = crop_edges
      )
    )

  # convert to raster stack
  r <- terra::rast(result)

  # remove any duplicated sample count layers
  r <- rm_duplicate_sample_count(r)

  return(r)
}

#' Helper function for mapping over stats
#' @noRd
window_gd_stats <- function(gen, coords, lyr, stat, wdim, fact,
                            rarify, rarify_n, rarify_nit, min_n,
                            fun, L, rarify_alleles, sig,
                            crop_edges, ...) {
  # check that the input file is a vcf or a path to a vcf object
  vcf <- vcf_check(gen)

  # check that coords and vcf align and reformat data, if necessary
  # note: list2env adds the new, corrected vcf and coords back to the environment
  list2env(check_data(vcf, coords), envir = environment())

  # convert vcf based on statistic being calculated
  x <- convert_vcf(vcf, stat)

  results <- window_general(
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
    sig = sig,
    crop_edges = crop_edges,
    ...
  )

  return(results)
}

#' General function for making moving window maps
#'
#' Generate a continuous raster map using moving windows. While \link[wingen]{window_gd} is
#' built specifically for making moving window maps of genetic diversity from vcfs,
#' `window_general` can be used to make moving window maps from different data inputs.
#' See details for how to format data inputs for different statistics.
#'
#' @param x data to be summarized by the moving window (*note:* order matters! `coords` should be in the same order, there are currently no checks for this). The class of `x` required depends on the statistic being calculated (see the `stat` argument and the function description for more details)
#' @param stat moving window statistic to calculate (can either be `"pi"` for nucleotide diversity (`x` must be a dosage matrix), `"Ho"` for average observed heterozygosity across all loci (`x` must be a heterozygosity matrix) , `"allelic_richness"` for average allelic richness across all loci (`x` must be a `genind` type object), `"biallelic_richness"` to get average allelic richness across all loci for a biallelic dataset (`x` must be a dosage matrix). `stat` can also be set to any function that will take `x`as input and return a single numeric value (for example, `x` can be a vector and `stat` can be set equal to a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `stat`, additional arguments to pass to the `stat` function (e.g. if `stat = mean`, users may want to set `na.rm = TRUE`)
#' @inheritParams window_gd
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#'
#' @details
#'
#' To calculate genetic diversity statistics with the built in wingen functions, data must be formatted as such:
#' - for `"pi"` or  `"biallelic_richness"`, `x` must be a dosage matrix with values of 0, 1, or 2
#' - for `"Ho"`, `x` must be a heterozygosity matrix where values of 0 = homozygosity and values of 1 = heterozygosity
#' - for `"allelic_richness"` or `"hwe`, `x` must be a `genind` type object
#' - for `"basic_stats"`, `x` must be a `hierfstat` type object
#'
#' Otherwise, `stat` can be any function that takes a matrix or data frame and outputs a
#' single numeric value (e.g., a function that produces a custom diversity index);
#' however, this should be attempted with caution since this functionality has
#' not have been tested extensively and may produce errors.
#'
#' @export
window_general <- function(x, coords, lyr, stat, wdim = 3, fact = 0,
                           rarify = FALSE, rarify_n = NULL, rarify_nit = 5, min_n = 2,
                           fun = mean, L = "nvariants", rarify_alleles = TRUE, sig = 0.05,
                           crop_edges = FALSE, ...) {
  # check and aggregate layer and coords (only lyr is returned)
  lyr <- layer_coords_check(lyr, coords, fact)

  # check wdim
  wdim <- wdim_check(wdim)

  # make neighbor matrix
  nmat <- wdim_to_mat(wdim)

  # get cell index for each coordinate
  coord_cells <- terra::extract(lyr, coords, cell = TRUE)[, "cell"]

  # run general moving window
  result <- run_general(
    x = x,
    lyr = lyr,
    coords = coords,
    coord_cells = coord_cells,
    nmat = nmat,
    distmat = NULL,
    maxdist = NULL,
    stat = stat,
    rarify = rarify,
    rarify_n = rarify_n,
    rarify_nit = rarify_nit,
    min_n = min_n,
    fun = fun,
    L = L,
    rarify_alleles = rarify_alleles,
    sig = sig,
    ...
  )

  # crop resulting raster
  if (crop_edges) result <- edge_crop(result, wdim)

  return(result)
}
