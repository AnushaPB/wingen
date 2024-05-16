#' Run distance based moving window calculations
#'
#' Function used by both circle_gd and resist_gd
#'
#' @noRd
dist_gd <- function(gen, coords, lyr, stat = "pi",
                    fact = fact,
                    maxdist, distmat,
                    rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                    fun = mean, L = "nvariants", rarify_alleles = TRUE, sig = 0.05,
                    ...) {
  # convert maxdist to SpatRaster
  if (!is.numeric(maxdist)) {
    if (!inherits(maxdist, "SpatRaster")) maxdist <- terra::rast(maxdist)
    if (fact > 0) maxdist <- terra::aggregate(maxdist, fact, fun = mean)
  }

  # run moving window
  result <-
    purrr::map(
      stat,
      ~ dist_gd_stats(
        gen = gen,
        coords = coords,
        lyr = lyr,
        stat = .x,
        maxdist = maxdist,
        distmat = distmat,
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
    )


  # convert to raster stack
  r <- terra::rast(result)

  # remove any duplicated sample count layers
  r <- rm_duplicate_sample_count(r)

  return(r)
}


#' Helper function for mapping over stats
#'
#' @noRd
dist_gd_stats <- function(gen, coords, lyr, stat, maxdist, distmat,
                          rarify, rarify_n, rarify_nit, min_n,
                          fun, L, rarify_alleles, sig, ...) {
  # check that the input file is a vcfR or a path to a vcf object
  vcf <- vcf_check(gen)

  # check that coords and gen align and reformat data, if necessary
  # note: list2env adds the new, corrected gen and coords back to the environment
  list2env(check_data(vcf, coords = coords, distmat = distmat), envir = environment())

  # convert gen based on statistic being calculated
  x <- convert_vcf(vcf, stat)

  results <- dist_general(
    x = x,
    coords = coords,
    lyr = lyr,
    stat = stat,
    maxdist = maxdist,
    distmat = distmat,
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

  return(results)
}

#' General function for dist_gd
#'
#' @param x data to be summarized by the moving window (*note:* order matters! `coords` should be in the same order, there are currently no checks for this). The class of `x` required depends on the statistic being calculated (see the `stat` argument and the function description for more details)
#' @param stat moving window statistic to calculate (can either be `pi` for nucleotide diversity (`x` must be a dosage matrix), `Ho` for average observed heterozygosity across all loci (`x` must be a heterozygosity matrix) , "allelic_richness" for average allelic richness across all loci (`x` must be a `genind` type object), "biallelic_richness" to get average allelic richness across all loci for a biallelic dataset (`x` must be a dosage matrix). `stat` can also be set to any function that will take `x`as input and return a single numeric value (for example, `x` can be a vector and `stat` can be set equal to a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `stat`, additional arguments to pass to the `stat` function (e.g. if `stat = mean`, users may want to set `na.rm = TRUE`)
#' @inheritParams dist_gd
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#'
#' @noRd
dist_general <- function(x, coords, lyr, stat, maxdist, distmat,
                         rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                         fun = mean, L = NULL, rarify_alleles = TRUE, sig = 0.05, ...) {
  # check lyr and distmat
  lyr <- layer_coords_check(lyr, coords)
  if (terra::ncell(lyr) != ncol(distmat)) stop("Number of cells in raster layer and number of columns of distmat do not match")

  # run general moving window
  result <- run_general(
    x = x,
    lyr = lyr,
    coords = coords,
    coord_cells = NULL,
    nmat = NULL,
    distmat = distmat,
    maxdist = maxdist,
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

  return(result)
}

get_dist_index <- function(i, distmat, maxdist) {
  # get maxdist
  if (!is.numeric(maxdist)) md <- as.numeric(maxdist[i]) else md <- maxdist

  # modify dist matrix
  distmat[distmat > md] <- NA

  return(which(!is.na(distmat[i, ])))
}
