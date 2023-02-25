
#' Run distance based moving window calculations
#'
#' Function used by both circle_gd and resist_gd
#'
#' @noRd
dist_gd <- function(gen, coords, lyr, stat = "pi", maxdist, distmat,
                    rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                    fun = mean, L = "nvariants", rarify_alleles = TRUE,
                    parallel = FALSE, ncores = NULL, crop_edges = FALSE) {

  # run moving window
  result <-
    purrr::map(stat,
               ~dist_gd_stats(
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
                 parallel = parallel,
                 ncores = ncores,
                 crop_edges = crop_edges
               ))


  # convert to raster stack
  r <- terra::rast(result)

  # remove any duplicated sample count layers
  r <- rm_duplicate_sample_count(r)

  return(result)
}


#' Helper function for mapping over stats
#'
#' @noRd
dist_gd_stats <- function(gen, coords, lyr, stat = "pi", maxdist, distmat,
                          rarify = FALSE, rarify_n = NULL, rarify_nit = 5, min_n = 2,
                          fun = mean, L = "nvariants", rarify_alleles = TRUE,
                          parallel = FALSE, ncores = NULL, crop_edges = FALSE){

  # check that the input file is a vcfR or a path to a vcf object
  vcf <- vcf_check(gen)

  # check that coords and gen align and reformat data, if necessary
  # note: list2env adds the new, corrected gen and coords back to the environment
  list2env(check_data(vcf, coords), envir = environment())

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
    parallel = parallel,
    ncores = ncores,
    crop_edges = crop_edges
  )

  return(results)

}

#' General function for dist_gd
#' TODO: CREATE GENERAL CIRCLE GD AND RESIST GD
#' @param x data to be summarized by the moving window (*note:* order matters! `coords` should be in the same order, there are currently no checks for this). The class of `x` required depends on the statistic being calculated (see the `stat` argument and the function description for more details)
#' @param stat moving window statistic to calculate (can either be `pi` for nucleotide diversity (`x` must be a dosage matrix), `Ho` for average observed heterozygosity across all loci (`x` must be a heterozygosity matrix) , "allelic_richness" for average allelic richness across all loci (`x` must be a `genind` type object), "biallelic_richness" to get average allelic richness across all loci for a biallelic dataset (`x` must be a dosage matrix). `stat` can also be set to any function that will take `x`as input and return a single numeric value (for example, `x` can be a vector and `stat` can be set equal to a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `stat`, additional arguments to pass to the `stat` function (e.g. if `stat = mean`, users may want to set `na.rm = TRUE`)
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#'
#' @export
dist_general <- function(x, coords, lyr, stat, maxdist, distmat,
                         rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                         fun = mean, L = "nvariants", rarify_alleles = TRUE,
                         parallel = FALSE, ncores = NULL, crop_edges = FALSE, ...) {

  # Modify dist matrix
  distmat[distmat > maxdist] <- NA

  # run sliding window calculations
  if (parallel) {

    if (is.null(ncores)) ncores <- future::availableCores() - 1

    future::plan(future::multisession, workers = ncores)

    # currently, terra uses a C++ pointer which means SpatRasters cannot be directly passed to nodes on a computer cluster
    # instead of saving the raster layer to a file, I am converting it to a RasterLayer temporarily (it will get switched back)
    lyr <- raster::raster(lyr)

    rast_vals <- furrr::future_map(1:terra::ncell(lyr), dist_helper,
                                   lyr = lyr, x = x, distmat = distmat,
                                   stat = stat_function, rarify = rarify, rarify_n = rarify_n, rarify_nit = rarify_nit,
                                   min_n = min_n, fun = fun, L = L, rarify_alleles = rarify_alleles,
                                   .options = furrr::furrr_options(seed = TRUE, packages = c("wingen", "terra", "raster", "adegenet"))
    )

    # convert back to SpatRast
    lyr <- terra::rast(lyr)
  } else {

    rast_vals <- purrr::map(1:terra::ncell(lyr), dist_helper,
                            lyr = lyr, x = x, distmat = distmat,
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
#' @param lyr SpatRaster
#' @param x genetic data
#' @param distmat distance matrix
#'
#' @inheritParams dist_general
#'
#' @return genetic diversity and counts for a single cell
#'
#' @noRd
dist_helper <- function(i, lyr, x, distmat, stat_function,
                        rarify, rarify_n, rarify_nit, min_n,
                        fun, L = NULL, rarify_alleles = TRUE){

  # if rarify = TRUE, min_n = rarify_n (i.e. minimum defaults to rarify_n)
  if (rarify) min_n <- rarify_n

  # skip if raster value is NA
  # note: need to provide ns to give some output so the cell is counted
  # don't need to provide gd as this is repaired later
  if (is.na(lyr[i])) {
    return(c(sample_count = NA))
  }

  # get sample indices in window
  sub <- get_coord_index(i, distmat)

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

  # if gd has no name give call it custom
  if (is.null(names(gd))) names(gd) <- "custom"

  # return vector
  if (all(is.na(gd))) return(c(sample_count = ns)) else return(c(gd, sample_count = ns))
}


get_coord_index <- function(i, distmat) which(!is.na(distmat[i,]))
