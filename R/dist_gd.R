
#' Run distance based moving window calculations
#'
#' Function used by both circle_gd and resist_gd
#'
#' @noRd
dist_gd <- function(vcf, coords, lyr, distmat = NULL, stat = "pi", fact = 0,
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
  result <- dist_general(
    x = x,
    coords = coords,
    lyr = lyr,
    distmat = distmat,
    stat = stat,
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
#' Generate a continuous raster map using moving windows. While \link[wingen]{dist_gd} is built specifically for making moving window maps of genetic diversity from vcfs,
#' `dist_general` can be used to make moving window maps from different data inputs. Unlike `dist_gd`, `dist_general` will not convert your data into
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
#' @inheritParams dist_gd
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#'
#' @export
dist_general <- function(x, coords, lyr, stat, maxdist, distmat, fact = 0,
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
  if (all(is.na(distmat[i,]))) {
    return(data.frame(gd = NA, ns = NA))
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

  return(data.frame(gd = gd, ns = ns))
}


get_coord_index <- function(i, distmat) which(!is.na(distmat[i,]))
