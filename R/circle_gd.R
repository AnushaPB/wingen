
#' Create a moving window map of genetic diversity
#'
#' Generate a continuous raster map of genetic diversity using moving windows
#'
#' @param vcf object of type vcf or a path to a vcf file (*note:* order matters! The coordinate and genetic data should be in the same order; there are currently no checks for this)
#' @param coords coordinates of samples as an sf object, a two-column matrix, or a data.frame representing x and y coordinates. Should be in a Euclidean system (i.e., not longitude latitude) or the window cell height and width will not be equal (see details).
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
#' As such, longitude-latitude systems should be transformed before using circle_gd. Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#' @export
#'
#' @examples
#'
#' load_mini_ex()
#' wpi <- circle_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE)
#' plot_gd(wpi, main = "Window pi")
#' plot_count(wpi)
#'
circle_gd <- function(vcf, coords, lyr, radius, distmat = NULL, stat = "pi", fact = 0,
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
  result <- circle_general(
    x = x,
    coords = coords,
    lyr = lyr,
    radius = radius,
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
#' Generate a continuous raster map using moving windows. While \link[wingen]{circle_gd} is built specifically for making moving window maps of genetic diversity from vcfs,
#' `circle_general` can be used to make moving window maps from different data inputs. Unlike `circle_gd`, `circle_general` will not convert your data into
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
#' @inheritParams circle_gd
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#'
#' @export
circle_general <- function(x, coords, lyr, stat, radius, distmat = NULL, fact = 0,
                           rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                           fun = mean, L = "nvariants", rarify_alleles = TRUE,
                           parallel = FALSE, ncores = NULL, crop_edges = FALSE, ...) {

  # check layers and coords (only lyr is modified and returned)
  lyr <- layer_coords_check(lyr, coords)

  # set L if pi is being calculated
  if (!is.null(L) & !is.numeric(L)) if (L == "nvariants") L <- ncol(x)

  # Get function to calculate the desired statistic
  stat_function <- return_stat(stat, ...)

  # check that coords and x align and reformat data, if necessary
  # note: list2env adds the new, corrected x and coords back to the environment
  list2env(check_data(x, coords), envir = environment())

  # make aggregated raster
  if (fact == 0) lyr <- lyr * 0 else lyr <- terra::aggregate(lyr, fact, fun = mean) * 0

  # make distmat if not provided
  if (is.null(distmat)) distmat <- make_distmat(coords, lyr, parallel = parallel, ncores = ncores)
  distmat[distmat > radius] <- NA

  # run sliding window calculations
  if (parallel) {

    if (is.null(ncores)) ncores <- future::availableCores() - 1

    future::plan(future::multisession, workers = ncores)

    # currently, terra uses a C++ pointer which means SpatRasters cannot be directly passed to nodes on a computer cluster
    # instead of saving the raster layer to a file, I am converting it to a RasterLayer temporarily (it will get switched back)
    lyr <- raster::raster(lyr)

    rast_vals <- furrr::future_map(1:terra::ncell(lyr), circle_helper,
                                       lyr = lyr, x = x, distmat = distmat,
                                       stat = stat_function, rarify = rarify, rarify_n = rarify_n, rarify_nit = rarify_nit,
                                       min_n = min_n, fun = fun, L = L, rarify_alleles = rarify_alleles,
                                       .options = furrr::furrr_options(seed = TRUE, packages = c("wingen", "terra", "raster", "adegenet"))
    )

    # convert back to SpatRast
    lyr <- terra::rast(lyr)
  } else {

    rast_vals <- purrr::map(1:terra::ncell(lyr), circle_helper,
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
#' @param coord_cells cell indices for each coordinate
#' @param nmat neighborhood matrix
#'
#' @inheritParams circle_general
#'
#' @return genetic diversity and counts for a single cell
#'
#' @noRd
circle_helper <- function(i, lyr, x, distmat, stat_function,
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

make_distmat <- function(coords, lyr, parallel = FALSE, ncores = NULL){
  lyr_df <- terra::as.data.frame(lyr, xy = TRUE, na.rm = FALSE)
  lyr_sf <- sf::st_as_sf(lyr_df, coords = c("x", "y"), crs = terra::crs(lyr))

  if (parallel){
    if (is.null(ncores)) ncores <- future::availableCores() - 1
    distls <- furrr::future_map(1:nrow(lyr_sf), ~ sf::st_distance(.y[.x,], coords), lyr_sf,
                                .options = furrr::furrr_options(seed = TRUE, packages = c("sf")))
  } else {
    distls <- purrr::map(1:nrow(lyr_sf), ~ sf::st_distance(.y[.x,], coords), lyr_sf)
  }


  distmat <- matrix(unlist(distls), ncol = length(distls[[1]]), byrow = TRUE)

  return(distmat)
}

get_coord_index <- function(i, distmat) which(!is.na(distmat[i,]))


