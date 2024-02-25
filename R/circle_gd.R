#' Create a moving window map of genetic diversity using a circle window
#'
#' Generate a continuous raster map of genetic diversity using circle moving windows
#'
#' @param maxdist maximum geographic distance used to define neighborhood; any samples further than this distance will not be included (this can be thought of as the neighborhood radius)
#' Can either be (1) a single numeric value or (2) a SpatRaster where each pixel is the maximum distance to be used for that cell on the landscape (must be the same spatial scale as `lyr`).
#' @param distmat distance matrix output from \link[wingen]{get_geodist} (optional; can be used to save time on distance calculations)
#' @inheritParams window_gd
#' @details Coordinates and rasters should be in a Euclidean coordinate system (i.e., UTM coordinates) such that raster cell width and height are equal distances.
#' As such, longitude-latitude systems should be transformed before using dist_gd. Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#' @export
#'
#' @examples
#' \donttest{
#' load_mini_ex()
#' cpi <- circle_gd(mini_vcf, mini_coords, mini_lyr, fact = 2, maxdist = 5)
#' }
circle_gd <- function(gen, coords, lyr, maxdist, distmat = NULL, stat = "pi", fact = 0,
                      rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                      fun = mean, L = "nvariants", rarify_alleles = TRUE, sig = 0.05) {
  # convert lyr to SpatRaster
  if (!inherits(lyr, "SpatRaster")) lyr <- terra::rast(lyr)

  # convert coords if not in sf
  if (!inherits(coords, "sf")) coords <- coords_to_sf(coords)

  # make aggregated raster
  if (fact > 0) lyr <- terra::aggregate(lyr, fact, fun = mean)

  # make distmat
  if (is.null(distmat)) suppressWarnings(distmat <- get_geodist(coords, lyr))

  # run dist_gd
  results <-
    dist_gd(
      gen = gen,
      coords = coords,
      lyr = lyr,
      maxdist = maxdist,
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
      sig = sig
    )

  return(results)
}


#' General function for making circular moving window maps
#'
#' Generate a continuous raster map using circular moving windows.
#' While \link[wingen]{resist_gd} is built specifically for making maps
#' of genetic diversity from vcfs,`circle_general` can be used to make maps
#' from different data inputs. Unlike `resist_gd`, `resist_general`
#' will not convert your data into the correct format for calculations of different
#' diversity metrics. See details for how to format data inputs for different statistics.
#'
#' @param x data to be summarized by the moving window (*note:* order matters! `coords` should be in the same order, there are currently no checks for this). The class of `x` required depends on the statistic being calculated (see the `stat` argument and the function description for more details)
#' @param stat moving window statistic to calculate (see details). `stat` can generally be set to any function that will take `x`as input and return a single numeric value (for example, `x` can be a vector and `stat` can be set equal to a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `stat`, additional arguments to pass to the `stat` function (e.g. if `stat = mean`, users may want to set `na.rm = TRUE`)
#' @inheritParams window_general
#' @inheritParams circle_gd
#'
#' @details
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
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#'
#' @export
circle_general <- function(x, coords, lyr, maxdist, distmat = NULL, stat, fact = 0,
                           rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                           fun = mean, L = NULL, rarify_alleles = TRUE, sig = 0.05, ...) {
  # check and aggregate layer and coords  (only lyr is returned)
  lyr <- layer_coords_check(lyr = lyr, coords = coords, fact = fact)

  # make distmat
  if (is.null(distmat)) suppressWarnings(distmat <- get_geodist(coords, lyr))

  # run general resist
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


#' Get a matrix of geographic distances for \link[wingen]{circle_gd}
#'
#' Create a distance matrix based on coordinates and a raster layer.
#' The output is a distance matrix where rows represent cells on the landscape
#' and columns represent individual locations on the landscape. Each value is
#' the geographic distance between each individual and each cell calculated
#' using \link[sf]{st_distance}. This matrix is used by \link[wingen]{circle_gd}.
#' If coords_only = TRUE, the result is a distance matrix for the sample coordinates
#' only.
#'
#' @param lyr SpatRaster or RasterLayer for generating distances (not required if coords_only = TRUE)
#' @param coords_only whether to return distances only for sample coordinates
#' @inheritParams circle_gd
#'
#' @return a distance matrix used by \link[wingen]{circle_gd}
#' @export
#'
#' @examples
#' load_mini_ex()
#' distmat <- get_geodist(mini_coords, mini_lyr)
#'
get_geodist <- function(coords, lyr = NULL, fact = 0, coords_only = FALSE) {
  # convert coords if not in sf
  if (!inherits(coords, "sf")) coords <- coords_to_sf(coords)

  # create distance matrix using only coordinates
  if (coords_only) {
    return(sf::st_distance(coords))
  }

  # check crs
  layer_coords_check(lyr = lyr, coords = coords)

  # apply CRS so that if one is NA and the other isn't you don't get an error
  # assumes CRS are the same (warning printed by layer_coords_check)
  if (is.na(sf::st_crs(lyr)) & !is.na(sf::st_crs(coords))) sf::st_crs(lyr) <- sf::st_crs(coords)
  if (!is.na(sf::st_crs(lyr)) & is.na(sf::st_crs(coords))) sf::st_crs(coords) <- sf::st_crs(lyr)

  # aggregate raster
  if (fact != 0) lyr <- terra::aggregate(lyr, fact, fun = mean)

  # make into df
  lyr_df <- terra::as.data.frame(lyr, xy = TRUE, na.rm = FALSE)
  lyr_sf <- sf::st_as_sf(lyr_df, coords = c("x", "y"), crs = terra::crs(lyr))

  # .y = lyr_sf
  # .x = index
  distls <-
    furrr::future_map(
      1:nrow(lyr_sf),
      ~ sf::st_distance(.y[.x, ], coords),
      lyr_sf,
      .options = furrr::furrr_options(seed = TRUE, packages = c("sf")),
      .progress = TRUE
    )

  # convert into matrix
  distmat <- matrix(unlist(distls), ncol = length(distls[[1]]), byrow = TRUE)

  # transpose distmat so that the rows are the samples
  distmat <- t(distmat)

  return(distmat)
}
