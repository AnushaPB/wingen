#' Create a moving window map of genetic diversity based on resistance
#'
#' Generate a continuous raster map of genetic diversity using resistance distances calculated with a conductivity surface
#'
#' @param maxdist maximum cost distance used to define neighborhood; any samples further than this cost distance will not be included (this can be thought of as the neighborhood radius, but in terms of cost distance).
#' Can either be (1) a single numeric value or (2) a SpatRaster where each pixel is the maximum distance to be used for that cell on the landscape (must be the same spatial scale as `lyr`).
#' @param lyr conductivity layer (higher values should mean greater conductivity) to move window across. Can be either a SpatRaster or RasterLayer.
#' @param distmat distance matrix output from \link[wingen]{get_resdist} (optional; can be used to save time on distance calculations)
#' @param transitionFunction function to calculate transition values from grid values (defaults to mean)
#' @param directions directions in which cells are connected (4, 8, 16, or other), see \link[raster]{adjacent} (defaults to 8)
#' @param geoCorrection whether to apply correction to account for local distances (defaults to TRUE). Geographic correction is necessary for all objects of the class Transition that are either: (1) based on a grid in a geographic (lonlat) projection and covering a large area; (2) made with directions > 4 (see \link[gdistance]{geoCorrection} for more details).
#' @inheritParams window_gd
#' @details Coordinates and rasters should be in a Euclidean coordinate system (i.e., UTM coordinates) such that raster cell width and height are equal distances.
#' As such, longitude-latitude systems should be transformed before using dist_gd. Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' load_mini_ex()
#' rpi <- resist_gd(mini_vcf, mini_coords, mini_lyr, maxdist = 500)
#' plot_gd(rpi, main = "Resist pi")
#' plot_count(rpi)
#' }
resist_gd <- function(gen, coords, lyr, maxdist, distmat = NULL, stat = "pi", fact = 0,
                      rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                      fun = mean, L = "nvariants", rarify_alleles = TRUE,
                      transitionFunction = mean, directions = 8, geoCorrection = TRUE,
                      parallel = FALSE, ncores = NULL) {

  # check and aggregate layer and coords  (only lyr is returned)
  lyr <- layer_coords_check(lyr = lyr, coords = coords, fact = fact)

  # make distmat
  if (is.null(distmat)) distmat <- get_resdist(coords, lyr = lyr, transitionFunction = transitionFunction, directions = directions, geoCorrection = geoCorrection, parallel = parallel, ncores = ncores)

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
      parallel = parallel,
      ncores = ncores
    )

  return(results)
}

#' General function for making resistance-based maps
#'
#' Generate a continuous raster map using resistance distances.
#' While \link[wingen]{resist_gd} is built specifically for making maps
#' of genetic diversity from vcfs,`resist_general` can be used to make maps
#' from different data inputs. Unlike `resist_gd`, `resist_general`
#' will not convert your data into the correct format for calculations of different
#' diversity metrics. See details for how to format data inputs for different statistics.
#'
#' @param x data to be summarized by the moving window (*note:* order matters! `coords` should be in the same order, there are currently no checks for this). The class of `x` required depends on the statistic being calculated (see the `stat` argument and the function description for more details)
#' @param stat moving window statistic to calculate (see details). `stat` can generally be set to any function that will take `x`as input and return a single numeric value (for example, `x` can be a vector and `stat` can be set equal to a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `stat`, additional arguments to pass to the `stat` function (e.g. if `stat = mean`, users may want to set `na.rm = TRUE`)
#' @inheritParams window_general
#' @inheritParams resist_gd
#'
#' @details
#' To calculate genetic diversity statistics with the built in wingen functions, data must be formatted as such:
#' - for `"pi"` or `"biallelic_richness"`, `x` must be a dosage matrix with values of 0, 1, or 2
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
resist_general <- function(x, coords, lyr, maxdist, distmat, stat, fact = 0,
                           rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                           fun = mean, L = NULL, rarify_alleles = TRUE,
                           transitionFunction = mean, directions = 8, geoCorrection = TRUE,
                           parallel = FALSE, ncores = NULL, ...) {
  # check and aggregate layer and coords  (only lyr is returned)
  lyr <- layer_coords_check(lyr = lyr, coords = coords, fact = fact)

  # make distmat
  if (is.null(distmat)) distmat <- get_resdist(coords, lyr = lyr, transitionFunction = transitionFunction, directions = directions, geoCorrection = geoCorrection, parallel = parallel, ncores = ncores)

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
    parallel = parallel,
    ncores = ncores,
  )

  return(results)
}


#' Get a matrix of resistance distances for \link[wingen]{resist_gd}
#'
#' Create a distance matrix based on coordinates and a connectivity layer.
#' The output is a distance matrix where rows represent cells on the landscape
#' and columns represent individual locations on the landscape. Each value is
#' the resistance distance between each individual and each cell calculated
#' using the gdistance package. This matrix is used by \link[wingen]{resist_gd}.
#'
#' @inheritParams resist_gd
#' @return a distance matrix used by \link[wingen]{resist_gd}
#' @export
#'
#' @examples
#' \dontrun{
#' load_mini_ex()
#' distmat <- get_resdist(mini_coords, mini_lyr)
#' }
get_resdist <- function(coords, lyr, fact = 0, transitionFunction = mean, directions = 8, geoCorrection = TRUE, ncores = NULL, parallel = TRUE) {
  # convert lyr to raster
  if (!inherits(lyr, "RasterLayer")) lyr <- raster::raster(lyr)

  # aggregate raster
  if (fact != 0) lyr <- terra::aggregate(lyr, fact, fun = mean)

  # convert coords to dataframe and rename
  coords_df <- coords_to_df(coords)

  # create transition surface
  trSurface <- gdistance::transition(lyr, transitionFunction = transitionFunction, directions = directions)
  # note: type = "c" is for least cost distances
  if (geoCorrection) trSurface <- gdistance::geoCorrection(trSurface, type = "c", scl = FALSE)

  # get layer coordinates
  lyr_coords <- terra::as.data.frame(lyr, xy = TRUE, na.rm = FALSE)[, 1:2]

  # get all combinations of layer and sample coordinate indices
  params <- expand.grid(list(lyr = 1:nrow(lyr_coords), coords = 1:nrow(coords_df)))

  # make vector of distances
  if (parallel) {
    if (is.null(ncores)) ncores <- future::availableCores() - 1

    future::plan(future::multisession, workers = ncores)

    suppressWarnings({
      distvec <- furrr::future_map2_dbl(params$lyr, params$coords, run_gdist, trSurface, lyr_coords, coords_df, .progress = TRUE)
    })

    future::plan("sequential")
  } else {
    suppressWarnings({
      distvec <- purrr::map2_dbl(params$lyr, params$coords, run_gdist, trSurface, lyr_coords, coords_df, .progress = TRUE)
    })
  }

  # convert from vector to matrix
  distmat <-
    data.frame(params, dist = distvec) %>%
    tidyr::pivot_wider(names_from = "coords", values_from = "dist") %>%
    dplyr::select(-1) %>%
    as.matrix()

  return(distmat)
}

#' Run gdistance
#'
#' @param x index for lyr_coords
#' @param y index for coords_df
#' @param trSurface transition surface
#' @param lyr_coords raster coordinates
#' @param coords_df individual coordinates
#'
#' @return single resistance distance value
#'
#' @noRd
run_gdist <- function(x, y, trSurface, lyr_coords, coords_df) {
  # make spatial points
  sp <- sp::SpatialPoints(rbind(lyr_coords[x, ], coords_df[y, ]))

  # calculate circuit distances
  d <- possible_gdist(trSurface, sp)

  # get distance
  if (!all(is.na(d))) d <- d[1, 2]

  return(d)
}

#' Possibly wrapper for gdistance to return NA instead of error
#'
#' @param sp pair of coordinates
#' @param trSurface transition surface
#'
#' @return distance matrix
#'
#' @noRd
possible_gdist <- purrr::possibly(function(trSurface, sp) as.matrix(gdistance::commuteDistance(trSurface, sp)), NA)

#' Convert coordinates to dataframe
#'
#' @param coords coordinates
#'
#' @return dataframe of coordinates
#'
#' @noRd
coords_to_df <- function(coords) {
  if (inherits(coords, "SpatVector")) coords <- sf::st_as_sf(coords)
  if (is.matrix(coords)) coords <- data.frame(coords)
  if (inherits(coords, "sf")) coords <- as.data.frame(sf::as_Spatial(coords))
  colnames(coords) <- c("x", "y")
  return(coords)
}
