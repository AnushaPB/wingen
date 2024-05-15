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
#' @details
#'
#' Coordinates and rasters should be in a projected (planar) coordinate system such that raster cells are of equal sizes.
#' Therefore, spherical systems (including latitute-longitude coordinate systems) should be projected prior to use.
#' Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' Current genetic diversity metrics that can be specified with `stat` include:
#' - `"pi"` for nucleotide diversity (default) calculated using `hierfstat` \link[hierfstat]{pi.dosage}. Use the `L` argument to set the sequence length (defaults to dividing by the number of variants).
#' - `"Ho"` for average observed heterozygosity across all sites
#' - `"allelic_richness"` for average number of alleles across all sites
#' - `"biallelic_richness"` for average allelic richness across all sites for a biallelic dataset (this option is faster than `"allelic_richness"`)
#' - `"hwe"` for the proportion of sites that are not in Hardy–Weinberg equilibrium, calculated using `pegas` \link[pegas]{hw.test} at the 0.05 level (other alpha levels  can be specified by adding the sig argument; e.g., `sig = 0.10`).
#' - `"basic_stats"` for a series of statistics produced by `hierfstat` \link[hierfstat]{basic.stats} including
#' mean observed heterozygosity (same as Ho), mean gene diversities within population (Hs),
#' Gene diversities overall (Ht), and Fis following Nei (1987). Population-based statistics (e.g., FST) normally reported by \link[hierfstat]{basic.stats}
#' are not included as they are not meaningful within the individual-based moving windows.
#'
#' @examples
#' load_mini_ex()
#' rpi <- resist_gd(mini_vcf, mini_coords, mini_lyr, maxdist = 50)
resist_gd <- function(gen, coords, lyr, maxdist, distmat = NULL, stat = "pi", fact = 0,
                      rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                      fun = mean, L = "nvariants", rarify_alleles = TRUE, sig = 0.05,
                      transitionFunction = mean, directions = 8, geoCorrection = TRUE) {
  # check and aggregate layer and coords  (only lyr is returned)
  lyr <- layer_coords_check(lyr = lyr, coords = coords, fact = fact)

  # make distmat
  if (is.null(distmat)) suppressWarnings(distmat <- get_resdist(coords, lyr = lyr, transitionFunction = transitionFunction, directions = directions, geoCorrection = geoCorrection))

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
resist_general <- function(x, coords, lyr, maxdist, distmat = NULL, stat, fact = 0,
                           rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                           fun = mean, L = "nvariants", rarify_alleles = TRUE, sig = 0.05,
                           transitionFunction = mean, directions = 8, geoCorrection = TRUE, ...) {
  # check and aggregate layer and coords  (only lyr is returned)
  lyr <- layer_coords_check(lyr = lyr, coords = coords, fact = fact)

  # make distmat
  if (is.null(distmat)) {
    suppressWarnings(
      distmat <-
        get_resdist(
          coords,
          lyr = lyr,
          transitionFunction = transitionFunction,
          directions = directions,
          geoCorrection = geoCorrection
        )
    )
  }

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


#' Get a matrix of resistance distances for \link[wingen]{resist_gd}
#'
#' Create a distance matrix based on coordinates and a connectivity layer.
#' The output is a distance matrix where rows represent cells on the landscape
#' and columns represent individual locations on the landscape. Each value is
#' the resistance distance between each sample and each cell calculated
#' using the gdistance package. This matrix is used by \link[wingen]{resist_gd}.
#' If coords_only = TRUE, the result is a distance matrix for the sample coordinates
#' only.
#' @param lyr conductivity layer (higher values should mean greater conductivity) for generating distances. Can be either a SpatRaster or RasterLayer.
#' @inheritParams resist_gd
#' @param coords_only whether to return distances only for sample coordinates
#' @return a distance matrix used by \link[wingen]{resist_gd}
#' @export
#'
#' @examples
#' load_mini_ex()
#' distmat <- get_resdist(mini_coords, mini_lyr)
get_resdist <- function(coords, lyr, fact = 0, transitionFunction = mean, directions = 8, geoCorrection = TRUE, coords_only = FALSE) {
  # check lyr and coords
  lyr <- layer_coords_check(lyr = lyr, coords = coords, fact = fact)

  # convert lyr to raster
  if (!inherits(lyr, "RasterLayer")) lyr <- raster::raster(lyr)

  # convert coords to matrix and rename
  coords_mat <- as.matrix(coords_to_df(coords))

  # create transition surface
  trSurface <- gdistance::transition(lyr, transitionFunction = transitionFunction, directions = directions)

  # note: type = "c" is for least cost distances
  if (geoCorrection) trSurface <- gdistance::geoCorrection(trSurface, type = "c", scl = FALSE)

  # create distance matrix using only coordinates
  if (coords_only) {
    return(as.matrix(gdistance::commuteDistance(trSurface, coords_mat)))
  }

  # make vector of distances
  distrasts <-
    furrr::future_map(
      1:nrow(coords_mat),
      ~ gdistance::accCost(trSurface, coords_mat[.x, ]),
      .options = furrr::furrr_options(seed = TRUE, packages = c("gdistance")),
      .progress = TRUE
    )

  # convert from raster to matrix
  diststack <- terra::rast(raster::stack(distrasts))
  distmat <- as.matrix(terra::as.data.frame(diststack, na.rm = FALSE))

  # transpose distmat so that the rows are the samples
  distmat <- t(distmat)

  return(distmat)
}

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
  if (inherits(coords, "sf")) coords <- data.frame(sf::st_coordinates(coords))
  colnames(coords) <- c("x", "y")
  return(coords)
}
