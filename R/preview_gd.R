#' Preview moving window and sample counts
#'
#' Generate preview of moving window size and sample counts based on the coordinates and parameters to be supplied to \link[wingen]{window_gd}
#' @param lyr SpatRaster or RasterLayer to slide the window across (see Details for important information about projections). For `method = "resist"` this should also be the conductivity layer (see \link[wingen]{resist_gd()})
#' @param method which method to use to create preview (`"window"` for \link[wingen]{window_gd()}, `"circle"` for \link[wingen]{circle_gd()}, or `"resist"` for \link[wingen]{resist_gd()}; defaults to `"window"`)
#' @param sample_count whether to create plot of sample counts for each cell (defaults to TRUE)
#' @param wdim if `method = "window"`, dimensions (height x width) of window; if only one value is provided, a square window is created (defaults to 3 x 3 window)
#' @param distmat if `method = "circle"` or `method = "resist"`, an optional distance matrix to be used output from either \link[wingen]{get_geodist()} or \link[wingen]{get_resistdist()}, respectively. If not provided, one will be automatically calculated.
#' @param maxdist if `method = "circle"` or `method = "resist`, the maximum geographic distance used to define the neighborhood; any samples further than this distance will not be included (see \link[wingen]{get_geodist()} or \link[wingen]{get_resistdist()})
#' @param min_n min number of samples to use in calculations (any focal cell with a window containing less than this number of samples will be assigned a value of NA)
#' @param plot whether to plot results (default = TRUE)
#' @param parallel whether to use parallelization for calculating the distance matrices for `method = "circle"` or `method = "resist` (defaults to `FALSE`).
#' @inheritParams window_gd
#' @details
#'
#' Coordinates and rasters should be in a projected (planar) coordinate system such that raster cells are of equal sizes.
#' Therefore, spherical systems (including latitute-longitude coordinate systems) should be projected prior to use.
#' Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return SpatRaster with sample counts layer (if sample_count = TRUE)
#' @export
#'
#' @examples
#' load_mini_ex()
#' preview_gd(mini_lyr, mini_coords, wdim = 3, fact = 3, sample_count = TRUE, min_n = 2)
preview_gd <- function(lyr, coords, method = "window", wdim = NULL, maxdist = NULL, distmat = NULL, fact = 0, sample_count = TRUE, min_n = 0, plot = TRUE, parallel = FALSE, ncores = NULL) {
  # convert to spat rast
  if (!inherits(lyr, "SpatRaster")) lyr <- terra::rast(lyr)
  if (fact != 0) lyr <- terra::aggregate(lyr, fact)

  if (method == "window"){
    # plot window preview
    if (plot) preview_window(lyr, wdim, coords)
  } else {
    # convert coords if not in sf
    if (!inherits(coords, "sf")) coords <- coords_to_sf(coords)

    # check distmat
    if (!is.null(distmat)) if (terra::ncell(lyr) != nrow(distmat)) stop("Number of cells in raster layer and number of columns of distmat do not match")

    # make distmat
    if (is.null(distmat) & method == "circle") distmat <- get_geodist(coords, lyr, parallel = parallel, ncores = ncores)
    if (is.null(distmat) & method == "resist") distmat <- get_resdist(coords, con_lyr = lyr, parallel = parallel, ncores = ncores)

    # Modify dist matrix
    distmat[distmat > maxdist] <- NA
  }

  # plot count preview and return count raster
  if (sample_count) {
    lyrc <- preview_count(lyr = lyr, coords = coords, wdim = wdim, distmat = distmat, min_n = min_n, plot = plot)
    return(lyrc)
  }
}

#' Plot preview of moving window
#'
#' @param lyr RasterLayer
#' @param wdim window dimensions
#' @param coords coordinates
#'
#' @noRd
preview_window <- function(lyr, wdim, coords = NULL) {
  # convert wdim to matrix
  nmat <- wdim_to_mat(wdim)

  # get center of raster
  center <- get_center(lyr)

  # get adjacent cells to center cell
  adjc <- terra::adjacent(lyr, center, directions = nmat)
  # get list of indices of coords in that set of cells
  adjci <- purrr::map_dbl(adjc, 1, ~ seq(.x[1], .x[2]))

  # fill in window
  lyrw <- lyr * 0
  lyrw[adjci] <- 1
  lyrw[center] <- 2

  # suppress irrelevant plot warnings
  suppressWarnings({
    terra::plot(lyrw, col = viridis::mako(3, direction = -1), legend = FALSE, axes = FALSE, box = FALSE)
    graphics::legend("bottomleft", c("raster layer", "window", "focal cell"), col = viridis::mako(3, direction = -1), pch = 15)
    if (!is.null(coords)) {
      if (is.matrix(coords)) coords <- data.frame(coords)
      # note: sf coords also inherit "data.frame" so second condition is needed
      if (inherits(coords, "data.frame") & !inherits(coords, "sf")) terra::points(coords, pch = 3, col = viridis::magma(1, begin = 0.7))
      if (inherits(coords, "sf") | inherits(coords, "SpatVector")) terra::plot(coords, pch = 3, col = viridis::magma(1, begin = 0.7), add = TRUE)
    }
  })
}

#' Get center cell of a raster
#'
#' @param x raster
#'
#' @noRd
get_center <- function(x) {
  e <- as.vector(terra::ext(x))
  c <- c(mean(e[c(1, 2)]), mean(e[c(3, 4)]))
  center <- terra::cellFromXY(x, xy = matrix(c, ncol = 2))
  return(center)
}

#' Plot preview of sample count layer
#'
#' @param lyr RasterLayer
#' @param coords coordinates
#' @param wdim window dimensions
#' @param plot whether to plot resuls
#'
#' @noRd
preview_count <- function(lyr, coords, wdim, distmat, min_n, plot = TRUE) {
  # get coord cells
  coord_cells <- terra::extract(lyr, coords, cell = TRUE)[, "cell"]

  # make copy of raster for counting
  lyrc <- lyr

  # get counts for each raster cell
  if (!is.null(wdim)) nc <- purrr::map_dbl(1:terra::ncell(lyr), sample_count_wdim, lyr, wdim_to_mat(wdim), coord_cells)
  if (!is.null(distmat)) nc <- purrr::map_dbl(1:terra::ncell(lyr), sample_count_dist, distmat)

  # assign values to raster
  lyrc <- terra::setValues(lyr, nc)

  # mask areas where counts are less than min value
  lyrc[lyrc < min_n] <- NA

  # set name
  names(lyrc) <- "sample_count"

  # plot results
  if (plot) {
    suppressWarnings(terra::plot(lyrc, col = viridis::mako(100), box = FALSE, axes = FALSE))
    graphics::title(main = list("Sample Count", font = 1), adj = 0, line = -0.5)
  }

  return(lyrc)
}

#' Count samples in rectangular window around focal cell
#'
#' @param x focal cell index
#' @param lyr raster layer
#' @param nmat neighborhood matrix
#' @param coord_cells cell indexes of coordinates
#'
#' @noRd
sample_count_wdim <- function(x, lyr, nmat, coord_cells) {
  sub <- get_adj(x, lyr, nmat, coord_cells)
  return(length(sub))
}

#' Count samples in circle or resistance window around focal cell
#'
#' @param x focal cell index
#' @param distmat distance matrix
#'
#' @noRd
sample_count_dist <- function(x, distmat) {
  sub <- get_dist_index(x, distmat)
  return(length(sub))
}
