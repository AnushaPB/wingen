
#' Preview moving window and sample counts
#'
#' Generate preview of moving window size and sample counts based on the coordinates and parameters to be supplied to \link[wingen]{window_gd}
#'
#' @param sample_count whether to create plot of sample counts for each cell (defaults to TRUE)
#' @param min_n min number of samples to use in calculations (any focal cell with a window containing less than this number of samples will be assigned a value of NA)
#' @inheritParams window_gd
#' @details
#'
#' Coordinates and rasters should be in a projected (planar) coordinate system such that raster cells are of equal sizes.
#' Therefore, spherical systems (including latitute-longitude coordinate systems) should be projected prior to use.
#' Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return RasterStack with example window layer and sample counts (if sample_count = TRUE)
#' @export
#'
#' @examples
#' load_mini_ex()
#' preview_gd(mini_lyr, mini_coords, wdim = 3, fact = 3, sample_count = TRUE, min_n = 2)
preview_gd <- function(lyr, coords, wdim, fact = 0, sample_count = TRUE, min_n = 0) {
  # convert to spat rast
  if (!inherits(lyr, "SpatRaster")) lyr <- terra::rast(lyr)
  if (fact != 0) lyr <- terra::aggregate(lyr, fact)

  # convert wdim to matrix
  nmat <- wdim_to_mat(wdim)

  # plot window preview
  preview_window(lyr, nmat, coords)

  # plot count preview and return count raster
  if (sample_count) {
    lyrc <- preview_count(lyr, coords, nmat, min_n)
    return(lyrc)
  }
}

#' Plot preview of moving window
#'
#' @param lyr RasterLayer
#' @param nmat neighbor matrix
#' @param coords coordinates
#'
#' @noRd
preview_window <- function(lyr, nmat, coords = NULL) {
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

  # suppress annoying and irrelevant plot warnings
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
#' @param nmat neighborhood matrix
#'
#' @noRd
preview_count <- function(lyr, coords, nmat, min_n) {
  # get coord cells
  coord_cells <- terra::extract(lyr, coords, cell = TRUE)[, "cell"]

  # make copy of raster for counting
  lyrc <- lyr

  # get counts for each raster cell
  nc <- purrr::map_dbl(1:terra::ncell(lyr), sample_count, lyr, nmat, coord_cells)

  # assign values to raster
  lyrc <- terra::setValues(lyr, nc)

  # mask areas where counts are less than min value
  lyrc[lyrc < min_n] <- NA

  # set name
  names(lyrc) <- "sample_count"

  # plot results
  suppressWarnings(terra::plot(lyrc, col = viridis::mako(100), box = FALSE, axes = FALSE))
  graphics::title(main = list("Sample Count", font = 1), adj = 0, line = -0.5)

  return(lyrc)
}

#' Count samples in window around focal cell
#'
#' @param x focal cell index
#' @param lyr raster layer
#' @param nmat neighborhood matrix
#' @param coord_cells cell indexes of coordinates
#'
#' @noRd
sample_count <- function(x, lyr, nmat, coord_cells) {
  sub <- get_adj(x, lyr, nmat, coord_cells)
  return(length(sub))
}
