#' Create a raster from coordinates
#'
#' Generate a raster layer from coordinates which can be used in \link[wingen]{window_gd} as the RasterLayer to move the window across
#'
#' @param coords coordinates of samples as sf points, a SpatVector, a two-column matrix, or a data.frame with x and y coordinates
#' @param buffer size of buffer to add to edge of raster (defaults to 0)
#' @param res desired resolution of raster (defaults to 1). Can be a single value for square cells or a vector with two values representing x and y resolutions
#' @param agg aggregation factor to apply to raster (defaults to NULL)
#' @param disagg disaggregation factor to apply to raster (defaults to NULL)
#' @param plot whether to plot resulting raster with coords (defaults to FALSE)
#'
#' @return RasterLayer
#' @export
#'
#' @examples
#' load_mini_ex()
#' coords_to_raster(mini_coords, buffer = 1, plot = TRUE)
coords_to_raster <- function(coords, buffer = 0, res = 1, agg = NULL, disagg = NULL, plot = FALSE) {
  # make a matrix
  r <- make_raster(coords, buffer = buffer, res = res)

  # aggregate or disaggregate
  if (!is.null(agg) & !is.null(disagg)) {
    warning("both agg and disagg were provided. Did you mean to do this? (if so, note that aggregation will occur first and then disaggregation second")
  }
  if (!is.null(agg)) r <- terra::aggregate(r, agg)
  if (!is.null(disagg)) r <- terra::disagg(r, disagg)

  # assign values to make it easier to visualize the resolution
  r <- terra::init(r, fun = 1:terra::ncell(r))

  # plot raster
  if (plot) {
    terra::plot(r, legend = FALSE, col = viridis::mako(terra::ncell(r)))
    if (is.matrix(coords)) coords <- data.frame(coords)
    terra::points(coords, col = viridis::magma(1, begin = 0.7), pch = 3, lwd = 2)
  }

  return(r)
}

#' coords to raster converter
#'
#' @inheritParams coords_to_raster
#'
#' @noRd
make_raster <- function(coords, buffer = 0, res = 1) {
  # format coords if not a spat vector
  if (inherits(coords, "sf")) coords <- terra::vect(coords)
  if (inherits(coords, "data.frame") | inherits(coords, "matrix")) coords <- terra::vect(as.matrix(coords), type = "points", atts = NULL)

  # turn into raster
  if (length(res) > 2) stop("invalid res provided")
  r <- terra::rast(coords, res = res)

  # extend based on buffer
  r <- terra::extend(r, buffer)

  return(r)
}
