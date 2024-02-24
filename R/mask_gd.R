#' Mask moving window maps
#'
#' Mask genetic diversity layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x Raster object to mask
#' @param y Raster object or Spatial object to use as mask
#' @param minval if y is a Raster object, value of y below which to mask
#' @param maxval if y is a Raster object, value of y above which to mask
#'
#' @return RasterLayer
#' @export
#'
#' @examples
#' data("mini_lyr")
#' mpi <- mask_gd(mini_lyr, mini_lyr, minval = 0.01)
#'
mask_gd <- function(x, y, minval = NULL, maxval = NULL) {
  # make sure x is a SpatRaster
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)

  # match raster layers
  if (inherits(y, "RasterLayer") | inherits(y, "RasterStack") | inherits(y, "RasterBrick") | inherits(y, "SpatRaster")) {
    # convert raster
    if (!inherits(y, "SpatRaster")) y <- terra::rast(y)

    # y areas below min/max val if provided
    if (!is.null(minval)) {
      y[y < minval] <- NA
    }

    if (!is.null(maxval)) {
      y[y > maxval] <- NA
    }
  }

  # if not a raster convert to a spat vector
  if (!inherits(y, "SpatRaster") & !inherits(y, "SpatVector")) y <- terra::vect(y)

  # mask
  x <- terra::mask(x, y)

  return(x)
}
