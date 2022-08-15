
#' Mask moving window maps
#'
#' Mask genetic diversity layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x Raster object to mask
#' @param mask Raster object or Spatial object to use as mask
#' @param resample if x and mask are non matching rasters, which layer to resample to match them (defaults to "mask")
#' @param minval if mask is a Raster object, value of mask below which to mask
#' @param maxval if mask is a Raster object, value of mask above which to mask
#'
#' @return RasterLayer
#' @export
#'
#' @examples
#' library("raster")
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, L = 10, rarify = TRUE)
#' kpi <- krig_gd(wpi, mini_lyr)
#' mpi <- mask_gd(kpi, mini_lyr, minval = 0.01)
#' plot_gd(mpi, main = "Kriged and Masked Pi")
#'
mask_gd <- function(x, mask, resample = "mask", minval = NULL, maxval = NULL) {

  # match raster layers
  if (class(mask) == "RasterLayer" | class(mask) == "RasterStack" | class(mask) == "RasterBrick") {

    # mask areas below min/max val if provided
    if (!is.null(minval)) {
      mask[mask < minval] <- NA
    }

    if (!is.null(maxval)) {
      mask[mask > maxval] <- NA
    }

    if (!raster::compareRaster(x, mask, stopiffalse = FALSE)) {
      if (resample == "mask") mask <- raster::resample(mask, x)
      if (resample == "x") x <- raster::resample(x, mask)
      if (resample != "x" & resample != "mask") stop("invalid arugment provided for resample (must be \"x\" or \"mask\")")
    }
  }

  # mask
  x <- raster::mask(x, mask)

  return(x)
}
