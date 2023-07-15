#' Check CRS of coords and layer for window_gd
#'
#' @param lyr RasterLayer or SpatRaster
#' @param coords sf object, data frame, or matrix representing coordinates
#'
#' @return NULL
#'
#' @noRd
crs_check_window <- function(lyr, coords) {
  # change matrix to df
  if (is.matrix(coords)) coords <- data.frame(coords)

  # get CRS of each object
  coords_crs <- sf::st_crs(coords)
  lyr_crs <- sf::st_crs(lyr)

  if (is.na(coords_crs)) warning("No CRS found for the provided coordinates. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)\n")

  if (is.na(lyr_crs)) warning("No CRS found for the provided raster. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)\n")

  if (!is.na(lyr_crs) & !is.na(coords_crs)) {
    if (coords_crs != lyr_crs) stop("CRS of the provided coordinates and raster do not match")
  }
}

#' Check CRS of coords and layer for krig_gd
#'
#' @param r RasterLayer or SpatRaster
#' @param grd RasterLayer or SpatRaster
#' @param coords sf object, data frame, or matrix representing coordinates
#'
#' @return NULL
#'
#' @noRd
crs_check_krig <- function(r, grd = NULL, coords = NULL) {
  r_crs <- sf::st_crs(r)
  if (is.na(r_crs)) {
    if (!is.null(grd) & is.null(coords)) warning("No CRS found for the provided raster (r). Make sure that r and grd have the same projection.\n")
    if (!is.null(grd) & !is.null(coords)) warning("No CRS found for the provided raster (r). Make sure that r, grd, and coords have the same projection.\n")
    if (is.null(grd) & !is.null(coords)) warning("No CRS found for the provided raster (r). Make sure that r and coords have the same projection.\n")
    if (is.null(grd) & is.null(coords)) warning("No CRS found for the provided raster (r).\n")
  }

  if (!is.null(grd)) {
    grd_crs <- sf::st_crs(grd)

    if (is.na(grd_crs)) warning("No CRS found for the provided raster (grd). Make sure that r and grd have the same projection.\n")

    if (!is.na(r_crs) & !is.na(grd_crs)) {
      if (grd_crs != r_crs) stop("CRS of the two rasters (r and grd) do not match")
    }
  }

  if (!is.null(coords)) {
    coords_crs <- sf::st_crs(coords)

    if (is.na(coords_crs)) warning("No CRS found for the provided coordinates. Make sure the coordinates and the rasters have the same projection.\n")

    if (!is.na(r_crs) & !is.na(coords_crs)) {
      if (coords_crs != r_crs) stop("CRS of the provided coordinates and raster (r) do not match")
    }
  }

  if (!is.null(coords) & !is.null(grd)) {
    if (!is.na(grd_crs) & !is.na(coords_crs)) {
      if (coords_crs != grd_crs) stop("CRS of the provided coordinates and raster (grd) do not match")
    }
  }
}
