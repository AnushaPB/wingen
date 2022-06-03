
#' Krige map
#'
#' @param r raster for kriging
#' @param spdf SpatialPointsDataFrame object to create grid for kriging
#' @param xy whether to co-krige with x and y (~x+y)
#' @param n_cell number of cells to interpolate across
#'
#' @return Raster with kriged values
#' @export
#'
#' @examples
krig_div <- function(r, spdf, xy = TRUE, n_cell = 10000) {

  # convert raster to points
  pa_df <- raster::rasterToPoints(r)

  # create grid
  krig_grid <- spdf_to_grid(spdf, n_cell = n_cell)

  # TODO: FIX PROJECTION STUFF
  krig_df <- coord_proj(pa_df[, c("x", "y")], spdf)

  # Assign values to df
  krig_df$layer <- pa_df[, 3]

  # remove na values
  krig_df <- krig_df[!is.na(krig_df$layer), ]

  # Krige
  if (xy) {
    krig_res <- automap::autoKrige(layer ~ x + y, krig_df, krig_grid)
  }
  if (!xy) {
    krig_res <- automap::autoKrige(layer ~ 1, krig_df, krig_grid)
  }

  # Get kriged spdf
  krig_spdf <- krig_res$krige_output

  # turn SPDF into raster and stack
  krig_raster <- raster::rasterFromXYZ(krig_spdf, crs = raster::crs(krig_grid))

  # plot
  raster::plot(krig_raster, col = viridis::turbo(100), legend = FALSE, box = FALSE, axes = FALSE)
  graphics::lines(spdf)

  return(krig_raster)
}

#' Make grid for kriging
#'
#' @param xrange range of xvalues
#' @param yrange range of yvalues
#' @param len length provided to expand.grid
#'
#' @return
#' @export
#'
#' @examples
range_to_grid <- function(xrange, yrange, len) {
  grd <- expand.grid(
    x = seq(from = xrange[1], to = xrange[2], len = len),
    y = seq(from = yrange[1], to = yrange[2], len = len)
  )
  sp::coordinates(grd) <- ~ x + y
  sp::gridded(grd) <- TRUE
  return(grd)
}

#' Make grid from Spatial Points Data Frame
#'
#' @param spdf Spatial Points Dataframe to create grid for kriging
#' @param n_cell number of grid cells to use when kriging
#'
#' @note code from: https://stackoverflow.com/questions/43436466/create-grid-in-r-for-kriging-in-gstat
#' @return
#' @export
#'
#' @examples
spdf_to_grid <- function(spdf, n_cell = 1000) {
  # make grid from spdf
  grd <- sp::makegrid(spdf, n = n_cell)
  colnames(grd) <- c("x", "y")
  sp::coordinates(grd) <- ~ x + y

  # Next, convert the grid to `SpatialPoints` and subset these points by the polygon.
  grd_pts <- sp::SpatialPoints(
    coords      = grd,
    proj4string = raster::crs(spdf)
  )

  # subset all points in `grd_pts` that fall within `spdf`
  krig_grd <- grd_pts[spdf, ]

  return(krig_grd)
}


#' TEMP COORDS FUNCTION
#'
#' @param coords coordinates
#' @param spdf SpatialPointsDataFrame
#' @param crop_to_spdf whether to crop the coordinates to the spdf
#'
#' @return
#' @export
#'
#' @examples
coord_proj <- function(coords, spdf, crop_to_spdf = FALSE) {
  # make df
  coords_spdf <- data.frame(x = coords[, 1], y = coords[, 2])

  # make into SPFD
  sp::coordinates(coords_spdf) <- ~ x + y

  # Assign CRS
  raster::crs(coords_spdf) <- raster::crs("+proj=longlat +datum=WGS84 +no_defs")

  # Transform to CRS of spdf
  coords_spdf <- sp::spTransform(coords_spdf, raster::crs(spdf))

  # Crop points to SPDF
  if (crop_to_spdf) {
    coords_spdf <- raster::crop(coords_spdf, spdf)
  }

  return(coords_spdf)
}
