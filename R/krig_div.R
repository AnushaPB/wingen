
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
krig_div <- function(r, spdf, xy = TRUE, n_cell = 10000){

  # convert raster to points
  pa_df <- rasterToPoints(r)

  # create grid
  krig_grid <- spdf_to_grid(spdf, n_cell = n_cell)

  # TODO: FIX PROJECTION STUFF
  krig_df <- coord_proj(pa_df[,c("x","y")], spdf)

  # Assign values to df
  krig_df$layer <- pa_df[,3]

  # remove na values
  krig_df <- krig_df[!is.na(krig_df$layer),]

  # Krige
  if(xy){krig_res <- autoKrige(layer ~ x + y, krig_df, krig_grid)}
  if(!xy){krig_res <- autoKrige(layer ~ 1, krig_df, krig_grid)}

  # Get kriged spdf
  krig_spdf <- krig_res$krige_output

  # turn SPDF into raster and stack
  krig_raster <- raster::rasterFromXYZ(krig_spdf, crs = raster::crs(krig_grid))

  # plot
  raster::plot(krig_raster, col = turbo(100), legend = FALSE, box = FALSE, axes = FALSE)
  lines(spdf)

  return(krig_raster)
}


