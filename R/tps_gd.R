#' Thin-Plate Spline Interpolation of Windowed Diversity
#'
#' Perform interpolation of the raster(s) produced by \link[wingen]{window_gd} using \link[fields]{Tps} (thin-plate spline).
#'
#' @param r SpatRaster produced by \link[wingen]{window_gd}
#' @param grd object to create grid for interpolating; can be a SpatRaster or RasterLayer. If undefined, will use \code{r} to create a grid.
#' @param index integer index of the layer in the raster stack to interpolate (defaults to 1; i.e., the first layer)
#' @param coords if provided, kriging will occur based only on values at these coordinates. Can be provided as an sf points, a two-column matrix, or a data.frame representing x and y coordinates
#' @param agg_grd factor to use for aggregation of `grd`, if provided (this will decrease the resolution of the final interpolated raster; defaults to NULL)
#' @param disagg_grd factor to use for disaggregation of `grd`, if provided (this will increase the resolution of the final interpolated raster; defaults to NULL)
#' @param agg_r factor to use for aggregation of `r`, if provided (this will decrease the number of points used in the kriging model; defaults to NULL)
#' @param disagg_r factor to use for disaggregation, of `r` if provided (this will increase the number of points used in the kriging model; defaults to NULL)
#' @param m Numeric; order of derivative for the spline (default = 2 for thin-plate spline).
#' @param scale.type Character; how to scale x/y coordinates before fitting (passed to fields::Tps; default = "range").
#' @param lambda Optional numeric; smoothing parameter. If NULL (default), chosen by GCV.
#' @param lower_bound if TRUE (default), converts all values in the interpolated raster less than the minimum value of the input raster, to that minimum.
#' @param upper_bound if TRUE (default), converts all values in the interpolated raster greater than the maximum value of the input raster, to that maximum.
#' @param resample whether to resample `grd` or `r`. Set to `"r"` to resample `r` to `grd`. Set to `"grd"` to resample `grd` to `r` (defaults to FALSE for no resampling)
#' @param resample_first if aggregation or disaggregation is used in addition to resampling, specifies whether to resample before (resample_first = TRUE) or after (resample_first = FALSE) aggregation/disaggregation (defaults to TRUE)
#'
#' @return A SpatRaster of the interpolated surface, with the same CRS as `r`.
#'
#' @examples
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, L = 10, rarify = TRUE)
#' tpspi <- tps_gd(wpi, mini_lyr)
#' plot_gd(tpspi, main = "TPS Pi")
#'
#' @export
tps_gd <- function(r, grd = NULL, index = 1, coords = NULL,
                   agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL,
                   lower_bound = TRUE, upper_bound = TRUE,
                   m = 2, scale.type = "range", lambda = NULL,
                   resample = FALSE, resample_first = TRUE) {
  # Ensure inputs are SpatRasters
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  if (!is.null(grd) && inherits(grd, "RasterLayer")) grd <- terra::rast(grd)

  # Subset desired layer
  if (terra::nlyr(r) > 1) r <- r[[index]]

  # Apply aggregation/disaggregation/resampling
  stk <- raster_transform(
    r = r,
    grd = grd,
    agg_grd = agg_grd,
    disagg_grd = disagg_grd,
    agg_r = agg_r,
    disagg_r = disagg_r,
    resample = resample,
    resample_first = resample_first
  )
  r_t <- stk[[names(r)]]
  grd_t <- stk[["grd"]]

  # Build point data.frame for fitting
  if (!is.null(coords)) {
    pts_df <- terra::as.data.frame(
      terra::extract(r_t, coords, ID = FALSE, xy = TRUE),
      xy = TRUE
    )
    names(pts_df)[3] <- "layer"
  } else {
    pts_df <- terra::as.data.frame(r_t, xy = TRUE, na.rm = TRUE)
    names(pts_df)[3] <- "layer"
  }

  # Fit thin-plate spline
  tps_mod <- fields::Tps(
    x          = as.matrix(pts_df[, c("x", "y")]),
    Y          = pts_df$layer,
    m          = m,
    scale.type = scale.type,
    lambda     = lambda
  )

  # Build prediction grid
  grd_df <- terra::as.data.frame(grd_t, xy = TRUE, na.rm = FALSE)

  # Predict values
  pred_vals <- fields::predict.Krig(
    object = tps_mod,
    x      = as.matrix(grd_df[, c("x", "y")])
  )
  grd_df$layer <- pred_vals

  # Apply bounds
  if (is.numeric(lower_bound)) {
    grd_df$layer[grd_df$layer < lower_bound] <- lower_bound
  }
  if (is.numeric(upper_bound)) {
    grd_df$layer[grd_df$layer > upper_bound] <- upper_bound
  }
  if (isTRUE(lower_bound)) {
    mn <- min(pts_df$layer, na.rm = TRUE)
    grd_df$layer[grd_df$layer < mn] <- mn
  }
  if (isTRUE(upper_bound)) {
    mx <- max(pts_df$layer, na.rm = TRUE)
    grd_df$layer[grd_df$layer > mx] <- mx
  }

  # Convert back to SpatRaster
  out <- terra::rast(grd_df[, c("x", "y", "layer")], crs = terra::crs(r_t))
  names(out) <- names(r_t)

  return(out)
}
