#' GAM-Based Interpolation of Windowed Diversity
#'
#' Perform generalized additive model based interpolation of the raster(s) produced by \link[wingen]{window_gd} using \link[mgcv]{gam}.
#'
#' @param r SpatRaster produced by \link[wingen]{window_gd}
#' @param grd object to create grid for interpolation; can be a SpatRaster or RasterLayer. If undefined, will use \code{r} to create a grid.
#' @param index integer index of the layer in the raster stack to interpolate (defaults to 1; i.e., the first layer)
#' @param coords if provided, interpolation will occur based only on values at these coordinates. Can be provided as an sf points, a two-column matrix, or a data.frame representing x and y coordinates
#' @param agg_grd factor to use for aggregation of `grd`, if provided (this will decrease the resolution of the final interpolated raster; defaults to NULL)
#' @param disagg_grd factor to use for disaggregation of `grd`, if provided (this will increase the resolution of the final interpolated raster; defaults to NULL)
#' @param agg_r factor to use for aggregation of `r`, if provided (this will decrease the number of points used in the interpolation model; defaults to NULL)
#' @param disagg_r factor to use for disaggregation, of `r` if provided (this will increase the number of points used in the interpolation model; defaults to NULL)
#' @param lower_bound if TRUE (default), converts all values in the interpolated raster less than the minimum value of the input raster, to that minimum.
#' @param upper_bound if TRUE (default), converts all values in the interpolated raster greater than the maximum value of the input raster, to that maximum.
#' @param resample whether to resample `grd` or `r`. Set to `"r"` to resample `r` to `grd`. Set to `"grd"` to resample `grd` to `r` (defaults to FALSE for no resampling)
#' @param resample_first if aggregation or disaggregation is used in addition to resampling, specifies whether to resample before (resample_first = TRUE) or after (resample_first = FALSE) aggregation/disaggregation (defaults to TRUE)
#' @param bs Character; the basis type for the spatial smoother (passed to \code{mgcv::s()}, default = "tp").
#' @param k Integer; the basis dimension (maximum wiggliness + 1) for the smoother. Defaults to the closest integar value of the square root of the number of values in `r` (or `coords`, if provided). If the number of values is less than 10, defaults to k = 3.
#' @param method Character; smoothing-parameter estimation method for \code{mgcv::gam()} (default = "REML").
#' @param gamma Numeric; inflation factor for the smoothing selection criterion (default = 1).
#' @param select Logical; whether to apply shrinkage to smooth terms (default = FALSE).
#'
#' @return A SpatRaster of the interpolated surface, with the same CRS as `r`.
#'
#' @examples
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, L = 10, rarify = TRUE)
#' gampi <- gam_gd(wpi, mini_lyr)
#' plot_gd(gampi, main = "GAM Pi")
#'
#' @export
gam_gd <- function(r, grd = NULL, index = 1, coords = NULL,
                   agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL,
                   lower_bound = TRUE, upper_bound = TRUE,
                   resample = FALSE, resample_first = TRUE,
                   bs = "tp", k = NULL, method = "REML", gamma = 1, select = FALSE) {

  # Ensure terra objects
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  if (!is.null(grd) && inherits(grd, "RasterLayer")) grd <- terra::rast(grd)

  # Subset layer
  if (terra::nlyr(r) > 1) r <- r[[index]]

  # Transform rasters (agg/disagg/resample)
  stk <- raster_transform(
    r               = r,
    grd             = grd,
    agg_grd         = agg_grd,
    disagg_grd      = disagg_grd,
    agg_r           = agg_r,
    disagg_r        = disagg_r,
    resample        = resample,
    resample_first  = resample_first
  )
  r_t <- stk[[names(r)]]
  grd_t <- stk[["grd"]]

  # Build point data.frame
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

  # Set K if not provided
  if (is.null(k)) {
    n_obs <- nrow(pts_df)
    if (n_obs < 10) {
      k <- 3  # Minimum basis dimension
    } else {
      k <- floor(sqrt(n_obs))  # Default based on number of observations
    }
    message("Setting k = ", k)
  }

  # Fit GAM
  gam_mod <- mgcv::gam(
    formula = layer ~ s(x, y, bs = bs, k = k),
    data    = pts_df,
    method  = method,
    gamma   = gamma,
    select  = select
  )
  
  # Print gam check results
  message("GAM model check results:")
  gam.check(gam_mod)

  # Prediction grid
  grd_df <- terra::as.data.frame(grd_t, xy = TRUE, na.rm = FALSE)

  # Predict values
  pred_vals <- stats::predict(gam_mod, newdata = grd_df)
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

  # Reconstruct SpatRaster
  out <- terra::rast(grd_df[, c("x", "y", "layer")],
    type = "xyz",
    crs  = terra::crs(r_t)
  )
  names(out) <- names(r_t)

  return(out)
}
