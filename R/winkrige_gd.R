#' Krige moving window maps with variogram selection
#'
#' Perform ordinary kriging of the raster(s) produced by \link[wingen]{window_gd} using the gstat package to fit variograms and perform model selection. This function replaces the older \link[wingen]{krig_gd} function to provide more flexibility in variogram model selection. While I have not formally validated the default parameters, they have performed well in practice for kriging `wingen` outputs from both simulated and empirical datasets.
#'
#' The function fits multiple variogram models (Spherical, Exponential, Gaussian, Matern by default) and selects the best fit based on SSErr. It also includes optional weighting to account for sample count variation.
#' 
#' By default, starting values for variogram parameters are set heuristically:
#' - partial sill = ~80% of global variance
#' - nugget = ~20% of global variance
#' - range = 50% of maximum pairwise distance
#'
#' The variogram fitting method defaults to Ordinary Least Squares (\code{fit_method = 6}), which tends to produce more stable fits in noisy or irregular datasets. Weighted methods (1, 2, and 7) may offer improved accuracy in large, well-distributed datasets.
#'
#' For more fine-scale control over variogram fitting and kriging, consult the `gstat` package documentation.
#' 
#' @param r SpatRaster produced by \link[wingen]{window_gd}. Only the first layer is used if multiple layers are present.
#' @param grd Object to create grid for kriging; can be a SpatRaster or RasterLayer. If undefined, \code{r} is used to create a grid.
#' @param weight_r Optional \link[terra]{SpatRaster} with sample counts per cell, used to compute location-specific measurement variance and weights for kriging. If \code{NULL} (default), no weighting is applied.
#' @param nmax Integer. Maximum number of neighboring observations to use for kriging at each prediction location (default: \code{Inf}). Users are encouraged to experiment with \code{nmax} to balance smoothness and local detail; starting with a value of 30 is recommended to reduce computational cost while still capturing local variability.
#' @param maxdist Maximum distance to consider for neighboring observations (default: \code{Inf}). If set together with \code{nmax}, both parameters limit the number of neighbors.
#' @param models Character vector of variogram model names to try (default: \code{c("Sph", "Exp", "Gau", "Mat")}).
#' @param psill_start Optional starting value for partial sill. If \code{NULL} (default), a heuristic value is used (see Note).
#' @param nugget_start Optional starting value for nugget effect. If \code{NULL} (default), a heuristic value is used (see Note).
#' @param range_start Optional starting value for range parameter. If \code{NULL} (default), a heuristic value is used (see Note).
#' @param max_range_frac Numeric. Maximum fraction of the range parameter to consider for neighboring observations (default: 0.5). This can help to limit the influence of distant points.
#' @param fit_method Integer. Variogram fitting method passed to \link[gstat]{fit.variogram}: 1 = weights \(N_j\); 2 = weights \(N_j / \gamma(h_j)^2\); 6 = Ordinary Least Squares (unweighted); 7 = weights \(N_j / h_j^2\). The default (\code{6}) uses OLS, which is generally more robust for small or noisy datasets (see Note).
#' @param model_output Logical. If \code{TRUE}, returns a list with the prediction raster, variogram, and fitted variogram model. If \code{FALSE} (default), returns only the prediction raster.
#'
#' @note
#'
#' **Convergence warnings from `gstat::fit.variogram()`**
#'
#' During variogram fitting, you may see:
#' ```
#' Warning: No convergence after 200 iterations: try different initial values?
#' ```
#' This means the optimizer reached its iteration limit before fully minimizing the error. Even in these cases, `gstat` returns the best-fit model found so far. Warnings often occur with small datasets or noisy empirical variograms. You can experiment with different `psill_start`, `range_start`, and `nugget_start` values, or increase the iteration limit using \code{options(gstat.fit.maxiter)}.
#'
#' @return A \link[terra]{SpatRaster} object of kriged predictions if \code{model_output = FALSE}. If \code{model_output = TRUE}, returns a list with:
#' \describe{
#'   \item{prediction}{Kriged prediction raster (\link[terra]{SpatRaster})}
#'   \item{variogram}{Empirical variogram (\link[gstat]{variogram})}
#'   \item{model}{Best-fit variogram model (\link[gstat]{vgm})}
#' }
#'
#' @details
#' This function uses \link[gstat]{gstat} for variogram fitting and kriging. Weights are computed as the inverse of estimated location-specific variance (\eqn{\sigma^2 / n}) if sample counts are provided.
#'
#' @examples
#' # Note: this toy example uses a very small dataset. 
#' # Warnings may occur due to limited points for variogram fitting.
#' suppressWarnings({
#'   load_mini_ex()
#'   wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, L = 10, rarify = TRUE)
#'   kpi <- krig_gd2(wpi[["pi"]], grd = mini_lyr, nugget = 0)
#'   plot_gd(kpi, main = "Kriged Pi")
#' })
#' @export
winkrig_gd <- function(r, grd = NULL, weight_r = NULL,
                     models = c("Sph", "Exp", "Gau", "Mat"),
                     nmax = Inf, maxdist = Inf,
                     psill_start = NULL, nugget_start = NULL, range_start = NULL, max_range_frac = 0.5, fit_method = 6,
                     model_output = FALSE) {

  # Add check for r being a SpatRaster with one layer
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  if (terra::nlyr(r) > 1) {
    message("Input raster has multiple layers. Using only the first layer.")
    r <- terra::subset(r, 1)
  }

  # Convert raster to sf
  r_pts <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(r_pts)[3] <- "value"
  r_sf <- sf::st_as_sf(r_pts, coords = c("x", "y"), crs = terra::crs(r))

  # Handle weights (convert sample counts to variance estimates)
  # Approximate location-specific measurement variance as σ²_global / n:
  # - σ²_global = variance of all raster values (assumes homoskedasticity)
  # - n = sample count in each cell
  # - Var(mean) = σ² / n because averaging reduces variance by factor n
  # Gstat expects weights = 1/variance to reflect confidence in observations.
  # More samples → lower variance → higher weight.
  if (!is.null(weight_r)) {
    # Convert weight_r to SpatRaster if not already
    if (!inherits(weight_r, "SpatRaster")) weight_r <- terra::rast(weight_r)

    # Extract sample counts
    sample_counts <- terra::extract(weight_r, r_pts[, c("x", "y")])[,2]
    
    # Confirm no sample count values are less than 1 or NA
    stopifnot(all(sample_counts >= 1))
    stopifnot(all(!is.na(sample_counts))) 

    # Estimate global variance
    sigma2_global <- stats::var(r_pts$value, na.rm = TRUE)

    # Estimate location-specific variance
    r_pts$variance <- sigma2_global / sample_counts

    # Convert variance to weights
    r_pts$weight <- 1 / r_pts$variance  
  } else {
    r_pts$weight <- NULL  # No location-specific variance information
  }

  # Fit variogram model
  v <- gstat::variogram(value ~ 1, data = r_sf)

  # Approximate starting values
  if (is.null(psill_start)) psill_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.8
  if (is.null(nugget_start)) nugget_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.2
  if (is.null(range_start)) range_start <- max(v$dist, na.rm = TRUE) * max_range_frac

  # Fit models
  fitted_models <- purrr::map(models, purrr::safely(function(mod) {
    start_model <- gstat::vgm(psill = psill_start, model = mod, range = range_start, nugget = nugget_start)
    fit <- gstat::fit.variogram(v, model = start_model, fit.method = fit_method)
    attr(fit, "model_name") <- mod
    fit
  }))

  # Extract successful results
  fitted_models <- purrr::map(fitted_models, "result")

  # Check if any models were successfully fitted
  fitted_models <- Filter(Negate(is.null), fitted_models)
  if (length(fitted_models) == 0) stop("Variogram fitting failed for all models.")

  # Assign SSErr to each fitted model
  sse <- sapply(fitted_models, function(f) attr(f, "SSErr"))
  best_fit <- fitted_models[[which.min(sse)]]
  best_model <- attr(best_fit, "model_name")

  # Print best model information
  message("Best model:", best_model, "\n")
  message("SSErr:", min(sse), "\n")

  # Use best model for kriging
  if (is.null(weight_r)) {
    krig_model <- gstat::gstat(formula = value ~ 1, locations = r_sf, model = best_fit, nmax = nmax, maxdist = maxdist)
  } else {
    krig_model <- gstat::gstat(formula = value ~ 1, locations = r_sf, model = best_fit, weights = r_pts$weight, nmax = nmax, maxdist = maxdist)
  }

  # Create prediction grid
  if (is.null(grd)) grd <- r
  if (!inherits(grd, "SpatRaster")) grd <- terra::rast(grd)
  grd_df <- terra::as.data.frame(grd, xy = TRUE, na.rm = FALSE)
  grd_sf <- sf::st_as_sf(grd_df, coords = c("x", "y"), crs = terra::crs(r))
  krig_pred <- predict(krig_model, newdata = grd_sf)

  # Convert predictions to rasters
  r_pred <- terra::rasterize(krig_pred, grd, crs = terra::crs(grd), field = "var1.pred")
  names(r_pred) <- names(r)
 
  if (model_output) {
    # If model output is requested, return the prediction and variogram
    return(list(raster = r_pred, variogram = v, model = best_fit))
  } else {
    # If only prediction is requested, return the raster
    return(r_pred)
  }
}
