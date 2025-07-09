
# gstat vignette (Pebesma, 2004) uses nmax = 20.
# Diggle & Ribeiro (2007), Model-based Geostatistics, discuss ~30 neighbors as typical.
# ArcGIS kriging defaults to 12–32 neighbors.
krig_gd2 <- function(r, grd = NULL, weight_r = NULL,
                     candidate_models = c("Sph", "Exp", "Gau", "Mat"),
                     max_range_frac = 0.5, nmax = 30, verbose = TRUE,
                     psill_start = NULL, nugget_start = NULL, range_start = NULL,
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
  formula_str <- value ~ 1
  v <- gstat::variogram(formula_str, data = r_sf)

  # Approximate starting values
  if (is.null(psill_start)) psill_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.8
  if (is.null(nugget_start)) nugget_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.2
  if (is.null(range_start)) range_start <- max(v$dist, na.rm = TRUE) * max_range_frac

  # Fit multiple models and pick the best
  fitted_models <- purrr::map(candidate_models, purrr::safely(function(mod) {
    start_model <- gstat::vgm(psill = psill_start, model = mod, range = range_start, nugget = nugget_start)
    fit <- gstat::fit.variogram(v, model = start_model, fit.method = 6)
    attr(fit, "model_name") <- mod
    return(fit)
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

  if (verbose) {
    cat("Best model:", best_model, "\n")
    cat("SSErr:", min(sse), "\n")
  }

  # Krige results
  if (is.null(weight_r)) {
    krig_model <- gstat::gstat(formula = value ~ 1, locations = r_sf, model = best_fit, nmax = nmax)
  } else {
    krig_model <- gstat::gstat(formula = value ~ 1, locations = r_sf, model = best_fit, weights = r_pts$weight, nmax = nmax)
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
    return(list(prediction = r_pred, variogram = v, model = best_fit))
  } else {
    # If only prediction is requested, return the raster
    return(r_pred)
  }
}


