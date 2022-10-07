
#' Krige moving window maps
#'
#' Perform interpolation of the raster(s) produced by \link[wingen]{window_gd} using \link[automap]{autoKrige}
#'
#' @param r RasterLayer or RasterStack produced by \link[wingen]{window_gd}
#' @param index integer indices of layers in raster stack to krige (defaults to 1; i.e., the first layer)
#' @param grd object to create grid for kriging; can be RasterLayer, SpatialPointsDataFrame, or a gridded object as defined by 'sp'. If undefined, will use \code{r} to create a grid.
#' @param coords if provided, kriging will occur based only on values at these coordinates
#' @param agg_grd factor to use for aggregation of `grd`, if provided (this will decrease the resolution of the final kriged raster; defaults to NULL)
#' @param disagg_grd factor to use for disaggregation of `grd`, if provided (this will increase the resolution of the final kriged raster; defaults to NULL)
#' @param agg_r factor to use for aggregation of `r`, if provided (this will decrease the number of points used in the kriging model; defaults to NULL)
#' @param disagg_r factor to use for disaggregation, of `r` if provided (this will increase the number of points used in the kriging model; defaults to NULL)
#' @param autoKrige_output whether to return full output from \link[automap]{autoKrige} including uncertainty rasters (defaults to FALSE). If TRUE, returns a list with the kriged input raster layer ("raster"), kriged variance ("var"), kriged standard deviation ("stdev"), and full autoKrige output ("autoKrige_output").
#' @param zero_correction if TRUE (default), converts all values in the kriged raster less than zero, to zero (since genetic diversity and sample count values can't be negative)
#' @param krig_method method to use for kriging. If `ordinary`, ordinary/simple kriging is performed (formula: ~ 1). If `universal`,  universal kriging is performed (default, formula = ~ x + y).
#' @param resample whether to resample `grd` or `r`. Set to `"r"` to resample `r` to `grd`. Set to `"grd"` to resample `grd` to `r` (defaults to FALSE for no resampling)
#' @param resample_first if aggregation or disaggregation is used in addition to resampling, specifies whether to resample before (resample_first = TRUE) or after (resample_first = FALSE) aggregation/disaggregation (defaults to TRUE)
#'
#' @return a Raster* object or a list of \link[automap]{autoKrige} outputs (if autoKrige_output = TRUE)
#' @export
#'
#' @examples
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, L = 10, rarify = TRUE)
#' kpi <- krig_gd(wpi, mini_lyr)
#' plot_gd(kpi, main = "Kriged Pi")
#'
krig_gd <- function(r, grd = NULL, index = 1, coords = NULL,
                    agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL,
                    autoKrige_output = FALSE,
                    zero_correction = TRUE,
                    krig_method = "universal",
                    resample = FALSE, resample_first = TRUE) {

  # subset desired layers
  if (raster::nlayers(r) > 1) {
    r <- r[[index]]
  }

  # convert from stack to list
  rls <- raster::as.list(r)

  if (is.null(grd)) {
    grd <- r[[1]]
    warning("no grd provided, defaults to using first raster layer to create grd")
  }

  # krige
  rstk <- purrr::map(rls,
    krig_gd_lyr,
    grd = grd,
    coords = coords,
    agg_grd = agg_grd,
    disagg_grd = disagg_grd,
    agg_r = agg_r,
    disagg_r = disagg_r,
    autoKrige_output = autoKrige_output,
    zero_correction = zero_correction,
    krig_method = krig_method,
    resample = resample,
    resample_first = resample_first
  )

  # give names from stack to list
  names(rstk) <- names(r)

  if (length(rls) == 1) {
    rstk <- rstk[[1]]
  }

  if (!autoKrige_output) {
    # convert from list to stack
    rstk <- raster::stack(rstk)
  }

  return(rstk)
}

#' Krige RasterLayer
#'
#' Helper function for \code{\link{krig_gd}}
#'
#' @inheritParams krig_gd
#'
#' @noRd
krig_gd_lyr <- function(r, grd = NULL, coords = NULL,
                        agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL,
                        autoKrige_output = FALSE,
                        zero_correction = TRUE,
                        krig_method = "ordinary",
                        resample = FALSE, resample_first = TRUE) {

  # Transform raster layer
  if (inherits(grd, "RasterLayer")) {
    stk <- raster_transform(
      r = r, grd = grd,
      agg_grd = agg_grd, disagg_grd = disagg_grd, agg_r = agg_r, disagg_r = disagg_r,
      resample = resample, resample_first = resample_first
    )
    r <- stk[[names(r)]]
    grd <- stk[["grd"]]
  }

  # create df
  krig_df <- make_krig_df(r, coords)

  # create grid
  krig_grid <- make_krige_grid(r, grd)

  # remove crs values (automap doesn't like latlon CRS)
  if (!raster::compareCRS(krig_df, krig_grid)) warning("the provided raster and grid have different crs")
  raster::crs(krig_df) <- NA
  raster::crs(krig_grid) <- NA

  # krige using autoKrige
  krig_r <- krig(krig_df, krig_grid,
    autoKrige_output = autoKrige_output,
    krig_method = krig_method, zero_correction = zero_correction
  )

  return(krig_r)
}


#' Perform kriging with autoKrige
#'
#' @param krig_df dataframe for kriging
#' @param krig_grid grid for kriging
#' @inheritParams krig_gd
#'
#' @noRd
krig <- function(krig_df, krig_grid, autoKrige_output = FALSE, krig_method = "ordinary", zero_correction = TRUE) {
  # autokrige
  if (krig_method == "ordinary") {
    krig_res <- automap::autoKrige(layer ~ 1, input_data = krig_df, new_data = krig_grid)
  } else if (krig_method == "universal") {
    krig_res <- automap::autoKrige(layer ~ x + y, input_data = krig_df, new_data = krig_grid)
  } else stop("invalid krig_method specified")

  # Get kriged spdf
  krig_spdf <- krig_res$krige_output

  # turn spdf into raster (automatically just uses the first variable)
  krig_r <- raster::rasterFromXYZ(krig_spdf, crs = raster::crs(krig_grid))

  # replace negative values with zero
  if (zero_correction) krig_r[krig_r < 0] <- 0

  # create results
  if (autoKrige_output) {
    krig_var <- raster::rasterFromXYZ(krig_spdf[, "var1.var"], crs = raster::crs(krig_grid))
    krig_stdev <- raster::rasterFromXYZ(krig_spdf[, "var1.stdev"], crs = raster::crs(krig_grid))
    result <- list(raster = krig_r, var = krig_var, stdev = krig_stdev, autoKrige_output = krig_res)
  } else {
    result <- krig_r
  }

  return(result)
}

#' Create df for kriging
#'
#' @inheritParams krig_gd
#'
#' @noRd
make_krig_df <- function(r, coords = NULL) {

  # convert raster to df
  krig_df <- data.frame(raster::rasterToPoints(r))

  # use coords if provided
  if (!is.null(coords)) {
    coords <- data.frame(coords)
    rex <- raster::extract(r, coords)
    krig_df <- data.frame(coords, layer = rex)
  }

  # Define colnames
  colnames(krig_df) <- c("x", "y", "layer")

  # convert to spdf
  sp::coordinates(krig_df) <- ~ x + y

  # remove na values
  krig_df <- krig_df[!is.na(krig_df$layer), ]

  return(krig_df)
}

#' Create grid for kriging
#'
#' @inheritParams krig_gd
#'
#' @noRd
make_krige_grid <- function(r = NULL, grd = NULL) {
  if (is.null(grd)) {
    krig_grid <- raster_to_grid(r)
  } else if (inherits(grd, "RasterLayer") | inherits(grd, "RasterStack")) {
    krig_grid <- raster_to_grid(grd)
  } else if (sp::gridded(grd)) {
    krig_grid <- grd
  } else {
    stop(" unable to find an inherited method for type of grd provided")
  }
  return(krig_grid)
}

#' Convert a raster to a grid
#'
#' @param x RasterLayer
#'
#' @return gridded SpatialPixelsDataFrame
#'
#' @noRd
raster_to_grid <- function(x) {
  grd <- data.frame(raster::rasterToPoints(x))
  sp::coordinates(grd) <- ~ x + y
  sp::gridded(grd) <- TRUE
  return(grd)
}

#' Transform raster
#'
#' @inheritParams krig_gd
#'
#' @noRd
raster_transform <- function(r, grd, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE) {
  if (raster::nlayers(r) > 1) stop(">1 layer provided for r")
  if (raster::nlayers(grd) > 1) stop(">1 layer provided for grd")

  if (resample_first) {
    if (resample == "r") r <- raster::resample(r, grd)
    if (resample == "grd") grd <- raster::resample(grd, r)
  }

  if (!is.null(agg_grd) & !is.null(disagg_grd)) stop("Both agg_grd and disagg_grd provided, when only one should be provided")
  if (!is.null(agg_grd)) grd <- raster::aggregate(grd, agg_grd)
  if (!is.null(disagg_grd)) grd <- raster::disaggregate(grd, disagg_grd)

  if (!is.null(agg_r) & !is.null(disagg_r)) stop("Both agg_r and disagg_r provided, when only one should be provided")
  if (!is.null(agg_r)) r <- raster::aggregate(r, agg_r)
  if (!is.null(disagg_r)) r <- raster::disaggregate(r, disagg_r)

  if (!resample_first) {
    if (resample == "r") r <- raster::resample(r, grd)
    if (resample == "grd") grd <- raster::resample(grd, r)
  }

  s <- list(r, grd)
  names(s) <- c(names(r), "grd")

  return(s)
}

