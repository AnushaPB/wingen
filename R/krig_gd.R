#' Krige moving window maps
#'
#' Perform interpolation of the raster(s) produced by \link[wingen]{window_gd} using \link[automap]{autoKrige}
#'
#' @param r SpatRaster produced by \link[wingen]{window_gd}
#' @param index integer indices of layers in raster stack to krige (defaults to 1; i.e., the first layer)
#' @param grd object to create grid for kriging; can be a SpatRaster or RasterLayer. If undefined, will use \code{r} to create a grid.
#' @param coords if provided, kriging will occur based only on values at these coordinates. Can be provided as an sf points, a two-column matrix, or a data.frame representing x and y coordinates
#' @param agg_grd factor to use for aggregation of `grd`, if provided (this will decrease the resolution of the final kriged raster; defaults to NULL)
#' @param disagg_grd factor to use for disaggregation of `grd`, if provided (this will increase the resolution of the final kriged raster; defaults to NULL)
#' @param agg_r factor to use for aggregation of `r`, if provided (this will decrease the number of points used in the kriging model; defaults to NULL)
#' @param disagg_r factor to use for disaggregation, of `r` if provided (this will increase the number of points used in the kriging model; defaults to NULL)
#' @param autoKrige_output whether to return full output from \link[automap]{autoKrige} including uncertainty rasters (defaults to FALSE). If TRUE, returns a list with the kriged input raster layer ("raster"), kriged variance ("var"), kriged standard deviation ("stdev"), and full autoKrige output ("autoKrige_output").
#' @param lower_bound if TRUE (default), converts all values in the kriged raster less than the minimum value of the input raster, to that minimum
#' @param upper_bound if TRUE (default), converts all values in the kriged raster greater than the maximum value of the input raster, to that maximum
#' @param krig_method method to use for kriging. If `ordinary`, ordinary/simple kriging is performed (formula: ~ 1; default). If `universal`,  universal kriging is performed (formula = ~ x + y).
#' @param resample whether to resample `grd` or `r`. Set to `"r"` to resample `r` to `grd`. Set to `"grd"` to resample `grd` to `r` (defaults to FALSE for no resampling)
#' @param resample_first if aggregation or disaggregation is used in addition to resampling, specifies whether to resample before (resample_first = TRUE) or after (resample_first = FALSE) aggregation/disaggregation (defaults to TRUE)
#'
#' @return a SpatRaster object or a list of \link[automap]{autoKrige} outputs (if autoKrige_output = TRUE)
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
                    lower_bound = TRUE, upper_bound = TRUE,
                    krig_method = "ordinary",
                    resample = FALSE, resample_first = TRUE) {
  # check CRS
  crs_check_krig(r = r, grd = grd, coords = coords)

  # Make sure grid and raster layer are SpatRasters
  if (!inherits(grd, "SpatRaster") & !is.null(grd) & inherits(grd, "RasterLayer")) grd <- terra::rast(grd)
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)

  # subset desired layers
  if (terra::nlyr(r) > 1) {
    r <- r[[index]]
  }

  # convert from stack to list
  rls <- terra::as.list(r)

  if (is.null(grd)) {
    grd <- r[[1]]
    warning("No grd provided, defaults to using first raster layer to create grd")
  }

  # Temporary: warning if non-spatRaster grd is used
  if (!inherits(grd, "SpatRaster")) warning("/n Starting October 2023, the krig_gd() grd argument will only accept SpatRaster and Raster objects")

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
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    krig_method = krig_method,
    resample = resample,
    resample_first = resample_first
  )

  # give names from stack to list
  names(rstk) <- names(r)

  if (length(rstk) == 1) {
    rstk <- rstk[[1]]
  }

  if (!autoKrige_output) {
    # convert to stack if a list
    if (!inherits(rstk, "SpatRaster")) rstk <- terra::rast(rstk)
  }

  return(rstk)
}

#' Krige SpatRaster
#'
#' Helper function for \code{\link{krig_gd}}
#'
#' @inheritParams krig_gd
#'
#' @noRd
krig_gd_lyr <- function(r, grd = NULL, coords = NULL,
                        agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL,
                        autoKrige_output = FALSE,
                        lower_bound = TRUE,
                        upper_bound = TRUE,
                        krig_method = "ordinary",
                        resample = FALSE, resample_first = TRUE) {
  # Transform raster layer
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

  r <- stk[[names(r)]]
  grd <- stk[["grd"]]

  # create df
  krig_df <- make_krige_df(r, coords)

  # create grid
  krig_grid <- make_krige_grid(grd)

  # krige using autoKrige
  krig_r <- krig(
    krig_df = krig_df,
    krig_grid = krig_grid,
    grd = grd,
    autoKrige_output = autoKrige_output,
    krig_method = krig_method,
    lower_bound = lower_bound,
    upper_bound = upper_bound
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
krig <- function(krig_df, krig_grid, grd, autoKrige_output = FALSE, krig_method = "ordinary", lower_bound = TRUE, upper_bound = TRUE) {
  # autokrige
  if (krig_method == "ordinary") {
    krig_res <- automap::autoKrige(layer ~ 1, input_data = krig_df, new_data = krig_grid)
  } else if (krig_method == "universal") {
    krig_res <- automap::autoKrige(layer ~ x + y, input_data = krig_df, new_data = krig_grid)
  } else {
    stop("invalid krig_method specified")
  }

  # Convert results to vector
  krig_vect <- terra::vect(krig_res$krige_output)

  # turn results into raster
  krig_r <- terra::rasterize(krig_vect, grd, field = names(krig_vect))

  # rename autoKrige output
  names(krig_r) <- c("pred", "var", "stdev")

  # perform bounding on pred layer
  if (is.numeric(lower_bound)) krig_r$pred[krig_r$pred < lower_bound] <- lower_bound
  if (is.numeric(upper_bound)) krig_r$pred[krig_r$pred > upper_bound] <- upper_bound

  if (is.logical(lower_bound)) {
    if (lower_bound) krig_r$pred[krig_r$pred < min(krig_df$layer, na.rm = TRUE)] <- min(krig_df$layer, na.rm = TRUE)
  }

  if (is.logical(upper_bound)) {
    if (upper_bound) krig_r$pred[krig_r$pred > max(krig_df$layer, na.rm = TRUE)] <- max(krig_df$layer, na.rm = TRUE)
  }

  # create results
  if (autoKrige_output) {
    return(list(raster = krig_r[["pred"]], autoKrige_output = krig_res))
  } else {
    return(krig_r[["pred"]])
  }
}

#' Create df for kriging
#'
#' @inheritParams krig_gd
#'
#' @noRd
make_krige_df <- function(r, coords = NULL) {
  # use coords if provided
  if (!is.null(coords)) {
    if (inherits(coords, "sf")) coords <- terra::vect(coords)
    krig_df <- terra::extract(r, coords, ID = FALSE, xy = TRUE)
  } else {
    # convert raster to df
    krig_df <- terra::as.data.frame(r, xy = TRUE)
  }

  # convert to sf
  krig_sf <- make_krige_sf(krig_df)

  # reassign crs
  sf::st_crs(krig_sf) <- sf::st_crs(r)

  # remove na values
  krig_sf <- krig_sf[!is.na(krig_sf$layer), ]

  return(krig_sf)
}

#' Create grid for kriging
#'
#' @inheritParams krig_gd
#'
#' @noRd
make_krige_grid <- function(grd = NULL) {
  # make sf
  krig_sf <- make_krige_sf(grd)

  # reassign crs
  sf::st_crs(krig_sf) <- sf::st_crs(grd)

  return(krig_sf)
}

#' Convert a raster or data.frame of coordinates to a sf object for kriging
#'
#' @param x SpatRaster
#'
#' @return sf points
#'
#' @noRd
make_krige_sf <- function(x) {
  # convert from raster to coords if not already a data.frame
  if (inherits(x, "SpatRaster")) x <- terra::as.data.frame(x, xy = TRUE, na.rm = FALSE)

  # convert to sf
  krig_sf <- sf::st_as_sf(x, coords = c("x", "y"))

  # rename layer
  names(krig_sf)[1] <- "layer"

  # add columns for x and y
  # Note: this is necessary for universal kriging with the formula "layer ~ x + y"
  krig_sf$x <- x$x
  krig_sf$y <- x$y

  return(krig_sf)
}

#' Transform raster
#'
#' @inheritParams krig_gd
#'
#' @noRd
raster_transform <- function(r, grd, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE) {
  if (terra::nlyr(r) > 1) stop(">1 layer provided for r")
  if (terra::nlyr(grd) > 1) stop(">1 layer provided for grd")

  if (resample_first) {
    if (resample == "r") r <- terra::resample(r, grd)
    if (resample == "grd") grd <- terra::resample(grd, r)
  }

  if (!is.null(agg_grd) & !is.null(disagg_grd)) stop("Both agg_grd and disagg_grd provided, when only one should be provided")
  if (!is.null(agg_grd)) grd <- terra::aggregate(grd, agg_grd)
  if (!is.null(disagg_grd)) grd <- terra::disagg(grd, disagg_grd)

  if (!is.null(agg_r) & !is.null(disagg_r)) stop("Both agg_r and disagg_r provided, when only one should be provided")
  if (!is.null(agg_r)) r <- terra::aggregate(r, agg_r)
  if (!is.null(disagg_r)) r <- terra::disagg(r, disagg_r)

  if (!resample_first) {
    if (resample == "r") r <- terra::resample(r, grd)
    if (resample == "grd") grd <- terra::resample(grd, r)
  }

  s <- list(r, grd)
  names(s) <- c(names(r), "grd")

  return(s)
}
