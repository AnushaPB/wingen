
#' Krige moving window maps
#'
#' Perform interpolation of the raster(s) produced by \link[wingen]{window_gd} using \link[automap]{autoKrige}
#'
#' @param r SpatRaster produced by \link[wingen]{window_gd}
#' @param index integer indices of layers in raster stack to krige (defaults to 1; i.e., the first layer)
#' @param grd object to create grid for kriging; can be a SpatRaster, RasterLayer, SpatialPointsDataFrame, or a gridded object as defined by 'sp'. If undefined, will use \code{r} to create a grid.
#' @param coords if provided, kriging will occur based only on values at these coordinates. Can be provided as an sf object, a two-column matrix, or a data.frame representing x and y coordinates
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
    # convert to stack
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
  if (inherits(grd, "SpatRaster")) {
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

  # krige using autoKrige
  krig_r <- krig(krig_df, krig_grid,
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
krig <- function(krig_df, krig_grid, autoKrige_output = FALSE, krig_method = "ordinary", lower_bound = TRUE, upper_bound = TRUE) {
  # autokrige
  if (krig_method == "ordinary") {
    krig_res <- automap::autoKrige(layer ~ 1, input_data = krig_df, new_data = krig_grid)
  } else if (krig_method == "universal") {
    krig_res <- automap::autoKrige(layer ~ x + y, input_data = krig_df, new_data = krig_grid)
  } else {
    stop("invalid krig_method specified")
  }

  # Get kriged spdf
  krig_spdf <- krig_res$krige_output

  # turn spdf into raster (automatically just uses the first variable)
  krig_r <- terra::rast(krig_spdf, type = "xyz", crs = terra::crs(krig_grid))

  # perform bounding
  if (is.numeric(lower_bound)) krig_r[krig_r < lower_bound] <- lower_bound
  if (is.numeric(upper_bound)) krig_r[krig_r > upper_bound] <- upper_bound

  if (is.logical(lower_bound)) {
    if (lower_bound) krig_r[krig_r < min(krig_df$layer, na.rm = TRUE)] <- min(krig_df$layer, na.rm = TRUE)
  }

  if (is.logical(upper_bound)) {
    if (upper_bound) krig_r[krig_r > max(krig_df$layer, na.rm = TRUE)] <- max(krig_df$layer, na.rm = TRUE)
  }

  # rename autoKrige output
  names(krig_r) <- c("pred", "var", "stdev")

  # create results
  if (autoKrige_output) {
    return(list(raster = krig_r, autoKrige_output = krig_res))
  } else {
    return(krig_r[[1]])
  }

  return(result)
}

#' Create df for kriging
#'
#' @inheritParams krig_gd
#'
#' @noRd
make_krig_df <- function(r, coords = NULL) {
  # use coords if provided
  if (!is.null(coords)) {
    if (inherits(coords, "sf")) coords <- terra::vect(coords)
    krig_df <- terra::extract(r, coords, ID = FALSE, xy = TRUE)
  } else {
    # convert raster to df
    krig_df <- terra::as.data.frame(r, xy = TRUE)
  }

  # convert to sf
  krig_sf <- sf::st_as_sf(krig_df, coords = c("x", "y"))

  # reassign crs
  sf::st_crs(krig_sf) <- sf::st_crs(r)

  # convert to sp
  krig_sp <- sf::as_Spatial(krig_sf)

  # edit names
  names(krig_sp) <- "layer"
  colnames(krig_sp@coords) <- c("x", "y")

  # remove na values
  krig_sp <- krig_sp[!is.na(krig_sp$layer), ]

  return(krig_sp)
}

#' Create grid for kriging
#'
#' @inheritParams krig_gd
#'
#' @noRd
make_krige_grid <- function(r = NULL, grd = NULL) {
  if (is.null(grd)) {
    krig_grid <- raster_to_grid(r)
  } else if (inherits(grd, "SpatRaster")) {
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
#' @param x SpatRaster
#'
#' @return gridded SpatialPixelsDataFrame
#'
#' @noRd
raster_to_grid <- function(x) {
  grd <- terra::as.data.frame(x, xy = TRUE)
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
