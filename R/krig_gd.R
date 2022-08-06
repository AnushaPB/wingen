
#' Raster interpolation using 'autoKrige'
#'
#' @param r RasterLayer or RasterStack
#' @param index integer indices of layers in raster stack to krige (defaults to 1, i.e. the first layer)
#' @inheritParams krig_gd_lyr
#' @return RasterLayer or RasterStack
#' @export
#'
#' @examples
#' library("raster")
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, L = 10, rarify = TRUE)
#' kpi <- krig_gd(wpi, mini_lyr)
#' plot_gd(kpi, main = "Kriged Pi")
#'
krig_gd <- function(r, grd = NULL, index = 1, coords = NULL, xy = FALSE, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE, n_cell = 10000) {

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
  rstk <- purrr::map(rls, krig_gd_lyr, grd, coords, xy, resample, agg_grd, disagg_grd, agg_r, disagg_r, n_cell)
  # convert from list to stack
  rstk <- raster::stack(rstk)
  names(rstk) <- names(r)

  return(rstk)
}

#' Krige RasterLayer
#'
#' Helper function for \code{\link{krig_gd}}
#'
#' @param r raster for kriging
#' @param grd object to create grid for kriging, can be RasterLayer, SpatialPointsDataFrame, or a gridded object as defined by 'sp'. If undefined, will use \code{r} to create a grid.
#' @param coords if provided, kriging will occur based only on values at these coordinates
#' @param xy whether to co-krige with x and y (~x+y)
#' @param resample whether to resample grd or r. Set to "r" to resample r to grd Set to "grd" to resample grd to r (defaults to FALSE)
#' @param agg_grd factor to use for aggregation of grd, if provided (this will decrease the resolution of the final kriged raster; defaults to NULL)
#' @param disagg_grd factor to use for disaggregation of grd, if provided (this will increase the resolution of the final kriged raster; defaults to NULL)
#' @param agg_r factor to use for aggregation of r, if provided (this will decrease the number of points used in the kriging model; defaults to NULL)
#' @param disagg_r factor to use for disaggregation, of r if provided (this will increase the number of points used in the kriging model; defaults to NULL)
#' @param resample_first if aggregation or disaggregation is used in addition to resampling, whether to resample before (resample_first = TRUE) or after (resample_first = FALSE) aggregation/disaggregation (defaults to TRUE)
#' @param n_cell number of cells to interpolate across if SpatialPointsDataFrame is provided for \code{grd}
#'
#' @return RasterLayer
#' @export
#'
#' @keywords internal
#'
krig_gd_lyr <- function(r, grd = NULL, coords = NULL, xy = FALSE, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE, n_cell = 1000) {

  # TODO: FIX COORDSS
  # if grd is NULL use r
  if (is.null(grd)) grd <- r

  # Transform raster layer
  if (inherits(grd, "RasterLayer")) {
    stk <- raster_transform(r, grd, resample, agg_grd, disagg_grd, agg_r, disagg_r, resample_first)
    r <- stk[[names(r)]]
    grd <- stk[["grd"]]
  }

  # convert raster to df
  krig_df <- data.frame(raster::rasterToPoints(r))

  if (!is.null(coords)) {
    coords <- data.frame(coords)
    colnames(coords) <- c("x", "y")
    rex <- raster::extract(r, coords)
    krig_df <- data.frame(coords, layer = rex)
  }

  # create grid
  if (is.null(grd)) {
    krig_grid <- raster_to_grid(r)
  } else if (inherits(grd, "RasterLayer")) {
    krig_grid <- raster_to_grid(grd)
  } else if (sp::gridded(grd)) {
    krig_grid <- grd
  } else {
    stop(" unable to find an inherited method for type of grd provided")
  }

  # Assign values to df
  krig_df$layer <- krig_df[, 3]

  # convert to spdf
  sp::coordinates(krig_df) <- ~ x + y

  # remove na values
  krig_df <- krig_df[!is.na(krig_df$layer), ]

  # remove crs values (automap doesn't like latlon CRS)
  if (!raster::compareCRS(krig_df, krig_grid)) {
    warning("the provided raster and grid have different crs")
  }
  raster::crs(krig_df) <- NA
  raster::crs(krig_grid) <- NA

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
  krig_gder <- raster::rasterFromXYZ(krig_spdf, crs = raster::crs(krig_grid))
  return(krig_gder)
}


#' Conver a raster to a grid
#'
#' @param x RasterLayer
#'
#' @return gridded SpatialPixelsDataFrame
#' @export
#'
#' @keywords internal
#'
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
#' @return stack of transformed rasters
#' @export
#'
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
