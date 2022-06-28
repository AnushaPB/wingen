
#' Raster interpolation using 'autoKrige'
#'
#' @param r RasterLayer or RasterStack
#' @inheritParams krig_gd_lyr
#' @return RasterLayer or RasterStack
#' @export
#'
#' @examples
#' library("raster")
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, nloci = 10, rarify = TRUE)
#' kpi <- krig_gd(wpi, mini_lyr)
#' plot_gd(kpi, main = "Kriged Pi")
#' plot_count(kpi)
#'

krig_gd <- function(r, grd = NULL, coords = NULL, xy = FALSE, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE, n_cell = 10000){

  rls <- raster::as.list(r)

  if(is.null(grd)){
    grd <- r[[1]]
    warning("no grd provided, defaults to using first raster layer to create grd")
  }

  rstk <- purrr::map(rls, krig_gd_lyr, grd, coords, xy, resample, agg_grd, disagg_grd, agg_r, disagg_r, n_cell)
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
#' @param resample whether to resample grd or r. Set to "r" to resample r to grd Set to "grd" to resample grd to r. Defaults to FALSE but we *highly recommend setting this to either "grd" or "r"*
#' @param agg_grd factor to use for aggregation of grd, if provided
#' @param disagg_grd factor to use for disaggregation of grd, if provided
#' @param agg_r factor to use for aggregation of r, if provided
#' @param disagg_r factor to use for disaggregation, of r if provided
#' @param resample_first if aggregation or disaggregation is used in addition to resampling, whether to resample before (resample_first = TRUE) or after (resample_first = FALSE) aggregation/disaggregation (defaults to TRUE)
#' @param n_cell number of cells to interpolate across if SpatialPointsDataFrame is provided for \code{grd}
#'
#' @return RasterLayer
#' @export
#'
#' @keywords internal
#'
#' @examples
krig_gd_lyr <- function(r, grd = NULL, coords = NULL, xy = FALSE, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE, n_cell = 1000) {

  # Transform raster layer
  if(class(grd) == "RasterLayer"){
    stk <- raster_transform(r, grd)
    r <- stk[[names(r)]]
    grd <- stk[["grd"]]
  }

  # convert raster to df
  krig_df <- data.frame(raster::rasterToPoints(r))

  if(!is.null(coords)){
    rex <- raster::extract(r, coords)
    krig_df <- data.frame(coords, layer = rex)
  }

  # create grid
  if(is.null(grd)){
    krig_grid <- raster_to_grid(r)
  } else if(class(grd) == "SpatialPointsDataFrame") {
    #TODO: FIX THIS NOT WORKING
    krig_grid <- spdf_to_grid(grd, n_cell = n_cell)
  } else if(class(grd) == "RasterLayer"){
    krig_grid <- raster_to_grid(grd)
  } else if(sp::gridded(grd)){
    krig_grid <- grd
  } else { stop(" unable to find an inherited method for type of grd provided")}

  # Assign values to df
  krig_df$layer <- krig_df[, 3]

  # convert to spdf
  sp::coordinates(krig_df) <- ~ x + y

  # remove na values
  krig_df <- krig_df[!is.na(krig_df$layer), ]

  # remove crs values (automap doesn't like latlon CRS)
  if(!raster::compareCRS(krig_df, krig_grid)){warning("the provided raster and grid have different crs")}
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
#' @inheritParams krig_gd
#'
#' @return gridded SpatialPixelsDataFrame
#' @export
#'
#' @keywords internal
#'
#' @examples
raster_to_grid <- function(x){
  grd <- data.frame(raster::rasterToPoints(x))
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
#' @return gridded SpatialPixelsDataFrame
#' @export
#'
#' @keywords internal
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


#' Transform raster
#'
#' @inheritParams krig_gd
#'
#' @return stack of transformed rasters
#' @export
#'
#' @examples
raster_transform <- function(r, grd, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE){
  if(nlayers(r) > 1) stop(">1 layer provided for r")
  if(nlayers(grd) > 1) stop(">1 layer provided for grd")

  if(resample_first){
    if(resample == "r") r <- raster::resample(r, grd)
    if(resample == "grd") grd <- raster::resample(grd, r)
  }

  if(!is.null(agg_grd) & !is.null(disagg_grd)) stop("Both agg_grd and disagg_grd provided, when only one should be provided")
  if(!is.null(agg_grd)) grd <- raster::aggregate(grd, agg_grd)
  if(!is.null(disagg_grd)) grd <- raster::disaggregate(grd, disagg_grd)

  if(!is.null(agg_r) & !is.null(disagg_r)) stop("Both agg_r and disagg_r provided, when only one should be provided")
  if(!is.null(agg_r)) r <- raster::aggregate(r, agg_r)
  if(!is.null(disagg_r)) r <- raster::disaggregate(r, disagg_r)

  if(!resample_first){
    if(resample == "r") r <- raster::resample(r, grd)
    if(resample == "grd") grd <- raster::resample(grd, r)
  }

  s <- list(r, grd)
  names(s) <- c(names(r), "grd")

  return(s)
}
