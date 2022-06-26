
#' Raster interpolation using 'autoKrige'
#'
#' @param r RasterLayer or RasterStack
#' @inheritParams krig_gd_lyr
#' @return RasterLayer or RasterStack
#' @export
#'
#' @examples
#' load_mini_ex()
#' wpi <- window_gd(vcf, coords, lyr, stat = "pi", nloci = 10, rarify_n = 4, rarify_nit = 5, rarify = TRUE)
#' kpi <- krig_gd(wpi, lyr)
#' plot_gd(kpi)
#' plot_count(kpi)
#'

krig_gd <- function(r, grd = NULL, coords = NULL, xy = FALSE, resample = FALSE, agg = NULL, disagg = NULL, n_cell = 10000){

  rls <- raster::as.list(r)

  if(is.null(grd)){
    grd <- r[[1]]
    warning("no grd provided, defaults to using first raster layer to create grd")
  }

  rstk <- purrr::map(rls, krig_gd_lyr, grd, coords, xy, resample, agg, disagg, n_cell)
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
#' @param resample if grd is a raster, whether to resample r to the resolution of grd (defaults to FALSE)
#' @param agg factor used for aggregation of grd if provided
#' @param disagg factor used for disaggregation of grd if provided
#' @param n_cell number of cells to interpolate across if SpatialPointsDataFrame is provided for \code{grd}
#'
#' @return RasterLayer
#' @export
#'
#' @keywords internal
#'
#' @examples
krig_gd_lyr <- function(r, grd = NULL, coords = NULL, xy = FALSE, resample = FALSE, agg = NULL, disagg = NULL, n_cell = 1000) {

  # resample raster layer to grd resolution if grd is a raster and resample = TRUE
  if(class(grd) == "RasterLayer" & resample){
    r <- raster::resample(r, grd)
  }

  # convert raster to df
  krig_df <- data.frame(raster::rasterToPoints(r))

  if(!is.null(coords)){
    rex <- raster::extract(r, coords)
    krig_df <- data.frame(coords, layer = rex)
  }

  # create grid
  if(is.null(grd)){
    krig_grid <- raster_to_grid(r, agg = agg, disagg = disagg)
  } else if(class(grd) == "SpatialPointsDataFrame") {
    #TODO: FIX THIS NOT WORKING
    krig_grid <- spdf_to_grid(grd, n_cell = n_cell)
  } else if(class(grd) == "RasterLayer"){
    krig_grid <- raster_to_grid(grd, agg = agg, disagg = disagg)
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
#' @examples
raster_to_grid <- function(x, agg = NULL, disagg = NULL){
  if(!is.null(agg)){x <- raster::aggregate(x, agg)}
  if(!is.null(disagg)){x <- raster::disaggregate(x, disagg)}

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

