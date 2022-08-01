#' Create RasterLayer from coordinates
#'
#' @param coords coordinates (two columns, the first should be x and the second should be y)
#' @param buffer buffer to add to edge of raster (defaults to 0)
#' @param res desired resolution of raster (defaults to NULL). Can be a single value for square cells or a vector with two values representing x and y resolutions.
#' @param agg aggregation factor to apply to raster (defaults to NULL)
#' @param disagg disaggregation factor to apply to raster (defaults to NULL)
#' @param plot whether to plot resulting raster with coords (defaults to FALSE)
#'
#' @return RasterLayer
#' @export
#'
#' @examples
#' load_mini_ex()
#' coords_to_raster(mini_coords, buffer = 1, plot = TRUE)
coords_to_raster <- function(coords, buffer = 0, res = NULL, agg = NULL, disagg = NULL, plot = FALSE){

  # format coords
  coords <- as.data.frame(coords)
  colnames(coords) <- c("x", "y")

  # coords
  xmin <- min(coords$x, na.rm = TRUE) - buffer
  xmax <- max(coords$x, na.rm = TRUE) + buffer
  ymin <- min(coords$y, na.rm = TRUE) - buffer
  ymax <- max(coords$y, na.rm = TRUE) + buffer

  # make a matrix
  if(is.null(res)){
    nrow = (ymax - ymin)
    ncol = (xmax - xmin)
  } else {
    if(length(res) == 1){
      nrow <- (ymax - ymin)/res
      ncol <- (xmax - xmin)/res
    } else if (length(res) == 2) {
      nrow <- (ymax - ymin)/(res[1])
      ncol <- (xmax - xmin)/(res[2])
    } else {
      stop("invalid res provided")
    }
  }

  m <- matrix(nrow = nrow, ncol = ncol)

  # turn into raster
  r <- raster::raster(m)

  # set extent
  raster::extent(r) <- c(xmin, xmax, ymin, ymax)

  # aggregate or disaggregate
  if(!is.null(agg) & !is.null(disagg)){warning("both agg and disagg were provided. Did you mean to do this? (if so, note that aggregation will occur first and then disaggregation second")}
  if(!is.null(agg)) r <- raster::aggregate(r, agg)
  if(!is.null(disagg)) r <- raster::disaggregate(r, disagg)

  # assign values to make it easier to visualize the resolution
  r <- raster::init(r)
  r[] <- 1:raster::ncell(r)

  # plot raster
  if(plot){
    raster::plot(r, legend = FALSE, col = mako(ncell(r)))
    points(coords, col = viridis::magma(1, begin = 0.7), pch = 3, lwd = 2)
  }

  return(r)
}
