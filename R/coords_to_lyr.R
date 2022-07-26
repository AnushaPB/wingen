coords_to_lyr <- function(coords, buffer = 1, agg = NULL, disagg = NULL, plot = TRUE){
  colnames(coords) <- c("x", "y")
  # coords
  xmin <- min(coords$x, na.rm = TRUE) - buffer
  xmax <- max(coords$x, na.rm = TRUE) + buffer
  ymin <- min(coords$y, na.rm = TRUE) - buffer
  ymax <- max(coords$y, na.rm = TRUE) + buffer

  # make a matrix
  m <- matrix(nrow = (ymax - ymin), ncol = (xmax - xmin))

  # turn into raster
  r <- raster(m)

  # set extent
  extent(r) <- c(xmin, xmax, ymin, ymax)

  # aggregate or disaggregate
  if(!is.null(agg)) r <- aggregate(r, agg)
  if(!is.null(disagg)) r <- disaggregate(r, disagg)

  # assign values to make it easier to visualize the resolution
  r <- init(r)
  r[] <- 1:ncell(r)

  # plot raster
  if(plot){
    raster::plot(r, legend = FALSE, col = mako(ncell(r)))
    points(coords, col = viridis::magma(1, begin = 0.7), pch = 3, lwd = 2)
  }

  return(r)
}
