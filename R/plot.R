#' Plot genetic diversity results
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where first layer is genetic diversity)
#' @inheritParams raster::plot
#'
#' @return plot of genetic diversity
#' @export
#'
#' @examples
plot_gd <- function(x, col = viridis::magma(100), zlim = NULL, main = NULL){
  raster::plot(x[[1]], col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE)
}

#' Plot sample counts
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where second layer is sample counts)
#' @inheritParams raster::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
plot_count <- function(x, col = viridis::mako(100), zlim = NULL, main = NULL){
  raster::plot(x[[2]], col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE)
}
