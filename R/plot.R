#' Plot genetic diversity results
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where first layer is genetic diversity)
#' @param bkg RasterLayer or other spatial object that will be plotted as the "background" in gray
#' @param zlim limits of the color scale values
#' @inheritParams raster::plot
#'
#' @return plot of genetic diversity
#' @export
#'
#' @examples
plot_gd <- function(x, bkg = NULL, col = viridis::magma(100), zlim = NULL, main = NULL, legend = TRUE){

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({
    if(raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, using first layer")

    if(!is.null(bkg)) {
      raster::plot(x[[1]], col = "white", legend = FALSE, main = main, axes = FALSE, box = FALSE)
      raster::plot(bkg, col = "lightgray", border = "white", axes = FALSE, box = FALSE, add = TRUE, legend = FALSE)
      raster::plot(x[[1]], col = col, zlim = zlim, add = TRUE, axes = FALSE, box = FALSE, legend = legend)
    } else {
      raster::plot(x[[1]], col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE, legend = legend)
    }
  }
  )

}

#' Plot sample counts
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where second layer is sample counts)
#' @param zlim limits of the color scale values
#' @inheritParams raster::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
plot_count <- function(x, col = viridis::mako(100), zlim = NULL, main = NULL, legend = legend){

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({

  if(nlayers(x) > 2) warning("More than two raster layers in stack provided, using second layer")
  if(nlayers(x) == 2) raster::plot(x[[2]], col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE, legend = legend)
  if(nlayers(x) == 1) raster::plot(x, col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE, legend = legend)

  })
}
