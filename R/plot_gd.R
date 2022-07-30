#' Plot genetic diversity results
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where first layer is genetic diversity)
#' @param col color pallete to use for plotting (defaults to viridis::magma pallete)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param bkg RasterLayer or other spatial object that will be plotted as the "background" in gray
#' @param zlim limits of the color scale values
#' @inheritParams raster::plot
#'
#' @return plot of genetic diversity
#' @export
#'
#' @examples
plot_gd <- function(x, bkg = NULL, col = viridis::magma(breaks), breaks = 10, zlim = NULL, main = NULL, legend = TRUE, legend.width = 1, axis.args = list(cex.axis = 1)){

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({
    if(raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, using first layer")

    if(!is.null(bkg)) {
      raster::plot(x[[1]], col = "white", legend = FALSE, main = main, axes = FALSE, box = FALSE)
      raster::plot(bkg, col = "lightgray", border = "white", axes = FALSE, box = FALSE, add = TRUE, legend = FALSE)
      raster::plot(x[[1]], col = col, zlim = zlim, add = TRUE, axes = FALSE, box = FALSE, legend = legend, legend.width = legend.width, axis.args = axis.args)
    } else {
      raster::plot(x[[1]], col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE, legend = legend, legend.width = legend.width,  axis.args = axis.args)
    }
  }
  )

}

#' Plot sample counts
#'
#' @param x RasterLayer of counts or RasterStack where indexed layer is sample counts
#' @param index if RasterStack is provided, index of the sample count layer to plot (defaults to 2)
#' @param col color pallete to use for plotting (defaults to viridis::magma pallete)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param zlim limits of the color scale values
#' @inheritParams raster::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
plot_count <- function(x, index = 2, breaks = 10, col = viridis::mako(breaks), zlim = NULL, main = NULL, legend = TRUE, legend.width = 1, axis.args = list(cex.axis = 1)){

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({

  if(raster::nlayers(x) > 1) raster::plot(x[[index]], col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE, legend = legend, legend.width = legend.width, axis.args = axis.args)
  if(raster::nlayers(x) == 1) raster::plot(x, col = col, zlim = zlim, main = main, axes = FALSE, box = FALSE, legend = legend, legend.width = legend.width, axis.args = axis.args)

  })
}
