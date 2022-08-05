#' Plot genetic diversity results
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where first layer is genetic diversity)
#' @param bkg RasterLayer or other spatial object that will be plotted as the "background" in gray
#' @param col color pallete to use for plotting (defaults to viridis::magma pallete)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param index if RasterStack is provided, index of the sample count layer to plot (defaults to plotting first layer)
#' @param zlim limits of the color scale values
#' @param legend if FALSE, no legend bar is displayed.
#' @param legend.width Width in characters of the legend bar
#' @param axis.args A list of additional arguments for the legend axis
#' @inheritParams raster::plot
#'
#' @return plot of genetic diversity
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_gd(mini_lyr)
#'
plot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(breaks), breaks = 10, zlim = NULL, main = NULL, legend = TRUE, legend.width = 1, axis.args = list(cex.axis = 1)) {
  if (is.null(index) & raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, plotting first layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 1

  # suppress irrelevant plot warnings
  suppressWarnings({
    if (!is.null(bkg)) {
      purrr::map(index, plot_gd_bkg,
        x = x, bkg = bkg, col = col, breaks = breaks, zlim = zlim,
        main = main, legend = legend, legend.width = legend.width, axis.args = axis.args
      )
    } else {
      raster::plot(x[[index]],
        col = col,
        zlim = zlim,
        main = main,
        axes = FALSE,
        box = FALSE,
        legend = legend,
        legend.width = legend.width,
        axis.args = axis.args
      )
    }
  })

  return()
}

#' Helper function for plot_gd
#'
#' @inheritParams plot_gd
#'
#' @keywords internal
#'
#' @export
#'
plot_gd_bkg <- function(index, x, bkg = NULL, col = viridis::magma(breaks), breaks = 10, zlim = NULL,
                        main = NULL, legend = TRUE, legend.width = 1, axis.args = list(cex.axis = 1)) {

  # suppress irrelevant plot warnings
  suppressWarnings({

    # calculate extent
    xmin <- min(raster::extent(bkg)@xmin, raster::extent(x)@xmin)
    xmax <- max(raster::extent(bkg)@xmax, raster::extent(x)@xmax)
    ymin <- min(raster::extent(bkg)@ymin, raster::extent(x)@ymin)
    ymax <- max(raster::extent(bkg)@ymax, raster::extent(x)@ymax)

    raster::plot(bkg,
      col = "lightgray",
      border = "white",
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      axes = FALSE,
      box = FALSE,
      legend = FALSE,
      main = main
    )

    raster::plot(x[[index]],
      col = col,
      zlim = zlim,
      add = TRUE,
      axes = FALSE,
      box = FALSE,
      legend = legend,
      legend.width = legend.width,
      axis.args = axis.args
    )
  })


  return()
}

#' Plot sample counts
#'
#' @param x RasterLayer of counts or RasterStack where indexed layer is sample counts
#' @param index if RasterStack is provided, index of the sample count layer to plot (assumes this is a stacked output from window_gd and defaults to plotting second layer)
#' @param col color pallete to use for plotting (defaults to viridis::magma pallete)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param zlim limits of the color scale values
#' @inheritParams plot_gd
#' @inheritParams raster::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_count(mini_lyr)
plot_count <- function(x, index = NULL, breaks = 10, col = viridis::mako(breaks), zlim = NULL, main = NULL, legend = TRUE, legend.width = 1, axis.args = list(cex.axis = 1)) {
  if (is.null(index) & raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, plotting second layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 2

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({
    if (raster::nlayers(x) > 1) {
      raster::plot(x[[index]],
        col = col,
        zlim = zlim,
        main = main,
        axes = FALSE,
        box = FALSE,
        legend = legend,
        legend.width = legend.width,
        axis.args = axis.args
      )
    }

    if (raster::nlayers(x) == 1) {
      raster::plot(x,
        col = col,
        zlim = zlim,
        main = main,
        axes = FALSE,
        box = FALSE,
        legend = legend,
        legend.width = legend.width,
        axis.args = axis.args
      )
    }
  })

  return()
}
