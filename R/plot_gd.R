#' Plot moving window map of genetic diversity
#'
#' Plot genetic diversity layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (SpatRaster where first layer is genetic diversity)
#' @param index if a raster stack is provided, index of the layer to plot (defaults to plotting all layers except layers named "sample_count")
#' @param bkg optional SpatRaster or other spatial object that will be plotted as the "background" in gray
#' @param col color palette to use for plotting (defaults to \link[viridis]{magma} palette)
#' @param breaks number of breaks to use in color scale (defaults to 100)
#' @param box whether to include a box around the Raster plot (defaults to FALSE)
#' @param range numeric. minimum and maximum values to be used for the continuous legend
#' @param legend whether to include legend
#'
#' @inheritParams terra::plot
#'
#' @return plot of genetic diversity
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_gd(mini_lyr)
#'
plot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(breaks), breaks = 100, main = NULL, box = FALSE, range = NULL, legend = TRUE, ...) {
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)
  if (inherits(bkg, "Raster")) bkg <- terra::rast(bkg)

  # plot all layers except sample counts
  if (is.null(index)) {
    # drop sample count layer
    if (any(names(x) == "sample_count")) x <- x[[names(x) != "sample_count"]]
    index <- 1:terra::nlyr(x)
  }

  plot_general(x, index = index, bkg = bkg, breaks = breaks, col = col, main = main, box = box, range = range, legend = legend, ...)
}

#' Plot moving window map of sample counts
#'
#' Plot sample counts layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x single SpatRaster of counts or SpatRaster where indexed layer is sample counts
#' @param index if a raster stack is provided, index of the sample count layer to plot
#' (defaults to plotting the layer named "sample_count" or the last layer of the stack)
#' @param col color palette to use for plotting (defaults to viridis::magma palette)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param box whether to include a box around the raster plot (defaults to FALSE)
#' @param range numeric. minimum and maximum values to be used for the continuous legend
#' @inheritParams plot_gd
#' @inheritParams terra::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_count(mini_lyr)
plot_count <- function(x, index = NULL, breaks = 100, col = viridis::mako(breaks), main = NULL, box = FALSE, range = NULL, legend = TRUE, ...) {
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)

  # plot sample count layer
  if (is.null(index)) {
    if (any(names(x) == "sample_count")) index <- which(names(x) == "sample_count") else index <- terra::nlyr(x)
  }

  plot_general(x, index = index, breaks = breaks, col = col, main = main, box = box, range = range, legend = legend, ...)
}

#' Helper function for plot_gd and plot_count
#'
#' @inheritParams plot_gd
#'
#' @noRd
plot_general <- function(x, index = 1, bkg = NULL, col = viridis::magma(breaks), breaks = 100, main = NULL, box = FALSE, range = NULL, legend = TRUE, ...) {
  # suppress irrelevant plot warnings
  suppressWarnings({
    if (!is.null(bkg)) {
      purrr::walk(index, plot_bkg, x = x, bkg = bkg, col = col, breaks = breaks, main = main, box = box, range = range, legend = legend, ...)
    } else {
      purrr::walk(index, \(index) {
        terra::plot(x[[index]], col = col, axes = FALSE, box = box, range = range, legend = legend, ...)
        if (is.null(main)) main <- names(x[[index]])
        graphics::title(main = list(main, font = 1), adj = 0)
      })
    }
  })
}


#' Helper function for plot_general
#'
#' @inheritParams plot_gd
#'
#' @noRd
plot_bkg <- function(index, x, bkg, col = viridis::magma(breaks), breaks = 100, main = NULL, box = FALSE, range = NULL, legend = TRUE, ...) {
  if (is.null(main)) main <- names(x[[index]])

  # suppress irrelevant plot warnings
  suppressWarnings({
    # calculate extent
    extx <- terra::ext(x)
    extb <- terra::ext(bkg)
    xmin <- min(min(extx)[1], min(extb)[1])
    xmax <- max(max(extx)[1], max(extb)[1])
    ymin <- min(min(extx)[2], min(extb)[2])
    ymax <- max(max(extx)[2], max(extb)[2])

    terra::plot(x[[index]],
      col = col,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      axes = FALSE,
      box = box,
      range = range,
      legend = legend
    )

    terra::plot(bkg,
      col = "lightgray",
      border = "white",
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      axes = FALSE,
      box = FALSE,
      legend = FALSE,
      add = TRUE
    )

    terra::plot(x[[index]],
      col = col,
      add = TRUE,
      axes = FALSE,
      box = FALSE,
      range = range,
      legend = FALSE,
      ...
    )
  })

  graphics::title(main = list(main, font = 1), adj = 0)
}
