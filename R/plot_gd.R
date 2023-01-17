#' Plot moving window map of genetic diversity
#'
#' Plot genetic diversity layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (SpatRaster where first layer is genetic diversity)
#' @param index if a raster stack is provided, index of the sample count layer to plot (defaults to plotting first layer)
#' @param bkg optional SpatRaster or other spatial object that will be plotted as the "background" in gray
#' @param col color palette to use for plotting (defaults to viridis::magma palette)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param box whether to include a box around the Raster plot (defaults to FALSE)
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
plot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)
  if (!inherits(x, "SpatRaster")) bkg <- terra::rast(bkg)

  if (is.null(index) & terra::nlyr(x) > 2) warning("More than two raster layers in stack provided, plotting first layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 1

  # suppress irrelevant plot warnings
  suppressWarnings({
    if (!is.null(bkg)) {
      plt <- purrr::map(index, plot_gd_bkg, x = x, bkg = bkg, col = col, breaks = breaks, main = main, box = box, ...)
    } else {
      plt <- terra::plot(x[[index]],
        col = col,
        axes = FALSE,
        box = box,
        ...
      )
      graphics::title(main = list(main, font = 1), adj = 0)
    }
  })

  return(invisible(plt))
}

#' Helper function for plot_gd
#'
#' @inheritParams plot_gd
#'
#' @noRd
plot_gd_bkg <- function(index, x, bkg = NULL, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {
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
      col = "white",
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      axes = FALSE,
      box = box,
      legend = FALSE
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
      ...
    )
  })

  graphics::title(main = list(main, font = 1), adj = 0)

  return()
}

#' Plot moving window map of sample counts
#'
#' Plot sample counts layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x single SpatRaster of counts or SpatRaster where indexed layer is sample counts
#' @param index if a raster stack is provided, index of the sample count layer to plot (assumes this is a stacked output from window_gd and defaults to plotting second layer)
#' @param col color palette to use for plotting (defaults to viridis::magma palette)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param box whether to include a box around the raster plot (defaults to FALSE)
#' @inheritParams plot_gd
#' @inheritParams terra::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_count(mini_lyr)
plot_count <- function(x, index = NULL, breaks = 20, col = viridis::mako(breaks), main = NULL, box = FALSE, ...) {
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)

  if (is.null(index) & terra::nlyr(x) > 2) warning("More than two raster layers in stack provided, plotting second layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 2

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({
    if (terra::nlyr(x) > 1) {
      plt <- terra::plot(x[[index]],
        col = col,
        axes = FALSE,
        box = box,
        ...
      )
      graphics::title(main = list(main, font = 1), adj = 0)
    }

    if (terra::nlyr(x) == 1) {
      plt <- terra::plot(x,
        col = col,
        axes = FALSE,
        box = box,
        ...
      )
      graphics::title(main = list(main, font = 1), adj = 0)
    }
  })

  return(invisible(plt))
}
