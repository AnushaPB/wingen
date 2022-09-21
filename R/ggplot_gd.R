#' Plot moving window map of genetic diversity
#'
#' Plot genetic diversity layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where first layer is genetic diversity)
#' @param bkg optional RasterLayer or other spatial object that will be plotted as the "background" in gray
#' @param col color pallete to use for plotting (defaults to viridis::magma pallete)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param index if RasterStack is provided, index of the sample count layer to plot (defaults to plotting first layer)
#' @inheritParams raster::plot
#'
#' @return plot of genetic diversity
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_gd(mini_lyr)
#'
ggplot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {
  if (is.null(index) & raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, plotting first layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 1

  # suppress irrelevant plot warnings
  suppressWarnings({
    if (!is.null(bkg)) {
      plt <- purrr::map(index, plot_gd_bkg, x = x, bkg = bkg, col = col, breaks = breaks, main = main, box = box, ...)
    } else {
      plt <- raster::plot(x[[index]],
        col = col,
        axes = FALSE,
        box = box,
        ...
      )
      title(main = list(main, font = 1), adj = 0)
    }
  })

  return(invisible(plt))
}

#' Helper function for plot_gd
#'
#' @inheritParams plot_gd
#'
#' @export
#' @noRd
ggplot_gd_bkg <- function(index, x, bkg = NULL, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {

  x_df <- x[[index]] %>%
    raster::rasterToPoints() %>%
    tidyr::as_tibble()
  colnames(x_df) <- c("x", "y", "layer")

  # calculate extent
  xmin <- min(raster::extent(bkg)@xmin, raster::extent(x)@xmin)
  xmax <- max(raster::extent(bkg)@xmax, raster::extent(x)@xmax)
  ymin <- min(raster::extent(bkg)@ymin, raster::extent(x)@ymin)
  ymax <- max(raster::extent(bkg)@ymax, raster::extent(x)@ymax)

  # TODO: change so it is more cusotmizable
  # TODO: change so legend isn't cramped
  gg <- ggplot2::ggplot() +
    # change so works if bkg is a raster
    ggplot2::geom_polygon(data = bkg, ggplot2::aes(x = long, y = lat, group = group), fill = "lightgray") +
    ggplot2::geom_tile(data = x_df, ggplot2::aes(x = x, y = y, fill = layer)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_stepsn(n.breaks = breaks, colours = col) +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())

  return()
}

#' Plot moving window map of sample counts
#'
#' Plot sample counts layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x RasterLayer of counts or RasterStack where indexed layer is sample counts
#' @param index if RasterStack is provided, index of the sample count layer to plot (assumes this is a stacked output from window_gd and defaults to plotting second layer)
#' @param col color pallete to use for plotting (defaults to viridis::magma pallete)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @inheritParams plot_gd
#' @inheritParams raster::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_count(mini_lyr)
ggplot_count <- function(x, index = NULL, breaks = 20, col = viridis::mako(breaks), main = NULL, box = FALSE, ...) {
  if (is.null(index) & raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, plotting second layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 2

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({
    if (raster::nlayers(x) > 1) {
      plt <- raster::plot(x[[index]],
        col = col,
        axes = FALSE,
        box = box,
        ...
      )
      title(main = list(main, font = 1), adj = 0)
    }

    if (raster::nlayers(x) == 1) {
      plt <- raster::plot(x,
        col = col,
        axes = FALSE,
        box = box,
        ...
      )
      title(main = list(main, font = 1), adj = 0)
    }
  })

  return(invisible(plt))
}
