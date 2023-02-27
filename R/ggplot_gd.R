#' Plot moving window map of genetic diversity
#'
#' Plot genetic diversity layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where first layer is genetic diversity)
#' @param bkg optional RasterLayer or other spatial object that will be plotted as the "background" in gray
#' @param col color pallete to use for plotting (defaults to viridis::magma pallete)
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
ggplot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(100), main = NULL, box = FALSE, ...) {

  if (is.null(index) & terra::nlyr(x) > 2) warning("More than two raster layers in stack provided, plotting first layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 1

  # make df
  x_df <- x[[index]] %>%
    terra::as.data.frame(xy = TRUE) %>%
    tidyr::as_tibble()

  # plot results
  plts <-
    x_df  %>%
    dplyr::select(-x, -y) %>%
    purrr::map(\(var) ggplot_helper(var = var, x_df = x_df, bkg = bkg))

  purrr::walk(plts, print)

  return(plts)
}

ggplot_helper <- function(var, x_df, bkg = NULL){
  gg <- ggplot2::ggplot()

  if (!is.null(bkg)) {
    if (inherits(bkg, "sf")) gg <- gg + ggplot2::geom_sf(data = bkg, fill = "lightgray")
    if (inherits(bkg, "SpatRaster")) {
      bkg_df <- bkg %>%
        terra::as.data.frame(xy = TRUE) %>%
        tidyr::as_tibble()
      gg <- gg + ggplot2::geom_tile(data = bkg_df, ggplot2::aes(x = x, y = y), fill = "lightgray")
    }

  }

  gg <- gg +
    # change so works if bkg is a raster
    ggplot2::geom_tile(data = x_df, ggplot2::aes(x = x, y = y, fill = {{var}})) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradientn(colours = col, na.value = rgb(0,0,0,0)) +
    ggplot2::labs(fill = deparse(substitute(var))) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())

  return(gg)
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
