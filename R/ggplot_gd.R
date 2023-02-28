#' Plot moving window map of genetic diversity
#'
#' Plot genetic diversity layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x output from \link[wingen]{window_gd} or \link[wingen]{krig_gd} (RasterStack where first layer is genetic diversity)
#' @param bkg optional RasterLayer or sf polygon
#' @param col color pallete to use for plotting (defaults to \link[viridis]{magma} pallete)
#' @param index if RasterStack is provided, index of the sample count layer to plot (defaults to plotting first layer)
#'
#' @return list of ggplots
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_gd(mini_lyr)
#'
ggplot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(100), ...) {

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


#' Plot moving window map of sample counts
#'
#' Plot sample counts layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x single SpatRaster of counts or SpatRaster where indexed layer is sample counts
#' @param index if a raster stack is provided, index of the sample count layer to plot
#' (assumes this is a stacked output from window_gd and defaults to plotting the last layer (which should be sample counts))
#' @param col color palette to use for plotting (defaults to viridis::magma palette)
#'
#' @return list of ggplots
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_count(mini_lyr)
ggplot_count <- function(x, index = NULL, breaks = 100, col = viridis::mako(breaks), ...) {
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)
  if (is.null(index) & terra::nlyr(x) > 2) warning("More than two raster layers in stack provided, plotting second layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- terra::nlyr(x)

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

#' Helper function for ggplot_gd and ggplot_count
#'
#' @param var variable to plot
#' @param x_df ggplot dataframe
#' @param bkg background raster layer or sf objects
#' @return ggplot
#' @noRd
ggplot_helper <- function(var, x_df, bkg = NULL){
  # create ggplot
  gg <- ggplot2::ggplot()

  # add background
  if (!is.null(bkg)) {
    if (inherits(bkg, "sf")) gg <- gg + ggplot2::geom_sf(data = bkg, fill = "lightgray")
    if (inherits(bkg, "SpatRaster")) {
      bkg_df <- bkg %>%
        terra::as.data.frame(xy = TRUE) %>%
        tidyr::as_tibble()
      gg <- gg + ggplot2::geom_tile(data = bkg_df, ggplot2::aes(x = x, y = y), fill = "lightgray")
    }

  }

  # plot result
  gg <- gg +
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


