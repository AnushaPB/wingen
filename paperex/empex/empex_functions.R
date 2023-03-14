

#' Helper function to test different parameter combinations
#'
#' @param params vector of rarify_n, wdim, and disagg values
#' @inheritParams window_gd
#'
#'
#'
test_params_empex <- function(x, vcf, coords){
  rarify_n <- as.numeric(x["rarify_n"])
  wdim <- as.numeric(x["wdim"])
  res <- as.numeric(x["res"])

  lyr <- coords_to_raster(coords, res = res, plot = FALSE)

  pg <- window_gd(vcf,
                  coords,
                  lyr,
                  stat = "pi",
                  wdim = wdim,
                  fact = 0,
                  rarify = TRUE,
                  rarify_n = rarify_n,
                  rarify_nit = 5,
                  parallel = TRUE,
                  ncores = 4)
  return(pg)
}

#' Create plots from default time test raster results
#'
#' @param r raster
#' @param bkg background plot
#' @param zlim zlimit for plotting
#'
#' @return
#'
#'
#' @examples
test_empex_plot <- function(r, bkg, zlim = c(0.02, 0.11)){
  r <- raster::stack(r)
  raster_plot_gd(r, bkg = bkg, zlim = zlim, legend = FALSE, breaks = 100)
  return(NULL)
}

#' Convert dataframe to list of vectors
#'
#' @param x dataframe
#'
#'
#'
df_to_ls <- function(x){
  x <- split(x, seq(nrow(x)))
  return(x)
}

#' Get minimimum and maximum of two RasterLayers
#'
#' @param x RasterLayer
#' @param y RasterLayer
#'
#'
#'
get_minmax <- function(x, y){
  if(inherits(x, "SpatRaster")) x <- raster(x)
  if(inherits(y, "SpatRaster")) y <- raster(y)
  x <- x[[1]]
  y <- y[[1]]
  mn <- min(cellStats(x, na.rm = TRUE, min), cellStats(y, na.rm = TRUE, min))
  mx <- max(cellStats(x, na.rm = TRUE, max), cellStats(y, na.rm = TRUE, max))
  res <- c(mn, mx)
  names(res) <- c("min", "max")
  return(res)
}

#' Krige and mask with appropriate projections
#'
#' @param x RasterLayer
#' @param lyr layer for kriging
#' @param mask longlat layer for masking
empex_krig_mask <- function(x, lyr, mask){
  kx <- krig_gd(x, index = 1, lyr, disagg_grd = 4)

  mx <-
    terra::project(kx, terra::crs(mask)) %>%
    terra::mask(mask) %>%
    trim()

  return(mx)
}

#' Plot_gd rewritten for raster
#'
#'
raster_plot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {
  if (inherits(x, "SpatRaster")) x <- raster::raster(x)
  if (inherits(x, "SpatRaster")) bkg <- raster::raster(bkg)

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
plot_gd_bkg <- function(index, x, bkg, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {
  # suppress irrelevant plot warnings
  suppressWarnings({
    # calculate extent
    extx <- raster::extent(x)
    extb <- raster::extent(bkg)
    xmin <- min(extx[1], extb[1])
    xmax <- max(extx[2], extb[2])
    ymin <- min(extx[3], extb[3])
    ymax <- max(extx[4], extb[4])

    raster::plot(x[[index]],
                 col = "white",
                 xlim = c(xmin, xmax),
                 ylim = c(ymin, ymax),
                 axes = FALSE,
                 box = box,
                 legend = FALSE
    )

    raster::plot(bkg,
                 col = "lightgray",
                 border = "white",
                 xlim = c(xmin, xmax),
                 ylim = c(ymin, ymax),
                 axes = FALSE,
                 box = FALSE,
                 legend = FALSE,
                 add = TRUE
    )

    raster::plot(x[[index]],
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
raster_plot_count <- function(x, index = NULL, breaks = 20, col = viridis::mako(breaks), main = NULL, box = FALSE, ...) {
  if (inherits(x, "SpatRaster")) x <- raster::raster(x)

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
      graphics::title(main = list(main, font = 1), adj = 0)
    }

    if (raster::nlayers(x) == 1) {
      plt <- raster::plot(x,
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
