
#' Generate preview of moving window and sample counts
#'
#' @param coords coordinates (two columns, the first should be x and the second should be y)
#' @param sample_count whether to create plot of sample counts for each cell (defaults to TRUE)
#' @param min_n min number of samples to use in calculations (any focal cell with a window containing less than this number of samples will be assigned a value of NA)
#' @inheritParams window_gd
#'
#' @return RasterStack with example window layer and sample counts (if sample_count = TRUE)
#' @export
#'
#' @examples
#' load_mini_ex()
#' preview_gd(mini_lyr, mini_coords, wdim = 3, fact = 3, sample_count = TRUE, min_n = 2)
preview_gd <- function(lyr, coords, wdim, fact = 0, sample_count = TRUE, min_n = 0) {
  if (fact != 0) lyr <- raster::aggregate(lyr, fact)

  # convert wdim to matrix
  nmat <- wdim_to_mat(wdim)

  # get center of raster
  e <- as.vector(raster::extent(lyr))
  c <- c(mean(e[c(1, 2)]), mean(e[c(3, 4)]))
  center <- raster::cellFromXY(lyr, c)

  # get adjacent cells to center cell
  adjc <- raster::adjacent(lyr, center, directions = nmat)
  # get list of indices of coords in that set of cells
  adjci <- purrr::map_dbl(adjc, 1, function(x) {
    seq(x[1], x[2])
  })
  # fill in window
  lyrw <- lyr * 0
  lyrw[adjci] <- 1
  lyrw[center] <- 2

  raster::plot(lyrw, col = viridis::mako(3, direction = -1), legend = FALSE, axes = FALSE, box = FALSE)
  graphics::legend("bottomleft", c("raster layer", "window", "focal cell"), col = viridis::mako(3, direction = -1), pch = 15)

  if (!is.null(coords)) graphics::points(coords, pch = 3, col = viridis::magma(1, begin = 0.7))

  if (sample_count) {

    # get coord cells
    coord_cells <- raster::extract(lyr, coords, cell = TRUE)[, "cells"]

    # count
    lyrc <- lyr
    nc <- purrr::map_dbl(1:raster::ncell(lyr), function(x, lyr, nmat, coord_cells) {
      sub <- get_adj(x, lyr, nmat, coord_cells)
      return(length(sub))
    }, lyr, nmat, coord_cells)
    lyrc <- raster::setValues(lyr, nc)

    lyrc[lyrc < min_n] <- NA
    raster::plot(lyrc, col = viridis::mako(100), box = FALSE, axes = FALSE, main = "sample count")

    # make stack
    lyrw <- raster::stack(lyrw, lyrc)
    names(lyrw) <- c("window", "sample_count")
  }

  return(lyrw)
}
