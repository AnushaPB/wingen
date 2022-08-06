

#' Helper function to test different parameter combinations
#'
#' @param params vector of rarify_n, wdim, and disagg values
#' @inheritParams window_gd
#'
#' @export
#'
test_params_empex <- function(params, vcf, coords){
  rarify_n <- as.numeric(params["rarify_n"])
  wdim <- as.numeric(params["wdim"])
  disagg <- as.numeric(params["disagg"])

  lyr <- coords_to_raster(coords, disagg = disagg, plot = FALSE)

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
#' @export
#'
#' @examples
test_empex_plot <- function(r, bkg, zlim = c(0.02, 0.11)){
  r <- raster::stack(r)
  plot_gd(r, bkg = bkg, zlim = zlim, legend = FALSE, breaks = 100)
  return(NULL)
}

#' Convert dataframe to list of vectors
#'
#' @param x dataframe
#'
#' @export
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
#' @export
#'
get_minmax <- function(x, y){
  x <- x[[1]]
  y <- y[[1]]
  mn <- min(cellStats(x, na.rm = TRUE, min), cellStats(y, na.rm = TRUE, min))
  mx <- max(cellStats(x, na.rm = TRUE, max), cellStats(y, na.rm = TRUE, max))
  res <- c(mn, mx)
  names(res) <- c("min", "max")
  return(res)
}


