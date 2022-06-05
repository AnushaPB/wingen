
#' Sliding window allelic richness
#'
#' @param ar_df dataframe of allelic richness (*note:* rows should be loci, columns should be individuals)
#' @param coords coordinates (two columns, the first should be x and the second should be y and the order should be the same as ar_df),
#' @param lyr RasterLayer to slide window across
#' @param fact aggregation factor to apply to the RasterLayer (reduces computational time)
#' @param wdim dimensions (height x width) of window, if one value is provided a square window is created
#' @param rarify whether to use rarefaction in calculations
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction
#' @param min_n min number of samples to calculate allelic richness for (equal to rarify_n if provided, otherwise defaults to 2)
#' @param fun function to use to summarize data in window (defaults to mean)
#' @param plot_its whether to produce a plot every plot_its_steps iterations to show progress
#' @param plot_its_steps when to plot an iteration if plot_its is true
#'
#' @return RasterStack that includes a raster of allelic richness and a raster of the number of samples used in the calculation of allelic richness for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' window_ar(ar_df, coords, lyr)
#' }
#'
window_ar <- function(genind, coords, lyr, fact = 0, wdim = 10, rarify = FALSE, rarify_n = 4, rarify_nit = 10, min_n = 2, fun = mean, plot_its = FALSE, plot_its_steps = 1) {

  #genind pop
  genind$pop <- rep(factor(1), nrow(genind$tab))

  # TODO: ADD FUNCTIONALITY SO RARIFY CAN EQUAL 1
  # check to make sure coords and ar_df align
  if (nrow(genind$tab) != nrow(coords)) {
    stop("nrow of the coords data and ncol genind data are not equal, make sure rows are individuals for the coords data and cols are individuals for the allelic richness data")
  }

  # make neighbor matrix for window
  nmat <- wdim_to_mat(wdim)

  # make aggregated raster
  ragg <- raster::aggregate(lyr, fact) * 0

  # get cell index for each coordinate
  coord_cells <- raster::extract(ragg, coords, cell = TRUE)[,"cells"]

  rast_vals <- foreach(i = 1:raster::ncell(ragg), .combine = rbind, .packages = c("raster", "purrr", "hierfstat", "stats", "adegenet")) %dopar% {

    result <- window_helper(i, ragg, genind, coord_cells, nmat, min_n, fun)

    return(result)

  }

  # make copies of rasters
  aragg <- ragg
  nragg <- ragg
  # assign values to rasters
  aragg[] <- rast_vals[,"gd"]
  nragg[] <- rast_vals[,"nr"]

  # name rasters
  names(aragg) <- "allelic_richness"
  names(nragg) <- "sample_count"

  results <- raster::stack(aragg, nragg)

  return(results)
}

#' Helper function for window calculations
#'
#' @param i
#' @param ragg
#' @param genind
#' @param coord_cells
#' @param nmat
#' @param min_n
#' @param fun
#'
#' @return
#' @export
#'
#' @examples
window_helper <- function(i, ragg, genind, coord_cells, nmat, min_n, fun){
  # skip if raster value is NA
  if (is.na(ragg[i])) {
    return(data.frame(gd = NA, nr = NA))
  }

  sub <- get_adj(i, ragg, nmat, coord_cells)

  # if there are too few samples in that window assign the cell value NA
  if (length(sub) < min_n) {
    gd <- NA
    nr <- length(sub)
  } else {
    # get allelic richness
    gd <- hierfstat::allelic.richness(genind[sub,])$Ar[,1]
    gd <- fun(stats::na.omit(gd))

    # count the number of points used in the calculation
    nr <- length(sub)
  }

  return(data.frame(gd = gd, nr = nr))
}

#' Rarefaction function
#'
#' @param ar_df dataframe of allelic richness (*note:* cols should be individuals, rows should be loci)
#' @param sub subset of row indices
#' @param rarify_nit number of iterations (i.e. random samples of size rarify_n to draw)
#' @param rarify_n number of points to use for each rarefaction iteration
#' @param fun function to use to calculate the mean (defaults to base R mean)
#'
#' @return rarified mean allelic richness for a subsample
#' @export
#'
#' @keywords internal
#'
#' @examples
rarify_ar <- function(ar_df, sub, rarify_nit = 10, rarify_n = 4, fun = mean) {

  # check to make sure sub is greater than 4
  if (!(length(sub) > rarify_n)) {
    stop("rarify_n is less than the number of samples provided")
  }

  # get all possible combos (transpose so rows are unique combos)
  cmb <- t(utils::combn(sub, rarify_n))

  # define subsample to rarify
  # (note: this is done so when the number of unique combos < rarify_nit, extra calcs aren't performed)
  if (nrow(cmb) < rarify_nit) {
    rarify_sub <- nrow(cmb)
  } else {
    rarify_sub <- rarify_nit
  }

  # for each of the possible combos get AR
  ar <- apply(cmb[1:rarify_sub, ], 1, sample_ar, ar_df = ar_df, fun = fun)

  return(ar)
}


#' Calculate mean allelic richness of a sample
#'
#' @param ar_df dataframe of allelic richness (*note:* cols should be individuals, rows should be loci)
#' @param sub row indices of subsample
#' @param fun function to use to calculate the mean (defaults to base R mean)
#'
#' @return mean allelic richness of a subsample
#' @export
#'
#' @keywords internal
#'
#' @examples
sample_ar <- function(ar_df, sub, fun = mean) {
  ar <- ar_df[, sub] %>%
    rowMeans(na.rm = TRUE) %>%
    stats::na.omit() %>%
    fun()
  return(ar)
}
