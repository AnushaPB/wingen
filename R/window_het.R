
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
window_het <- function(ar_df, coords, lyr, fact = 0, wdim = 10, rarify = FALSE, rarify_n = 4, rarify_nit = 10, min_n = 2, fun = mean, plot_its = FALSE, plot_its_steps = 1) {

  # TODO: ADD FUNCTIONALITY SO RARIFY CAN EQUAL 1
  # check to make sure coords and ar_df align
  if (ncol(ar_df) != nrow(coords)) {
    stop("nrow of the coords data and ncol allelic richness data are not equal, make sure rows are individuals for the coords data and cols are individuals for the allelic richness data")
  }

  # make ar_df into dataframe
  ar_df <- data.frame(ar_df)

  # make aggregated raster
  ragg <- raster::aggregate(lyr, fact) * 0

  # make neighbor matrix for window
  n <- wdim_to_mat(wdim)

  # make copy of raster for ar surface
  gdagg <- ragg
  # make copy of raster to count number of individuals used in the calculation
  nragg <- ragg

  # get cell index for each coordinate
  coord_cells <- raster::extract(ragg, coords, cell = TRUE)[,"cells"]

  # for every cell in aragg, calculate allelic richness
  pb <- progress::progress_bar$new(total = raster::ncell(aragg))

  for (i in 1:raster::ncell(ragg)) {
    # skip if raster value is NA
    if (is.na(ragg[i])) {
      # skip to next
      next
    }

    # get adjacent cells to cell i
    adjc <- raster::adjacent(ragg, i, directions = n, include = TRUE, sorted = TRUE)
    # get indices of adjacent cells
    adjci <- purrr::map_dbl(adjc, 1, function(x) {
      seq(x[1], x[2])
    })
    # get list of indices of coords that fall in cells
    sub <- which(coord_cells %in% adjci)

    # if there are too few samples in that window assign the cell value NA and skip to the next iteration
    if (length(sub) < min_n) {
      gdagg[i] <- NA
      nragg[i] <- length(sub)

      next
    }

    if (rarify) {

      # if number of samples is less than rarify_n, assign the value NA
      if (length(sub) < rarify_n) {
        gdagg[i] <- NA
        nragg[i] <- length(sub)

        next
      }

      # if number of samples is greater than rarify_n, rarify
      if (length(sub) > rarify_n) {
        gdagg[i] <- rarify_ar(ar_df, sub, rarify_nit = rarify_nit, rarify_n = rarify_n, fun = fun) %>%
          stats::na.omit() %>%
          fun()
      }

      # if the number of samples is equal to rarify_n, calculate the raw mean
      if (length(sub) == rarify_n) {
        gdagg[i] <- ar_df[, sub] %>%
          rowMeans(na.rm = TRUE) %>%
          stats::na.omit() %>%
          fun()
      }
    } else {

      # get allelic richness, first averages by loci and then by individual
      gdagg[i] <- ar_df[, sub] %>%
        rowMeans(na.rm = TRUE) %>%
        stats::na.omit() %>%
        fun()
    }

    # every set number of steps plot aragg
    if (plot_its & i %% plot_its_steps == 0) {
      raster::plot(gdagg, col = viridis::turbo(100), legend = FALSE, axes = FALSE, box = FALSE, main = "Value")
    }

    # count the number of points used in the calculation
    nragg[i] <- length(sub)

  }

  # name rasters
  names(gdagg) <- "heterozygosity"
  names(nragg) <- "sample_count"

  results <- raster::stack(gdagg, nragg)

  return(results)
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
rarify_het <- function(ar_df, sub, rarify_nit = 10, rarify_n = 4, fun = mean) {

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

  # for each of the possible combos get heterozygosity
  gd <- apply(cmb[1:rarify_sub, ], 1, sample_ar, ar_df = ar_df, fun = fun)

  return(gd)
}


#' Calculate mean heterozygosity of a sample
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
  het <- ar_df[, sub] %>%
    rowMeans(na.rm = TRUE) %>%
    stats::na.omit() %>%
    fun()
  return(het)
}

