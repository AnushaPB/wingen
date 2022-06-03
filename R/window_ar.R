
#' Sliding window allelic richness
#'
#' @param ar_df dataframe of allelic richness (*note:* rows should be loci, columns should be individuals)
#' @param coords coordinates (two columns, the first should be x and the second should be y and the order should be the same as ar_df),
#' @param lyr raster layer to slide window across
#' @param fact factor of aggregation to apply to the raster layer (reduces computational time)
#' @param wdim dimensions (height x width) of window, if one value is provided a square window is created
#' @param rarify whether to use rarefaction in calculations
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction
#' @param min_n min number of samples to calculate allelic richness for (equal to rarify_n if provided, otherwise defaults to 2)
#' @param fast if TRUE uses Rcpp mean to calculate mean, which is slightly faster than base R (*note:* compromises on numerical accuracy)
#' @param plot_its whether to produce a plot every plot_its_steps iterations to show progress
#' @param plot_its_steps when to plot an iteration if plot_its is true
#'
#' @return RasterStack that includes a raster of allelic richness and a raster of the number of samples used in the calculation of allelic richness for each cell
#' @export
#'
#' @examples
window_ar <- function(ar_df, coords, lyr, fact = 0, wdim = 10, rarify = TRUE, rarify_n = 4, rarify_nit = 10, min_n = 2, fast = FALSE, plot_its = FALSE, plot_its_steps = 1){

  # check to make sure coords and ar_df align
  if(ncol(ar_df) != nrow(coords)){stop("nrow of the coords data and ncol allelic richness data are not equal,
                                       make sure rows are individuals for the coords data and cols are individuals
                                       for the allelic richness data")}

  # define which mean function to use
  if(fast){meanf <- meanC} else {meanf <- mean}

  # make ar_df into dataframe
  ar_df <- data.frame(ar_df)

  # make aggregated raster
  agg <- aggregate(lyr, fact)*0
  # make a copy that will later be used for window calcs
  ragg <- agg

  # wdim has to be odd
  if(length(wdim) == 1 & wdim %% 2 == 0){wdim <- wdim + 1; warning(paste("wdim must be odd, using wdim =", wdim, "instead"))}
  if(length(wdim) == 2 & wdim[1] %% 2 == 0){wdim[1] <- wdim[1] + 1; warning(paste("wdim must be odd, using wdim[1] =", wdim, "instead"))}
  if(length(wdim) == 2 & wdim[2] %% 2 == 0){wdim[2] <- wdim[2] + 1; warning(paste("wdim must be odd, using wdim[2] =", wdim, "instead"))}

  # make neighbor matrix for window
  if(length(wdim) == 2){n <- matrix(1, wdim[1], wdim[2])}
  if(length(wdim) == 1){n <- matrix(1, wdim, wdim)}

  # focal cell (center of matrix) has to be zero
  n[wdim/2 + 0.5, wdim/2 + 0.5] <- 0

  # make copy of raster for ar surface
  aragg <- ragg
  # make copy of raster to count number of individuals used in the calculation
  nragg <- ragg

  # for every cell in aragg, calculate allelic richness
  pb <- progress_bar$new(total = ncell(aragg))
  for(i in 1:ncell(aragg)){
    # skip if raster value is NA
    if(is.na(aragg[i])) next

    # reset raster every iteration
    ragg <- agg

    # get adjacent cells to cell i
    adjc <- adjacent(ragg, i, directions = n, include = TRUE, sorted = TRUE)
    # get indices of adjacent cells
    adjci <- map_dbl(adjc, 1, function(x){seq(x[1],x[2])})
    # assign adjacent cells a value of 1
    ragg[adjci] <- 1

    # extract coordinates to determine which cells are in that cluster of cells
    coords$group <- raster::extract(ragg, coords[,c("x","y")])
    # get list of indices (COORDS AND AR_DF must be in same order)
    sub <- which(coords$group == 1)

    # if there are too few samples in that window assign the cell value NA and skip to the next iteration
    if(length(sub) < min_n){aragg[i] <- NA; nragg[i] <- length(sub); next}

    if(rarify){

      # if number of samples is less than rarify_n, assign the value NA
      if(length(sub) < rarify_n){aragg[i] <- NA; nragg[i] <- length(sub); next}

      # if number of samples is greater than rarify_n, rarify
      if(length(sub) > rarify_n){aragg[i] <- rarify_ar(ar_df, sub, rarify_nit = rarify_nit, rarify_n = rarify_n, meanf = meanf) %>% na.omit() %>% meanf()}

      # if the number of samples is equal to rarify_n, calculate the raw mean
      if(length(sub) == rarify_n){aragg[i] <- ar_df[,sub] %>% rowMeans(na.rm = TRUE) %>% na.omit() %>% meanf()}

    } else {

      # get allelic richness, first averages by loci and then by individual
      aragg[i] <- ar_df[,sub] %>% rowMeans(na.rm = TRUE) %>% na.omit() %>% meanf()

    }

    # every set number of steps plot aragg
    if(plot_its & i %% plot_its_steps == 0){
      raster::plot(aragg, col = turbo(100), legend = FALSE, axes = FALSE, box = FALSE, main = "Value")

    }

    # count the number of points used in the calculation
    nragg[i] <- length(sub)

    # progress bar
    pb$tick()

  }

  # name rasters
  names(aragg) <- "allelic_richness"
  names(nragg) <- "sample_count"

  results <- stack(aragg, nragg)

  return(results)

}



#' Rarefaction function
#'
#' @param ar_df dataframe of allelic richness (*note:* cols should be individuals, rows should be loci)
#' @param sub subset of row indices
#' @param rarify_nit number of iterations (i.e. random samples of size rarify_n to draw)
#' @param rarify_n number of points to use for each rarefaction iteration
#' @param meanf function to use to calculate the mean (defaults to base R mean)
#'
#' @return rarified mean allelic richness for a subsample
#' @export
#'
#' @examples
rarify_ar <- function(ar_df, sub, rarify_nit = 10, rarify_n = 4, meanf = mean){

  # check to make sure sub is greater than 4
  if(!(length(sub) > rarify_n)){stop("rarify_n is less than the number of samples provided")}

  # get all possible combos (transpose so rows are unique combos)
  cmb <- t(combn(sub, rarify_n))

  # define subsample to rarify (note: this is done so when the number of unique combos < rarify_nit, extra calcs aren't performed)
  if(nrow(cmb) < rarify_nit){rarify_sub <- nrow(cmb)} else {rarify_sub <- rarify_nit}

  # for each of the possible combos get AR
  ar <- apply(cmb[1:rarify_sub,], 1, sample_ar, ar_df = ar_df, meanf = meanf)

  return(ar)
}


#' Calculate mean allelic richness of a sample
#'
#' @param ar_df dataframe of allelic richness (*note:* cols should be individuals, rows should be loci)
#' @param sub row indices of subsample
#' @param meanf function to use to calculate the mean (defaults to base R mean)
#'
#' @return mean allelic richness of a subsample
#' @export
#'
#' @examples
sample_ar <- function(ar_df, sub, meanf = mean){
  ar <- ar_df[,sub] %>% rowMeans(na.rm = TRUE) %>% na.omit() %>% meanf()
  return(ar)
}

