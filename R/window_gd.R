

window_gd <- function(gen, coords, lyr, stat = "het", fact = 0, wdim = 10, rarify = FALSE, rarify_n = 4, rarify_nit = 10, min_n = 2, fun = mean){
  # check to make sure coords and gen align
  check_data(gen, coords)

  # calc stats
  if(stat == "allelic.richness"){
    if(class(gen) == "vcfR"){
      gen <- vcfR::vcfR2genind(gen)
    } else if (class(gen) != "genind"){
      stop("object must be of type vcfR or genind to calculate allelic richness")
    }

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_ar, fact, wdim, rarify, rarify_n, rarify_nit, min_n, fun)

    names(results[[1]]) <- "allelic_richness"

  }

  if(stat == "het"){
    if (class(gen) != "vcfR"){stop("object must be of type vcfR to calculate heterozygosity")}
    gen <- is.het(extract.gt(gen))

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_het, fact, wdim, rarify, rarify_n, rarify_nit, min_n, fun)

    names(results[[1]]) <- "heterozygosity"
  }

  if(stat == "pi"){
    if (class(gen) != "vcfR"){stop("object must be of type vcfR or genind to calculate pi")}
    gen <- vcf_to_dosage(gen)

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_pi, fact, wdim, rarify, rarify_n, rarify_nit, min_n, fun)

    names(results[[1]]) <- "pi"
   }



  names(results[[2]]) <- "sample_count"

  return(results)

}
#' Sliding window genetic diversity
#'
#' @param gen (*note:* order matters! the coordinate and genetic data should be in the same order, there are currently no checks for this.)
#' @param coords coordinates (two columns, the first should be x and the second should be y and the order should be the same as gen),
#' @param lyr RasterLayer to slide window across
#' @param fact aggregation factor to apply to the RasterLayer (reduces computational time)
#' @param wdim dimensions (height x width) of window, if one value is provided a square window is created
#' @param rarify whether to use rarefaction in calculations
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction
#' @param min_n min number of samples to calculate allelic richness for (equal to rarify_n if provided, otherwise defaults to 2)
#' @param fun function to use to summarize data in window (defaults to mean)
#'
#' @return RasterStack that includes a raster of allelic richness and a raster of the number of samples used in the calculation of allelic richness for each cell
#' @export
#'
#' @details ADD WARNINGS
#'
#' @examples
#' \dontrun{
#' window_gd(gen, coords, lyr)
#' }
#'
window_gd_general <- function(gen, coords, lyr, stat = calc_mean_ar, fact = 0, wdim = 10, rarify = FALSE, rarify_n = 4, rarify_nit = 10, min_n = 2, fun = mean) {

  # TODO: ADD FUNCTIONALITY SO RARIFY CAN EQUAL 1

  # make neighbor matrix for window
  nmat <- wdim_to_mat(wdim)

  # make aggregated raster
  ragg <- raster::aggregate(lyr, fact) * 0

  # get cell index for each coordinate
  coord_cells <- raster::extract(ragg, coords, cell = TRUE)[,"cells"]

  rast_vals <- foreach(i = 1:raster::ncell(ragg), .combine = rbind, .packages = c("raster", "purrr", "hierfstat", "stats", "adegenet")) %dopar% {

    result <- window_helper(i, ragg, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun)

    return(result)

  }

  # make copies of rasters
  aragg <- ragg
  nsagg <- ragg
  # assign values to rasters
  aragg[] <- rast_vals[,"gd"]
  nsagg[] <- rast_vals[,"ns"]

  results <- raster::stack(aragg, nsagg)

  return(results)
}

#' Helper function for window calculations
#'
#' @param i
#' @param ragg
#' @param gen
#' @param coord_cells
#' @param nmat
#' @param min_n
#' @param fun
#'
#' @keywords internal
#'
#' @return
#' @export
#'
#' @examples
window_helper <- function(i, ragg, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun){
  # skip if raster value is NA
  if (is.na(ragg[i])) {
    return(data.frame(gd = NA, ns = NA))
  }

  # get sample indices in window
  sub <- get_adj(i, ragg, nmat, coord_cells)

  # if there are too few samples in that window assign the cell value NA
  if (length(sub) < min_n) {
    gd <- NA
  } else if (rarify){
    gd <- rarify_helper(gen, sub, rarify_n, rarify_nit, stat, fun)
  } else {
    gd <- sample_gd(gen, sub, stat)
  }

  # count the number of samples in the window
  ns <- length(sub)

  return(data.frame(gd = gd, ns = ns))
}



rarify_helper <- function(gen, sub, rarify_n, rarify_nit, stat, fun = mean){
    # if number of samples is less than rarify_n, assign the value NA
    if (length(sub) < rarify_n) {
      gd <- NA
    }

    # if number of samples is greater than rarify_n, rarify
    if (length(sub) > rarify_n) {
      gd <- rarify_gd(gen, sub, rarify_nit = rarify_nit, rarify_n = rarify_n, stat = stat, fun = fun)
    }

    # if the number of samples is equal to rarify_n, calculate stat
    if (length(sub) == rarify_n) {
      gd <- sample_gd(gen, sub, stat)
    }

  return(gd)
}


#' Rarify subsample and calculate genetic diversity
#'
#' @param gen dataframe of allelic richness (*note:* cols should be individuals, rows should be loci)
#' @param sub subset of row indices
#' @param rarify_nit number of iterations (i.e. random samples of size rarify_n to draw)
#' @param rarify_n number of points to use for each rarefaction iteration
#' @param stat function to calculate genetic diversity stat
#'
#' @return rarified mean allelic richness for a subsample
#' @export
#'
#' @keywords internal
#'
#' @examples
rarify_gd <- function(gen, sub, rarify_nit = 10, rarify_n = 4, stat = stat, fun = fun) {

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

  # for each of the possible combos get gendiv stat
  gdrar <- apply(cmb[1:rarify_sub, ], 1, sample_gd, gen = gen, stat = stat)

  # summarize rarefaction results
  gd <- stats::na.omit(fun(gdrar))

  return(gd)
}


#' Calculate genetic diversity of a sample
#'
#' @param gen
#' @param sub row indices of subsample
#' @param fun function to use
#'
#' @return mean allelic richness of a subsample
#' @export
#'
#' @keywords internal
#'
#' @examples
sample_gd <- function(gen, sub, stat) {
  gd <- stat(gen[sub,])
  return(gd)
}

calc_mean_ar <- function(gen){
  gen$pop <- rep(factor(1), nrow(gen$tab))
  ar_byloci <- hierfstat::allelic.richness(gen)$Ar[,1]
  gd <- mean(stats::na.omit(ar_byloci))
  return(gd)
}

calc_mean_het <- function(gen){
  gd_byloci <- colMeans(gen, na.rm = TRUE)
  gd <- stats::na.omit(mean(gd_byloci))
  return(gd)
}

calc_mean_pi <- function(gen){
  gd <- pi.dosage(gen, L = ncol(gen))
  return(gd)
}

check_data <- function(gen, coords){

  # check number of samples
  if(class(gen) == "genind"){
    nind <- nrow(gen$tab)
    if (nrow(gen$tab) != nrow(coords)) {
      stop("number of samples in coords data and number of samples in gen data are not equal")
    }
  }

  if(class(gen) == "vcfR"){nind <- ncol(ex_vcf)}
  # check to make sure coords and gen align


}
