

#' Sliding window map of genetic diversity
#'
#' @param vcf object of type vcf ( (*note:* order matters! the coordinate and genetic data should be in the same order, there are currently no checks for this.))
#' @param stat genetic diversity stat to calculate (can either be "pi", "allelic.richness", or "het")
#' @param parallel whether to use parallelization
#'
#' @inheritParams window_gd_general
#' @return RasterStack that includes a raster of genetic diversity and a raster of the number of samples within the window for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' window_gd(vcf, coords, lyr, stat = "het")
#' window_gd(vcf, coords, lyr, stat = "pi")
#' window_gd(vcf, coords, lyr, stat = "allelic.richness")
#' }
window_gd <- function(vcf, coords, lyr, stat = "het", fact = 0, wdim = 10, rarify = FALSE, rarify_n = 4, rarify_nit = 10, min_n = 2, fun = mean, parallel = FALSE){
  # check that the input file is a vcf or a path to a vcf object
  if(class(vcf) != "vcfR" & is.character(vcf)){
    vcf <- read.vcfR(vcf)
  } else if (class(vcf) != "vcfR" & !is.character(vcf)){
    stop("gen object must be of type vcfR or a path to a .vcf files")
  }

  # check to make sure coords and gen align
  check_data(vcf, coords)

  # calc stats
  if(stat == "allelic.richness"){
    gen <- vcfR::vcfR2genind(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_ar, fact, wdim, rarify, rarify_n, rarify_nit, min_n, fun, parallel)

    names(results[[1]]) <- "allelic_richness"

  }

  if(stat == "het"){
    gen <- vcfR::is.het(vcfR::extract.gt(vcf), na_is_false = FALSE)

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_het, fact, wdim, rarify, rarify_n, rarify_nit, min_n, fun, parallel)

    names(results[[1]]) <- "heterozygosity"
  }

  if(stat == "pi"){
    gen <- vcf_to_dosage(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = calc_pi, fact, wdim, rarify, rarify_n, rarify_nit, min_n, fun, parallel)

    names(results[[1]]) <- "pi"
   }



  names(results[[2]]) <- "sample_count"

  return(results)

}
#' Helper function for window_gd
#'
#' @param gen genetic data (*note:* order matters! the coordinate and genetic data should be in the same order, there are currently no checks for this.)
#' @param coords coordinates (two columns, the first should be x and the second should be y and the order should be the same as the genetic data),
#' @param lyr RasterLayer to slide window across
#' @param fact aggregation factor to apply to the RasterLayer (*note:* increasing this value reduces computational time)
#' @param wdim dimensions (height x width) of window, if only one value is provided a square window is created
#' @param rarify if rarify = TRUE, rarefaction is performed
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction
#' @param min_n min number of samples to calculate allelic richness for (equal to rarify_n if provided, otherwise defaults to 2)
#' @param fun function to use to summarize data in window (defaults to base R mean)
#'
#' @return RasterStack that includes a raster of genetic diversity and a raster of the number of samples within the window for each cell
#' @export
#'
#' @importFrom foreach %dopar%
#'
#' @keywords internal
#'
#' @examples
#'
window_gd_general <- function(gen, coords, lyr, stat = calc_mean_ar, fact = 0, wdim = 10, rarify = FALSE, rarify_n = 4, rarify_nit = 10, min_n = 2, fun = mean, parallel = FALSE) {

  # TODO: ADD FUNCTIONALITY SO RARIFY CAN EQUAL 1

  # make neighborhood matrix for window
  nmat <- wdim_to_mat(wdim)

  # make aggregated raster
  if(fact == 0){lyr <- lyr * 0} else {lyr <- raster::aggregate(lyr, fact) * 0}


  # get cell index for each coordinate
  coord_cells <- raster::extract(lyr, coords, cell = TRUE)[,"cells"]

  if(parallel){
    rast_vals <- foreach::foreach(i = 1:raster::ncell(lyr), .combine = rbind, .packages = c("raster", "purrr", "hierfstat", "stats", "adegenet")) %dopar% {

      result <- window_helper(i, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun)

      return(result)

    }
  } else {

    rast_vals <- purrr::map_dfr(1:raster::ncell(lyr), window_helper, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun)

  }


  # make copies of rasters
  alyr <- lyr
  nsagg <- lyr
  # assign values to rasters
  alyr[] <- rast_vals[,"gd"]
  nsagg[] <- rast_vals[,"ns"]

  results <- raster::stack(alyr, nsagg)

  return(results)
}

#' Helper function for window calculations
#'
#' @param i cell index
#' @param coord_cells cell indices for each coordinate
#' @param nmat neighborhood matrix
#'
#' @inheritParams window_gd_general
#'
#' @keywords internal
#'
#' @return
#' @export
#'
#' @examples
window_helper <- function(i, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun){
  # skip if raster value is NA
  if (is.na(lyr[i])) {
    return(data.frame(gd = NA, ns = NA))
  }

  # get sample indices in window
  sub <- get_adj(i, lyr, nmat, coord_cells)

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



#' Rarefaction helper function
#'
#' @inheritParams window_gd_general
#'
#' @keywords internal
#'
#' @return
#' @export
#'
#' @examples
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


#' Helper function to rarify subsample and calculate genetic diversity
#'
#' @inheritParams window_gd_general
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
rarify_gd <- function(gen, sub, rarify_nit = 10, rarify_n = 4, stat = stat, fun = fun) {

  # check to make sure sub is greater than rarify_n
  if (!(length(sub) > rarify_n)) {
    stop("rarify_n is less than the number of samples provided")
  }

  # define subsample to rarify
  # (note: this combo step is done so when the number of unique combos < rarify_nit, extra calcs aren't performed)
  if (choose(length(sub), rarify_n) < rarify_nit) {
    # get all possible combos (transpose so rows are unique combos)
    cmb <- t(utils::combn(sub, rarify_n))

  } else {
    # random sample subsets of size rarify_nit (transpose so rows are unique combos)
    cmb <- t(replicate(rarify_nit, sample(sub, rarify_n), simplify = TRUE))
  }

  # for each of the possible combos get gendiv stat
  gdrar <- apply(cmb, 1, sample_gd, gen = gen, stat = stat)

  # summarize rarefaction results
  gd <- stats::na.omit(fun(gdrar))

  return(gd)
}


#' Helper function to calculate genetic diversity of a sample
#'
#' @inheritParams window_gd_general
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

#' Calculate mean allelic richness
#'
#' @param genind genind object
#'
#' @return allelic richness averaged across all loci
#' @export
#'
#' @examples
calc_mean_ar <- function(genind){
  genind$pop <- rep(factor(1), nrow(genind$tab))
  ar <- hierfstat::allelic.richness(genind)$Ar[,1]
  gd <- mean(stats::na.omit(ar))
  return(gd)
}

#' Calculate mean heterozygosity
#'
#' @param hetmat matrix of heterozygosity (0 = homozygote, 1 = heterozygote)
#'
#' @return heterozygosity averaged across all individuals and then all loci
#' @export
#'
#' @examples
calc_mean_het <- function(hetmat){
  gd_byloci <- colMeans(hetmat, na.rm = TRUE)
  gd <- stats::na.omit(mean(gd_byloci))
  return(gd)
}

#' Calculate nucleotide diversity (pi) from dosage data
#'
#' Wrapper for \link[hierfstat]{pi.dosage} function
#'
#' @param dos a ni X nl dosage matrix containing the number of derived/alternate alleles each individual carries at each SNP
#' @param L length of the sequence (*note:* defaults to number of loci in the provided dosage matrix; TODO: COME BACK AND FIX THIS)
#'
#' @return
#' @export
#'
#' @examples
calc_pi <- function(dos, L = ncol(dos)){
  gd <- hierfstat::pi.dosage(dos)
  return(gd)
}

#' Check coordinate and genetic data
#'
#' Check that the number of individuals in each data set align
#'
#' @param gen genetic data
#' @param coords coordinates
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
check_data <- function(gen, coords){

  # check number of samples
  if(class(gen) == "genind"){
    nind <- nrow(gen$tab)
  }

  if(class(gen) == "vcfR"){
    nind <- (ncol(gen@gt) - 1)
  }

  # check to make sure coords and gen align
  if (nind != nrow(coords)) {
    stop("number of samples in coords data and number of samples in gen data are not equal")
  }

}

#' Helper function to get adjacent cells to a given cell index
#'
#' @param i cell index
#' @param r RasterLayer
#' @param n neighborhood matrix
#' @param coord_cells cell numbers of coordinates
#'
#' @return indices of coordinates that are adjacent to the given cell
#' @export
#'
#' @keywords internal
#'
#' @examples
get_adj <- function(i, r, n, coord_cells){
  # get adjacent cells to cell i
  adjc <- raster::adjacent(r, i, directions = n, include = TRUE, sorted = TRUE)
  # get indices of adjacent cells
  adjci <- purrr::map_dbl(adjc, 1, function(x) {seq(x[1], x[2])})
  # get list of indices of coords in that set of cells
  sub <- which(coord_cells %in% adjci)

  return(sub)
}



