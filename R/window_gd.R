

#' Sliding window map of genetic diversity
#'
#' @param vcf object of type vcf ( (*note:* order matters! the coordinate and genetic data should be in the same order, there are currently no checks for this.))
#' @param stat genetic diversity stat to calculate (can either be "pi" for nucleotide diversity, "het" for average heterozygosity across all loci, "allelic.richness" for average allelic richness across all loci, or "biallelic.richness" to get average allelic richness across all loci for a biallelic dataset (this option faster than "allelic.richness"))
#'
#' @inheritParams window_gd_general
#' @return RasterStack that includes a raster of genetic diversity and a raster of the number of samples within the window for each cell
#' @export
#'
#' @examples
#' library("raster")
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, nloci = 10, rarify = TRUE)
#' plot_gd(wpi, main = "Window pi")
#' plot_count(wpi)
#'
window_gd <- function(vcf, coords, lyr, stat = "pi", wdim = 5, fact = 0, rarify = FALSE, rarify_n = 4, rarify_nit = 5, min_n = 2, fun = mean, parallel = FALSE, nloci = NULL){
  # check that the input file is a vcf or a path to a vcf object
  if(class(vcf) != "vcfR" & is.character(vcf)){
    vcf <- vcfR::read.vcfR(vcf)
  } else if (class(vcf) != "vcfR" & !is.character(vcf)){
    stop("gen object must be of type vcfR or a path to a .vcf files")
  }

  # check to make sure coords and gen align
  check_data(vcf, coords)

  # calc stats
  if(stat == "allelic.richness"){
    # convert from vcf to genind
    gen <- vcfR::vcfR2genind(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_ar, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel)

    names(results[[1]]) <- "allelic_richness"

  }

  if(stat == "het" | stat == "heterozygosity"){
    # convert from vcf to heterozygosity matrix
    gen <- vcfR::is.het(vcfR::extract.gt(vcf), na_is_false = FALSE)

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_het, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel)

    names(results[[1]]) <- "heterozygosity"
  }

  if(stat == "pi"){
    # convert from vcf to dosage matrix
    gen <- vcf_to_dosage(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = calc_pi, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel, nloci)

    names(results[[1]]) <- "pi"
   }

  if(stat == "biallelic.richness"){
    #convert vcf to dosage matrix
    gen <- vcf_to_dosage(vcf)

    results <- window_gd_general(gen, coords, lyr, stat = calc_mean_biar, wdim, fact, rarify, rarify_n, rarify_nit, min_n, fun, parallel)

    names(results[[1]]) <- "biallelic_richness"

  }

  names(results[[2]]) <- "sample_count"

  return(results)

}
#' Helper function for window_gd
#'
#' @param gen genetic data (*note:* order matters! the coordinate and genetic data should be in the same order, there are currently no checks for this.)
#' @param coords coordinates (two columns, the first should be x and the second should be y and the order should be the same as the genetic data),
#' @param lyr RasterLayer to slide window across
#' @param stat function to calculate genetic diversity (can either be calc_mean_arcalc_pi, calc_mean_biar, or calc_mean_het)
#' @param wdim dimensions (height x width) of window, if only one value is provided a square window is created
#' @param fact aggregation factor to apply to the RasterLayer (*note:* increasing this value reduces computational time)
#' @param rarify if rarify = TRUE, rarefaction is performed
#' @param rarify_n if rarify = TRUE, number of points to use for rarefaction
#' @param rarify_nit if rarify = TRUE, number of iterations to use for rarefaction
#' @param min_n min number of samples to calculate allelic richness for (equal to rarify_n if provided, otherwise defaults to 2)
#' @param fun function to use to summarize data in window (defaults to base R mean)
#' @param parallel whether to parallelize the function (see vignette for setting up a cluster to do so)
#' @param nloci for calculating pi, if nloci=NULL (default), returns the sum over SNPs of nucleotide diversity; otherwise return the average nucleotide diversity per nucleotide given the length nloci of the sequence (L argument in \link[hierfstat]{pi.dosage} function)
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
window_gd_general <- function(gen, coords, lyr, stat = calc_mean_ar, wdim = 5, fact = 0, rarify = FALSE, rarify_n = 4, rarify_nit = 10, min_n = 2, fun = mean, parallel = FALSE, nloci = NULL) {

  # TODO: ADD FUNCTIONALITY SO RARIFY CAN EQUAL 1

  # make neighborhood matrix for window
  nmat <- wdim_to_mat(wdim)

  # make aggregated raster
  if(fact == 0){lyr <- lyr * 0} else {lyr <- raster::aggregate(lyr, fact) * 0}

  # get cell index for each coordinate
  coord_cells <- raster::extract(lyr, coords, cell = TRUE)[,"cells"]

  # ignore this: (need to assign i something so that R CMD Check recognizes it as a defined global variable - also this is useful for testing)
  i <- 1

  if(parallel){

    rast_vals <- foreach::foreach(i = 1:raster::ncell(lyr), .combine = rbind, .packages = c("raster", "purrr", "hierfstat", "stats", "adegenet")) %dopar% {

      result <- window_helper(i, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, nloci)

      return(result)

    }

  } else {

    rast_vals <- purrr::map_dfr(1:raster::ncell(lyr), window_helper, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, nloci)

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
window_helper <- function(i, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, nloci = NULL){
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
    gd <- rarify_helper(gen, sub, rarify_n, rarify_nit, stat, fun, nloci)
  } else {
    gd <- sample_gd(gen, sub, stat, nloci)
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
rarify_helper <- function(gen, sub, rarify_n, rarify_nit, stat, fun = mean, nloci = NULL){
    # if number of samples is less than rarify_n, assign the value NA
    if (length(sub) < rarify_n) {
      gd <- NA
    }

    # if number of samples is greater than rarify_n, rarify
    if (length(sub) > rarify_n) {
      gd <- rarify_gd(gen, sub, rarify_nit = rarify_nit, rarify_n = rarify_n, stat = stat, fun = fun, nloci = nloci)
    }

    # if the number of samples is equal to rarify_n, calculate stat
    if (length(sub) == rarify_n) {
      gd <- sample_gd(gen, sub, stat, nloci)
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
rarify_gd <- function(gen, sub, rarify_nit = 10, rarify_n = 4, stat, fun, nloci = NULL) {

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
  gdrar <- apply(cmb, 1, sample_gd, gen = gen, stat = stat, nloci = nloci)

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
sample_gd <- function(gen, sub, stat, nloci = NULL) {
  if(is.null(nloci)){gd <- stat(gen[sub,])} else {gd <- stat(gen[sub,], nloci)}
  return(gd)
}


#' Calculate mean allelic richness
#'
#' @param genind genind
#'
#' @return allelic richness averaged across all loci
#' @export
#'
#' @keywords internal
#'
#' @examples
calc_mean_ar <- function(genind){
  genind$pop <- rep(factor(1), nrow(genind$tab))
  #note [,1] references the first column which is AR for each locus across all inds (nrow(AR) == nloci)
  ar <- hierfstat::allelic.richness(genind)$Ar[,1]
  gd <- mean(stats::na.omit(ar))
  return(gd)
}

#' Calculate mean heterozygosity
#'
#' @param hetmat matrix of heterozygosity (0/FALSE = homozygote, 1/TRUE = heterozygote)
#'
#' @return heterozygosity averaged across all individuals and then all loci
#' @export
#'
#' @keywords internal
#'
#' @examples
calc_mean_het <- function(hetmat){
  # if hetmat is not a matrix (e.g. has NULL dimensions/is a vector or has one value, calculate and return the mean)
  if(is.null(dim(hetmat)) & length(hetmat) > 0){
    return(stats::na.omit(mean(hetmat)))
  }

  gd <- mean(hetmat, na.rm = TRUE)

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
#' @keywords internal
#'
#' @examples
calc_pi <- function(dos, nloci = NULL){
  gd <- hierfstat::pi.dosage(dos, L = nloci)
  return(gd)
}

#' Calculate mean allelic richness for biallelic data
#'
#' @param dos dosage matrix
#'
#' @return allelic richness averaged across all loci
#' @export
#'
#' @keywords internal
#'
#' @examples
calc_mean_biar <- function(dos){
  if(!all(dos %in% c(0,1,2))){stop("to calculate biallelic richness, all values in genetic matrix must be 0, 1 or 2")}
  ar_by_locus <- apply(dos, 2, helper_calc_biar)
  mean_ar <- mean(stats::na.omit(ar_by_locus))
  return(mean_ar)
}

#' Helper function to calculate allelic richness for a biallelic locus
#'
#' @param loc genotypes at a biallelic locus (must have values of 0, 1, or 2)
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
helper_calc_biar <- function(loc){
  uq <- stats::na.omit(unique(loc))
  if (1 %in% uq){
    return(2)
  } else if (0 %in% uq & 2 %in% uq){
    return(2)
  } else {
    return(1)
  }
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

preview_window <- function(lyr, coords, wdim, fact = 0, sample_count = TRUE, min_n = 0){
  if(fact != 0) lyr <- aggregate(lyr, fact)

  # convert wdim to matrix
  nmat <- wdim_to_mat(wdim)

  # get center of raster
  e <- as.vector(extent(lyr))
  c <- c(mean(e[c(1,2)]),mean(e[c(3,4)]))
  center <- cellFromXY(lyr, c)

  # get adjacent cells to center cell
  adjc <- raster::adjacent(lyr, center, directions = nmat)
  # get list of indices of coords in that set of cells
  adjci <- purrr::map_dbl(adjc, 1, function(x) {seq(x[1], x[2])})
  # fill in window
  lyrw <- lyr*0
  lyrw[adjci] <- 1
  lyrw[center] <- 2

  raster::plot(lyrw, col = mako(3, direction = -1), legend = FALSE, axes = FALSE, box = FALSE)
  legend("bottomleft", c("raster layer", "window", "focal cell"), col = mako(3, direction = -1), pch = 15)
  if(!is.null(coords)) points(coords, pch = 3)

  if(sample_count){

    # get coord cells
    coord_cells <- raster::extract(lyr, coords, cell = TRUE)[,"cells"]

    # count
    lyrc <- lyr
    nc <- map_dbl(1:ncell(lyr), function(x, lyr, nmat, coord_cells){sub <- get_adj(x, lyr, nmat, coord_cells); return(length(sub))}, lyr, nmat, coord_cells)
    lyrc <- setValues(lyr, nc)

    lyrc[lyrc < min_n] <- NA
    raster::plot(lyrc, col = mako(100), box = FALSE, axes = FALSE, main = "sample count")
  }
}
