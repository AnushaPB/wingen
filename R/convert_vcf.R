#' Convert a vcf to a dosage matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return dosage matrix
#'
#' @export
vcf_to_dosage <- function(x) {
  # convert to genlight
  genlight <- vcfR::vcfR2genlight(x)

  # convert to dosage matrix
  gen <- as.matrix(genlight)

  return(gen)
}


#' Convert vcf to heterozygosity matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return heterozygosity matrix
#'
#' @noRd
vcf_to_het <- function(x) {
  het <- vcfR::is.het(vcfR::extract.gt(x), na_is_false = FALSE)

  # IMPORTANT: transform matrix so that rows are individuals and cols are loci
  het <- t(het)

  # if gen is a vector of only one site, turn into matrix with one column
  if (nrow(x@gt) == 1) {
    het <- matrix(het, ncol = 1)
  }

  return(het)
}

#' Convert vcf to hierfstat
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#' @param pop population assignments (defaults to 1 for moving window functions)
#' @return hierfstat object
#'
#' @noRd
vcf_to_hf <- function(x, pop = 1) {
  genind <- vcfR::vcfR2genind(x)
  hf <- hierfstat::genind2hierfstat(genind, pop = pop)
  return(hf)
}

#' Check if an object is a vcf or a path to a vcf
#'
#' @param x vcfR object or path to vcf
#'
#' @return vcf object
#'
#' @noRd
vcf_check <- function(x) {
  if (inherits(x, "vcfR")) {
    vcf <- x
  } else if (is.character(x)) {
    if (file.exists(x)) {
      vcf <- vcfR::read.vcfR(x)
    } else {
      stop("Cannot open file: No such file or directory")
    }
  } else {
    stop("Input is expected to be an object of class 'vcfR' or a path to a .vcf file")
  }

  return(vcf)
}
