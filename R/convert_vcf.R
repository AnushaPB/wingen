
#' Convert a vcf to a dosage matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return dosage matrix
#'
#' @export
vcf_to_dosage <- function(x) {
  # check vcf
  vcf <- vcf_check(x)

  # convert to genlight
  genlight <- vcfR::vcfR2genlight(vcf)

  # convert to dosage matrix
  gen <- as.matrix(genlight)

  return(gen)
}

#' Wrapper for vcfR2genind function that assigns pops
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#' @param pops if NULL (default), and there are no pops detected from the vcf, each individual is assigned its own pop. If FALSE then genind$pop is left NULL. Alternatively, a vector of population assignments for each individual can be provided
#' @param warning whether to retunr warning if new pops are assigned
#'
#' @return returns genind object
#'
#' @noRd
vcf_to_genind <- function(x, pops = NULL, warning = FALSE) {

  # check vcf
  vcf <- vcf_check(x)

  # convert to genind
  genind <- vcfR::vcfR2genind(vcf)

  # leave pops NULL if pops is FALSE
  if (is.logical(pops)) if (!pops) {
    return(genind)
  }

  # assign pops if null or pop vector provided
  if (is.null(pops)) genind$pop <- as.factor(1:nrow(genind@tab))
  if (is.vector(pops)) {
    if (length(pops) != nrow(genind@tab)) stop("length of pops does not match number of individuals in genind")
  }
  genind$pop <- as.factor(pops)

  return(genind)
}

#' Convert vcf to heterozygosity matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#' #'
#' @return heterozygosity matrix
#'
#' @noRd
vcf_to_het <- function(x) {
  # check vcf
  vcf <- vcf_check(x)

  het <- vcfR::is.het(vcfR::extract.gt(vcf), na_is_false = FALSE)

  # IMPORTANT: transform matrix so that rows are individuals and cols are loci
  het <- t(het)

  # if gen is a vector of only one locus, turn into matrix with one column
  if (nrow(vcf@gt) == 1) {
    het <- matrix(het, ncol = 1)
  }

  return(het)
}

#' Check if an object is a vcf or a path to a vcf
#'
#' @param x vcfR object or path to vcf
#'
#' @return vcf object
#'
#' @noRd
vcf_check <- function(x) {
  if (class(x)[1] == "vcfR") {
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
