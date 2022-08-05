
#' Convert vcf to dosage matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return returns dosage matrix
#' @export
#'
#' @keywords internal
#'
#' @examples
#' data("mini_vcf")
#' vcf_to_dosage(mini_vcf)
#'
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
#'
#' @return returns genind object
#' @export
#'
#' @keywords internal
#'
#' @examples
#' data("mini_vcf")
#' vcf_to_genind(mini_vcf)
#'
vcf_to_genind <- function(x, pops = NULL) {

  # check vcf
  vcf <- vcf_check(x)

  # convert to genind
  genind <- vcfR::vcfR2genind(vcf)

  # TODO: Clean this up - Check if pops is false
  if (is.logical(pops)) if (!pops) {
    return(genind)
  }

  # assign pops if null or pop vector provided
  if (is.null(genind$pop) | is.vector(pops)) {
    if (is.null(pops)) {
      genind$pop <- as.factor(1:nrow(genind@tab))
      warning("no pops were provided, assigning a pop to each individual (to stop this, set pops = FALSE)")
    }

    if (is.vector(pops)) {
      if (length(pops) != nrow(genind@tab)) {
        stop("length of pops does not match number of individuals in genind")
      } else {
        genind$pop <- as.factor(pops)
      }
    }
  }

  return(genind)
}

#' Check if an object is a vcf or a path to a vcf
#'
#' @param x vcfR object or path to vcf
#'
#' @return vcf object
#' @export
#'
#' @keywords internal
#'
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
