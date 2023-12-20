#' Calculate mean allelic richness
#'
#' @param genind genind
#'
#' @return allelic richness averaged across all loci
#'
#' @noRd
calc_mean_ar <- function(genind) {
  ar <- helper_calc_ar(genind)
  gd <- mean(ar, na.rm = TRUE)
  names(gd) <- "allelic_richness"
  return(gd)
}

#' Helper function to calculate allelic richness
#'
#' @param genind genind object
#'
#' @return allelic richness
#'
#' @noRd
helper_calc_ar <- function(genind) {
  # get number of individuals
  nind <- nrow(genind@tab)

  # assign pops so that the whole sample is treated as one pop
  genind$pop <- rep(factor(1), nind)

  # note: [,1] references the first column which is AR for each site across all inds (nrow(AR) == L)
  # note: these are rarified allele counts (to make them not rarified set min.n = nind*2)
  # min.n is the The number of alleles down to which the number of alleles should be rarefied.
  # The default is the minimum number of individuals genotyped (times 2 for diploids). However, if there
  # are NA values then it doesn't count those as genotypes.
  ar <- hierfstat::allelic.richness(genind)$Ar[, 1]
  return(ar)
}


#' Calculate mean heterozygosity
#'
#' @param hetmat matrix of heterozygosity (0/FALSE = homozygote, 1/TRUE = heterozygote)
#'
#' @return heterozygosity averaged across all individuals then all loci
#'
#' @noRd
calc_mean_het <- function(hetmat) {
  if (is.null(dim(hetmat))) {
    gd <- mean(hetmat, na.rm = TRUE)
  } else {
    het_by_locus <- colMeans(hetmat, na.rm = TRUE)
    gd <- mean(het_by_locus, na.rm = TRUE)
  }
  names(gd) <- "Ho"
  return(gd)
}

#' Calculate nucleotide diversity (pi) from dosage data
#'
#' Wrapper for \link[hierfstat]{pi.dosage} function
#'
#' @param dos a ni X nl dosage matrix containing the number of derived/alternate alleles each individual carries at each SNP
#' @param L length of the sequence (*note:* defaults to NULL which returns the sum over SNPs of nucleotide diversity)
#'
#' @return nucleotide diversity (pi)
#'
#' @noRd
calc_pi <- function(dos, L = NULL) {
  gd <- hierfstat::pi.dosage(dos, L = L)
  names(gd) <- "pi"
  return(gd)
}

#' Calculate mean allelic richness for biallelic data
#'
#' @param dos dosage matrix
#'
#' @return allelic richness averaged across all loci
#'
#' @noRd
calc_mean_biar <- function(dos, rarify_alleles = TRUE) {
  if (!all(dos %in% c(0, 1, 2, NA))) {
    stop("to calculate biallelic richness, all values in genetic matrix must be NA, 0, 1 or 2")
  }

  # if rarify_alleles = TRUE, get min.n for rarefaction
  # min.n is set to the number of non-NA genotypes * 2 (for diploid)
  if (rarify_alleles) {
    min.n <- get_minn(dos)
  } else {
    min.n <- NULL
  }

  # if null dimensions (e.g., only one value provided), use helper_calc_biar directly
  # otherwise apply across columns (i.e., loci)
  if (is.null(dim(dos))) {
    ar_by_locus <- sapply(dos, helper_calc_biar, rarify_alleles, min.n)
  } else {
    ar_by_locus <- apply(dos, 2, helper_calc_biar, rarify_alleles, min.n)
  }

  gd <- mean(ar_by_locus, na.rm = TRUE)
  names(gd) <- "biallelic_richness"
  return(gd)
}

#' Helper function to calculate allelic richness for a biallelic locus
#'
#' @param loc genotypes at a biallelic locus (must have values of 0, 1, or 2)
#'
#' @return biallelic richness value
#'
#' @noRd
helper_calc_biar <- function(loc, rarify_alleles = TRUE, min.n = NULL) {
  # check if all NA and if so return NA
  if (all(is.na(loc))) {
    return(NA)
  }

  if (rarify_alleles) {
    # omit NAs before counting alleles
    loc_NArm <- stats::na.omit(loc)

    # make df of counts of reference and alternate
    counts <- c(R = sum(loc_NArm), A = 2 * length(loc_NArm) - sum(loc_NArm, na.rm = TRUE))

    # rarefied counts calculation taken from hierfstat::allelic.richness function
    AR <- raref(counts, min.n = min.n)
  } else {
    # calculate number of unique alleles
    # note: has to be na.omit (na.rm is not an argument for unique)
    uq <- unique(stats::na.omit(loc))
    if (1 %in% uq) {
      AR <- 2
    } else if (0 %in% uq & 2 %in% uq) {
      AR <- 2
    } else {
      AR <- 1
    }
  }

  return(AR)
}

#' Rarify allele counts (based on \link[hierfstat]{allelic.richness}code)
#'
#' @param x allele counts
#' @param min.n the number of alleles down to which the number of alleles should be rarefied.
#'
#' @noRd
raref <- function(x, min.n) {
  nn <- sum(x)
  dum <- exp(lchoose(nn - x, min.n) - lchoose(nn, min.n))
  dum[is.na(dum)] <- 0
  return(sum(1 - dum))
}


#' Calculate min.n
#'
#' @param dos dosage matrix
#'
#' @noRd
get_minn <- function(dos) {
  if (is.null(nrow(dos))) {
    # if nrow is NULL then only one sample is included so min.n must be 2
    min.n <- 2
  } else {
    min.n <- 2 * min(apply(dos, 2, countgen), na.rm = TRUE)
  }
  return(min.n)
}

#' Count not NA genotypes and return NA if all are NA
#'
#' @param x dosage matrix
#'
#' @noRd
countgen <- function(x) {
  notNA <- sum(!is.na(x))
  if (notNA == 0) notNA <- NA
  return(notNA)
}


#' Calculate proportion of snps out of HWE using \link[pegas]{hw.test}
#'
#' @param genind genind object
#' @param sig significance level for HWE test
#'
#' @return proportion of snps out of HWE
#'
#' @noRd
calc_prop_hwe <- function(genind, sig = 0.05) {
  hwe <- pegas::hw.test(genind)
  prop <- mean(hwe[, "Pr.exact"] < sig, na.rm = TRUE)
  names(prop) <- "hwe"
  return(prop)
}

#' Calculate basic stats using \link[hierfstat]{basic.stats}
#'
#' @param hf hierfstat object
#' @param ... additional arguments to pass to basic stats
#'
#' @return vector of overall stats produced by \link[hierfstat]{basic.stats}
#'
#' @noRd
calc_mean_basic_stats <- function(hf) {
  # reassign pop so that if the levels present in the original factor change you don't get an error
  # e.g., if original pop was 1 & 2, but a subset only has 2, you will get an error
  hf$pop <- as.numeric(as.character(hf$pop))
  hfstat <- hierfstat::basic.stats(hf)
  mean_stats <- hfstat$overall

  # drop stats that aren't meaningful if there is only one populatoin
  mean_stats <- mean_stats[c("Ho", "Hs", "Ht", "Fis")]

  names(mean_stats) <- paste0(names(mean_stats), "_hierfstat")

  return(mean_stats)
}

#' Helper function to get genetic diversity functions
#'
#' @param stat moving window statistic to calculate (can either be `pi` for nucleotide diversity, `Ho` for average observed heterozygosity across all loci, "allelic_richness" for average number of alleles across all loci, "biallelic_richness" to get average number of alleles across all loci for a biallelic dataset. `stat` can also be set to functions that will return a single numeric value from the input data (for example a summary statistic like `mean`, `sum`, or `sd`)
#' @param ... if a function is provided for `x`, additional arguments to pass to the `x` function (e.g. if `x = mean`, users may want to set `na.rm = TRUE`)
#'
#' @return function corresponding with desired statistic
#'
#' @noRd
return_stat <- function(stat, ...) {
  if (inherits(stat, "function")) {
    return(purrr::partial(stat, ...))
  }

  if (stat == "pi") {
    return(calc_pi)
  }

  if (stat == "biallelic_richness") {
    return(calc_mean_biar)
  }

  if (stat == "allelic_richness") {
    return(calc_mean_ar)
  }

  if (stat == "Ho") {
    return(calc_mean_het)
  }

  if (stat == "hwe") {
    return(calc_prop_hwe)
  }

  if (stat == "basic_stats") {
    return(calc_mean_basic_stats)
  }

  stop(paste(stat, "is an invalid argument for stat"))
}
