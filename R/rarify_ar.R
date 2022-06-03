
#' Rarefaction function
#'
#' @param ar_df dataframe of allelic richness (*note:* rows should be individuals, columns should be loci)
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
