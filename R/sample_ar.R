
#' Calculate allelic richness of a sample
#'
#' @param ar_df dataframe of allelic richness (*note:* rows should be individuals, columns should be loci)
#' @param sub row indices of subsample
#' @param meanf function to use to calculate the mean (defaults to base R mean)
#'
#' @return mean allelic richness of a subsample
#' @export
#'
#' @examples
sample_ar <- function(ar_df, sub, meanf = mean){
  ar <- ar_df[sub,] %>% colMeans(na.rm = TRUE) %>% na.omit() %>% meanf()
  return(ar)
}
