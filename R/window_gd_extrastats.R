
calc_prop_hwe <- function(genind, sig = 0.05){
  hwe <- pegas::hw.test(genind)
  prop <- mean(hwe[, "Pr.exact"] < sig, na.rm = TRUE)
  return(prop)
}

calc_mean_basic_stats <- function(hf){
  hfstat <- hierfstat::basic.stats(hf)
  mean_stats <- hfstat$overall
  names(mean_stats) <- paste0(names(mean_stats), "_hierfstat")
  return(mean_stats)
}
