
calc_prop_hwe <- function(genind, sig = 0.05){
  hwe <- pegas::hw.test(genind)
  prop <- mean(hwe[, "Pr.exact"] < sig, na.rm = TRUE)
  return(prop)
}

calc_mean_fis <- function(hf){
  hfstat <- hierfstat::basic.stats(hf)
  # check why this is different from overall when NA values are present
  #Fis <- mean(hfstat$Fis[,1], na.rm = TRUE)
  Fis <- hfstat$overall["Fis"]
  return(Fis)
}
