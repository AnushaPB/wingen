vario_gd <- function(x, lotr_coords, plot = TRUE){

  if (inherits(x, "vcfR")) x <- vcf_to_dosage(x)
  colnames(x) <- purrr::map_chr(colnames(x), ~ifelse(regexpr("[0-9]", .x) == 1, paste0("X", .x), .x))

  df <- data.frame(lotr_coords, x)
  sp::coordinates(df) <- ~ x + y

  f <- as.formula(paste(paste0(colnames(x), collapse = "+"), "~ 1"))
  variogram <- automap::autofitVariogram(f, df)

  if (plot) plot(variogram)

  return(list(range = variogram$var_model[2, "range"], variogram = variogram))
}
