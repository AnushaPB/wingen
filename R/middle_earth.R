#' Load middle earth example
#'
#' Loads middle earth example data and assigns to simple names
#'
#' @return three objects are assigned in the GlobalEnv (vcf, coords, and lyr)
#' @export
#'
#' @examples
load_middle_earth_ex <- function(){

  # load all data
  data("middle_earth_vcf")
  data("middle_earth_lyr")
  data("middle_earth_coords")

  # assign data to simpler object names
  assign("vcf", middle_earth_vcf, envir = .GlobalEnv)
  assign("coords", middle_earth_coords, envir = .GlobalEnv)
  assign("lyr", middle_earth_lyr, envir = .GlobalEnv)

  return(message("middle earth loaded..."))
}
