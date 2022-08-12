#' Load middle earth example
#'
#' Loads middle earth example data and assigns to simple names
#'
#' @param quiet whether to hide message (defaults to FALSE)
#'
#' @return three objects are loaded (lotr_vcf, lotr_coords, and lotr_lyrs)
#' @export
#'
#' @examples
#' load_middle_earth_ex()
load_middle_earth_ex <- function(quiet = FALSE) {

  # load all data
  utils::data(list = c("lotr_vcf", "lotr_lyr", "lotr_coords", "lotr_range"))

  if (!quiet) {
    # give message with information about objects
    return(message(cat(
      crayon::cyan(crayon::bold("\n-------------- middle earth example --------------\n")),
      crayon::silver("\nObjects loaded:"),
      crayon::yellow(crayon::bold("\n*lotr_vcf*")),
      crayon::yellow(paste0("vcfR object (100 loci x 100 samples)")),
      crayon::green(crayon::bold("\n*lotr_coords*")), crayon::green("dataframe with x and y coordinates"),
      crayon::magenta(crayon::bold("\n*lotr_lyr*")), crayon::magenta("middle earth RasterLayer (100 x 100)"),
      crayon::blue(crayon::bold("\n*lotr_range*")), crayon::blue("SpatialPolygonsDataFrame of spp range"),
      crayon::cyan(crayon::bold("\n\n--------------------------------------------------\n"))
    )))
  }
}

#' Load mini middle earth example
#'
#' Loads mini middle earth example data and assigns to simple names
#'
#' @param quiet whether to hide message (defaults to FALSE)
#'
#' @return three objects are assigned in the GlobalEnv (vcf, coords, and lyr)
#' @export
#'
#' @examples
#' load_mini_ex()
load_mini_ex <- function(quiet = FALSE) {

  # load all data
  utils::data(list = c("mini_vcf", "mini_coords", "mini_lyr"))

  # give message with information about objects
  if (!quiet) {
    message(cat(
      crayon::cyan(crayon::bold("\n---------- mini middle earth example ----------\n")),
      crayon::blue("\nObjects loaded:"),
      crayon::yellow(crayon::bold("\n*mini_vcf*")), crayon::yellow("vcfR object (10 loci x 10 samples)"),
      crayon::green(crayon::bold("\n*mini_coords*")), crayon::green("dataframe with x and y coordinates"),
      crayon::magenta(crayon::bold("\n*mini_lyr*")), crayon::magenta("middle earth RasterLayer (10 x 10)"),
      crayon::cyan(crayon::bold("\n\n-----------------------------------------------"))
    ))
  }
}
