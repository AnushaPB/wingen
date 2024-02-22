#' Middle earth example
#'
#' Loads middle earth example data
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
    return(message(
      crayon::cyan(crayon::bold("\n-------------- middle earth example --------------\n")),
      crayon::silver("\nObjects loaded:"),
      crayon::yellow(crayon::bold("\nlotr_vcf: ")),
      crayon::yellow(paste0("vcfR object (100 variants x 100 samples)")),
      crayon::green(crayon::bold("\nlotr_coords: ")), crayon::green("dataframe with x and y coordinates"),
      crayon::magenta(crayon::bold("\nlotr_lyr: ")), crayon::magenta("middle earth RasterLayer (100 x 100)"),
      crayon::blue(crayon::bold("\nlotr_range: ")), crayon::blue("SpatialPolygonsDataFrame of spp range"),
      crayon::cyan(crayon::bold("\n\n--------------------------------------------------\n"))
    ))
  }
}

#' Mini middle earth example
#'
#' Loads mini middle earth example data
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
  utils::data(list = c("mini_vcf", "mini_vcf_NA", "mini_coords", "mini_lyr"))

  # give message with information about objects
  if (!quiet) {
    message(
      crayon::cyan(crayon::bold("\n-------------------- mini middle earth example --------------------\n")),
      crayon::blue("\nObjects loaded:"),
      crayon::yellow(crayon::bold("\nmini_vcf: ")), crayon::yellow("vcfR object (10 variants x 10 samples) and no missing data"),
      crayon::yellow(crayon::bold("\nmini_vcf_NA: ")), crayon::yellow("vcfR object (10 variants x 10 samples) and missing data"),
      crayon::green(crayon::bold("\nmini_coords: ")), crayon::green("dataframe with x and y coordinates"),
      crayon::magenta(crayon::bold("\nmini_lyr: ")), crayon::magenta("middle earth RasterLayer (10 x 10)"),
      crayon::cyan(crayon::bold("\n\n-------------------------------------------------------------------"))
    )
  }
}
