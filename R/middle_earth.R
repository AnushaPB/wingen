#' Load middle earth example
#'
#' Loads middle earth example data and assigns to simple names
#'
#' @return three objects are loaded (lotr_vcf, lotr_coords, and lotr_lyrs)
#' @export
#'
#' @examples
load_middle_earth_ex <- function() {

  # load all data
  utils::data("lotr_vcf")
  utils::data("lotr_lyr")
  utils::data("lotr_coords")
  utils::data("lotr_range")

  # give message with information about objects
  return(message(cat(
    crayon::cyan(crayon::bold("\n-------------- middle earth example --------------\n")),
    crayon::silver("\nObjects loaded:"),
    crayon::yellow(crayon::bold("\n*lotr_vcf*")),
    crayon::yellow(paste0("vcfR object (", nrow(lotr_vcf@gt), " loci x ", ncol(lotr_vcf@gt) - 1, " samples)")),
    crayon::green(crayon::bold("\n*lotr_coords*")), crayon::green("dataframe with x and y coordinates"),
    crayon::magenta(crayon::bold("\n*lotr_lyr*")), crayon::magenta("middle earth RasterLayer (100 x 100)"),
    crayon::blue(crayon::bold("\n*lotr_range*")), crayon::blue("SpatialPolygonsDataFrame of spp range"),
    crayon::cyan(crayon::bold("\n\n--------------------------------------------------\n"))
  )))
}

#' Load mini middle earth example
#'
#' Loads mini middle earth example data and assigns to simple names
#'
#' @return three objects are assigned in the GlobalEnv (vcf, coords, and lyr)
#' @export
#'
#' @examples
load_mini_ex <- function() {

  # load all data
  utils::data("mini_vcf", envir = environment())
  utils::data("mini_coords", envir = environment())
  utils::data("mini_lyr", envir = environment())

  # give message with information about objects
  message(cat(
    crayon::cyan(crayon::bold("\n---------- mini middle earth example ----------\n")),
    crayon::blue("\nObjects loaded:"),
    crayon::yellow(crayon::bold("\n*mini_vcf*")), crayon::yellow("vcfR object (10 loci x 10 samples)"),
    crayon::green(crayon::bold("\n*mini_coords*")), crayon::green("dataframe with x and y coordinates"),
    crayon::magenta(crayon::bold("\n*mini_lyr*")), crayon::magenta("middle earth RasterLayer (10 x 10)"),
    crayon::cyan(crayon::bold("\n\n-----------------------------------------------"))
  ))
}
