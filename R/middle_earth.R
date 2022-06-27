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

  # give message with information about objects
  message(cat(crayon::cyan(crayon::bold("\n------------ middle earth loaded -------------\n")),
              crayon::blue("\nAdded to GlobalEnv:"),
              crayon::yellow(crayon::bold("\nvcf")), crayon::yellow("vcfR object with (1000 loci x 200 samples)"),
              crayon::green(crayon::bold("\ncoords")), crayon::green("dataframe with x and y coordinates"),
              crayon::magenta(crayon::bold("\nlyr")), crayon::magenta("middle earth RasterLayer (100 x 100)"),
              crayon::cyan(crayon::bold("\n\n----------------------------------------------"))
              )
          )
}

#' Load mini middle earth example
#'
#' Loads mini middle earth example data and assigns to simple names
#'
#' @return three objects are assigned in the GlobalEnv (vcf, coords, and lyr)
#' @export
#'
#' @examples
load_mini_ex <- function(){

  # load all data
  data("mini_vcf")
  data("mini_coords")
  data("mini_lyr")

  # assign data to simpler object names
  assign("vcf", mini_vcf, envir = .GlobalEnv)
  assign("coords", mini_coords, envir = .GlobalEnv)
  assign("lyr", mini_lyr, envir = .GlobalEnv)

  # give message with information about objects
  message(cat(crayon::cyan(crayon::bold("\n-------- mini middle earth loaded ---------\n")),
              crayon::blue("\nAdded to GlobalEnv:"),
              crayon::yellow(crayon::bold("\nvcf")), crayon::yellow("vcfR object with (10 loci x 10 samples)"),
              crayon::green(crayon::bold("\ncoords")), crayon::green("dataframe with x and y coordinates"),
              crayon::magenta(crayon::bold("\nlyr")), crayon::magenta("middle earth RasterLayer (10 x 10)"),
              crayon::cyan(crayon::bold("\n\n-------------------------------------------"))
              ))
}
