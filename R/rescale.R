
#' Rescale values from zero to one
#'
#' @param x vector of numbers
#' @param x.min current min
#' @param x.max current max
#' @param new.min new min
#' @param new.max new max
#'
#' @return values rescales to new min and max
#' @export
#'
#' @examples
rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x)
  if(is.null(x.max)) x.max = max(x)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}
