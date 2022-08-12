
#' Runtime evaluation
#'
#' @param val parameter value
#' @param var parameter
#' @inheritParams window_gd
#'
#' @return
#' @export
#'
#' @examples
time_eval <- function(val, var, vcf, coords, lyr, stat = "pi", wdim = 3, fact = 0, fun = mean, parallel = FALSE){
  # reassign argument
  assign(var, val)

  # start timing
  ptm <- Sys.time()

  # run moving window
  # note: min_n needs to be 1 so that landscape/sampling doesn't affect the results (e.g. no NA cells)
  window_gd(vcf,
            coords,
            lyr,
            stat,
            wdim,
            fact,
            rarify = FALSE,
            min_n = 1,
            fun,
            parallel)

  # create aggregated layer for calculating ncell
  if(fact != 0) lyr <- aggregate(lyr, fact)

  # make df of results
  df <- data.frame(time = as.numeric(Sys.time() - ptm, units = "secs"),
                   wsize = wdim*wdim,
                   fact = fact,
                   wdim = wdim,
                   ncell = ncell(lyr))

  return(df)
}

#' Run multiple iterations of time_eval
#'
#' @param x iteration number
#' @param vals values to test
#' @inheritParams time_eval
#' @return
#' @export
#'
#' @examples
time_eval_its <- function(x, vals, var, vcf, coords, lyr){
  df <- purrr::map_dfr(vals,
                       time_eval,
                       var,
                       vcf,
                       coords,
                       lyr)
  df$it <- x
  return(df)
}
