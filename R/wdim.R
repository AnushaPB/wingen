#' Helper function to convert wdim object to neighbor matrix
#'
#' @param wdim dimensions (height x width) of window, if only one value is provided a square window is created
#'
#' @return neighborhood matrix
#'
#' @noRd
wdim_to_mat <- function(wdim) {
  if (length(wdim) == 2) {
    n <- matrix(1, wdim[1], wdim[2])
    # focal cell (center of matrix) has to be zero
    n[wdim[1] / 2 + 0.5, wdim[2] / 2 + 0.5] <- 0
  } else if (length(wdim) == 1) {
    n <- matrix(1, wdim, wdim)
    # focal cell (center of matrix) has to be zero
    n[wdim / 2 + 0.5, wdim / 2 + 0.5] <- 0
  } else {
    stop("wdim must be a single integer or a vector two integers")
  }


  return(n)
}

#' Helper function to check that wdim object is correctly assigned
#'
#' @param wdim dimensions (height x width) of window, if only one value is provided a square window is created
#'
#' @return corrected wdim
#'
#' @noRd
wdim_check <- function(wdim) {
  if (any(wdim < 3)) {
    stop("wdim cannot be less than 3")
  }

  if (length(wdim) == 1) {
    if (wdim %% 2 == 0) {
      wdim <- wdim + 1
      warning(paste("wdim must be odd, using wdim =", wdim, "instead"))
    }
  }
  if (length(wdim) == 2) {
    if (wdim[1] %% 2 == 0) {
      wdim[1] <- wdim[1] + 1
      warning(paste("wdim must be odd, using wdim[1] =", wdim[1], "instead"))
    }

    if (wdim[2] %% 2 == 0) {
      wdim[2] <- wdim[2] + 1
      warning(paste("wdim must be odd, using wdim[2] =", wdim[2], "instead"))
    }
  }

  return(wdim)
}
