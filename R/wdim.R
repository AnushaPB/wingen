
wdim_to_mat <- function(wdim) {
  if (any(wdim < 3)) {
    stop("wdim cannot be less than 3")
  }

  if (length(wdim) == 2) {
    n <- matrix(1, wdim[1], wdim[2])
  } else if (length(wdim) == 1) {
    n <- matrix(1, wdim, wdim)
  } else {
    stop("wdim must be a number or a vector of length 2")
  }

  # focal cell (center of matrix) has to be zero
  n[wdim / 2 + 0.5, wdim / 2 + 0.5] <- 0

  return(n)
}

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
