#' Core function used by \link[wingen]{window_gd}, \link[wingen]{circle_gd}, and \link[wingen]{resist_gd}
#'
#' @return SpatRaster of genetic diversity and sample counts
#'
#' @noRd
run_general <- function(x, lyr, coords,
                        coord_cells = NULL, nmat = NULL,
                        distmat = NULL, maxdist = NULL,
                        stat, rarify, rarify_n, rarify_nit, min_n,
                        fun, L, rarify_alleles, sig,
                        ...) {

  # check that any stats will be calculated
  counts <- preview_count(lyr = lyr, coords = coords, distmat = distmat, nmat = nmat, min_n = min_n, plot = FALSE)
  if (all(is.na(terra::values(counts)))) stop("Minimum sample size (min_n) is not met for any window across this raster")

  # set L if pi is being calculated
  if (is.character(stat) & !is.null(L)) if (stat == "pi" & L == "nvariants") L <- ncol(x)

  # Get function to calculate the desired statistic
  stat_function <- return_stat(stat = stat, ...)

  # check that coords and x align and reformat data, if necessary
  corrected_data <- check_data(x, coords = coords, distmat = distmat)
  x <- corrected_data$x
  coords <- corrected_data$coords
  distmat <- corrected_data$distmat

  # transpose distmat so that the rows are the landscape cells
  if (!is.null(distmat)) distmat <- t(distmat)

  # run sliding window calculations

  # wrapping SpatRaster so it can be passed to future_map
  wlyr <- terra::wrap(lyr)

  rast_vals <-
    furrr::future_map(
      1:terra::ncell(lyr),
      ~ window_helper(
        i = .x,
        wlyr = wlyr,
        x = x,
        coord_cells = coord_cells,
        nmat = nmat,
        distmat = distmat,
        maxdist = maxdist,
        stat_function = stat_function,
        rarify = rarify,
        rarify_n = rarify_n,
        rarify_nit = rarify_nit,
        min_n = min_n,
        fun = fun,
        L = L,
        rarify_alleles = rarify_alleles,
        sig = sig
      ),
      .options = furrr::furrr_options(
        seed = TRUE,
        packages = c("wingen", "terra")
      ),
      .progress = TRUE
    )

  # format resulting raster values
  result <- vals_to_lyr(lyr, rast_vals, stat)

  return(result)
}

#' Helper function for window calculations
#'
#' Provide nmat for `window_gd()` and distmat for `circle_gd()`/`resist_gd()`/`dist_gd()`
#' @return genetic diversity and counts for a single cell
#'
#' @noRd
window_helper <- function(i, x, wlyr,
                          coord_cells = NULL, nmat = NULL,
                          distmat = NULL, maxdist = NULL,
                          stat_function,
                          rarify, rarify_n, rarify_nit, min_n,
                          fun, L, rarify_alleles, sig) {
  # unwrap layer
  lyr <- terra::unwrap(wlyr)

  # if rarify = TRUE and rarify_n isn't specified, rarify_n = min_n (i.e. rarify_n defaults to min_n)
  if (is.null(rarify_n)) rarify_n <- min_n

  # skip if raster value is NA
  # note: need to provide ns to give some output so the cell is counted
  # don't need to provide gd as this is repaired later
  if (is.na(lyr[i])) {
    return(c(sample_count = NA))
  }

  # get sample indices in neighborhood rectangle if nmat is provided (window_gd)
  if (!is.null(nmat)) sub <- get_adj(i, lyr, nmat, coord_cells)

  # get sample indices based on distance (circle_gd/resist_gd/dist_gd)
  if (!is.null(distmat)) sub <- get_dist_index(i, distmat, maxdist)

  # if there are too few samples in that window assign the cell value NA
  if (length(sub) < min_n) {
    gd <- NA
  } else if (rarify) {
    gd <- rarify_helper(x, sub, rarify_n, rarify_nit, stat_function, fun, L = L, rarify_alleles = rarify_alleles, sig = sig)
  } else {
    gd <- sample_gd(x, sub, stat_function, L = L, rarify_alleles = rarify_alleles, sig = sig)
  }

  # count the number of samples in the window
  ns <- length(sub)

  # if gd has no name give call it custom
  if (is.null(names(gd))) names(gd) <- "custom"

  # return vector
  if (all(is.na(gd))) {
    return(c(sample_count = ns))
  } else {
    return(c(gd, sample_count = ns))
  }
}

#' Rarefaction helper function
#'
#' @inheritParams window_general
#'
#' @noRd
#'
#' @return genetic diversity statistic for a rarified subsample
#'
#' @noRd
rarify_helper <- function(x, sub, rarify_n, rarify_nit, stat_function,
                          fun, L, rarify_alleles, sig) {
  # if number of samples is less than rarify_n, assign the value NA
  if (length(sub) < rarify_n) {
    gd <- NA
  }

  # if number of samples is greater than rarify_n, rarify
  if (length(sub) > rarify_n) {
    gd <- rarify_gd(x, sub, rarify_nit = rarify_nit, rarify_n = rarify_n, stat_function = stat_function, fun = fun, L = L, rarify_alleles = rarify_alleles, sig = sig)
  }

  # if the number of samples is equal to rarify_n, calculate stat
  if (length(sub) == rarify_n) {
    gd <- sample_gd(x, sub, stat_function, L = L, rarify_alleles = rarify_alleles, sig = sig)
  }

  return(gd)
}


#' Helper function to rarify subsample and calculate genetic diversity
#'
#' @inheritParams window_general
#'
#' @return rarified genetic diversity statistic
#'
#' @noRd
rarify_gd <- function(x, sub, rarify_nit, rarify_n, stat_function,
                      fun, L, rarify_alleles, sig) {
  # check to make sure sub is greater than rarify_n
  if (!(length(sub) > rarify_n)) {
    stop("rarify_n is less than the number of samples provided")
  }

  # define subsample to rarify
  if (rarify_nit == "all") {
    cmb <- utils::combn(sub, rarify_n)
  } else if (choose(length(sub), rarify_n) < rarify_nit) {
    # (note: this combo step is done so when the number of unique combos < rarify_nit, extra calcs aren't performed)
    # get all possible combos (transpose so rows are unique combos)
    cmb <- utils::combn(sub, rarify_n)
  } else {
    # randomly sample subsets of size rarify_nit (transpose so rows are unique combos)
    # note: replace is set to FALSE so the same individual cannot be drawn multiple times within the same sample
    # however, individuals can be drawn multiple times across different samples
    cmb <- replicate(rarify_nit, sample(sub, rarify_n, replace = FALSE), simplify = TRUE)
  }

  # get all possible combos (rows are unique combos)
  # for each of the possible combos get gendiv stat
  cmb_ls <- as.list(data.frame(cmb))
  gdrar <- purrr::map(cmb_ls, \(sub) sample_gd(x = x, sub = sub, stat_function = stat_function, L = L, rarify_alleles = rarify_alleles, sig = sig))

  # summarize rarefaction results
  gd <- purrr::list_transpose(gdrar) %>% purrr::map_dbl(fun, na.rm = TRUE)

  return(gd)
}


#' Helper function to calculate genetic diversity of a sample
#'
#' @return mean allelic richness of a subsample
#'
#' @noRd
sample_gd <- function(x, sub, stat_function, L, rarify_alleles, sig) {
  if (isTRUE(all.equal(stat_function, calc_mean_biar))) {
    return(stat_function(x[sub, ], rarify_alleles = rarify_alleles))
  }

  if (isTRUE(all.equal(stat_function, calc_pi))) {
    return(stat_function(x[sub, ], L = L))
  }

  if (isTRUE(all.equal(stat_function, calc_prop_hwe))) {
    return(stat_function(x[sub, ], sig = sig))
  }

  return(stat_function(x[sub, ]))
}


#' Helper function to get adjacent cells to a given cell index
#'
#' @param i cell index
#' @param r SpatRast
#' @param n neighborhood matrix
#' @param coord_cells cell numbers of coordinates
#'
#' @return indices of coordinates that are adjacent to the given cell
#'
#' @noRd
get_adj <- function(i, r, n, coord_cells) {
  # get adjacent cells to cell i
  adjc <- terra::adjacent(r, i, directions = n, include = TRUE)

  # get indices of adjacent cells
  adjci <- purrr::map_dbl(adjc, 1, ~ seq(.x[1], .x[2]))

  # remove NA values from indices
  ## note: if NA values are not removed, coord_cells with NA values will be included
  adjci_nona <- adjci[!is.na(adjci)]

  # get vector of indices of coords in that set of cells
  sub <- which(coord_cells %in% adjci_nona)

  return(sub)
}

#' Check coordinate and genetic data
#'
#' Check that the number of individuals in each data set align
#'
#' @param x moving window data
#' @param coords coordinates
#' @param distmat distance matrix
#'
#' @noRd
#'
check_data <- function(x, coords = NULL, distmat = NULL) {
  # if x is a vector, convert to a dataframe
  if (is.vector(x)) x <- data.frame(x)

  # check and format coords
  if (!is.null(coords)) {
    if (nrow(coords) == 1) stop("cannot run window_gd with only one individual")
  }

  # check number of samples
  nind <- NULL

  if (inherits(x, "genind")) nind <- nrow(x$tab)

  if (inherits(x, "vcfR")) nind <- (ncol(x@gt) - 1)

  if (inherits(x, "data.frame") | inherits(x, "matrix") | inherits(x, "sf")) nind <- nrow(x)

  # if no other type matches, try and calculate based on nrow() or length()
  if (is.null(nind)) {
    nind_row <- nrow(x)
    nind_length <- length(x)

    if ((!is.null(nind_row) & !is.null(nind_length)) | (is.null(nind_row) & is.null(nind_length))) stop("Unable to determine length or numeber of rows from the provided x")
    if (is.null(nind_row) & !is.null(nind_length)) nind <- nind_length
    if (!is.null(nind_row) & is.null(nind_length)) nind <- nind_row
  }

  # check coords
  if (!is.null(coords)) {
    if (nind != nrow(coords)) {
      stop("number of samples in coords data and number of samples in gen data are not equal")
    }
  }

  # check distmat
  if (!is.null(distmat)) {
    if (nind != nrow(distmat)) {
      stop("number of samples in distmat data and number of samples in gen data are not equal")
    }
  }

  # check for rows or columns with missing data in a vcf and give warning if there are invariant sites
  if (inherits(x, "vcfR")) {
    return(check_vcf_NA(vcf = x, coords = coords, distmat = distmat))
  } else {
    return(list(x = x, coords = coords, distmat = distmat))
  }
}

#' Check vcf for loci and individuals with all NAs and return corrected vcf and coords
#'
#' @param vcf vcfR
#' @param coords coordinates
#' @param distmat distance matrix
#'
#' @noRd
check_vcf_NA <- function(vcf, coords = NULL, distmat = NULL) {
  # check for mismatch before indexing
  if (!is.null(coords)) {
    if ((ncol(vcf@gt) - 1) != nrow(coords)) {
      stop("number of samples in coords data and number of samples in vcf are not equal")
    }
  }

  if (nrow(vcf@fix) == 1) {
    NA_ind <- get_allNA(vcf@gt[-1])
    NA_row <- FALSE
  } else {
    NA_ind <- get_allNA(vcf@gt[, -1], MARGIN = 2)
    NA_row <- get_allNA(vcf@gt[, -1], MARGIN = 1)
  }

  if (any(NA_row)) {
    warning("Markers with no scored alleles have been removed")
    vcf <- vcf[!NA_row, ]
  }

  if (any(NA_ind)) {
    warning("Individuals with no scored loci have been removed")
    vcf <- vcf[, c(TRUE, !NA_ind)]
  }

  # check for invariant sites
  if (any(!vcfR::is.polymorphic(vcf, na.omit = TRUE))) warning("invariant sites found in vcf")

  # make results
  if (!is.null(coords)) coords <- coords[!NA_ind, ]
  if (!is.null(distmat)) distmat <- distmat[!NA_ind, ]

  results <- list(vcf = vcf, coords = coords, distmat = distmat) %>% purrr::discard(is.null)
  return(results)
}

#' Helper function to get NA values
#'
#' @inheritParams base::array
#'
#' @noRd
get_allNA <- function(x, MARGIN = NULL) {
  if (is.null(dim(x))) allNA <- is.na(x)
  if (!is.null(dim(x))) {
    allNA <- apply(x, MARGIN, function(x) {
      all(is.na(x))
    })
  }
  return(allNA)
}

#' Convert vcf to correct format based on stat
#'
#' @param vcf vcfR
#' @param stat genetic diversity statistic
#'
#' @noRd
convert_vcf <- function(vcf, stat) {
  if (stat == "allelic_richness" | stat == "hwe") {
    return(vcfR::vcfR2genind(vcf))
  }

  if (stat == "Ho") {
    return(vcf_to_het(vcf))
  }

  if (stat == "pi" | stat == "biallelic_richness") {
    return(vcf_to_dosage(vcf))
  }

  if (stat == "basic_stats") {
    return(vcf_to_hf(vcf))
  }

  stop(paste0(stat, " is an invalid arugment for stat"))
}


#' Helper function to check and convert lyr and coords
#'
#' @param lyr RasterLayer or SpatRaster
#' @param coords sf points, data frame, or matrix representing coordinates
#' @param fact factor of aggregation
#'
#' @return SpatRaster
#'
#' @noRd
layer_coords_check <- function(lyr, coords, fact = 0) {
  # check coords and lyr
  crs_check_window(lyr, coords)

  # convert to terra
  if (!inherits(lyr, "SpatRaster")) lyr <- terra::rast(lyr)

  # check number of layers
  nlayers <- terra::nlyr(lyr)
  if (nlayers > 1) {
    warning(paste0(nlayers, " provided, but only one is need. Defaults to using the first layer."))
    lyr <- lyr[[1]]
  }
  # make aggregated raster
  if (fact != 0) lyr <- terra::aggregate(lyr, fact, fun = mean)

  return(lyr)
}

#' Convert values into new raster layers
#'
#' @param lyr SpatRaster
#' @param rast_vals dataframe of gd and ns
#'
#' @return SpatRaster
#'
#' @noRd
vals_to_lyr <- function(lyr, rast_vals, stat) {
  # bind rows to get vector per layer
  # note: an important repair will happen here to fill in NAs where there are no genetic diversity values
  ls <-
    rast_vals %>%
    dplyr::bind_rows() %>%
    dplyr::relocate("sample_count", .after = tidyselect::last_col()) %>%
    as.list()

  # assign vector values to rasters
  rast_list <- purrr::map(ls, ~ terra::setValues(lyr, .x))

  # convert from list to raster stack
  rast_stack <- terra::rast(rast_list)

  # assign back names (necessary if only one layer (sample_count) is returned)
  names(rast_stack) <- names(rast_list)

  # give warning if only sample_count is returned
  if (terra::nlyr(rast_stack) == 1) warning("All values NA, only a sample count layer will be returned")

  return(rast_stack)
}


#' Crop edge off raster
#'
#' @param x SpatRaster
#' @param wdim window dimensions
#'
#' @return SpatRaster
#'
#' @noRd
edge_crop <- function(x, wdim) {
  if (length(wdim) == 1) wdim <- c(wdim, wdim)

  # get extent
  x_ext <- terra::ext(x)

  # calculate x edge buffer
  x_edge_size <- terra::res(x)[1] * ((wdim[1] - 1) / 2)
  xmin <- x_ext$xmin + x_edge_size
  xmax <- x_ext$xmax - x_edge_size

  # calculate y edge buffer
  y_edge_size <- terra::res(x)[2] * ((wdim[2] - 1) / 2)
  ymin <- x_ext$ymin + y_edge_size
  ymax <- x_ext$ymax - y_edge_size

  # crop raster
  x_crop <- terra::crop(x, terra::ext(xmin, xmax, ymin, ymax))

  return(x_crop)
}

#' Remove duplicate sample count layers
#'
#' @param r SpatRaster
#'
#' @return SpatRaster
#'
#' @noRd
rm_duplicate_sample_count <- function(r) {
  # subset one sample_count layer
  sample_count <- r[[which(names(r) == "sample_count")]][[1]]

  # subset genetic diversity layers
  gd <- r[[which(names(r) != "sample_count")]]

  # recombine
  r <- c(gd, sample_count)
  return(r)
}
