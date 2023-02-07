#' Create a moving window map of genetic diversity based on resistance
#'
#' Generate a continuous raster map of genetic diversity using resistance distances calculated with a conductivity surface
#'
#' @param maxdist maximum cost distance used to define neighborhood; any samples further than this cost distance will not be included (this can be thought of as the neighborhood radius, but in terms of cost distance)
#' @param cond_lyr conductivity layer (higher values should mean greater conductivity). Can be either a SpatRaster or RasterLayer. If not provided, `lyr` will be used.
#' @param distmat distance matrix output from \link[wingen]{get_resdist} (optional; can be used to save time on distance calculations)
#' @details Coordinates and rasters should be in a Euclidean coordinate system (i.e., UTM coordinates) such that raster cell width and height are equal distances.
#' As such, longitude-latitude systems should be transformed before using dist_gd. Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#' @export
#'
#' @examples
#'
#' load_mini_ex()
#' wpi <- resist_gd(mini_vcf, mini_coords, mini_lyr)
#' plot_gd(wpi, main = "Window pi")
#' plot_count(wpi)
#'
resist_gd <- function(vcf, coords, lyr, maxdist, cond_lyr = NULL, distmat = NULL, stat = "pi", fact = 0,
                   rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                   fun = mean, L = "nvariants", rarify_alleles = TRUE,
                   parallel = FALSE, ncores = NULL){

  # convert coords if not in sf
  if (!inherits(coords, "sf")) coords <- coords_to_sf(coords)

  # make distmat
  if(!is.null(con_lyr)) con_lyr <- lyr
  if(!is.null(distmat)) distmat <- get_resdist(coords, cond.r = cond_layer, parallel = parallel, ncores = ncores)
  distmat[distmat > maxdist] <- NA

  # run dist_gd
  results <-
    dist_gd(vcf = vcf,
            coords = coords,
            lyr = lyr,
            maxdist = maxdist,
            distmat = distmat,
            stat = stat,
            fact = fact,
            rarify = rarify,
            rarify_n = rarify_n,
            rarify_nit = rarify_nit,
            min_n = min_n,
            fun = fun,
            L = L,
            rarify_alleles = rarify_alleles,
            parallel = parallel,
            ncores = ncores)

  return(results)
}



get_resdist <- function(coords, cond.r, ncores = 1, parallel = TRUE, progress = TRUE){
  # rename coords
  colnames(coords) <- c("x", "y")

  # Create transition surface
  trSurface <- gdistance::transition(cond.r, transitionFunction = mean, directions = 8)
  trSurface <- gdistance::geoCorrection(trSurface, type = "c", scl = FALSE)

  # get layer coordinates
  lyr_coords <- terra::as.data.frame(cond.r, xy = TRUE, na.rm = FALSE)[,1:2]

  # get all combinations of layer and sample coordinate indices
  params <- expand.grid(list(lyr = 1:nrow(lyr_coords), coords = 1:nrow(coords)))

  # make vector of distances
  future::plan(future::multisession, workers = ncores)

  suppressWarnings({
    distvec <- furrr::future_map2_dbl(params$lyr, params$coords, run_gdist, trSurface, lyr_coords, coords, .progress = progress)
  })

  # convert from vector to matrix
  distmat <-
    data.frame(params, dist = distvec) %>%
    tidyr::pivot_wider(names_from = coords, values_from = dist) %>%
    dplyr::select(-1) %>%
    as.matrix()

  return(distmat)
}


run_gdist <- function(x, y, trSurface, lyr_coords, coords){
  # Make spatial points
  sp <- sp::SpatialPoints(rbind(lyr_coords[x,], coords[y,]))

  # Calculate circuit distances
  distmat <- possible_gdist(trSurface, sp)

  # Get distance
  if(!all(is.na(distmat))) distmat <- distmat[1,2]

  return(distmat)
}

possible_gdist <- purrr::possibly(function(trSurface, sp) as.matrix(gdistance::commuteDistance(trSurface, sp)), NA)
