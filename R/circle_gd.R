
#' Create a moving window map of genetic diversity using a circle window
#'
#' Generate a continuous raster map of genetic diversity using circle moving windows
#'
#' @param maxdist maximum geographic distance used to define neighborhood; any samples further than this distance will not be included (this can be thought of as the neighborhood radius)
#' @param distmat distance matrix output from \link[wingen]{get_geodist} (optional; can be used to save time on distance calculations)
#' @details Coordinates and rasters should be in a Euclidean coordinate system (i.e., UTM coordinates) such that raster cell width and height are equal distances.
#' As such, longitude-latitude systems should be transformed before using dist_gd. Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).
#'
#' @return SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
#' @export
#'
#' @examples
#'
#' load_mini_ex()
#' wpi <- circle_gd(mini_vcf, mini_coords, mini_lyr, rarify = TRUE)
#' plot_gd(wpi, main = "Window pi")
#' plot_count(wpi)
#'
circle_gd <- function(vcf, coords, lyr, maxdist, distmat = NULL, stat = "pi", fact = 0,
                      rarify = FALSE, rarify_n = 2, rarify_nit = 5, min_n = 2,
                      fun = mean, L = "nvariants", rarify_alleles = TRUE,
                      parallel = FALSE, ncores = NULL){

  # convert coords if not in sf
  if (!inherits(coords, "sf")) coords <- coords_to_sf(coords)

  # make distmat
  if(!is.null(distmat)) distmat <- get_geodist(coords, lyr, parallel = parallel, ncores = ncores)
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


get_geodist <- function(coords, lyr, parallel = FALSE, ncores = NULL){
  lyr_df <- terra::as.data.frame(lyr, xy = TRUE, na.rm = FALSE)
  lyr_sf <- sf::st_as_sf(lyr_df, coords = c("x", "y"), crs = terra::crs(lyr))

  if (parallel){
    if (is.null(ncores)) ncores <- future::availableCores() - 1
    distls <- furrr::future_map(1:nrow(lyr_sf), ~ sf::st_distance(.y[.x,], coords), lyr_sf,
                                .options = furrr::furrr_options(seed = TRUE, packages = c("sf")))
  } else {
    distls <- purrr::map(1:nrow(lyr_sf), ~ sf::st_distance(.y[.x,], coords), lyr_sf)
  }


  distmat <- matrix(unlist(distls), ncol = length(distls[[1]]), byrow = TRUE)

  return(distmat)
}
