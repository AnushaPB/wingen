
get_resdist <- function(coords, cond.r, ncores = 1, progress = TRUE){
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
