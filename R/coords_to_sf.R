coords_to_sf <- function(coords) {
  if (inherits(coords, "sf")) {
    return(coords)
  }
  if (inherits(coords, "SpatVector")) {
    return(sf::st_as_sf(coords))
  }
  if (is.matrix(coords)) coords <- data.frame(coords)
  if (is.data.frame(coords)) colnames(coords) <- c("x", "y")
  return(sf::st_as_sf(coords, coords = c("x", "y")))
}
