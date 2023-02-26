coords_to_sf <- function(coords) {
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  sf_coords <- sf::st_as_sf(coords, coords = c("x", "y"))
  return(sf_coords)
}
