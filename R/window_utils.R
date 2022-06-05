
get_adj <- function(i, r, n, coord_cells){
  # get adjacent cells to cell i
  adjc <- raster::adjacent(r, i, directions = n, include = TRUE, sorted = TRUE)
  # get indices of adjacent cells
  adjci <- purrr::map_dbl(adjc, 1, function(x) {seq(x[1], x[2])})
  # get list of indices of coords in that set of cells
  sub <- which(coord_cells %in% adjci)

  return(sub)
}


