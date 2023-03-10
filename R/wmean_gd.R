
wmean_gd <- function(x, lyr = NULL, parallel = TRUE, ncores = 1){
  # make df from raster
  x_df <- terra::as.data.frame(x, xy = TRUE, na.rm = FALSE)

  # make distance matrix
  if (is.null(lyr)) {
    distmat <- as.matrix(Rfast::Dist(x_df))
  } else {
    lyr_df <- terra::as.data.frame(lyr, xy = TRUE)
    distmat <- get_resdist(lyr_df[, c("x", "y")], lyr, parallel = parallel, ncores = ncores)
  }

  # make inverse distance weighted matrix
  idw <- 1/distmat
  # TODO: think about this:
  idw[is.infinite(idw)] <- max(idw[!is.infinite(idw)], na.rm = TRUE)
  idw[is.na(idw)] <- 0

  # get raster values
  vals <- x_df[,3]
  # get weighted mean of raster values based on cell
  wmean <- purrr::map_dbl(1:nrow(x_df), ~get_wmean(.x, vals, idw))
  wmean_rast <- terra::setValues(x, wmean)

  return(wmean_rast)
}

get_wmean <- function(i, vals, idw){
  x_dist <- idw[i, ]
  wmean <- weighted.mean(x = vals, w = x_dist, na.rm = TRUE)
  return(wmean)
}
