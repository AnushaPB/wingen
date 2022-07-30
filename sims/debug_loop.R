rast_vals <- purrr::map_dfr(1:raster::ncell(lyr), window_helper, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, nloci)

for(i in 1:raster::ncell(lyr)){
  window_helper(i, lyr, gen, coord_cells, nmat, stat, rarify, rarify_n, rarify_nit, min_n, fun, nloci)
}
