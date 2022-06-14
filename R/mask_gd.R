
#' Mask diversity map based on sample number
#'
#' @param x RasterStack where the first layer is genetic diversity and the second layer is sample count
#' @param min_n minimum number of samples (everything less than this number is masked)
#' @param plot whether to plot results
#' @param bkg.col background color for plotting map
#' @param col.pal color palette to use for plotting
#'
#' @return RasterLayer
#' @export
#'
#' @examples
mask_gd <- function(x, min_n, plot = FALSE, bkg.col = "white", col.pal = viridis::magma(100)){

  sc_index <- grepl("^sample_count", names(x))
  ar <- x[[which(!sc_index)]]
  counts <- x[[which(sc_index)]]

  for(i in 1:raster::nlayers(ar)){
    counts[[i]][counts[[i]] < min_n] <- NA
    ar[[i]] <- raster::mask(ar[[i]], counts[[i]])
  }


  if(plot){
    for(i in 1:raster::nlayers(ar)){
      raster::plot(x[[which(!sc_index)]][[i]], col = bkg.col,  box = FALSE, axes = FALSE, legend = FALSE,  main = names(ar[[i]]))
      raster::plot(ar[[i]], col = col.pal, box = FALSE, axes = FALSE, add = TRUE)
    }
  }

  return(ar)
}
