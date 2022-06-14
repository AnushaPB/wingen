grid_samp <- function(pts, npts, ldim, full = FALSE){
  inc <- ldim/sqrt(npts)
  xgrid <- ygrid <- seq(0, ldim, inc)
  subs <- c()
  #first round of sampling: entire grid
  for(i in 1:(length(xgrid)-1)){
    for(j in 1:(length(ygrid)-1)){
      gridsq = subset(pts, y > ygrid[j] & y < ygrid[j+1] & x > xgrid[i] & x < xgrid[i+1])
      if(dim(gridsq)[1]>0){ subs = rbind(subs, gridsq[sample(1:dim(gridsq)[1],1 ), ]) }
    }
  }

  #second round of sampling: cycle through gridcells again until number of desired samples is reached
  if(full){
    while(nrow(subs) != npts){
      #reset grid indices
      i=1
      j=1
      while(nrow(subs) != npts & i < (length(xgrid)-1) & j < (length(ygrid)-1)){
        i = i+1
        j = j+1
        gridsq = subset(pts, y > ygrid[j] & y < ygrid[j+1] & x > xgrid[i] & x < xgrid[i+1])
        if(dim(gridsq)[1]>0){subs = rbind(subs, gridsq[sample(1:dim(gridsq)[1],1 ), ])}
      }
    }
  }
  return(subs$idx)
}


sim <- function(vcf, coords, lyr, stat, fact = 0, wdim = 10, min_n = 2, rarify = FALSE, rarify_n = 4, rarify_nit = 10, parallel = FALSE){

  if(parallel){
    cores <- 6
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  }

  # Start the clock!
  ptm <- Sys.time()

  res <- window_gd(vcf, coords, lyr, stat, fact, wdim, rarify = rarify, rarify_n, rarify_nit, min_n, fun = mean, parallel)

  # Stop the clock
  pt <- (Sys.time() - ptm)

  plot(res,  col = magma(100))

  if(parallel){
    stopCluster(cl)
  }

  return(list(pt = pt, res = res))
}

time_test <- function(val, var, vcf, coords, lyr, stat = "pi", fact = 2, wdim = 5, rarify = TRUE, rarify_n = 10, rarify_nit = 10, min_n = 2){

  # reassign argument
  assign(var, val)

  ptm <- Sys.time()
  gdmapr <- window_gd(vcf, coords, lyr, stat, fact, wdim, rarify, rarify_n, rarify_nit,  min_n)
  #plot(gdmapr, col = magma(100))

  msk_count <- mask(gdmapr[[2]], gdmapr[[1]])
  if(rarify){
    no_rarify <- cellStats(msk_count == rarify_n, "sum", na.rm = TRUE)*rarify_n
    yes_rarify <- cellStats(msk_count > rarify_n, "sum", na.rm = TRUE)*rarify_n*rarify_nit
    total_count = yes_rarify + no_rarify
  } else {
    total_count <- cellStats(msk_count, "sum", na.rm = TRUE)
  }


  df <- data.frame(time = (Sys.time() - ptm),
                   total_count = total_count,
                   ncell = ncell(aggregate(lyr, fact)),
                   wsize = wdim*wdim,
                   wprop = wdim*wdim/ncell(aggregate(lyr, fact)),
                   fact = fact)

  return(list(df, gdmapr))
}

unlist_test <- function(res){
  df <- purrr::map_dfr(res, function(x){x[[1]]})
  r <- purrr::map(res, function(x){x[[2]]})
  return(list(df = df, raster = r))
}


plot_time_test <- function(res, stat = "all"){
  res <- unlist_test(res)$df
  if(stat == "all"){
    par(pty = "s", mfrow = c(1,5))
    plot(res$total_count, res$time, type = "b")
    plot(res$ncell, res$time, type = "b")
    plot(res$wsize, res$time, type = "b")
    plot(res$wprop, res$time, type = "b")
    plot(res$fact, res$time, type = "b")
  } else {
    par(pty = "s", mfrow = c(1,1))
    plot(res[,stat], res$time, type = "b", xlab = stat, ylab = "System Time (seconds)")
    return(glm(res[,stat] ~ res$time))
  }

}

plot_rast_test <- function(res, zlim1 = NULL, zlim2 = NULL){
  res <- unlist_test(res)$raster
  par(pty = "s", mfrow = c(2, length(res)), mar = rep(2,4), oma = rep(2,4))
  lapply(res, function(x){plot(x[[1]], axes = FALSE, box = FALSE, col = viridis::magma(100), zlim = zlim1)})
  lapply(res, function(x){plot(x[[2]], axes = FALSE, box = FALSE, col = viridis::mako(100), zlim = zlim2)})
}
