get_exdir <- function(){
  return(here::here("paperex", "simex"))
}

load_middle_earth <- function(){
  # get wdir
  wdir <- get_exdir()

  # load genetic data
  vcf <- vcfR::read.vcfR(here::here(wdir, "data", "mod-sim_params_it-0_t-1000_spp-spp_0.vcf"))
  assign("vcf", vcf, envir = .GlobalEnv)

  # load coords
  geo <- read.csv(here::here(wdir, "data", "mod-sim_params_it-0_t-1000_spp-spp_0.csv"))
  coords <- geo[,c("idx", "x", "y")]
  coords$y <- -coords$y
  assign("coords", coords, envir = .GlobalEnv)

  # load rasters
  lyr <- read.csv(here::here(wdir, "data", "middle_earth.csv"), header = FALSE)
  lyr <- raster::raster(as.matrix(lyr))
  raster::extent(lyr) <- raster::extent(0,100,-100,0)
  assign("lyr", lyr, envir = .GlobalEnv)

  return()
}

default_time_test <- function(stat, vcf, coords, lyr, wdim = 3, fact = 3, rarify, rarify_n = 2, rarify_nit = 5,
                              min_n = 2, fun = mean, parallel, ncores, file.name){

  # get wdir
  wdir <- get_exdir()

  ptm <- Sys.time()
  gdmapr <- window_gd(vcf, coords, lyr, stat, wdim = wdim, fact = fact, rarify = rarify, rarify_n = rarify_n, rarify_nit = rarify_nit,
                      min_n = min_n, fun = mean, parallel = parallel, ncores = ncores)

  df <- data.frame(time = as.numeric(Sys.time() - ptm, units = "secs"),
                   fact = fact,
                   wdim = wdim)

  if(rarify){
    df$rarify_n <- rarify_n
    df$rarify_nit <- rarify_nit
  } else {
    df$min_n <- min_n
  }

  # make ls of results
  results <- list(df, gdmapr)

  write_rast_test(results, here(wdir, "outputs", paste0(file.name,"_rarify", rarify, "_nsamp", nrow(coords), "_nsnps", nrow(vcf@gt))))

  message("calculation of ", stat, " complete...")

  return(results)
}

run_default_time_test <- function(vcf, coords, lyr, rarify, parallel, ncores, file.name,
                                  stats =  c("pi", "het", "biallelic.richness")){
  # get wdir
  wdir <- get_exdir()

  results <- purrr::map(stats, default_time_test, vcf = vcf, coords = coords, lyr = lyr, rarify = rarify,
                        parallel = parallel, ncores = ncores, file.name = file.name)
  write_time_test(results, here(wdir, "outputs", paste0(file.name,"_rarify", rarify, "_nsamp", nrow(coords), "_nsnps", nrow(vcf@gt), "_parallel", parallel, "_time_results.csv")))
  purrr::map(results, write_rast_test, here(wdir, "outputs", paste0(file.name,"_rarify", rarify, "_nsamp", nrow(coords), "_nsnps", nrow(vcf@gt))))
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
  purrr::map(res, function(x){plot(x[[1]], axes = FALSE, box = FALSE, col = viridis::magma(100), zlim = zlim1)})
  purrr::map(res, function(x){plot(x[[2]], axes = FALSE, box = FALSE, col = viridis::mako(100), zlim = zlim2)})
}


write_time_test <- function(res, file.name){
  df <- purrr::map_dfr(res, function(x){x[[1]]})
  write.csv(df, file.name)
}

write_rast_test <- function(res, file.name){
  if(class(res[[2]]) == "RasterStack"){
    resl <- res[[2]]
    write_rast_helper(resl, file.name)
  } else {
    resl <- purrr::map(res, function(x){x[[2]]})
    purrr::map(resl, write_rast_helper, file.name)
  }
}

write_rast_helper <- function(resl, file.name){
  lyrname <- names(resl)[1]
  terra::writeRaster(terra::rast(resl), paste0(file.name,"_", lyrname, ".tif"), overwrite = TRUE)
}

get_divout <- function(file.name, rarify = NULL, measure = NULL, nsamp = NULL, file.type = ".tif", rootPath = here(get_exdir(), "outputs")){
  # Searches for file in directory
  listFiles <- list.files(rootPath, recursive = FALSE)
  presentFile <- grepl(file.name, listFiles) & grepl(file.type, listFiles)
  if(!is.null(rarify)){ presentFile <- presentFile & grepl(paste0("rarify", rarify), listFiles)}
  if(!is.null(measure)){ presentFile <- presentFile & grepl(measure, listFiles)}
  if(file.name != "FULL" & !is.null(nsamp)){ presentFile <- presentFile & grepl(paste0("nsamp", nsamp), listFiles)}
  locFile <- listFiles[presentFile]

  if(!any(presentFile)){warning(paste("File does not exist for", file.name, rarify, measure, "- returning NULL")); return(NULL)}
  file <- here(rootPath, locFile)
  if(length(file) > 1) {print(file); stop("more than one file provided")}
  r <- raster::stack(file)
  r <- r[[seq(1, nlayers(r),2)]]

  return(r)
}

get_timeout <- function(file.name, rarify = NULL, parallel = NULL, nsamp = NULL, file.type = ".csv", rootPath = here(get_exdir(), "outputs")){
  # Searches for file in directory
  listFiles <- list.files(rootPath, recursive = FALSE)
  presentFile <- grepl(paste0(file.name, "_rarify"), listFiles) & grepl(file.type, listFiles)
  if(!is.null(rarify)){ presentFile <- presentFile & grepl(paste0("rarify", rarify), listFiles)}
  if(!is.null(parallel)){ presentFile <- presentFile & grepl(paste0("parallel", parallel), listFiles)}
  if(file.name != "FULL" & !is.null(nsamp)){ presentFile <- presentFile & grepl(paste0("nsamp", nsamp), listFiles)}
  locFile <- listFiles[presentFile]

  if(!any(presentFile)){warning(paste("File does not exist for", file.name, rarify, parallel, "- skipping")); return()}
  file <- here(rootPath, locFile)
  r <- purrr::map_dfr(file, read.csv)

  r$stat <- c("pi", "heterozygosity", "allelic richness")
  r$parallel <- parallel
  r$rarify <- rarify
  r$nsamp <- nsamp
  r$dataset <- file.name

  return(r)
}

