#' Get directory for files
#'
#' @export
#'
get_exdir <- function(){
  return(here::here("paperex", "simex"))
}


#' Load middle earth data
#' @param subset if TRUE, subsets data
#' @param quiet if TRUE, no message is printed
#'
#' @export
#'
load_middle_earth <- function(subset = FALSE, quiet = FALSE){
  # get wdir
  wdir <- get_exdir()

  # load genetic data
  # check if file exists locally and if not download it
  file <- here::here(wdir, "data", "mod-sim_params_it-0_t-1000_spp-spp_0.vcf")
  if(!file.exists(file)){
    message("downloading vcf and storing locally...this will take some time, but only has to be done once")
    download.file("https://zenodo.org/record/7112468/files/mod-sim_params_it-0_t-1000_spp-spp_0.vcf?download=1", file)
  }

  vcf <- vcfR::read.vcfR(file, verbose = FALSE)
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

  # Create background layer for plotting
  bkg <- lyr
  bkg[bkg < 0.01] <- NA
  assign("bkg", bkg, envir = .GlobalEnv)

  # subset
  if (subset){
    # get subsampled dataset
    subdata <- subset_data(vcf, coords, nsamples = 200, nvariants = 10000)
    subvcf <- subdata[["vcf"]]
    subcoords <- subdata[["coords"]]
    assign("subcoords", subcoords, envir = .GlobalEnv)
    assign("subvcf", subvcf, envir = .GlobalEnv)


    if (!quiet) {
      # give message with information about objects
      return(message(cat(
        crayon::cyan(crayon::bold("\n--------------------- middle earth data ---------------------\n")),
        crayon::silver("\nObjects loaded:"),
        crayon::yellow(crayon::bold("\n*vcf*")),
        crayon::yellow(paste0("vcfR object (142129 loci x 1697 samples)")),
        crayon::green(crayon::bold("\n*coords*")), crayon::green("dataframe with x and y coordinates"),
        crayon::yellow(crayon::bold("\n*subvcf*")),
        crayon::yellow(paste0("vcfR object (10000 loci x 200 samples)")),
        crayon::green(crayon::bold("\n*subcoords*")), crayon::green("dataframe with x and y coordinates for 200 samples"),
        crayon::magenta(crayon::bold("\n*lyr*")), crayon::magenta("middle earth RasterLayer (100 x 100)"),
        crayon::blue(crayon::bold("\n*bkg*")), crayon::blue("background layer"),
        crayon::cyan(crayon::bold("\n\n-------------------------------------------------------------\n"))
      )))
    }

  } else {
    if (!quiet) {
      # give message with information about objects
      return(message(cat(
        crayon::cyan(crayon::bold("\n---------------- middle earth data ----------------\n")),
        crayon::silver("\nObjects loaded:"),
        crayon::yellow(crayon::bold("\n*lotr_vcf*")),
        crayon::yellow(paste0("vcfR object (100 loci x 100 samples)")),
        crayon::green(crayon::bold("\n*lotr_coords*")), crayon::green("dataframe with x and y coordinates"),
        crayon::magenta(crayon::bold("\n*lotr_lyr*")), crayon::magenta("middle earth RasterLayer (100 x 100)"),
        crayon::blue(crayon::bold("\n*lotr_range*")), crayon::blue("SpatialPolygonsDataFrame of spp range"),
        crayon::cyan(crayon::bold("\n\n---------------------------------------------------\n"))
      )))
    }
  }


  return()
}

#' Subset data using defined sample sets
#'
#' @param vcf full vcf
#' @param coords full coordinates
#' @param nsamples number of subsamples (either 100 or 200)
#' @param nvariants number of variants (either 10000 or 100000)
#'
#' @return
#' @export
#'
#' @examples
subset_data <- function(vcf, coords, nsamples, nvariants){

  # make file names
  variant.file <- paste0("variants_seed42_", as.integer(nvariants), ".csv")
  sample.file <- paste0("samples_seed42_", as.integer(nsamples), ".csv")

  # get variants
  variants <- read.csv(here("paperex", "simex", "data", "samples", variant.file))[,1]
  # get ids of inds to sample
  samples <- read.csv(here("paperex", "simex", "data", "samples", sample.file))[,1]

  # subset coodinates
  subcoords <- coords[samples,]
  # subset vcf
  subvcf <- vcf[variants, c(1, samples + 1)]
  # check match
  stopifnot(colnames(subvcf@gt)[-1] == as.character(subcoords$idx))
  # confirm that correct set is being used
  message(paste("nvariants", nrow(subvcf@gt), "/ nind", nrow(subcoords)))

  return(list(vcf = subvcf, coords = subcoords[,c("x","y")]))
}

#' Run default time tests and produce raster outputs
#'
#' @inheritParams window_gd
#' @param file.name file name to append to beginning of outputs
#'
#' @export
default_time_test <- function(stat, vcf, coords, lyr, wdim = 7, fact = 3, rarify, rarify_n = 2, rarify_nit = 5,
                              min_n = 2, fun = mean, rarify_alleles = TRUE, parallel = FALSE, ncores = 10, file.name){

  # get wdir
  wdir <- get_exdir()

  ptm <- Sys.time()
  gdmapr <- window_gd(vcf = vcf,
                      coords = coords,
                      lyr = lyr,
                      stat = stat,
                      wdim = wdim,
                      fact = fact,
                      rarify = rarify,
                      rarify_n = rarify_n,
                      rarify_nit = rarify_nit,
                      min_n = min_n,
                      fun = mean,
                      rarify_alleles = rarify_alleles,
                      parallel = parallel, ncores = ncores)

  df <- data.frame(time = as.numeric(Sys.time() - ptm, units = "secs"),
                   stat = stat,
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

#' Helper function for default_time_test
#'
#' @inheritParams default_time_test
#'
#' @export
run_default_time_test <- function(vcf, coords, lyr, rarify, rarify_alleles = TRUE, parallel = TRUE, ncores = 10, file.name,
                                  stats =  c("pi", "Ho", "biallelic_richness")){
  # get wdir
  wdir <- get_exdir()

  results <- purrr::map(stats,
                        default_time_test,
                        vcf = vcf,
                        coords = coords,
                        lyr = lyr,
                        rarify = rarify,
                        rarify_alleles = rarify_alleles,
                        parallel = parallel,
                        ncores = ncores,
                        file.name = file.name)

  write_time_test(results, here(wdir, "outputs", paste0(file.name,"_rarify", rarify, "_nsamp", nrow(coords), "_nsnps", nrow(vcf@gt), "_parallel", parallel, "_time_results.csv")))
  purrr::map(results, write_rast_test, here(wdir, "outputs", paste0(file.name,"_rarify", rarify, "_nsamp", nrow(coords), "_nsnps", nrow(vcf@gt))))
}


#' Write out default time test results
#'
#' @param res results
#' @inheritParams default_time_test
#'
#' @export
#'
write_time_test <- function(res, file.name){
  df <- purrr::map_dfr(res, function(x){x[[1]]})
  write.csv(df, file.name)
}

#' Write out default time test raster results
#'
#' @param res results
#' @inheritParams default_time_test
#'
#' @export
#'
write_rast_test <- function(res, file.name){
  if(class(res[[2]]) == "RasterStack"){
    resl <- res[[2]]
    write_rast_helper(resl, file.name)
  } else {
    resl <- purrr::map(res, function(x){x[[2]]})
    purrr::map(resl, write_rast_helper, file.name)
  }
}

#' Helper function for write_rast_test
#'
#' @inheritParams write_rast_test
#'
#' @export
#'
write_rast_helper <- function(resl, file.name){
  lyrname <- names(resl)[1]
  terra::writeRaster(terra::rast(resl), paste0(file.name,"_", lyrname, ".tif"), overwrite = TRUE)
}

#' Get raster outputs from default time test files
#'
#' @inheritParams default_time_test
#' @param file.type type of file (defaults to tif)
#' @param rootPath directory root
#'
#' @export
get_divout <- function(file.name, rarify = NULL, stat = NULL, nsamp = NULL, file.type = ".tif", rootPath = here(get_exdir(), "outputs")){
  # Searches for file in directory
  listFiles <- list.files(rootPath, recursive = FALSE)
  presentFile <- grepl(file.name, listFiles) & grepl(file.type, listFiles)
  if(!is.null(rarify)){ presentFile <- presentFile & grepl(paste0("rarify", rarify), listFiles)}
  if(!is.null(stat)){ presentFile <- presentFile & grepl(stat, listFiles)}
  if(file.name != "FULL" & !is.null(nsamp)){ presentFile <- presentFile & grepl(paste0("nsamp", nsamp), listFiles)}
  locFile <- listFiles[presentFile]

  if(!any(presentFile)){warning(paste("File does not exist for", file.name, rarify, stat, "- returning NULL")); return(NULL)}
  file <- here(rootPath, locFile)
  if(length(file) > 1) {print(file); stop("more than one file provided")}
  r <- raster::stack(file)
  r <- r[[seq(1, nlayers(r),2)]]

  return(r)
}

#' Get time test outputs from default time test files
#'
#' @inheritParams default_time_test
#' @param file.type type of file (defaults to csv)
#' @param rootPath directory root
#'
#' @export
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

  r$stat <- c("pi", "Ho", "allelic richness")
  r$parallel <- parallel
  r$rarify <- rarify
  r$nsamp <- nsamp
  r$dataset <- file.name

  return(r)
}

#' Helper function to test different parameter combinations
#'
#' @param params vector of wdim and fact values
#' @inheritParams window_gd
#'
#' @return
#' @export
#'
#' @examples
test_params_simex <- function(params, vcf, coords, lyr, stat = "pi"){
  wdim <- as.numeric(params["wdim"])
  fact <- as.numeric(params["fact"])
  res <- window_gd(vcf,
                   coords,
                   lyr,
                   stat = stat,
                   wdim = wdim,
                   fact = fact,
                   rarify = TRUE,
                   rarify_n = 2,
                   rarify_nit = 5,
                   parallel = TRUE,
                   ncores = 4)
  return(res)
}

#' Convert dataframe to list of vectors
#'
#' @param x dataframe
#'
#' @export
#'
df_to_ls <- function(x){
  x <- split(x, seq(nrow(x)))
  return(x)
}

#' Helper function to get rasters from default time tests and mask them
#'
#' @param params vector with dataset type, rarify value, and stat
#' @param nsamp number of samples
#' @param msk_lyr mask layer
#'
#' @return
#' @export
#'
#' @examples
test_datasets_simex <- function(params, nsamp, msk_lyr){
  file.name <- as.character(params[["datasets"]])
  rarify <- as.character(params[["rarify"]])
  stat <- as.character(params[["stat"]])

  r <- get_divout(file.name = file.name,
                  rarify = rarify,
                  stat = stat,
                  nsamp = nsamp)

  if(is.null(r)) return(NULL)
  r <- mask(r, msk_lyr)
  names(r) <- stat
  return(r)
}

#' Create plots from default time test raster results
#'
#' @param r raster
#' @param bkg background plot
#' @param legend whether to plot legend
#' @param ... Graphical parameters. Any argument that can be passed to image.plot and to base plot, such as axes=FALSE, main='title', ylab='latitude'
#'
#' @return
#' @export
#'
#' @examples
test_simex_plot <- function(r, bkg = NULL, legend = FALSE, ...){
  stat <- names(r)[1]

  if(stat == "pi") zlim <- c(0, 0.31)
  if(stat == "biallelic_richness") zlim <- c(1, 1.95)
  if(stat == "Ho") zlim <- c(0, 0.29)
  if(!exists("zlim")) stop(paste(names(r)[1], "is not a valid stat"))

  plot_gd(r, bkg = bkg, zlim = zlim, legend = legend, breaks = 100, box = TRUE, ...)
  return(NULL)
}

