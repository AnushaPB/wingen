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
  vcf <- vcfR::read.vcfR(here::here(wdir, "data", "mod-sim_params_it-0_t-1000_spp-spp_0.vcf"), verbose = FALSE)
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
    # Create file of individual coordinates if it doesn't exist already
    set.seed(42)
    if(!file.exists(here(wdir, "data", "samples_seed42.csv"))){
      message("creating new file")
      si <- sample(nrow(coords), 200)
      write.csv( data.frame(inds = si), here(wdir, "data", "samples_seed42.csv"), row.names = FALSE)
    } else {
      message("loading existing file")
      df <- read.csv(here(wdir, "data", "samples_seed42.csv"))
      si <- df$inds
    }

    # Subset coords
    subcoords <- coords[si,]
    subcoords <- subcoords[,c("x","y")]
    assign("subcoords", subcoords, envir = .GlobalEnv)

    # Sample 10k loci
    l <- sample(nrow(vcf@gt), 10000)
    # Subset VCF (note: first column of VCF are sample IDs)
    subvcf <- vcf[l, c(1, si + 1)]
    # Check match between VCF and coords
    stopifnot(all(colnames(subvcf@gt)[-1] == as.character(geo$idx[si])))
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

#' Run default time tests and produce raster outputs
#'
#' @inheritParams window_gd
#' @param file.name file name to append to beginning of outputs
#'
#' @export
default_time_test <- function(stat, vcf, coords, lyr, wdim = 3, fact = 3, rarify, rarify_n = 2, rarify_nit = 5,
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
                                  stats =  c("pi", "het", "biallelic.richness")){
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

  r$stat <- c("pi", "heterozygosity", "allelic richness")
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
#'
#' @return
#' @export
#'
#' @examples
test_simex_plot <- function(r, bkg, legend = FALSE){
  stat <- names(r)[1]

  if(stat == "pi"){zlim <- c(0, 0.30)}
  if(stat == "biallelic_richness"){zlim <- c(1, 1.86)}
  if(stat == "heterozygosity"){zlim <- c(0, 0.30)}

  plot_gd(r, bkg = bkg, zlim = zlim, legend = legend, breaks = 100, box = TRUE)
  return(NULL)
}

