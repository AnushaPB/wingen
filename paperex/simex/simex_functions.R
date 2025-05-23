#' Get directory for files
#'
#'
#'
get_exdir <- function(){
  return(here::here("paperex", "simex"))
}

#' Load middle earth data
#' @param subset if TRUE, subsets data
#' @param quiet if TRUE, no message is printed
#'
#'
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
        crayon::yellow(crayon::bold("\n*vcf*")),
        crayon::yellow(paste0("vcfR object (142129 loci x 1697 samples)")),
        crayon::green(crayon::bold("\n*coords*")), crayon::green("dataframe with x and y coordinates"),
        crayon::magenta(crayon::bold("\n*lyr*")), crayon::magenta("middle earth RasterLayer (100 x 100)"),
        crayon::blue(crayon::bold("\n*bkg*")), crayon::blue("background layer"),
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
#'
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
#'
default_time_test <- function(stat, vcf, coords, lyr, wdim = 7, fact = 3, rarify, rarify_n = 2, rarify_nit = 5,
                              min_n = 2, fun = mean, rarify_alleles = TRUE, parallel = FALSE, ncores = 10, file.name){

  # get wdir
  wdir <- get_exdir()

  ptm <- Sys.time()
  gdmapr <- window_gd(gen = vcf,
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
#'
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
#'
write_time_test <- function(res, file.name){
  df <- purrr::map_dfr(res, function(x){x[[1]]})
  write.csv(df, file.name)
}

#' Write out default time test raster results
#'
#' @param res results
#' @inheritParams default_time_test
write_rast_test <- function(res, file.name){
  if(inherits(res[[2]], "RasterStack")) res[[2]] <- terra::rast(res[[2]])
  if(inherits(res[[2]], "SpatRaster")){
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
write_rast_helper <- function(resl, file.name){
  if(!inherits(resl, "SpatRaster")) resl <- terra::rast(resl)
  lyrname <- names(resl)[1]
  terra::writeRaster(resl, paste0(file.name,"_", lyrname, ".tif"), overwrite = TRUE)
}

#' Get raster outputs from default time test files
#'
#' @inheritParams default_time_test
#' @param file.type type of file (defaults to tif)
#' @param rootPath directory root
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

  # remove AR calculation and rename biallelic richness for plotting
  r <-
    r %>%
    filter(stat != "allelic_richness") %>%
    mutate(stat = case_when(stat == "biallelic_richness" ~ "allelic richness",
                            .default = stat))

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
                   rarify_nit = 5)
  return(res)
}

#' Convert dataframe to list of vectors
#'
#' @param x dataframe
df_to_ls <- function(x){
  x <- split(x, seq(nrow(x)))
  return(x)
}

#' Helper function to get rasters from default time tests and mask them
#'
#' @param params vector with dataset type, rarify value, and stat
#' @param nsamp number of samples
#' @param msk_lyr mask layer
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
test_simex_plot <- function(r, bkg = NULL, legend = FALSE, zlim = NULL, polyx = 99, polyy = -99,...){
  stat <- names(r)[1]

  if(is.null(zlim)){
    if(stat == "pi") zlim <- c(0, 0.31)
    if(stat == "biallelic_richness") zlim <- c(1, 1.95)
    if(stat == "Ho") zlim <- c(0, 0.29)
  }

  if(inherits(r, "SpatRaster")) r <- raster::raster(r)

  par(pty = "s")
  raster_plot_gd(r, bkg = bkg, zlim = zlim, legend = legend, breaks = 100, box = FALSE, ...)
  polygon(x = c(0, 0, polyx, polyx), y = c(polyy, 0, 0, polyy), border = "black")

  return(NULL)
}

#' Create difference plots from default time test raster results
#'
#' @param r raster
#' @param bkg background plot
#' @param legend whether to plot legend
#' @param ... Graphical parameters. Any argument that can be passed to image.plot and to base plot, such as axes=FALSE, main='title', ylab='latitude'
test_simex_dif_plot <- function(r, bkg = NULL, legend = FALSE, polyx = 99, polyy = -99, ...){

  coldiv <- colorRampPalette(c("#0056A4", "#008DBC", "#00CCD3", 'gray96', "#FAAB36", "#F78104", "#FD5901"))

  raster_plot_gd(r, bkg = bkg, col = coldiv(100), zlim = c(-0.48, 0.48), legend = legend, breaks = 100, box = FALSE, ...)
  polygon(x = c(0, 0, polyx, polyx), y = c(polyy, 0, 0, polyy), border = "black")

  return(NULL)
}

#' Calculate difference between full and subsampled rasters
#' @param dataset dataset type (i.e., WGS or rr)
#' @param params parameters for filtering
#' @param nsamp number of samples (i.e., 100 or 200)
#' @param msk_lyr layer for masking
#'
simex_get_dif <- function(dataset, params, nsamp){

  # Get example layers for masking (doesn't matter which parameters other than nsamp)
  msk_lyr <- get_divout(file.name = "rr", rarify = TRUE, stat = "pi", nsamp = nsamp)

  subsample <-
    keep(params, ~.x$datasets == dataset) %>%
    purrr::map(test_datasets_simex, nsamp = nsamp, msk_lyr = msk_lyr)

  full <-
    keep(params, ~.x$datasets == "FULL") %>%
    purrr::map(test_datasets_simex, nsamp = nsamp, msk_lyr = msk_lyr)

  dif <-
    purrr::map2(subsample, full, ~(.x - .y)) %>%
    purrr::map2(subsample, function(x, y) {names(x) <- names(y); return(x)})

  return(dif)
}

#' raster_plot_gd rewritten for raster
#'
#'
raster_plot_gd <- function(x, bkg = NULL, index = NULL, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {
  if (inherits(x, "SpatRaster")) x <- raster::raster(x)
  if (inherits(x, "SpatRaster")) bkg <- raster::raster(bkg)

  if (is.null(index) & raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, plotting first layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 1

  # suppress irrelevant plot warnings
  suppressWarnings({
    if (!is.null(bkg)) {
      plt <- purrr::map(index, raster_plot_gd_bkg, x = x, bkg = bkg, col = col, breaks = breaks, main = main, box = box, ...)
    } else {
      plt <- raster::plot(x[[index]],
                          col = col,
                          axes = FALSE,
                          box = box,
                          ...
      )
      graphics::title(main = list(main, font = 1), adj = 0)
    }
  })

  return(invisible(plt))
}

#' Helper function for raster_plot_gd
#'
#' @inheritParams raster_plot_gd
#'
#' @noRd
raster_plot_gd_bkg <- function(index, x, bkg, col = viridis::magma(breaks), breaks = 20, main = NULL, box = FALSE, ...) {
  # suppress irrelevant plot warnings
  suppressWarnings({
    # calculate extent
    extx <- raster::extent(x)
    extb <- raster::extent(bkg)
    xmin <- min(extx[1], extb[1])
    xmax <- max(extx[2], extb[2])
    ymin <- min(extx[3], extb[3])
    ymax <- max(extx[4], extb[4])

    raster::plot(x[[index]],
                 col = "white",
                 xlim = c(xmin, xmax),
                 ylim = c(ymin, ymax),
                 axes = FALSE,
                 box = box,
                 legend = FALSE
    )

    raster::plot(bkg,
                 col = "lightgray",
                 border = "white",
                 xlim = c(xmin, xmax),
                 ylim = c(ymin, ymax),
                 axes = FALSE,
                 box = FALSE,
                 legend = FALSE,
                 add = TRUE
    )

    raster::plot(x[[index]],
                 col = col,
                 add = TRUE,
                 axes = FALSE,
                 box = FALSE,
                 ...
    )
  })

  graphics::title(main = list(main, font = 1), adj = 0)

  return()
}

#' Plot moving window map of sample counts
#'
#' Plot sample counts layer produced by \link[wingen]{window_gd} or \link[wingen]{krig_gd}
#'
#' @param x single SpatRaster of counts or SpatRaster where indexed layer is sample counts
#' @param index if a raster stack is provided, index of the sample count layer to plot (assumes this is a stacked output from window_gd and defaults to plotting second layer)
#' @param col color palette to use for plotting (defaults to viridis::magma palette)
#' @param breaks number of breaks to use in color scale (defaults to 10)
#' @param box whether to include a box around the raster plot (defaults to FALSE)
#' @inheritParams raster_plot_gd
#' @inheritParams terra::plot
#'
#' @return plot of sample counts
#' @export
#'
#' @examples
#' data("mini_lyr")
#' plot_count(mini_lyr)
raster_plot_count <- function(x, index = NULL, breaks = 100, col = viridis::mako(breaks), main = NULL, box = FALSE, ...) {
  if (inherits(x, "SpatRaster")) x <- raster::raster(x)

  if (is.null(index) & raster::nlayers(x) > 2) warning("More than two raster layers in stack provided, plotting second layer (to change this behavior use the index argument)")
  if (is.null(index)) index <- 2

  # suppress annoying and irrelevant plot warnings
  suppressWarnings({
    if (raster::nlayers(x) > 1) {
      plt <- raster::plot(x[[index]],
                         col = col,
                         axes = FALSE,
                         box = box,
                         ...
      )
      graphics::title(main = list(main, font = 1), adj = 0)
    }

    if (raster::nlayers(x) == 1) {
      plt <- raster::plot(x,
                          col = col,
                          axes = FALSE,
                          box = box,
                          ...
      )
      graphics::title(main = list(main, font = 1), adj = 0)
    }
  })

  return(invisible(plt))
}

