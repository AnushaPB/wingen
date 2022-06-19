library(wingen)
library(raster)
library(vcfR)
library(viridis)
library(foreach)
library(doParallel)
library(adegenet)
library(hierfstat)
library(here)
source(here("sims/sim_functions.R"))

set.seed(42)

vcf <- read.vcfR("data/mod-sim_params_it-0_t-500_spp-spp_0.vcf")
geo <- read.csv("data/mod-sim_params_it-0_t-500_spp-spp_0.csv")
coords <- geo[,c("x","y")]
coords$y <- -coords$y

lyr <- read.csv("data/middle_earth.csv", header = FALSE)
lyr <- raster(as.matrix(lyr))
extent(lyr) <- extent(0,100,-100,0)


cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

stats <- c("pi", "het", "allelic.richness")
for(i in 1:length(stats)){
  ptm <- Sys.time()
  gdmapr <- window_gd(vcf, coords, lyr, stat = stats[i], fact = 3, wdim = 5, rarify = TRUE, rarify_n = 4, rarify_nit = 5, min_n = 4, fun = mean, parallel = TRUE, nloci = nrow(vcf@gt))

  df <- data.frame(time = (Sys.time() - ptm),
                   fact = 3,
                   wdim = 5,
                   rarify_n = 4,
                   rarify_nit = 5)

  results[[i]] <- list(df, gdmapr)

  # temp: see results as they get output
  write_rast_test(results, here("sims/outputs/FULL"))
}

save(results, file = here("sims/outputs/FULL_rarify.RData"))

resls <- unlist_test(results)
write.csv(resls$df, here("sims/outputs/time_results_FULL.csv"))
write_rast_test(results, here("sims/outputs/FULL"))


stopCluster(cl)
