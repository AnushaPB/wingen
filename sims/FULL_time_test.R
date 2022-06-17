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

#results <- purrr::map(c("het", "pi", "allelic.richness"), time_test, "stat", vcf, coords, lyr, fact = 3, wdim = 5,  rarify = FALSE, min_n = 4, parallel = TRUE, nloci = nrow(vcf@gt))

#save(results, file = "results_FULL2.RData")

results <- purrr::map(c("het", "pi", "allelic.richness"), time_test, "stat", vcf, coords, lyr, fact = 3, wdim = 5, rarify_n = 4, rarify_nit = 5, parallel = TRUE, nloci = nrow(vcf@gt))

save(results, file = "results_FULL_rarify.RData")


stopCluster(cl)
