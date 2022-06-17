library(wingen)
library(raster)
library(vcfR)
library(viridis)
library(foreach)
library(doParallel)
library(adegenet)
library(hierfstat)

source("sim_functions.R")

vcf <- read.vcfR("data/mod-sim_params_it-0_t-500_spp-spp_0.vcf")
geo <- read.csv("data/mod-sim_params_it-0_t-500_spp-spp_0.csv")
coords <- geo[,c("idx","x","y")]

set.seed(42)

coords$y <- -coords$y


lyr <- read.csv("data/middle_earth.csv", header = FALSE)
lyr <- raster(as.matrix(lyr))
extent(lyr) <- extent(0,100,-100,0)

cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

full_norarify <-time_test("pi", "stat", vcf, coords[,c("x","y")], lyr, fact = 3, wdim = 5, min_n = 2, rarify = FALSE, nloci = nrow(vcf@gt))

save(full_norarify, file = "results_FULL_noRarify.RData")

full_rarify <- time_test("pi", "stat", vcf, coords[,c("x","y")], lyr, fact = 3, wdim = 5, min_n = 4, rarify_n = 4, rarify_nit = 5, nloci = nrow(vcf@gt))

save(full_rarify, file = "results_FULL_Rarify.RData")

stopCluster(cl)
