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

vcf <- read.vcfR(here("sims/data/mod-sim_params_it-0_t-500_spp-spp_0.vcf"))
geo <- read.csv(here("sims/data/mod-sim_params_it-0_t-500_spp-spp_0.csv"))
coords <- geo[,c("idx","x","y")]
IDS <- read.csv(here("sims/data/samples_seed42.csv"))

# get indexes of individuals
si <- IDS$inds
# sample loci
l <- sample(nrow(vcf@gt), 100000)

subvcf <- vcf[l,c(1, si+1)]

subcoords <- coords[si,]
subcoords <- subcoords[,c("x","y")]
subcoords$y <- -subcoords$y
coords$y <- -coords$y

# check match
all(colnames(subvcf@gt)[-1] == geo$idx[si])

lyr <- read.csv("data/middle_earth.csv", header = FALSE)
lyr <- raster(as.matrix(lyr))
extent(lyr) <- extent(0,100,-100,0)

cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

results <- purrr::map(c("het", "pi", "allelic.richness"), time_test, "stat", subvcf, subcoords, lyr, fact = 3, wdim = 5, rarify_n = 4, rarify_nit = 5, parallel = TRUE, nloci = nrow(subvcf@gt))

save(results, file = "results_WGS_10p_200ind.RData")

resls <- unlist_test(results)
write.csv(resls$df, "time_results_WGS_10p_200ind.csv")
write_rast_test(results, "WGS_200ind")


stopCluster(cl)
