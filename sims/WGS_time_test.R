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

vcf <- read.vcfR("data/mod-sim_params_it-0_t-500_spp-spp_0.vcf")
geo <- read.csv("data/mod-sim_params_it-0_t-500_spp-spp_0.csv")
coords <- geo[,c("idx","x","y")]

set.seed(42)
l <- sample(nrow(vcf@gt), 100000)
#s <- grid_samp(coords, npts = 350, ldim = 100)
#s <- s[1:200]
#si <- which(coords$idx %in% s)
si <- sample(nrow(coords), 100)

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

stopCluster(cl)

save(results, file = "results_WGS_10p.RData")
