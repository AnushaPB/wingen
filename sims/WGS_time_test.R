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

lyr <- read.csv(here("sims/data/middle_earth.csv"), header = FALSE)
lyr <- raster(as.matrix(lyr))
extent(lyr) <- extent(0,100,-100,0)

cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

results <- list()

stats <- c("pi", "het", "allelic.richness")
for(i in 1:length(stats)){
  ptm <- Sys.time()
  gdmapr <- window_gd(subvcf, subcoords, lyr, stat = stats[[i]], fact = 3, wdim = 5, rarify = TRUE, rarify_n = 4, rarify_nit = 5, min_n = 4, fun = mean, parallel = TRUE, nloci = nrow(subvcf@gt))

  df <- data.frame(time = (Sys.time() - ptm),
                   fact = 3,
                   wdim = 5,
                   rarify_n = 4,
                   rarify_nit = 5)

  results[[i]] <- list(df, gdmapr)

  # temp: see results as they get output
  write_rast_test(results, here("sims/outputs/WGS_200ind"))
}

save(results, file = here("sims/outputs/WGS_200ind.RData"))

resls <- unlist_test(results)
write.csv(resls$df, here("sims/outputs/time_results_WGS_10p_200ind.csv"))
write_rast_test(results, here("sims/outputs/WGS_200ind"))


stopCluster(cl)
