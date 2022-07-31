library(wingen)
library(foreach)
library(doParallel)
library(here)
source(here("sims/sims_functions.R"))


## RR TEST -----------------------------------------------------------------------------------------------
# load vcf, coords, and lyr
load_middle_earth()

# get ids of inds to sample
#IDS <- read.csv(here("sims/data/samples_seed42.csv"))
# get indexes of individuals
#si <- IDS$inds
set.seed(42)
si <- sample(nrow(coords), 200, replace = FALSE)
# sample loci
l <- sample(nrow(vcf@gt), 10000)
# subset coodinates
coords <- coords[si,]
# subset vcf
vcf <- vcf[l, c(1, si + 1)]

# confirm that correct set is being used
message(paste("nloci", nrow(vcf@gt), "/ nind", nrow(coords)))

# check match
stopifnot(colnames(vcf@gt)[-1] == as.character(coords$idx))

cores <- 20
cl <- makeCluster(cores)
registerDoParallel(cl)

# run test
run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = TRUE, file.name = "rr", stats = "allelic.richness")

stopCluster(cl)

## WGS TEST -----------------------------------------------------------------------------------------------
# load vcf, coords, and lyr
load_middle_earth()

# get ids of inds to sample
#IDS <- read.csv(here("sims/data/samples_seed42.csv"))
# get indexes of individuals
#si <- IDS$inds
set.seed(42)
si <- sample(nrow(coords), 200, replace = FALSE)
# sample loci
l <- sample(nrow(vcf@gt), 100000)
# subset coodinates
coords <- coords[si,]
# subset vcf
vcf <- vcf[l, c(1, si + 1)]

# confirm that correct set is being used
message(paste("nloci", nrow(vcf@gt), "/ nind", nrow(coords)))

# check match
stopifnot(colnames(vcf@gt)[-1] == as.character(coords$idx))

cores <- 20
cl <- makeCluster(cores)
registerDoParallel(cl)

# run test
run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = TRUE, file.name = "WGS", stats = "allelic.richness")

stopCluster(cl)
