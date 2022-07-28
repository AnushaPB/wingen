library(wingen)
library(foreach)
library(doParallel)
library(here)
source(here("sims/sims_functions.R"))

set.seed(42)

# load vcf, coords, and lyr
load_middle_earth()

IDS <- read.csv(here("sims/data/samples_seed42.csv"))
# get indexes of individuals
si <- IDS$inds
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

cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

# run test
run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = TRUE, file.name = "WGS")

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = FALSE, parallel = TRUE, file.name = "WGS")

stopCluster(cl)


# run test 100

# subset coordinates (using first 100 of already subsampled dataframe)
coords <- coords[1:100,]
# subset vcf (loci have already been subset previously)
vcf <- vcf[, c(1, 1:100 + 1)]

# check match
stopifnot(colnames(vcf@gt)[-1] == as.character(coords$idx))

# confirm that correct set is being used
message(paste("nloci", nrow(vcf@gt), "/ nind", nrow(coords)))

cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = TRUE, file.name = "WGS")

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = FALSE, parallel = TRUE, file.name = "WGS")

stopCluster(cl)
