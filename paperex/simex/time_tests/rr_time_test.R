library(wingen)
library(foreach)
library(doParallel)
library(here)
source(here("sims/sims_functions.R"))


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

# check match
stopifnot(colnames(vcf@gt)[-1] == as.character(coords$idx))

# confirm that correct set is being used
message(paste("nsnps", nrow(vcf@gt), "/ nind", nrow(coords)))

# run test 200
run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = FALSE, file.name = "rr")

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = FALSE, parallel = FALSE, file.name = "rr")

# run test 100

# subset coordinates (using first 100 of already subsampled dataframe)
coords <- coords[1:100,]
# subset vcf (loci have already been subset previously)
vcf <- vcf[, c(1, 1:100 + 1)]

# check match
stopifnot(colnames(vcf@gt)[-1] == as.character(coords$idx))

# confirm that correct set is being used
message(paste("nsnps", nrow(vcf@gt), "/ nind", nrow(coords)))

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = FALSE, file.name = "rr")

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = FALSE, parallel = FALSE, file.name = "rr")

