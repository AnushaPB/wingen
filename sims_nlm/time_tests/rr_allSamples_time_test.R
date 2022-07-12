library(wingen)
library(foreach)
library(doParallel)
library(here)
source(here("sims/sim_functions.R"))

set.seed(42)

# load vcf, coords, and lyr
load_middle_earth()

# sample loci
l <- sample(nrow(vcf@gt), 100000)
# subset vcf
vcf <- vcf[l, ]

# confirm that correct set is being used
message(paste("nloci", nrow(vcf@gt), "/ nind", nrow(coords)))

# check match
all(colnames(vcf@gt)[-1] == as.character(coords$idx))

# run test
run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = FALSE, file.name = "rr_allSamples")

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = FALSE, parallel = FALSE, file.name = "rr_allSamples")



