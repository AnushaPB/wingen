library(wingen)
library(foreach)
library(doParallel)
library(here)
source(here("sims/sim_functions.R"))

set.seed(42)

load_middle_earth()

IDS <- read.csv(here("sims/data/samples_seed42.csv"))
# get indexes of individuals
si <- IDS$inds
# sample loci
l <- sample(nrow(vcf@gt), 100000)
# subset coodinates
coords <- coords[si,]
# subset vcf
vcf <- vcf[l,c(1, si + 1)]

# check match
all(colnames(vcf@gt)[-1] == as.character(coords$idx))

# run test
run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = FALSE, file.name = "WGS")

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = FALSE, parallel = FALSE, file.name = "WGS")



