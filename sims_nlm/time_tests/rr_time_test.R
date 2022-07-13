library(wingen)
library(foreach)
library(doParallel)
library(here)
source(here("sims/sim_functions.R"))

set.seed(42)

# load vcf, coords, and lyr
load_middle_earth()

# sample inds
si <- sample(nrow(coords), 200)
# sample loci
l <- sample(nrow(vcf@gt), 10000)
# subset coodinates
coords <- coords[si,]
# subset vcf
vcf <- vcf[l, c(1, si + 1)]

# check match
stopifnot(colnames(vcf@gt)[-1] == as.character(coords$idx))

# confirm that correct set is being used
message(paste("nloci", nrow(vcf@gt), "/ nind", nrow(coords)))

# run test

cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = TRUE, parallel = TRUE, file.name = "rr")

run_default_time_test(vcf, coords[,c("x","y")], lyr, rarify = FALSE, parallel = TRUE, file.name = "rr")

stopCluster(cl)
