library(wingen)
library(foreach)
library(doParallel)
library(here)
source(here("sims/sim_functions.R"))

set.seed(42)

# load vcf, coords, and lyr
load_middle_earth()
coords <- coords[,c("x","y")]

# confirm that correct set is being used
message(paste("nloci", nrow(vcf@gt), "/ nind", nrow(coords)))

cores <- 10
cl <- makeCluster(cores)
registerDoParallel(cl)

run_default_time_test(vcf, coords, lyr, rarify = TRUE, parallel = TRUE, file.name = "FULL")

run_default_time_test(vcf, coords, lyr, rarify = FALSE, parallel = TRUE, file.name = "FULL")

stopCluster(cl)
