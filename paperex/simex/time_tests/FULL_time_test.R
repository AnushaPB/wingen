library(wingen)
library(here)
source(here("paperex", "simex", "simex_functions.R"))

set.seed(42)

# load vcf, coords, and lyr
load_middle_earth()
coords <- coords[,c("x","y")]

# confirm that correct set is being used
message(paste("nsnps", nrow(vcf@gt), "/ nind", nrow(coords)))

run_default_time_test(vcf, coords, lyr, rarify = TRUE, parallel = TRUE, ncores = 25, file.name = "FULL")

run_default_time_test(vcf, coords, lyr, rarify = FALSE, parallel = TRUE, ncores = 25, file.name = "FULL")
