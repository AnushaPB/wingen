library(wingen)
library(here)
source(here("paperex", "simex", "simex_functions.R"))

# load vcf, coords, and lyr
load_middle_earth()

# get subsampled dataset
subdata <- subset_data(vcf, coords, nsamples = 200, nvariants = 10000)
subvcf <- subdata[["vcf"]]
subcoords <- subdata[["coords"]]

# run test 200
run_default_time_test(subvcf, subcoords, lyr, rarify = TRUE, parallel = FALSE, file.name = "rr")

run_default_time_test(subvcf, subcoords, lyr, rarify = FALSE, parallel = FALSE, file.name = "rr")

# run test 100

# get subsampled dataset
subdata <- subset_data(vcf, coords, nsamples = 100, nvariants = 10000)
subvcf <- subdata[["vcf"]]
subcoords <- subdata[["coords"]]

run_default_time_test(subvcf, subcoords, lyr, rarify = TRUE, parallel = FALSE, file.name = "rr")

run_default_time_test(subvcf, subcoords, lyr, rarify = FALSE, parallel = FALSE, file.name = "rr")

