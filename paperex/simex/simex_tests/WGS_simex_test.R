library(wingen)
library(here)
source(here("paperex", "simex", "simex_functions.R"))

set.seed(42)

# load vcf, coords, and lyr
load_middle_earth()

# get subsampled dataset
subdata <- subset_data(vcf, coords, nsamples = 200, nvariants = 100000)
subvcf <- subdata[["vcf"]]
subcoords <- subdata[["coords"]]

# run test
run_default_time_test(subvcf, subcoords, lyr, rarify = TRUE, parallel = TRUE, ncores = 10, file.name = "WGS")

run_default_time_test(subvcf, subcoords, lyr, rarify = FALSE, parallel = TRUE, ncores = 10, file.name = "WGS")


# run test 100

# get subsampled dataset
subdata <- subset_data(vcf, coords, nsamples = 100, nvariants = 100000)
subvcf <- subdata[["vcf"]]
subcoords <- subdata[["coords"]]

run_default_time_test(subvcf, subcoords, lyr, rarify = TRUE, parallel = TRUE, ncores = 10, file.name = "WGS")

run_default_time_test(subvcf, subcoords, lyr, rarify = FALSE, parallel = TRUE, ncores = 10, file.name = "WGS")

