library(wingen)
library(here)
source(here("paperex", "simex", "simex_functions.R"))

# load vcf, coords, and lyr
load_middle_earth()

## RR TEST -----------------------------------------------------------------------------------------------

# get subsampled dataset
subdata <- subset_data(vcf, coords, nsamples = 200, nvariants = 10000)
subvcf <- subdata[["vcf"]]
subcoords <- subdata[["coords"]]

# run test
run_default_time_test(subvcf, subcoords, lyr, rarify = TRUE, parallel = TRUE, ncores = 20, file.name = "AR", stats = "allelic_richness")

## WGS TEST -----------------------------------------------------------------------------------------------

# get subsampled dataset
subdata <- subset_data(vcf, coords, nsamples = 200, nvariants = 100000)
subvcf <- subdata[["vcf"]]
subcoords <- subdata[["coords"]]

# run test
run_default_time_test(subvcf, subcoords, lyr, rarify = TRUE, parallel = TRUE, ncores = 30, file.name = "AR", stats = "allelic.richness")
