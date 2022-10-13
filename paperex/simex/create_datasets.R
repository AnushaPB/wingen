
# load genetic data ---------------------------------------------------------------------------------------------------

# check if the vcf file from the simulations exists locally and if not download it from Zenodo
file <- here::here("paperex", "simex", "data", "mod-sim_params_it-0_t-1000_spp-spp_0.vcf")
if(!file.exists(file)){
  message("downloading vcf and storing locally...this will take some time, but only has to be done once")
  download.file("https://zenodo.org/record/7112468/files/mod-sim_params_it-0_t-1000_spp-spp_0.vcf?download=1", file)
}

vcf <- vcfR::read.vcfR(file, verbose = FALSE)
assign("vcf", vcf, envir = .GlobalEnv)

# Create samples of individuals for datasets --------------------------------------------------------------------------

set.seed(42)
# (subtract 1 because vcf first col is FORMAT)
s200 <- sample(ncol(vcf@gt) - 1, 200, replace = FALSE)
s100 <- s200[1:100]

write.csv(s200, here::here("paperex", "simex", "data", "samples_seed42_200.csv"), row.names = FALSE)
write.csv(s100, here::here("paperex", "simex", "data", "samples_seed42_100.csv"), row.names = FALSE)

# Create samples of variants for datasets -----------------------------------------------------------------------------

subvcf100 <- vcf[, c(1, s100 + 1)]
# only retain polymorphic sites in smallest subsample
polysites <- vcfR::is.polymorphic(subvcf100)

set.seed(42)
var100000 <- sample(which(polysites), 100000, replace = FALSE)
var10000 <- var100000[1:10000]
write.csv(var10000, here::here("paperex", "simex", "data", "variants_seed42_10000.csv"), row.names = FALSE)
write.csv(var100000, here::here("paperex", "simex", "data", "variants_seed42_100000.csv"), row.names = FALSE)
