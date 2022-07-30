
lyr <- read.csv("data/middle_earth.csv", header = FALSE)
lyr <- raster(as.matrix(lyr))
extent(lyr) <- extent(0,100,-100,0)
bkg <- lyr
bkg[bkg < 0.01] <- NA

for (n in c(1000, 5000, 10000)){
  vcf <- read.vcfR(paste0("sims/data/mod-sim_params_it-0_t-", n, "_spp-spp_0.vcf"))
  geo <- read.csv(paste0("sims/data/mod-sim_params_it-0_t-", n, "_spp-spp_0.csv"))
  coords <- geo[,c("idx","x","y")]
  coords$y <- -coords$y

  set.seed(42)

  message("creating new file")
  si <- sample(nrow(coords), 200)

  l <- sample(nrow(vcf@gt), 10000)

  subvcf <- vcf[l, c(1, si + 1)]

  subcoords <- coords[si,]
  subcoords <- subcoords[,c("x","y")]

  # check match
  all(colnames(subvcf@gt)[-1] == as.character(geo$idx[si]))

  for(m in c("pi", "het", "biallelic.richness")){
    tr <- sim(subvcf, subcoords, lyr, stat = m, wdim = 3, fact = 3, rarify_n = 2, rarify_nit = 5, rarify = FALSE, parallel = TRUE, nloci = nrow(subvcf@gt))

    tr_k <- krig_gd(tr[["res"]][[1]], lyr, disagg_grd = 2)

    writeRaster(tr_k, here("sims", "outputs", paste0(m,"_",n,"krig.tif")), overwrite = TRUE)
  }
}


for (n in c(1000, 5000, 10000)){

  geo <- read.csv(paste0("sims/data/mod-sim_params_it-0_t-", n, "_spp-spp_0.csv"))
  coords <- geo[,c("idx","x","y")]
  coords$y <- -coords$y
  kde <- raster(MASS::kde2d(coords$x, coords$y, h = c(10,10), n = 100, lims = c(0,100,-100,0)))

  plot_gd(kde, breaks = 10)
  writeRaster(kde, here("sims", "outputs", paste0(n,"_kde.tif")), overwrite = TRUE)
}







