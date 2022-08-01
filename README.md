
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wingen <img src="man/figures/logo.png" align="right" height=150 />

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/AnushaPB/wingen/actions/workflows/check-release.yaml/badge.svg)](https://github.com/AnushaPB/wingen/actions/workflows/check-release.yaml) -->

[![codecov](https://codecov.io/gh/AnushaPB/wingen/branch/main/graph/badge.svg?token=P4Z35HFR4Y)](https://codecov.io/gh/AnushaPB/wingen)
[![build](https://github.com/AnushaPB/wingen/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/AnushaPB/wingen/actions/workflows/test-coverage.yaml)
[![license:
MIT](https://img.shields.io/badge/license-MIT-blue)](https://img.shields.io/badge/license-MIT-blue)
<!-- badges: end -->

Create maps of genetic diversity using a moving window approach.

## Installation

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AnushaPB/wingen")
```

## Example

``` r
library(wingen)

# load example data
load_middle_earth_ex()
```


    ------------- middle earth example -------------
     
    Objects loaded: 
    *lotr_vcf* vcfR object (100 loci x 100 samples) 
    *lotr_coords* dataframe with x and y coordinates 
    *lotr_lyr* middle earth RasterLayer (100 x 100) 

    ------------------------------------------------

``` r
# Run moving window calculations of pi with rarefaction
wgd <- window_gd(lotr_vcf,
          lotr_coords,
          lotr_lyr,
          stat = "pi",
          wdim = 3,
          fact = 5,
          rarify = TRUE,
          nloci = 1000)

# Krige results
kgd <- krig_gd(wgd, lotr_lyr)

# Mask results
mgd <- mask_gd(kgd, lotr_lyr, minval = 0.01)

# Plot results
par(mfrow = c(1,4), oma = rep(2,4), mar = rep(2,4))
plot_gd(wgd, bkg = mgd,  main = "Window pi")
plot_gd(kgd, main = "Kriged pi")
plot_gd(mgd, main = "Kriged & masked pi")
plot_count(wgd, main = "Window counts")
```

<img src="man/figures/README-example-1.png" width="100%" />

For an extended example check out the package vignette:

``` r
vignette("wingen-vignette")
```
