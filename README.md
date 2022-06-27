
<!-- README.md is generated from README.Rmd. Please edit that file -->
<STYLE type='text/css' scoped>
PRE.fansi SPAN {padding-top: .25em; padding-bottom: .25em};
</STYLE>

# wingen <img src="man/figures/logo.png" align="right" height=150 />

<!-- badges: start -->

[![R-CMD-check](https://github.com/AnushaPB/wingen/actions/workflows/check-release.yaml/badge.svg)](https://github.com/AnushaPB/wingen/actions/workflows/check-release.yaml)
<!-- badges: end -->

Create maps of genetic diversity using a sliding window approach.

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

<PRE class="fansi fansi-output"><CODE><span style='color: #00BBBB; font-weight: bold;'>
------------ middle earth example ------------
</span> <span style='color: #0000BB;'>
Objects loaded:</span> <span style='color: #BBBB00; font-weight: bold;'>
*lotr_vcf*</span> <span style='color: #BBBB00;'>vcfR object (1000 loci x 200 samples)</span> <span style='color: #00BB00; font-weight: bold;'>
*lotr_coords*</span> <span style='color: #00BB00;'>dataframe with x and y coordinates</span> <span style='color: #BB00BB; font-weight: bold;'>
*lotr_lyr*</span> <span style='color: #BB00BB;'>middle earth RasterLayer (100 x 100)</span> <span style='color: #00BBBB; font-weight: bold;'>

----------------------------------------------
</span>
</CODE></PRE>

``` r
# Run sliding window calculations of pi with rarefaction
wgd <- window_gd(lotr_vcf,
          lotr_coords,
          lotr_lyr,
          stat = "pi",
          fact = 3,
          wdim = 5,
          rarify = TRUE,
          nloci = 1000)

# Krige results
kgd <- krig_gd(wgd, lotr_lyr)

# Mask results
mgd <- mask_gd(kgd, min_n = 3)

# Plot results
par(mfrow = c(1,3), oma = rep(2,4), mar = rep(2,4))
plot_gd(wgd, main = "Window pi")
plot_gd(kgd, main = "Kriged pi")
plot_gd(mgd, main = "Kriged & masked pi")
```

<img src="man/figures/README-example-1.png" width="100%" />

For an extended example check out the package vignette:

``` r
vignette("wingen-vignette")
```
