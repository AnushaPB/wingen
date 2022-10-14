Simulation Example from Bishop et al. (submitted)
================

``` r
library(wingen)
library(raster)
library(vcfR)
library(viridis)
library(here)
library(ggplot2)
library(dplyr)
library(purrr)
library(adegenet)

wdir <- here("paperex", "simex")
source(here(wdir, "simex_functions.R"))
```

## Load simulation results

The following function loads the simulated data (produced by
`run_simex.sh`) and subsets it for the example walkthrough.

``` r
load_middle_earth(subset = TRUE)
```

    ## nvariants 10000 / nind 200

    ## 
    ## --------------------- middle earth data ---------------------
    ##  
    ## Objects loaded: 
    ## *vcf* vcfR object (142129 loci x 1697 samples) 
    ## *coords* dataframe with x and y coordinates 
    ## *subvcf* vcfR object (10000 loci x 200 samples) 
    ## *subcoords* dataframe with x and y coordinates for 200 samples 
    ## *lyr* middle earth RasterLayer (100 x 100) 
    ## *bkg* background layer 
    ## 
    ## -------------------------------------------------------------

    ## 

## **Figure 2:** Simulation Example

### Simulation setup plots

``` r
# Make kernel density raster from FULL coordinates
kde <- raster(MASS::kde2d(coords$x, coords$y, h = c(10,10), n = 100, lims = c(0,100,-100,0)))

# Plot simulation setup (original lyr, population density, and sample distribution)
par(mar = rep(2,4))
plot(lyr, col = viridis::magma(100), main = "Carrying Capacity/Conductance", axes = FALSE, box = FALSE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plot(kde, col = viridis::magma(100), main = "Population Density", axes = FALSE, box = FALSE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
plot(lyr, col = mako(1, begin = 0.1, alpha = 0.1), main = "Distribution of Individuals", axes = FALSE, box = FALSE, zlim = c(0.01,1), legend = FALSE)
points(coords$x, coords$y, pch = 16, cex = 1, col = mako(1, begin = 0.7, alpha = 0.4), xlab = "", ylab = "")
points(subcoords$x, subcoords$y, pch = 16, cex = 1, col = mako(1, begin = 0.2, alpha = 0.9), xlab = "", ylab = "")
legend(0,-80, 
       c("All Individuals", "Sampled Individuals"), 
       pch = c(16,16), 
       col = c(mako(1, begin = 0.7), mako(1, begin = 0.2)),
       bty = "n",
       text.col = "black",
       cex = 1.5)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

**Example moving window analysis:**

``` r
set.seed(42)

# run moving window
wg <- window_gd(subvcf, 
                subcoords, 
                lyr, 
                stat = "pi", 
                rarify = TRUE, 
                wdim = 7, 
                fact = 3, 
                rarify_n = 2, 
                rarify_nit = 5, 
                parallel = FALSE)

# krig results
kg <- krig_gd(wg, lyr, index = 1, disagg_grd = 2)
```

    ## [using ordinary kriging]

``` r
kc <- krig_gd(wg, lyr, index = 2, agg_r = 2, disagg_grd = 2)
```

    ## [using ordinary kriging]

``` r
# mask areas with less than one count
mg <- mask_gd(kg, kc, minval = 2)
mg <- mask_gd(mg, bkg)

par(mar = rep(0,4))
plot_gd(wg, bkg, zlim = c(0, 0.31))
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plot_count(wg, zlim = c(0,30), horizontal = TRUE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
plot_gd(kg, zlim = c(0, 0.31))
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
plot_count(kc, zlim = c(0,30))
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
plot_gd(mg, bkg,zlim = c(0, 0.31), legend = FALSE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->

### **Figure S1:** Window vs Aggregation Factor

``` r
params <- df_to_ls(expand.grid(wdim = c(3, 5, 7), fact = c(2, 3, 4)))

stk <- purrr::map(params, test_params_simex, subvcf, subcoords, lyr)

par(mfrow = c(3, 3), mar = rep(0, 4), oma = rep(0, 4), pty = "s")
purrr::walk(stk, test_simex_plot, bkg = bkg)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## **Figure 3 & Figure S2:** Comparison of datasets, statistics, and sample sizes

``` r
params <- df_to_ls(expand.grid(datasets = c("rr", "WGS", "FULL"),
                               rarify = c("TRUE", "FALSE"),
                               stat = c("pi", "biallelic_richness", "Ho")))

# Get example layers for masking (doesn't matter which parameters other than nsamp)
msk_lyr100 <- get_divout(file.name = "rr", rarify = TRUE, stat = "pi", nsamp = 100)

msk_lyr200 <- get_divout(file.name = "rr", rarify = TRUE, stat = "pi", nsamp = 200)

# Get output raster layers
stk100 <- purrr::map(params, test_datasets_simex, nsamp = 100, msk_lyr = msk_lyr100)

stk200 <- purrr::map(params, test_datasets_simex, nsamp = 200, msk_lyr = msk_lyr200)

# Plot results (note: legends are fixed to the same scale)
par(mfrow = c(2, 3), mar = c(1, 0, 1, 0), oma = rep(0, 4))
purrr::walk(stk100, test_simex_plot, legend = FALSE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
par(mfrow = c(2, 3), mar = c(1, 0, 1, 0), oma = rep(0, 4))
purrr::walk(stk200, test_simex_plot, legend = TRUE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->

## **Figure S5:** Computational time for simulation example

``` r
# Loop reads in outputs from simex_tests functions
# To recreate these outputs use the run_sims.sh script (or run the scripts within simex_tests individually)

tdf <- dplyr::bind_rows(get_timeout("rr", rarify = "TRUE", parallel = "FALSE", nsamp = 100),
                        get_timeout("rr", rarify = "TRUE", parallel = "FALSE", nsamp = 200),
                        get_timeout("WGS", rarify = "TRUE", parallel = "TRUE", nsamp = 100),
                        get_timeout("WGS", rarify = "TRUE", parallel = "TRUE", nsamp = 200)
                        )

tdf[tdf$dataset == "rr", "dataset"] <- "(a)"
tdf[tdf$dataset == "WGS", "dataset"] <- "(b)"

ggplot(data = tdf, aes(x = factor(nsamp), y = time, fill = stat)) +
  geom_hline(yintercept = 60, linetype = "dashed", col = "grey80", lwd = 1) + 
  geom_text(aes(factor(200), 60, label = "1 min", vjust = -1, hjust = -0.7), col = "grey80", fontface = "italic") +
  geom_hline(yintercept = 60*5, linetype = "dashed", col = "grey70", lwd = 1) + 
  geom_text(aes(factor(200), 60*5, label = "5 min", vjust = -1, hjust = -0.7), col = "grey70", fontface = "italic") +
  geom_hline(yintercept = 60*10, linetype = "dashed", col = "grey60", lwd = 1) + 
  geom_text(aes(factor(200), 60*10, label = "10 min", vjust = -1, hjust = -0.4), col = "grey60", fontface = "italic") +
  geom_hline(yintercept = 60*15, linetype = "dashed", col = "grey50", lwd = 1) + 
  geom_text(aes(factor(200), 60*15, label = "15 min", vjust = -1, hjust = -0.4), col = "grey50", fontface = "italic") +
  geom_col(position=position_dodge()) +
  geom_text(aes(label = round(time, 0), col = stat), 
            vjust = -0.5, position=position_dodge(width = .9)) + 
  scale_color_manual(values=c("pi"= mako(3, begin = 0.3, end = 0.8)[1], 
                              "allelic richness"= mako(3, begin = 0.3, end = 0.8)[2],
                              "Ho" = mako(3, begin = 0.3, end = 0.8)[3]), 
                    labels = c("pi" = expression(pi),
                               "allelic richness" = bquote(A[R]),
                               "Ho" = bquote(H[O]))) +
  scale_fill_manual(values=c("pi"=mako(3, begin = 0.3, end = 0.8)[1], 
                              "allelic richness"= mako(3, begin = 0.3, end = 0.8)[2],
                              "Ho"= mako(3, begin = 0.3, end = 0.8)[3]), 
                    labels = c("pi" = expression(pi),
                               "allelic richness" = bquote(A[R]),
                               "Ho" = bquote(H[O]))) +
  guides(color = guide_legend(override.aes = list(color = rgb(0,0,0,0)))) +
  facet_grid(~dataset,  scales = "free_y") +
  xlab("Number of Samples") +
  ylab("Time (seconds)") +
  ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, hjust = 0),
        legend.title=element_blank(),
        legend.text = element_text(size=12))
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
