Simulation Example
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

The following function loads the results from the simulation and subsets
the data for the example walkthrough.

``` r
load_middle_earth(subset = TRUE)
```

## Figure 2: Simulation Example

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
plot(kde, col = viridis::magma(10), main = "Population Density", axes = FALSE, box = FALSE)
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

### Figure \#X: Window vs Aggregation Factor

``` r
params <- df_to_ls(expand.grid(wdim = c(3, 5, 7), fact = c(2, 3, 4)))

stk <- purrr::map(params, test_params_simex, subvcf, subcoords, lyr)

par(mfrow = c(3, 3), mar = rep(0, 4), oma = rep(0, 4), pty = "s")
purrr::walk(stk, test_simex_plot, bkg = bkg)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Figure 3 & Figure 1S: Comparison of datasets, statistics, and sample sizes

``` r
params <- df_to_ls(expand.grid(datasets = c("rr", "WGS", "FULL"),
                               rarify = c("TRUE", "FALSE"),
                               stat = c("pi", "biallelic_richness", "heterozygosity")))

# Get example layers for masking (doesn't matter which parameters other than nsamp)
msk_lyr100 <- get_divout(file.name = "rr", rarify = TRUE, stat = "pi", nsamp = 100)

msk_lyr200 <- get_divout(file.name = "rr", rarify = TRUE, stat = "pi", nsamp = 200)

# Get output raster layers
stk100 <- purrr::map(params, test_datasets_simex, nsamp = 100, msk_lyr = msk_lyr100)

stk200 <- purrr::map(params, test_datasets_simex, nsamp = 200, msk_lyr = msk_lyr200)

# Plot results (note: legends are fixed to the same scale)
par(mfrow = c(2, 3), mar = rep(0, 4), oma = rep(0, 4), pty = "s")
purrr::walk(stk100, test_simex_plot, bkg = bkg, legend = FALSE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
par(mfrow = c(2, 3), mar = c(1, 0, 1, 0), oma = rep(0, 4), pty = "s")
purrr::walk(stk200, test_simex_plot, bkg = bkg, legend = FALSE)
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->![](simex_notebook_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->

## Figure \#X: Timing

``` r
# Loop reads in outputs from time_tests functions
# To recreate these outputs use the run_sims.sh script (or run the scripts within time_tests individually)

tdf <- dplyr::bind_rows(get_timeout("rr", rarify = "TRUE", parallel = "FALSE", nsamp = 100),
                        get_timeout("rr", rarify = "TRUE", parallel = "FALSE", nsamp = 200),
                        get_timeout("WGS", rarify = "TRUE", parallel = "TRUE", nsamp = 100),
                        get_timeout("WGS", rarify = "TRUE", parallel = "TRUE", nsamp = 200)
                        )

tdf[tdf$dataset == "rr", "dataset"] <- "10,000 loci (w/o Parallelization)"
tdf[tdf$dataset == "WGS", "dataset"] <- "100,000 loci (w/ Parallelization)"

ggplot(data = tdf, aes(x = factor(nsamp), y = time, fill = stat)) +
  geom_hline(yintercept = 60, linetype = "dashed", col = "darkgray", lwd = 1) + 
  geom_text(aes(factor(200), 60, label = "1 min", vjust = -1, hjust = -0.7), col = "gray", fontface = "italic") +
  geom_hline(yintercept = 120, linetype = "dashed", col = "gray", lwd = 1) + 
  geom_text(aes(factor(200), 120, label = "2 min", vjust = -1, hjust = -0.7), col = "lightgray", fontface = "italic") +
  geom_hline(yintercept = 180, linetype = "dashed", col = "lightgray", lwd = 1) + 
  geom_text(aes(factor(200), 180, label = "3 min", vjust = -1, hjust = -0.7), col = "lightgray", fontface = "italic") +
  geom_col(position=position_dodge()) +
  geom_text(aes(label = round(time, 0), col = stat), 
            vjust = -0.5, position=position_dodge(width = .9)) + 
  scale_color_manual(values=c("pi"=mako(3, begin = 0.3, end = 0.8)[1], 
                              "allelic richness"=mako(3, begin = 0.3, end = 0.8)[2],
                              "heterozygosity"=mako(3, begin = 0.3, end = 0.8)[3])) +
  scale_fill_manual(values=c("pi"=mako(3, begin = 0.3, end = 0.8)[1], 
                              "allelic richness"=mako(3, begin = 0.3, end = 0.8)[2],
                              "heterozygosity"=mako(3, begin = 0.3, end = 0.8)[3])) +
  guides(color = guide_legend(override.aes = list(color = rgb(0,0,0,0)))) +
  facet_grid(~dataset,  scales = "free_y") +
  xlab("number of samples") +
  ylab("time (seconds)") +
  ylim(0,200) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank())
```

![](simex_notebook_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->