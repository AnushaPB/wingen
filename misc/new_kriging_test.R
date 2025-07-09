library(wingen)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(terra)
library(here)

load_middle_earth_ex()

setwd(here("misc"))

# ORIGINAL DATA ANALYSIS ---------------------------------------------------
lotr_vcf <- lotr_vcf # already loaded in load_middle_earth_ex

wgd <- window_gd(lotr_vcf,
  lotr_coords,
  lotr_lyr,
  stat = "pi",
  wdim = 3,
  fact = 3,
  rarify = TRUE
)

test_wgd <- function(lotr_vcf, lotr_coords, lotr_lyr, wdim = 3, fact = 2){
  combos <- expand.grid(
    wdim = wdim,
    fact = fact,
    stringsAsFactors = FALSE
  )

  pmap(combos, function(wdim, fact) {
    window_gd(
      lotr_vcf,
      lotr_coords,
      lotr_lyr,
      stat = "pi",
      wdim = wdim,
      fact = fact,
      rarify = TRUE
    )
  }, .progress = TRUE)
}

test_krig <- function(wgd, fact) {
  map(fact, ~test_krig_helper(wgd, fact = .x), .progress = TRUE)
}

test_krig_helper <- function(wgd, fact = 1){
  grd <- aggregate(lotr_lyr, fact = fact)
  original_kg <- krig_gd(wgd, grd = grd)
  new_kg <- krig_gd2(wgd, grd = grd, weight_r = wgd[[2]])
  new_kg_weighted <- krig_gd2(wgd, grd = grd, weight_r = wgd[[2]])
  kg_list <- list(original_kg, new_kg, new_kg_weighted)
  names(kg_list) <- c("Original\nkriging", "New\nkriging", "New\nkriging\nweighted")
  rast(kg_list) %>% mask(lotr_range)
}

# UNBIASED ANALYSIS --------------------------------------------------------
wdim_vals <- c("wdim3", "wdim5", "wdim7")
fact_vals <- c(1, 2, 3)

wgd_list <- test_wgd(lotr_vcf, lotr_coords, lotr_lyr, wdim = c(3, 5, 7), fact = 2)
names(wgd_list) <- wdim_vals

wgd_kgd_list <- map(wgd_list, ~test_krig(.x, fact = fact_vals))

# Plot
mega_plt <- map(wdim_vals, ~{
  wdim <- .x
  plt_group <- map(fact_vals, \(fact) {
    krig_plots <- plot_grid(
      plotlist = ggplot_gd(wgd_kgd_list[[wdim]][[fact]], bkg = lotr_range),
      nrow = 1
    )
    wgd_plot <- ggplot_gd(wgd_list[[wdim]], bkg = lotr_range)
    plot_grid(wgd_plot, krig_plots, nrow = 1, rel_widths = c(1, 3))
  })
  plot_grid(plotlist = plt_group, ncol = 1,
            labels = paste("Kriged resolution:", fact_vals),
            label_size = 10, align = "hv")
})

mega_plt_labeled <- map2(mega_plt, wdim_vals, ~plot_grid(
  ggdraw() + draw_label(.y, fontface = "bold", size = 16, vjust = 1),
  .x,
  ncol = 1,
  rel_heights = c(0.05, 1)
))

pdf(here("kriging_unbiased_results.pdf"), width = 11, height = 7)
walk(mega_plt_labeled, print)
dev.off()


# BIASED DATA ANALYSIS -----------------------------------------------------
load("data/lotr_biased_coords.rda")
load("data/lotr_biased_vcf.rda")

# Run analyses on biased data
wgd_list_biased <- test_wgd(lotr_biased_vcf, lotr_biased_coords, lotr_lyr, wdim = c(3, 5, 7), fact = 2)
names(wgd_list_biased) <- wdim_vals

wgd_kgd_list_biased <- map(wgd_list_biased, ~test_krig(.x, fact = fact_vals))

# Plot
mega_plt_biased <- map(wdim_vals, ~{
  wdim <- .x
  plt_group <- map(fact_vals, \(fact) {
    krig_plots <- plot_grid(
      plotlist = ggplot_gd(wgd_kgd_list_biased[[wdim]][[fact]], bkg = lotr_range),
      nrow = 1
    )
    wgd_plot <- ggplot_gd(wgd_list_biased[[wdim]], bkg = lotr_range)
    plot_grid(wgd_plot, krig_plots, nrow = 1, rel_widths = c(1, 3))
  })
  plot_grid(plotlist = plt_group, ncol = 1,
            labels = paste("Kriged resolution:", fact_vals),
            label_size = 10, align = "hv")
})

mega_plt_biased_labeled <- map2(mega_plt_biased, wdim_vals, ~plot_grid(
  ggdraw() + draw_label(.y, fontface = "bold", size = 16, vjust = 1),
  .x,
  ncol = 1,
  rel_heights = c(0.05, 1)
))

pdf(here("kriging_biased_results.pdf"), width = 11, height = 7)
walk(mega_plt_biased_labeled, print)
dev.off()

