library(wingen)
library(cowplot)
library(ggplot2)
library(here)
source("krig_gd2.R")

# Load example Middle Earth data
load_middle_earth_ex()


# Window-based genetic diversity with small windows/sparse sampling
set.seed(42)
wgd <- window_gd(
  lotr_vcf,
  lotr_coords,
  lotr_lyr,
  stat = "pi",
  wdim = 3,
  fact = 2,
  rarify = TRUE 
)

# Note: everything plotted together at the end

# --- Original kriging using krig_gd() ---

kgd1A <- krig_gd(wgd[["pi"]]) %>% 
  mask(lotr_range) %>%
  ggplot_gd(bkg = lotr_range) + 
  labs(fill = "Original\nKriging\n(res = 2)")

kgd1B <- krig_gd(wgd[["pi"]], grd = lotr_lyr) %>% 
  mask(lotr_range) %>% 
  ggplot_gd(bkg = lotr_range) + 
  labs(fill = "Original\nKriging\n(res = 1)")

# --- New kriging using krig_gd2() ---

kgd2A <- krig_gd2(wgd[["pi"]]) %>% 
  mask(lotr_range) %>% 
  ggplot_gd(bkg = lotr_range) + 
  labs(fill = "New\nKriging\n(res = 2)")

kgd2B <- krig_gd2(wgd[["pi"]], grd = lotr_lyr) %>% 
  mask(lotr_range) %>% 
  ggplot_gd(bkg = lotr_range) + 
  labs(fill = "New\nKriging\n(res = 1)")

# --- New kriging using krig_gd() with weights based on sample count ---

kgd3A <- krig_gd2(wgd[["pi"]], weight_r = wgd[["sample_count"]]) %>% 
  mask(lotr_range) %>% 
  ggplot_gd(bkg = lotr_range) + 
  labs(fill = "New\nKriging\nWeighted\n(res = 2)")

kgd3B <- krig_gd2( wgd[["pi"]], grd = lotr_lyr, weight_r = wgd[["sample_count"]]) %>% 
  mask(lotr_range) %>% 
  ggplot_gd(bkg = lotr_range) + 
  labs(fill = "New\nKriging\nWeighted\n(res = 1)")

# Plotting
wgd_plot <- ggplot_gd(wgd, bkg = lotr_range) + 
  labs(fill = "Original\nwindow\n(res = 2)")
row1 <- plot_grid(kgd1A, kgd2A, kgd3A, nrow = 1, align = "hv")
row2 <- plot_grid(kgd1B, kgd2B, kgd3B, nrow = 1, align = "hv")
group1 <- plot_grid(row1, row2, nrow = 2, align = "hv")

plt <- plot_grid(wgd_plot, group1, nrow = 1, rel_widths = c(1, 3), align = "hv")
print(plt)

# Export PNG
png(
  here("krig_gd2_example.png"),
  width = 15 * 300, height = 5 * 300, res = 300
)
print(plt)
dev.off()
