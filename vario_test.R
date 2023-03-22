library(ggplot2)
library(dplyr)
library(tidyr)

distmat <- get_geodist(lotr_coords, terra::aggregate(lotr_lyr, 5))
res <- purrr::map(c(5, 10, 15, 20, 25, 30, 35, 40),
                  ~ circle_gd(lotr_vcf[1:100,],
                              lotr_coords,
                              lotr_lyr,
                              stat = "hwe",
                              maxdist = .x,
                              fact = 5,
                              rarify = FALSE,
                              distmat = distmat,
                              parallel = TRUE,
                              ncores = 30)[[1]])

res2 <- purrr::map(c(5, 10, 15, 20, 25, 30, 35, 40),
                  ~ circle_gd(lotr_vcf[1:100,],
                              lotr_coords,
                              lotr_lyr,
                              stat = "basic_stats",
                              maxdist = .x,
                              fact = 5,
                              rarify = FALSE,
                              distmat = distmat,
                              parallel = TRUE,
                              ncores = 30)[["Fis_hierfstat"]])

r <- terra::rast(res)
r2 <- terra::rast(res2)
percNA <- terra::global(r, fun = function(x) mean(is.na(x)))
avgHWE <- terra::global(r, fun = mean, na.rm = TRUE)
avgFIS <- terra::global(r2, fun = mean, na.rm = TRUE)

df <- data.frame(radius = c(5, 10, 15, 20, 25, 30, 35, 40), percNA = percNA, avgHWE = avgHWE, avgFIS = avgFIS)
colnames(df) <- c("radius", "percNA", "avgHWE", "avgFIS")

df <-
  df %>%
  pivot_longer(-radius, names_to = "stat", values_to = "value") %>%
  filter(stat != "percNA")

ggplot(data = df, aes(x = radius, y = value, col = stat)) +
  geom_line() +
  geom_point() +
  theme_classic()

v <- vario_gd(lotr_vcf, lotr_coords)
