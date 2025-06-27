library(wingen)
library(purrr)
library(tidyverse)
library(terra)

load_middle_earth_ex()

r <- window_gd(lotr_vcf, lotr_coords, lotr_lyr, fact = 3, wdim = 5, stat = "Ho", na.rm = TRUE)
plot_gd(r)
r <- r[[1]]

kr <- krig_gd(r, lotr_lyr)
mr <- mask_gd(kr, lotr_range)
plot_gd(mr)

library(terra)
library(sf)
library(gstat)
library(dplyr)

# — 1. Your input: a SpatRaster of window‐based diversity (“r”)
#    (e.g. the output from wingen::window_gd)

# 1.1 Turn that into a data.frame with x,y,layer
pts_df <- terra::as.data.frame(r, xy=TRUE, na.rm=TRUE)

# 1.2 Convert to sf points
pts_sf <- pts_df %>%
  st_as_sf(coords = c("x","y"), crs = crs(r)) 

# — 2. Build your prediction grid
#    If you already have a “grd” SpatRaster, substitute it here.  
#    Otherwise just reuse r’s extent & resolution:

grd <- r
grd_df <- terra::as.data.frame(grd, xy=TRUE, na.rm=FALSE)

grd_sf <- grd_df %>%
  st_as_sf(coords = c("x","y"), crs = crs(grd))

# — 3. Run IDW with sf inputs
#    “layer” is the column in pts_sf with your diversity value
idw_sf <- gstat::idw(
  formula  = Ho ~ 1,
  locations = pts_sf,
  newdata    = grd_sf,
  idp        = 2
)

# — 4. Convert the result back to a terra raster
#    idw_sf has geometry + a “var1.pred” column

# 4.1 pull out coordinates + predictions
out_df <- idw_sf %>%
  st_coordinates() %>% 
  as.data.frame() %>%
  bind_cols(pred = idw_sf$var1.pred)

# 4.2 make a raster from XYZ
r_idw <- rast(out_df, crs = crs(r))
r_idw <- mask(r_idw, lotr_range)  # Mask to the range if needed

library(fields)
pts_df <- terra::as.data.frame(r, xy=TRUE, na.rm=TRUE)

tps_mod <- Tps(pts_df[,c("x","y")], pts_df$Ho)
grd_df  <- terra::as.data.frame(lotr_lyr, xy=TRUE)
preds   <- predict.Krig(tps_mod, grd_df[,c("x","y")])

r_tps <- rast(cbind(grd_df[,1:2], layer=preds))
r_tps <- mask(r_tps, lotr_range)  # Mask to the range if needed

library(mgcv)
pts_df <- terra::as.data.frame(r, xy=TRUE, na.rm=TRUE)

gam_mod <- gam(Ho ~ s(x,y), data = pts_df)
grd_df  <- terra::as.data.frame(lotr_lyr, xy=TRUE)
grd_df$pred <- predict(gam_mod, newdata = grd_df)

r_gam <- rast(cbind(grd_df[,1:2], layer=grd_df$pred))
r_gam <- mask(r_gam, lotr_range)  # Mask to the range if needed


# — 1. Your windowed raster (response) and new predictor raster
r_response  <- r           # your diversity SpatRaster
r_predictor <- aggregate(lotr_lyr, 3)  # your environmental covariate

# — 2. Turn response → data.frame with coords + layer
pts_df <- as.data.frame(r_response, xy=TRUE, na.rm=TRUE)
names(pts_df) <- c("x","y","layer")

# — 3. Extract predictor at those same coords
#    This gives you one value of your env covariate per point
env_vals     <- terra::extract(r_predictor, pts_df[,c("x","y")], ID=FALSE)
pts_df$env   <- env_vals   # rename column as you like

# — 4. Fit GAM with a smooth on space + effect of your env surface
#    You can model env linearly or with its own smooth:
gam_mod <- gam(
  layer ~ s(x, y, k=100)     # spatial smoother
        + env,      # smooth effect of env
  data = pts_df
)

# — 5. Build your prediction grid (same extent & res as you used before)
grd_df <- as.data.frame(lotr_lyr, xy=TRUE, na.rm=FALSE)
# extract env on grid
grd_df$env <- terra::extract(lotr_lyr, grd_df[,c("x","y")], ID=FALSE)

# — 6. Predict across the landscape
grd_df$pred <- predict(gam_mod, newdata = grd_df)

# — 7. Turn back into a raster and plot
r_gam2 <- rast(cbind(grd_df[,1:2], layer=grd_df$pred))
crs(r_gam2) <- crs(r_response)
r_gam2 <- mask(r_gam2, lotr_range)  # Mask to the range if needed
plot(r_gam2, main="GAM w/ spatial + env")

# — 5. Visualize
par(mfrow = c(1, 6))
plot(r)
plot(r_idw, main = "IDW Smoothed Diversity")
plot(mr, main = "Kriged Diversity")
plot(r_tps, main="Thin‐plate spline")
plot(r_gam, main="GAM smooth")
plot(r_gam2, main="GAM smooth + env")



