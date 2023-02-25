
# IGNORE ---
```{r}
con_lyr <- mini_lyr*0+0.5
con_lyr[1:50] <- 1
# Create transition surface
trSurface <- gdistance::transition(con_lyr, transitionFunction = mean, directions = 8)
trSurface <- gdistance::geoCorrection(trSurface, type = "c", scl = TRUE)
#Make spatial points
# TODO: Figure out if CRS is needed here
sp <- sp::SpatialPoints(rbind(c(50,-10), c(52,-90)))
# Calculate circuit distances
distmat <- possible_gdist(trSurface, sp)

trSurface2 <- gdistance::transition(disaggregate(raster(trSurface), 2),  transitionFunction = mean, directions = 8)
#Make spatial points
# Calculate circuit distances
distmat2 <- possible_gdist(trSurface2, sp)
distmat/distmat2

trSurface2*0+1
plot(raster(trSurface2))
points(sp)
plot(raster(trSurface))
points(sp)
```

