# wingen 2.1.1

## Changes to krig_gd following removal of sp dependency
The legacy packages maptools, rgdal, and rgeos, underpinning the sp package will retire in October 2023, so we have made steps to remove sp as a required package to avoid any associated issues since sp is no longer being actively developed.
- the grd argument now only accepts either Raster or SpatRasters



# wingen 2.0.1

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8125610.svg)](https://doi.org/10.5281/zenodo.8125610)

### Deprecation of `parallel` and `ncores` arguments

The `parallel` and `ncores` arguments have been deprecated for all functions. 

Instead, `future::plan()` should be used to setup parallelization (see the package vignette for more details).

# wingen 2.0.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8125610.svg)](https://doi.org/10.5281/zenodo.8125610)

### Major updates
- Circular and resistance-based windows are now possible using the `circle_gd()`/`circle_general()` and `resist_gd()`/`resist_general()`( functions
- `preview_gd()` takes a new argument (`method`) to specify which method to generate a preview of (`"window"`, `"circle"`, or `"resist"`)
- Window size across the landscape can now be varied by providing a raster for `maxdist` for  `circle_gd()`/`resist_gd()`
- New population genetic statistic options (specified using `stat`)
- Create ggplots using `ggplot_gd()` and `ggplot_count()`

### Minor updates
- All moving window functions can now accept more than one statistic at a time in the form of a vector
- `plot_gd()` will now generate plots for all layers except the sample count layer and vise versa for `plot_count()`

# wingen 1.1.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7637712.svg)](https://doi.org/10.5281/zenodo.7637712)

### Update to terra and sf
- The guts of wingen have been updated from raster to terra and now incorporate sf
- Rasters can now be provided as terra SpatRasters or raster RasterLayers
- Coordinates can now be provided as sf points, a two-column matrix, or a data.frame representing x and y coordinates.

### Breaking changes
- Wingen is no longer agnostic to projection. Checks have been added that will give warnings for missing CRS and errors for mismatched CRS. The vignette and docs have also been updated to discuss the effects of projection and provide guidance on transformations.
- The first argument of `window_gd()` is now called `gen` instead of `vcf`. The argument is otherwise the same (i.e., takes the same input), however the name was generalized such that future versions could potentially take different types of genetic data without confusion.
- The second argument of `mask_gd()` is now called `y` instead of `mask` to avoid conflation with the function `mask()`
- `mask_gd()` will longer provide the option to resample rasters, rasters must be transformed before the function is used

# wingen 1.0.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7199558.svg)](https://doi.org/10.5281/zenodo.7199558)

Initial release of the wingen R package
