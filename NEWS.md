# wingen 2.1.2

Patch to fix `L` argument defaults and improve documentation. 

The wingen default for all functions is now `L = "nvariants"`, which sets `L` to the number of variants in the VCF when calculating pi. This was already the default in the original `window_gd()` function, but `circle_gd()` and `resist_gd()` had `L = NULL` as the default, which returns the sum over SNPs of nucleotide diversity. We have also improved documentation for `L`.

# wingen 2.1.1

Minor changes for CRAN

### Switch to `ggplot()` in examples and `preview_gd()`
- Switch from using base R plotting with `plot_gd()`/`plot_count()` to ggplot plotting with `ggplot_gd()`/`ggplot_count()` in the README, vignette, and` preview_gd()` functions
- This was done to avoid modifying the users environment with `par()`
- `plot_gd()` and `plot_count()` can still be used

### Simplification of vignette and examples
- To decrease code test time for CRAN the vignettes and the examples were reduced, however the core information remains the same

# wingen 2.1.0

### Fixed handling of arguments passed to functions via `...`

There was a bug that resulted in arguments not being passed to`*_general()`/`*_gd()` functions; this has been fixed for `*_general()` and deprecated for `*gd()` (see below)

### Deprecation of `...` for `window_gd()`/`circle_gd()`/`resist_gd()` and addition of `sig` argument

-   The use of `...` has been deprecated for these functions because it wasn't actually passing the `...` arguments to anywhere in the first place, so this should hopefully not affect any users
-   The only argument that may have been passed via `...` is `sig` which is used to set the alpha threshold when `stat = "hwe"`. We have addressed this issue by adding `sig` as an argument to all of the functions.

### Removal of `parallel` and `ncores` arguments

parallelization must now be set up outside of functions, as described in the vignette

### Removal of sp

The sp package is no longer imported by wingen

### Changes to grd argument

The `grd` argument in `krig_gd()` now only accepts SpatRaster and Raster objects

# wingen 2.0.1

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8125610.svg)](https://doi.org/10.5281/zenodo.8125610)

### Deprecation of `parallel` and `ncores` arguments

The `parallel` and `ncores` arguments have been deprecated for all functions.

Instead, `future::plan()` should be used to setup parallelization (see the package vignette for more details).

# wingen 2.0.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8125610.svg)](https://doi.org/10.5281/zenodo.8125610)

### Major updates

-   Circular and resistance-based windows are now possible using the `circle_gd()`/`circle_general()` and `resist_gd()`/`resist_general()`( functions
-   `preview_gd()` takes a new argument (`method`) to specify which method to generate a preview of (`"window"`, `"circle"`, or `"resist"`)
-   Window size across the landscape can now be varied by providing a raster for `maxdist` for `circle_gd()`/`resist_gd()`
-   New population genetic statistic options (specified using `stat`)
-   Create ggplots using `ggplot_gd()` and `ggplot_count()`

### Minor updates

-   All moving window functions can now accept more than one statistic at a time in the form of a vector
-   `plot_gd()` will now generate plots for all layers except the sample count layer and vise versa for `plot_count()`

# wingen 1.1.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7637712.svg)](https://doi.org/10.5281/zenodo.7637712)

### Update to terra and sf

-   The guts of wingen have been updated from raster to terra and now incorporate sf
-   Rasters can now be provided as terra SpatRasters or raster RasterLayers
-   Coordinates can now be provided as sf points, a two-column matrix, or a data.frame representing x and y coordinates.

### Breaking changes

-   wingen is no longer agnostic to projection. Checks have been added that will give warnings for missing CRS and errors for mismatched CRS. The vignette and docs have also been updated to discuss the effects of projection and provide guidance on transformations.
-   The first argument of `window_gd()` is now called `gen` instead of `vcf`. The argument is otherwise the same (i.e., takes the same input), however the name was generalized such that future versions could potentially take different types of genetic data without confusion.
-   The second argument of `mask_gd()` is now called `y` instead of `mask` to avoid conflation with the function `mask()`
-   `mask_gd()` will longer provide the option to resample rasters, rasters must be transformed before the function is used

# wingen 1.0.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7199558.svg)](https://doi.org/10.5281/zenodo.7199558)

Initial release of the wingen R package
