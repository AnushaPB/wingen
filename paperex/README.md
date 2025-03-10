Examples from Bishop et al. (2023)
================

## Simulation Example

**Basics:**

The simulation example was executed using the `run_simex.sh` script
which runs the Geonomics simulation and then the time tests to generate
all of the rasters and time data for the different dataset types.

The results and figures from the paper were created using the
`simex_notebook.Rmd` and `runtime_notebook.Rmd`

These analyses were run using wingen v1.0.0 which can be downloaded from
Zenodo:

    zenodo_path <- "https://zenodo.org/records/7199558/files/AnushaPB/wingen-v1.0.0.zip"
    install.packages(zenodo_path, repos = NULL, type = "source")

**Directory structure:**

    [simex]
    |   run_simex.sh* - main script to run simulations and generate results
    |   simex_notebook.Rmd* - main notebook for simulation analyses
    |   simex_notebook.md - knitted simex_notebook.Rmd
    │   simex_functions.R* - functions used for simulation example analysis
    |   create_middle_earth.R - used to create middle_earth.csv from DEM layers
    |   gnxsim.py - script to run geonomics simulations
    |   create_datasets - script to create subsampled datasets
    |   runtime_notebook.Rmd - main notebook for runtime analyses and figures
    |   runtime_notebook.md - knitted runtime_notebook.Rmd 
    |   runtime_functions.R - functions used for runtime analysis
    |   simex_figures.Rmd* - notebook used to create plots for figures
    |
    └───[data]
    │   │   middle_earth.csv - matrix used to create raster layer for simulations
    │   │   mod-sim_params_it-0_t-1000_spp-spp_0.csv - geospatial data from simulations
    |   |   mod-sim_params_it-0_t-1000_spp-spp_0.vcf - genomic data from simulations
    |   |   
    │   └───[middle_earth] - directory with original DEM raster layers from Robert 2020 (https://doi.org/10.21220/RKEZ-X707)
    │   │
    │   └───[samples] - directory with samples used in analysis
    |
    └───[simex_tests]
    │   │   AllelicRichness_simex_test.R* - generates results for allelic_richness (WGS and rr dataset)
    │   │   FULL_simex_test.R* - generates results for FULL dataset
    │   │   rr_simex_test.R* - generates results for rr dataset
    │   │   WGS_simex_test.R* - generates results for WGS dataset
    │   
    └───[outputs] - directory for results from time tests

## Empirical Example

**Basics:**

The results and figures from the paper were created using the
`empex_notebook.Rmd`

**Directory structure:**

    [paperex]
    |   empex_notebook.Rmd* - main notebook for empirical analyses and figures
    |   empex_notebook.md - knitted empex_notebook.Rmd
    |   empex_functions.R* - functions used for empirical example analysis
    |   empex_figures.Rmd* - notebook used to create plots for figures
    │
    └───[data]
        │   populations_r20.haplotypes.filtered_m70_randomSNP.ped - original ped file from Bouzid et al. 2022
        │   populations_r20.haplotypes.filtered_m70_randomSNP.map - original map file from Bouzid et al. 2022
        |   populations_r20.haplotypes.filtered_m70_randomSNP.vcf - vcf file generated from ped and map files
        │   Scelop.coord - sample coordinates from Bouzid et al. 2022

## Figures

PDF versions of figures from Bishop et al. (2023) are provided in the
`figures` directory.

*\*Note: files with parallelized tasks are marked with an asterix*
