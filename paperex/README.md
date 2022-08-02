Examples from Bishop et al. 202X
================

## Simulation Example

**Basics:**

The simulation example was executed using the `run_simex.sh` script
which runs the Geonomics simulation and then the time tests to generate
all of the rasters and time data for the different dataset types.

*Note: some of these are parallelized tasks set up to operate on a 32
processor system, this can be modified by going into the scripts
described below*

The results and figures from the paper were created using the
`simex_notebook.Rmd`

**Directory structure:**

    simex
    |   *run_simex.sh* - main script to rerun simulations
    |   *simex_notebook.Rmd* - main notebook to redo paper analyses and generate figures
    |   create_middle_earth.R - used to create middle_earth.csv from DEM layers
    |   gnxsim.py - script to run geonomics simulations
    │   sims_functions.R - functions used for simulation example analysis
    │
    └───data
    │   │   middle_earth.csv - matrix used to create raster layer for simulations
    │   │   mod-sim_params_it-0_t-1000_spp-spp_0.csv - geospatial data from simulations
    |   |   mod-sim_params_it-0_t-1000_spp-spp_0.vcf - genomic data from simulations
    |   |   samples_seed42.csv - samples used in analysis
    │   │
    │   └───middle_earth - directory with original DEM raster layers
    |
    └───time_tests
    │   │   AR_time_test.R - generates results for allelic.richness (WGS and rr dataset)
    │   │   FULL_time_test.R - generates results for FULL dataset
    │   │   rr_time_test.R - generates results for rr dataset
    │   │   WGS_time_test.R - generates results for WGS dataset
    │   
    └───outputs - directory for results from time tests

## Empirical Example

**Basics:**

The results and figures from the paper were created using the
`empex_notebook.Rmd`

**Directory structure:**

    paperex
    |   *empex_notebook.Rmd* - main notebook to redo paper analyses and generate figures
    │
    └───data
        │   populations_r20.haplotypes.filtered_m70_randomSNP.ped - original ped file from Bouzid et al. 2022
        │   populations_r20.haplotypes.filtered_m70_randomSNP.map - original map file from Bouzid et al. 2022
        |   populations_r20.haplotypes.filtered_m70_randomSNP.vcf - vcf file generated from ped and map files
        │   Scelop.coord - sample coordinates from Bouzid et al. 2022
