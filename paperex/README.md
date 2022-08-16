Examples from Bishop et al. 202X
================

## Simulation Example

**Basics:**

The simulation example was executed using the `run_simex.sh` script
which runs the Geonomics simulation and then the time tests to generate
all of the rasters and computational time results for the different dataset types.

*Note: some of these are parallelized tasks set up to operate on a 32
processor system, this can be modified by going into any of the files
marked with an asterix described below*

The results and figures from the paper were created using the
`simex_notebook.Rmd` and `runtime_notebook.Rmd`

**Directory structure:**

    [simex]
    |   run_simex.sh* - main script to run simulations and generate results
    |   creat_datasets.R - script to create example datasets
    |   simex_notebook.Rmd* - main notebook for simulation analyses and figures
    |   simex_notebook.md - knitted simex_notebook.Rmd
    │   simex_functions.R* - functions used for simulation example analysis
    |   create_middle_earth.R - used to create middle_earth.csv from DEM layers
    |   gnxsim.py - script to run geonomics simulations
    |   runtime_notebook.Rmd - main notebook for runtime analyses and figures
    |   runtime_notebook.md - knitted runtime_notebook.Rmd 
    |   runtime_functions.R - functions used for runtime analysis
    |
    └───[data]
    │   │   middle_earth.csv - matrix used to create raster layer for simulations
    │   │   mod-sim_params_it-0_t-1000_spp-spp_0.csv - geospatial data from simulations
    |   |   mod-sim_params_it-0_t-1000_spp-spp_0.vcf - genomic data from simulations
    │   │
    │   └───[middle_earth] - directory with original DEM raster layers
    │   │
    │   └───[samples] - directory with sample files used to create example datasets
    |
    └───[time_tests]
    │   │   AllelicRichness_time_test.R* - generates results for allelic.richness (WGS and rr dataset)
    │   │   FULL_time_test.R* - generates results for FULL dataset
    │   │   rr_time_test.R* - generates results for rr dataset
    │   │   WGS_time_test.R* - generates results for WGS dataset
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
    │
    └───[data]
        │   populations_r20.haplotypes.filtered_m70_randomSNP.ped - original ped file from Bouzid et al. 2022
        │   populations_r20.haplotypes.filtered_m70_randomSNP.map - original map file from Bouzid et al. 2022
        |   populations_r20.haplotypes.filtered_m70_randomSNP.vcf - vcf file generated from ped and map files
        │   Scelop.coord - sample coordinates from Bouzid et al. 2022
