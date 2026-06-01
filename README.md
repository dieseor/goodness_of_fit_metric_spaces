# Goodness-of-fit in metric spaces

This project explores goodness of fit tests and methods in metric space contexts.

Current supported directional families in the bootstrap pipeline include vMF, HvMF, Jones-Pewsey, spherical cardioid, and spherical Cauchy / Poisson-kernel on S^2.

## Project structure

```
goodness_of_fit_metric_spaces/
├── .Rprofile                 # R startup configuration
├── README.md                 # This file
├── main.R                    # Main analysis script
├── R/                        # R functions and modules
├── data/                     # Data files
├── output/                   # Generated outputs (plots, tables, etc.)
└── tests/                    # Test scripts
```

## Dependencies

The project automatically installs these packages if needed:
- ggplot2 (plotting)
- dplyr (data manipulation)
- readr (data reading)
- here (path management)
- ncdf4 (NetCDF reading for `wind/preprocess_risoe_modern_hvmf.R`)
- GIGrvg (GIG sampler used by `rhvmf_h2_gig()`)

## Reproducible development environment (renv)

This project uses `renv` to capture package dependencies so others can reproduce the environment.

To restore the project environment locally:

```bash
# From project root
# Install renv if needed and restore the environment
Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"
Rscript -e "remotes::install_github('rstudio/renv')"
Rscript -e "renv::restore()"
```

## Spherical Cauchy on S^2

The spherical Cauchy / Poisson-kernel family is implemented on S^2 with parameters `(mu, rho)`, where `mu` is a unit vector in `R^3` and `rho` lies in `[0, 1)`. The population distance profile uses a Legendre expansion as the primary path, with numerical safeguards for high `rho` and final clamping to `[0, 1]`.

The sequential calibration launcher for the default geodesic simple/composite study is:

```bash
Rscript scripts/run_spherical_cauchy_calibration_m500_b500_sequential.R
```
