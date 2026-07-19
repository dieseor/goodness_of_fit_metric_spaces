# Goodness-of-fit in metric spaces

This project explores goodness of fit tests and methods in metric space contexts.

Current supported directional families in the bootstrap pipeline include vMF, HvMF, Jones-Pewsey, spherical cardioid, and spherical Cauchy / Poisson-kernel on S^2.

## Project structure

```
goodness_of_fit_metric_spaces/
├── bootstrap/                # Model specs and bootstrap engines
├── distance_profiles/        # Distance-profile analyses and diagnostics
├── output/
│   ├── bahadur/              # Bahadur analyses by distribution
│   ├── calibration/          # Bootstrap calibration outputs
│   ├── catalog/              # Inventory and migration reports
│   ├── convergence/          # Empirical process convergence outputs
│   ├── diagnostics/          # Model diagnostics and quality checks
│   ├── distance_profiles/    # Distance-profile outputs by distribution
│   └── real_data/            # Outputs from real datasets (including comets)
├── scripts/                  # One-off runners and batch entry points
├── tests/                    # Test scripts
├── utils.R                   # Shared utilities and model helpers
└── wind/                     # Wind-related pipelines and outputs
```

Output aliases and the intermediate `output/structured/` layer have been removed.
The project now uses a single direct output tree under `output/`.

## Dependencies

The project automatically installs these packages if needed:
- ggplot2 (plotting)
- dplyr (data manipulation)
- readr (data reading)
- here (path management)
- ncdf4 (NetCDF reading for `wind/preprocess_risoe_modern_hvmf.R`)

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
