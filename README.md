# Goodness of Fit in Metric Spaces

This project explores goodness of fit tests and methods in metric space contexts.

## Project Structure

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