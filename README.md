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

## Getting Started

1. Open this folder in VS Code
2. Open an R script or create a new one (`.R` extension)
3. Run R code interactively using Ctrl+Enter (Cmd+Enter on Mac)

## Key Features

- **R Integration**: Full R language support with syntax highlighting and IntelliSense
- **Interactive Execution**: Run R code line by line or in blocks
- **Project Structure**: Organized directories for code, data, and outputs
- **Automatic Setup**: .Rprofile loads common packages automatically

## Usage

Create R scripts in the `R/` directory for your functions and analyses. Use `main.R` for your primary workflow.

## Dependencies

The project automatically installs these packages if needed:
- ggplot2 (plotting)
- dplyr (data manipulation)
- readr (data reading)
- here (path management)

Add additional packages to the `.Rprofile` file as needed.

## Reproducible development environment (renv)

This project uses `renv` to capture package dependencies so others can reproduce the environment.

To restore the project environment locally (recommended):

```bash
# From project root
# Install renv if needed and restore the environment
Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"
Rscript -e "remotes::install_github('rstudio/renv')"
Rscript -e "renv::restore()"
```

After restoring, you can run the test scripts under `tests/` by sourcing them:

```bash
Rscript -e "source('tests/test_vectorization_improvement.R')"
```

Tip: To create or update the lockfile after adding packages on your local machine, run:

```bash
Rscript -e "renv::snapshot(prompt = FALSE)"
git add renv.lock
git commit -m 'Update renv lockfile'
git push origin main
```

Notes:
- Some packages (e.g., `rgl`) may require system libraries or X11/OpenGL support when installed. On macOS this is usually available; on Linux (CI) you may need to install system packages. In CI, the workflow is configured to restore `renv` and source the test scripts but may skip tests that require interactive graphics.
- If you encounter package installation problems, try manually installing missing package system prerequisites then re-run `renv::restore()`.