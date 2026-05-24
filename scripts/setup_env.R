# Setup the project environment using renv
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("rstudio/renv")
# Project dependencies, including `ncdf4` for Risoe wind preprocessing, are restored from `renv.lock`.
renv::restore()
cat("renv restored.\n")
