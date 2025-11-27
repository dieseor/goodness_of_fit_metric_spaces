# Setup the project environment using renv
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("rstudio/renv")
renv::restore()
cat("renv restored.\n")
