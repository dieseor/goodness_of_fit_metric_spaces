# Loader for research scripts and tests that call source("utils.R") instead of
# installing gofmetric. It sources every file in R/ in package collation order;
# the compiled kernels are then built lazily from src/ with Rcpp::sourceCpp().
if (!exists("multiplier_bootstrap_gof", mode = "function")) {
  utils_roots <- c(".", "..", "../..")
  utils_root <- utils_roots[file.exists(file.path(utils_roots, "R", "model_specs.R"))][1L]
  if (is.na(utils_root)) stop("Could not locate the package R/ directory.")
  suppressPackageStartupMessages(
    for (pkg in c("parallel", "movMF", "sphunif", "pracma", "rotasym",
                  "mvtnorm", "gbutils", "ggplot2")) {
      library(pkg, character.only = TRUE)
    }
  )
  r_files <- sort(list.files(file.path(utils_root, "R"), "\\.R$"), method = "radix")
  for (r_file in setdiff(r_files, "RcppExports.R")) {
    source(file.path(utils_root, "R", r_file), local = environment())
  }
  source(file.path(utils_root, "scripts", "research_only_helpers.R"),
         local = environment())
  rm(utils_roots, utils_root, pkg, r_files, r_file)
  cat("Utility functions loaded successfully!\n")
}
