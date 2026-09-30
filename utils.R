# Compatibility loader for research scripts that call source("utils.R").
# Package users load these functions from dpgof instead.
utils_roots <- c(".", "..", "../..")
utils_root <- utils_roots[file.exists(file.path(utils_roots, "R", "utils-geometry.R"))][1L]
if (is.na(utils_root)) stop("Could not locate the package R/ directory.")
if (!exists("ensure_distance_profile_cpp_loaded", mode = "function")) {
  source(file.path(utils_root, "R", "aaa-backend-state.R"), local = environment())
  source(file.path(utils_root, "R", "distance_profile_backend.R"), local = environment())
}
suppressPackageStartupMessages({
  library(movMF)
  library(sphunif)
  library(pracma)
  library(rotasym)
})
for (module in c("geometry", "numerics", "normal", "directional", "simplex")) {
  source(file.path(utils_root, "R", paste0("utils-", module, ".R")),
         local = environment())
}
source(file.path(utils_root, "R", "aaa-numerical-state.R"), local = environment())
source(file.path(utils_root, "R", "zzz-backend-wrappers.R"), local = environment())
rm(utils_roots, utils_root, module)
cat("Utility functions loaded successfully!\n")
