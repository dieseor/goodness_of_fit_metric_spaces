# Compatibility loader for research runners that source the bootstrap engine.
# The installable implementation is in R/bootstrap-*.R.
bootstrap_roots <- c(".", "..", "../..")
bootstrap_root <- bootstrap_roots[file.exists(file.path(bootstrap_roots, "R", "bootstrap-preparation.R"))][1L]
if (is.na(bootstrap_root)) stop("Could not locate the package R/ directory.")
if (!exists("make_normal_spec", mode = "function") ||
    !exists("make_small_circle_weighted_mixture2_spec", mode = "function")) {
  source(file.path(bootstrap_root, "bootstrap", "model_specs.R"), local = environment())
}
for (module in c("preparation", "statistics", "fast-kernels", "execution",
                 "engine", "model-wrappers")) {
  source(file.path(bootstrap_root, "R", paste0("bootstrap-", module, ".R")),
         local = environment())
}
rm(bootstrap_roots, bootstrap_root, module)
