#!/usr/bin/env Rscript

# Add the missing beta = 0 designs to an existing final second-scenario run.
# Existing design identifiers are preserved, so the raw results remain safely
# resumable and the manifest records the null designs in the same campaign.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/run_second_scenarios_power.R")
source("scripts/run_hvmf_alt1_polar_power.R")

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
scenario <- cli$scenario %||% stop("Supply --scenario=normal|lg|vmf|hvmf_alt2|hvmf_alt1.")

specifications <- list(
  normal = list(output_dir = file.path(power_root, "final_M1000_B5000", "normal"), family = "normal", candidate = 3),
  lg = list(output_dir = file.path(power_root, "final_M1000_B5000", "lg"), family = "lg", candidate = 3),
  vmf = list(output_dir = file.path(power_root, "final_M1000_B5000", "vmf"), family = "vmf", candidate = 0.5),
  hvmf_alt2 = list(output_dir = file.path(power_root, "final_M1000_B5000", "hvmf"), family = "hvmf", candidate = list(kappa = 50, delta = 0.2)),
  hvmf_alt1 = list(output_dir = file.path(power_root, "final_hvmf_alt1_polar_kappa15_M1000_B5000"), family = "hvmf_alt1_polar", candidate = 1.5)
)
spec <- specifications[[scenario]] %||% stop("Unknown scenario.")
manifest_path <- file.path(spec$output_dir, "manifest.csv")
if (!file.exists(manifest_path)) stop("Final manifest not found: ", manifest_path)
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
if (any(manifest$M != 1000L) || any(manifest$B != 5000L)) stop("This launcher is only for M=1000, B=5000 final campaigns.")

if (identical(scenario, "hvmf_alt1")) {
  desired <- make_hvmf_alt1_design(spec$candidate, c(50L, 100L, 200L), 0)
} else {
  desired <- make_design(setNames(list(list(spec$candidate)), spec$family), c(50L, 100L, 200L), 0)
}

existing_null <- manifest[abs(manifest$beta) < 1e-12, , drop = FALSE]
if (nrow(existing_null)) {
  desired$design_id <- existing_null$design_id[match(desired$n, existing_null$n)]
  if (anyNA(desired$design_id)) stop("Existing beta=0 manifest rows do not cover n=50,100,200.")
} else {
  desired$design_id <- max(manifest$design_id) + seq_len(nrow(desired))
  extra <- manifest[rep(1L, nrow(desired)), , drop = FALSE]
  for (field in intersect(names(extra), names(desired))) extra[[field]] <- desired[[field]]
  write_atomic_csv(rbind(manifest, extra), manifest_path)
}

if (identical(scenario, "hvmf_alt1")) {
  results <- run_hvmf_alt1_jobs(desired, M = 1000L, B = 5000L, cores = 5L, output_dir = spec$output_dir,
                                stage = "final_alt1_polar_kappa15", base_seed = 20260716L)
  utils::write.csv(summarize_hvmf_alt1(results), file.path(spec$output_dir, "summary.csv"), row.names = FALSE)
} else {
  results <- run_power_jobs(desired, M = 1000L, B = 5000L, cores = 5L, output_dir = spec$output_dir,
                            stage = paste0("final_", spec$family), bootstrap_method = "fast_multiplier", hvmf_small_grid = FALSE)
  utils::write.csv(summarize_power(results), file.path(spec$output_dir, "summary.csv"), row.names = FALSE)
}
