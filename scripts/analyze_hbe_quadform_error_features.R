#!/usr/bin/env Rscript

# EN DUDA (2026-07-24): exploratory error-feature analysis only; it does not
# establish a production method-selection rule.
# Relates HBE's profile error to the fourth-cumulant mismatch that remains
# after its first three cumulants have been matched. Diagnostic only.

`%||%` <- function(x, y) if (is.null(x)) y else x

base_dir <- file.path(
  "real_data", "bootstrap_audit", "ks_cvm_real_data_20260723",
  "quadform_real_profile_focus_validation"
)
points_path <- file.path(base_dir, "quadform_exact_validation_points.csv")
results_dir <- file.path(
  "real_data", "logistic_gaussian", "screening", "fast",
  "paper_table_B5000_sampleks_fast_rerun_20260718"
)
if (!file.exists(points_path)) stop("Missing ", points_path)

source(file.path("bootstrap", "model_specs.R"))
points <- utils::read.csv(points_path, stringsAsFactors = FALSE)
paths <- list.files(results_dir, pattern = "_results[.]rds$", full.names = TRUE)

case_cache <- list()
load_case <- function(dataset) {
  if (!is.null(case_cache[[dataset]])) return(case_cache[[dataset]])
  idx <- vapply(paths, function(path) identical(readRDS(path)$dataset_name, dataset), logical(1))
  if (sum(idx) != 1L) stop("Could not find RDS for ", dataset)
  result <- readRDS(paths[[which(idx)]])
  raw <- result$bootstrap$raw_result
  theta <- raw$observed$theta_hat
  z <- normalize_logistic_gaussian_data(result$data_prep$X_closed)$ilr
  out <- list(
    theta = theta,
    z = z,
    lambda = as.numeric(theta$eigenvalues_full[theta$positive_idx]),
    eigenvectors = theta$eigenvectors_full[, theta$positive_idx, drop = FALSE]
  )
  case_cache[[dataset]] <<- out
  out
}

features <- lapply(split(points, seq_len(nrow(points))), function(row) {
  case <- load_case(row$dataset[[1L]])
  i <- as.integer(row$center_index[[1L]])
  shift <- case$theta$mu_ilr - case$z[i, ]
  delta <- as.numeric((as.numeric(shift %*% case$eigenvectors)^2) / case$lambda)
  kappa4 <- 48 * sum(case$lambda^4 * (1 + 4 * delta))
  kappa4_hbe <- 1.5 * row$kappa3[[1L]]^2 / row$kappa2[[1L]]
  data.frame(
    dataset = row$dataset[[1L]],
    center_index = i,
    threshold_index = row$threshold_index[[1L]],
    hbe_abs_error = abs(row$hbe_profile_stored[[1L]] - row$farebrother_cdf[[1L]]),
    hbe_signed_error = row$hbe_profile_stored[[1L]] - row$farebrother_cdf[[1L]],
    hbe_analytic_cutoff = row$hbe_analytic_cutoff[[1L]] && row$q[[1L]] > 1e-12,
    q_standardized_above_q0 = (row$q[[1L]] - row$hbe_q0[[1L]]) / sqrt(row$kappa2[[1L]]),
    kappa4_actual = kappa4,
    kappa4_hbe = kappa4_hbe,
    kappa4_actual_over_hbe = kappa4 / kappa4_hbe,
    stringsAsFactors = FALSE
  )
})
features <- do.call(rbind, features)

summary <- do.call(rbind, lapply(split(features, features$dataset), function(d) {
  noncut <- d[!d$hbe_analytic_cutoff, , drop = FALSE]
  data.frame(
    dataset = d$dataset[[1L]],
    all_cells = nrow(d),
    cutoff_cells = sum(d$hbe_analytic_cutoff),
    hbe_error_max = max(d$hbe_abs_error),
    hbe_error_median = stats::median(d$hbe_abs_error),
    hbe_error_median_cutoff = if (any(d$hbe_analytic_cutoff)) stats::median(d$hbe_abs_error[d$hbe_analytic_cutoff]) else NA_real_,
    hbe_error_median_noncutoff = if (nrow(noncut)) stats::median(noncut$hbe_abs_error) else NA_real_,
    hbe_error_max_noncutoff = if (nrow(noncut)) max(noncut$hbe_abs_error) else NA_real_,
    spearman_abs_error_vs_kappa4_ratio = if (stats::sd(d$kappa4_actual_over_hbe) > 0) stats::cor(d$hbe_abs_error, d$kappa4_actual_over_hbe, method = "spearman") else NA_real_,
    spearman_abs_error_vs_distance_above_cutoff = if (stats::sd(d$q_standardized_above_q0) > 0) stats::cor(d$hbe_abs_error, d$q_standardized_above_q0, method = "spearman") else NA_real_,
    kappa4_ratio_min = min(d$kappa4_actual_over_hbe),
    kappa4_ratio_median = stats::median(d$kappa4_actual_over_hbe),
    kappa4_ratio_max = max(d$kappa4_actual_over_hbe),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(features, file.path(base_dir, "hbe_error_feature_points.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(base_dir, "hbe_error_feature_summary.csv"), row.names = FALSE)
cat("Wrote HBE fourth-cumulant feature diagnostics.\n")
