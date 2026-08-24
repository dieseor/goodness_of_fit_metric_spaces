#!/usr/bin/env Rscript

# Extensive validation for the custom restricted-spiked Gaussian MLE.  This
# script deliberately does not run any GOF simulation.  It compares the
# profiled MLE with an independently specified OpenMx fit and records all
# diagnostics needed to audit the comparison.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  match <- args[startsWith(args, prefix)]
  if (!length(match)) return(default)
  sub(prefix, "", match[[1L]], fixed = TRUE)
}

output_dir <- arg_value("output_dir", "simulation_results/restricted_spiked_normal_mle_validation")
reps <- as.integer(arg_value("reps", "3"))
seed <- as.integer(arg_value("seed", "20260824"))
if (!is.finite(reps) || reps < 1L || !is.finite(seed)) {
  stop("`--reps` must be positive and `--seed` must be finite.")
}

source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "restricted_spiked_normal_openmx.R"))
if (!restricted_spiked_normal_openmx_available()) {
  stop("OpenMx is not installed; install it in this project's renv before running this validation.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

theta_design <- function(d) {
  e1 <- c(1, rep(0, d - 1L))
  diagonal <- rep(1 / sqrt(d), d)
  list(
    short_axis = 0.35 * e1,
    moderate_axis = 0.9 * e1,
    long_diagonal = 1.6 * diagonal
  )
}

dimensions <- c(2L, 5L)
n_values <- c(50L, 200L, 800L)
lambda_values <- c(0, 0.01, 0.5, 2, 8)
design <- do.call(rbind, lapply(dimensions, function(d) {
  theta_list <- theta_design(d)
  do.call(rbind, lapply(names(theta_list), function(theta_name) {
    theta <- theta_list[[theta_name]]
    expand.grid(
      d = d, n = n_values, lambda = lambda_values, replication = seq_len(reps),
      theta_name = theta_name, stringsAsFactors = FALSE
    )
  }))
}))
design$case_id <- seq_len(nrow(design))

run_case <- function(row) {
  d <- as.integer(row[["d"]])
  theta_true <- theta_design(d)[[row[["theta_name"]]]]
  lambda_true <- as.numeric(row[["lambda"]])
  x <- rrestricted_spiked_normal(as.integer(row[["n"]]), theta_true, lambda_true)
  custom <- fit_restricted_spiked_normal_theta(x, null = list(type = "composite"))

  # Deliberately separated reference starts: custom MLE, a fixed-scale
  # coordinate perturbation, and a rotated direction.  None is zero.
  random_direction <- rnorm(d)
  random_direction <- random_direction / sqrt(sum(random_direction^2))
  starts <- list(
    list(theta = custom$theta, lambda = custom$lambda),
    list(theta = theta_true * 1.4, lambda = max(0.05, lambda_true * 0.4)),
    list(theta = sqrt(sum(theta_true^2)) * random_direction, lambda = max(0.05, lambda_true * 2))
  )
  reference <- fit_restricted_spiked_normal_openmx_multistart(x, starts = starts)
  best <- reference$best
  c(
    as.list(row),
    list(
      theta_error = sqrt(sum((custom$theta - theta_true)^2)),
      lambda_error = custom$lambda - lambda_true,
      sigma_fro_error = sqrt(sum((custom$Sigma - restricted_spiked_normal_covariance(theta_true, lambda_true))^2)),
      custom_openmx_theta_error = sqrt(sum((custom$theta - best$theta)^2)),
      custom_openmx_lambda_error = custom$lambda - best$lambda,
      custom_openmx_sigma_error = sqrt(sum((custom$Sigma - best$Sigma)^2)),
      custom_openmx_loglik_error = custom$loglik * nrow(x) - best$loglik,
      openmx_manual_loglik_error = best$loglik - best$loglik_openmx,
      custom_lambda = custom$lambda,
      custom_profile_tau = custom$fit_diagnostics$profile_tau,
      custom_profile_tau_upper = custom$fit_diagnostics$profile_tau_upper,
      custom_profile_candidates = length(custom$fit_diagnostics$profile_local_candidate_indices),
      custom_profile_radius = custom$fit_diagnostics$profiled_radius,
      custom_radius_tolerance = custom$fit_diagnostics$radius_tolerance,
      openmx_all_converged = reference$all_converged,
      openmx_loglik_spread = diff(reference$loglik_range),
      openmx_successful_starts = length(reference$successful),
      openmx_status_code = best$fit_status_code
    )
  )
}

started <- proc.time()[["elapsed"]]
rows <- vector("list", nrow(design))
for (i in seq_len(nrow(design))) {
  rows[[i]] <- run_case(design[i, , drop = FALSE])
  if (i %% 10L == 0L || i == nrow(design)) {
    elapsed <- proc.time()[["elapsed"]] - started
    message(sprintf("restricted-spiked MLE validation: %d/%d cases, elapsed %.1fs", i, nrow(design), elapsed))
    utils::write.csv(do.call(rbind.data.frame, rows[seq_len(i)]),
                     file.path(output_dir, "case_results_partial.csv"), row.names = FALSE)
  }
}
results <- do.call(rbind.data.frame, rows)
utils::write.csv(results, file.path(output_dir, "case_results.csv"), row.names = FALSE)

summarize_group <- function(group) {
  data.frame(
    d = group$d[[1L]], n = group$n[[1L]], lambda = group$lambda[[1L]],
    theta_name = group$theta_name[[1L]],
    theta_bias = mean(group$theta_error), theta_rmse = sqrt(mean(group$theta_error^2)),
    lambda_bias = mean(group$lambda_error), lambda_rmse = sqrt(mean(group$lambda_error^2)),
    sigma_rmse = sqrt(mean(group$sigma_fro_error^2)),
    max_custom_openmx_theta_error = max(group$custom_openmx_theta_error),
    max_custom_openmx_lambda_error = max(abs(group$custom_openmx_lambda_error)),
    max_custom_openmx_sigma_error = max(group$custom_openmx_sigma_error),
    max_custom_openmx_loglik_error = max(abs(group$custom_openmx_loglik_error)),
    max_openmx_manual_loglik_error = max(abs(group$openmx_manual_loglik_error)),
    all_openmx_starts_converged = all(group$openmx_all_converged),
    max_openmx_start_loglik_spread = max(group$openmx_loglik_spread)
  )
}
groups <- split(results, interaction(results$d, results$n, results$lambda, results$theta_name, drop = TRUE))
summary_table <- do.call(rbind, lapply(groups, summarize_group))
utils::write.csv(summary_table, file.path(output_dir, "summary_by_design.csv"), row.names = FALSE)

# Independent large-sample moment, rotation, and profile-grid checks.
sanity_rows <- list()
for (d in dimensions) {
  theta <- theta_design(d)$long_diagonal
  lambda <- 2
  target_sigma <- restricted_spiked_normal_covariance(theta, lambda)
  x_large <- rrestricted_spiked_normal(100000L, theta, lambda)
  Q <- qr.Q(qr(matrix(rnorm(d^2), d, d)))
  fit <- fit_restricted_spiked_normal_theta(x_large[seq_len(1000L), , drop = FALSE], null = list(type = "composite"))
  fit_rotated <- fit_restricted_spiked_normal_theta(
    x_large[seq_len(1000L), , drop = FALSE] %*% Q, null = list(type = "composite")
  )
  xbar <- colMeans(x_large[seq_len(1000L), , drop = FALSE])
  S <- crossprod(sweep(x_large[seq_len(1000L), , drop = FALSE], 2L, xbar, FUN = "-")) / 1000L
  fine_profile <- restricted_spiked_normal_profile_maximize(xbar, S, grid_size = 20001L)
  sanity_rows[[as.character(d)]] <- data.frame(
    d = d,
    mean_error = sqrt(sum((colMeans(x_large) - theta)^2)),
    covariance_fro_error = sqrt(sum((stats::cov(x_large) * 99999 / 100000 - target_sigma)^2)),
    rotation_theta_error = sqrt(sum((fit_rotated$theta - drop(t(Q) %*% fit$theta))^2)),
    rotation_lambda_error = abs(fit_rotated$lambda - fit$lambda),
    rotation_sigma_error = sqrt(sum((fit_rotated$Sigma - t(Q) %*% fit$Sigma %*% Q)^2)),
    fitted_profile_vs_fine_grid = abs(fit$fit_diagnostics$profile_value - fine_profile$best$value),
    fitted_lambda_vs_fine_grid = abs(fit$lambda - fine_profile$best$lambda)
  )
}
sanity <- do.call(rbind, sanity_rows)
utils::write.csv(sanity, file.path(output_dir, "sanity_checks.csv"), row.names = FALSE)

manifest <- c(
  sprintf("seed: %d", seed),
  sprintf("reps: %d", reps),
  sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
  sprintf("n_values: %s", paste(n_values, collapse = ",")),
  sprintf("lambda_values: %s", paste(lambda_values, collapse = ",")),
  "theta configurations: short_axis=0.35e1; moderate_axis=0.9e1; long_diagonal=1.6(1,...,1)/sqrt(d)",
  "custom MLE: profiled one-dimensional tau=log(1+lambda), grid plus local refinement",
  "OpenMx: three independent starts per sample"
)
writeLines(manifest, file.path(output_dir, "manifest.txt"))
message(sprintf("Validation completed in %.1fs. Results: %s", proc.time()[["elapsed"]] - started, output_dir))
