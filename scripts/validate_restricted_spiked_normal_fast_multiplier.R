#!/usr/bin/env Rscript

# Diagnostic suite for the restricted mean-aligned spiked Gaussian adapter.
# It is intentionally opt-in: this file is not sourced by production runners.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  value <- args[startsWith(args, key)]
  if (!length(value)) return(default)
  sub(key, "", value[[1L]], fixed = TRUE)
}

cores <- as.integer(arg_value("cores", "1"))
B <- as.integer(arg_value("B", "99"))
mc_size <- as.integer(arg_value("mc_size", "100000"))
output_dir <- arg_value(
  "output_dir",
  "simulation_results/restricted_spiked_normal_fast_multiplier_validation"
)
seed <- as.integer(arg_value("seed", "20260825"))
if (any(!is.finite(c(cores, B, mc_size, seed))) || cores < 1L || B < 1L || mc_size < 1L) {
  stop("`cores`, `B`, `mc_size`, and `seed` must be positive finite integers.")
}

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

design <- do.call(rbind, lapply(c(2L, 5L), function(d) {
  e1 <- c(1, rep(0, d - 1L))
  diagonal <- rep(1 / sqrt(d), d)
  expand.grid(
    d = d,
    theta_name = c("axis", "diagonal"),
    lambda = c(1e-4, 0.05, 0.5, 2),
    stringsAsFactors = FALSE
  )
}))

theta_from_row <- function(row) {
  d <- as.integer(row$d)
  if (identical(as.character(row$theta_name), "axis")) {
    c(0.8, rep(0, d - 1L))
  } else {
    rep(0.8 / sqrt(d), d)
  }
}

numeric_score_error <- function(x, theta, lambda, step = 1e-6) {
  d <- length(theta)
  parameter <- c(theta, lambda)
  score <- as.numeric(restricted_spiked_normal_score_matrix(
    x, list(theta = theta, lambda = lambda)
  ))
  numerical <- vapply(seq_along(parameter), function(j) {
    left <- parameter; right <- parameter
    left[[j]] <- left[[j]] - step
    right[[j]] <- right[[j]] + step
    (
      restricted_spiked_normal_loglik(x, list(theta = right[seq_len(d)], lambda = right[[d + 1L]])) -
        restricted_spiked_normal_loglik(x, list(theta = left[seq_len(d)], lambda = left[[d + 1L]]))
    ) / (2 * step)
  }, numeric(1L))
  abs(score - numerical)
}

rows <- lapply(seq_len(nrow(design)), function(i) {
  row <- design[i, ]
  theta <- theta_from_row(row)
  lambda <- as.numeric(row$lambda)
  x_score <- rrestricted_spiked_normal(10L, theta, lambda)
  score_error <- unlist(lapply(seq_len(nrow(x_score)), function(j) {
    numeric_score_error(x_score[j, , drop = FALSE], theta, lambda)
  }))
  x_aux <- rrestricted_spiked_normal(mc_size, theta, lambda)
  score_aux <- restricted_spiked_normal_score_matrix(
    x_aux, list(theta = theta, lambda = lambda)
  )
  fisher <- restricted_spiked_normal_fisher_information(
    list(theta = theta, lambda = lambda)
  )
  empirical_fisher <- crossprod(score_aux) / nrow(score_aux)
  Q <- qr.Q(qr(matrix(rnorm(length(theta)^2), length(theta), length(theta))))
  rotated_theta <- drop(t(Q) %*% theta)
  Sigma <- restricted_spiked_normal_covariance(theta, lambda)
  Sigma_rotated <- restricted_spiked_normal_covariance(rotated_theta, lambda)
  Sigma_scaled_theta <- restricted_spiked_normal_covariance(2.75 * theta, lambda)
  sample_covariance <- stats::cov(x_aux) * (nrow(x_aux) - 1L) / nrow(x_aux)
  data.frame(
    d = row$d,
    theta_name = row$theta_name,
    lambda = lambda,
    score_error_max = max(score_error),
    score_error_mean = mean(score_error),
    fisher_mc_fro_error = sqrt(sum((fisher - empirical_fisher)^2)),
    fisher_mc_relative_fro_error = sqrt(sum((fisher - empirical_fisher)^2)) / sqrt(sum(fisher^2)),
    fisher_rcond = rcond(fisher),
    V_rcond = rcond(-fisher),
    fisher_min_eigenvalue = min(eigen(fisher, symmetric = TRUE, only.values = TRUE)$values),
    rotation_sigma_error = sqrt(sum((Sigma_rotated - t(Q) %*% Sigma %*% Q)^2)),
    scale_direction_sigma_error = sqrt(sum((Sigma_scaled_theta - Sigma)^2)),
    sample_covariance_fro_error = sqrt(sum((sample_covariance - Sigma)^2))
  )
})
diagnostics <- do.call(rbind, rows)
utils::write.csv(diagnostics, file.path(output_dir, "score_fisher_invariance.csv"), row.names = FALSE)

# Two deliberately small paired GOF checks.  The fast side exercises the
# sample KS grid and the fused C++ contiguous-double kernel; the reestimated
# side is retained only as a diagnostic comparator.
comparison_rows <- lapply(c(2L, 5L), function(d) {
  theta <- c(0.8, rep(0, d - 1L))
  lambda <- 0.5
  x <- rrestricted_spiked_normal(30L, theta, lambda)
  shared <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    n_cores = cores,
    seed = seed + d,
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE,
                bootstrap_thetas = FALSE),
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 1000L,
      derivative_mc_seed = seed + 100L + d,
      fast_multiplier_cvm_block_size = 50L,
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_multiplier_cache_corrections = "auto"
    ),
    distance_profile_backend = "r"
  )
  fast <- do.call(multiplier_bootstrap_restricted_spiked_normal, c(
    shared, list(bootstrap_method = "fast_multiplier")
  ))
  slow <- do.call(multiplier_bootstrap_restricted_spiked_normal, c(
    shared, list(bootstrap_method = "reestimated")
  ))
  data.frame(
    d = d,
    B = B,
    fast_ks_pvalue = fast$inference$ks$p_value,
    fast_cvm_pvalue = fast$inference$cvm$p_value,
    reestimated_ks_pvalue = slow$inference$ks$p_value,
    reestimated_cvm_pvalue = slow$inference$cvm$p_value,
    fast_backend = fast$diagnostics$fast_multiplier_backend_effective,
    fast_kernel = fast$diagnostics$fast_multiplier_cpp_kernel_effective,
    fast_fused = fast$diagnostics$fast_multiplier_fuse_ks_cvm_effective,
    fast_derivative = fast$diagnostics$derivative_method_effective,
    fast_vhat_method = fast$diagnostics$vhat_method,
    stringsAsFactors = FALSE
  )
})
comparison <- do.call(rbind, comparison_rows)
utils::write.csv(comparison, file.path(output_dir, "fast_vs_reestimated.csv"), row.names = FALSE)

writeLines(c(
  sprintf("seed: %d", seed),
  sprintf("cores: %d", cores),
  sprintf("B: %d", B),
  sprintf("score/Fisher Monte Carlo size: %d", mc_size),
  "fast multiplier: C++ contiguous_double, fused sample KS/CvM, score_mc derivative",
  "paper V=-Fisher; generic engine Vhat=Fisher is the algebraically equivalent internal representation"
), file.path(output_dir, "manifest.txt"))
