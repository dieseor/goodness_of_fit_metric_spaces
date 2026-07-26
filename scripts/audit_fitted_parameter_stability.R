#!/usr/bin/env Rscript

# EN DUDA (2026-07-24): exploratory diagnostic script.  It does not establish
# a production bootstrap or quadrature rule, and its conclusions require
# mathematical review and explicit approval before any use in the paper.

# Numerical and leave-one-out diagnostics for the fitted models in the
# selected KS/CvM audit cases.  This does not run any bootstrap simulation.

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else file.path(
  "real_data", "bootstrap_audit", "ks_cvm_real_data_20260723"
)

source("utils.R")
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "logistic_gaussian_model_spec.R"))
source(file.path("bootstrap", "hvmf_model_spec.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

coxite_path <- file.path(
  "real_data", "logistic_gaussian", "screening", "fast",
  "paper_table_B5000_sampleks_fast_rerun_20260718", "coxite_results.rds"
)
wind_run_path <- file.path(output_dir, "risoe_nov_dec_125m_start4_fast_result.rds")
matched_path <- file.path(output_dir, "matched_fast_vs_reestimated_coxite_wind.rds")

stopifnot(file.exists(coxite_path), file.exists(wind_run_path), file.exists(matched_path))

coxite_result <- readRDS(coxite_path)
coxite_x <- as.matrix(coxite_result$data_prep$X_closed)
coxite_theta <- fit_logistic_gaussian_theta(
  coxite_x,
  null = list(type = "composite"),
  unknown_param = "both"
)
coxite_eigen <- eigen(coxite_theta$Sigma_ilr, symmetric = TRUE, only.values = TRUE)$values

coxite_loo <- do.call(rbind, lapply(seq_len(nrow(coxite_x)), function(i) {
  theta_i <- fit_logistic_gaussian_theta(
    coxite_x[-i, , drop = FALSE],
    null = list(type = "composite"),
    unknown_param = "both"
  )
  eig_i <- eigen(theta_i$Sigma_ilr, symmetric = TRUE, only.values = TRUE)$values
  data.frame(
    index = i,
    min_eigenvalue = min(eig_i),
    max_eigenvalue = max(eig_i),
    covariance_condition_number = max(eig_i) / min(eig_i),
    mu_ilr_shift_l2 = sqrt(sum((theta_i$mu_ilr - coxite_theta$mu_ilr)^2)),
    covariance_relative_frobenius_shift = sqrt(sum((theta_i$Sigma_ilr - coxite_theta$Sigma_ilr)^2)) /
      sqrt(sum(coxite_theta$Sigma_ilr^2)),
    stringsAsFactors = FALSE
  )
}))
coxite_loo <- cbind(coxite_loo, as.data.frame(coxite_x))

wind_run <- readRDS(wind_run_path)
wind_data <- utils::read.csv(wind_run$data_path, stringsAsFactors = FALSE)
wind_x <- as.matrix(wind_data[, c("x0", "x1", "x2")])
wind_theta <- hvmf_mle_h2(wind_x)

wind_distance <- function(x, y) hyperbolic_geodesic_distance_h2(x, y)
wind_loo <- do.call(rbind, lapply(seq_len(nrow(wind_x)), function(i) {
  theta_i <- hvmf_mle_h2(wind_x[-i, , drop = FALSE])
  data.frame(
    index = i,
    kappa = theta_i$kappa,
    chi = theta_i$chi,
    theta_deg = theta_i$theta_deg,
    resultant_ratio = theta_i$R / theta_i$W,
    mu_hyperbolic_distance_to_full_fit = wind_distance(theta_i$mu, wind_theta$mu),
    stringsAsFactors = FALSE
  )
}))
wind_loo <- cbind(wind_loo, wind_data[, c("datetime", "speed", "speed_scaled", "direction_deg")])

# The score evaluated at the closed-form MLE should be zero (to floating
# precision).  Its observed outer-product matrix is a direct conditioning
# diagnostic, independent of the auxiliary Monte Carlo derivative sample.
dmu_dchi <- c(
  sinh(wind_theta$chi),
  cosh(wind_theta$chi) * cos(wind_theta$theta),
  cosh(wind_theta$chi) * sin(wind_theta$theta)
)
dmu_dtheta <- c(
  0,
  -sinh(wind_theta$chi) * sin(wind_theta$theta),
  sinh(wind_theta$chi) * cos(wind_theta$theta)
)
mink_mu <- apply(wind_x, 1L, minkowski_inner_product, y = wind_theta$mu)
wind_scores <- cbind(
  chi = wind_theta$kappa * apply(wind_x, 1L, minkowski_inner_product, y = dmu_dchi),
  theta = wind_theta$kappa * apply(wind_x, 1L, minkowski_inner_product, y = dmu_dtheta),
  kappa = 1 / wind_theta$kappa + 1 + mink_mu
)
wind_opg <- crossprod(wind_scores) / nrow(wind_scores)
wind_opg_eigen <- eigen(wind_opg, symmetric = TRUE, only.values = TRUE)$values

matched <- readRDS(matched_path)
wind_fast_diagnostics <- matched$risoe_nov_dec_125m_start4$fast$diagnostics

jackknife_se <- function(x) {
  n <- nrow(x)
  sqrt((n - 1) / n * colSums((x - matrix(colMeans(x), n, ncol(x), byrow = TRUE))^2))
}
wind_jackknife <- jackknife_se(as.matrix(wind_loo[, c("kappa", "chi", "theta_deg", "resultant_ratio")]))

summary <- data.frame(
  case = c("Coxite", "Risoe Nov-Dec 125m"),
  n = c(nrow(coxite_x), nrow(wind_x)),
  parameter_dimension = c(14L, 3L),
  model = c("Logistic Gaussian (4-dimensional ilr)", "HvMF on H^2"),
  primary_condition_number = c(max(coxite_eigen) / min(coxite_eigen), wind_fast_diagnostics$Vhat_condition_number),
  primary_rcond = c(min(coxite_eigen) / max(coxite_eigen), wind_fast_diagnostics$Vhat_rcond),
  stringsAsFactors = FALSE
)

coxite_summary <- data.frame(
  metric = c(
    "minimum_simplex_component", "minimum_covariance_eigenvalue",
    "maximum_covariance_eigenvalue", "covariance_condition_number",
    "smallest_eigenvalue_variance_share", "minimum_loo_covariance_eigenvalue",
    "maximum_loo_covariance_condition_number", "maximum_loo_mu_ilr_shift_l2",
    "maximum_loo_covariance_relative_frobenius_shift"
  ),
  value = c(
    min(coxite_x), min(coxite_eigen), max(coxite_eigen), max(coxite_eigen) / min(coxite_eigen),
    min(coxite_eigen) / sum(coxite_eigen), min(coxite_loo$min_eigenvalue),
    max(coxite_loo$covariance_condition_number), max(coxite_loo$mu_ilr_shift_l2),
    max(coxite_loo$covariance_relative_frobenius_shift)
  )
)
wind_summary <- data.frame(
  metric = c(
    "kappa_hat", "chi_hat", "theta_deg_hat", "resultant_ratio_R_over_W",
    "distance_from_mle_boundary_R_over_W_minus_1", "observed_score_max_abs_mean",
    "observed_score_opg_condition_number", "fast_multiplier_Vhat_condition_number",
    "minimum_loo_chi", "maximum_loo_mu_hyperbolic_distance",
    "jackknife_se_kappa", "jackknife_se_chi", "jackknife_se_theta_deg"
  ),
  value = c(
    wind_theta$kappa, wind_theta$chi, wind_theta$theta_deg, wind_theta$R / wind_theta$W,
    wind_theta$R / wind_theta$W - 1, max(abs(colMeans(wind_scores))),
    max(wind_opg_eigen) / min(wind_opg_eigen), wind_fast_diagnostics$Vhat_condition_number,
    min(wind_loo$chi), max(wind_loo$mu_hyperbolic_distance_to_full_fit),
    wind_jackknife[["kappa"]], wind_jackknife[["chi"]], wind_jackknife[["theta_deg"]]
  )
)

utils::write.csv(summary, file.path(output_dir, "fitted_parameter_stability_summary.csv"), row.names = FALSE)
utils::write.csv(coxite_summary, file.path(output_dir, "coxite_fitted_parameter_diagnostics.csv"), row.names = FALSE)
utils::write.csv(wind_summary, file.path(output_dir, "risoe_hvmf_fitted_parameter_diagnostics.csv"), row.names = FALSE)
utils::write.csv(coxite_loo, file.path(output_dir, "coxite_leave_one_out_fit_diagnostics.csv"), row.names = FALSE)
utils::write.csv(wind_loo, file.path(output_dir, "risoe_hvmf_leave_one_out_fit_diagnostics.csv"), row.names = FALSE)
saveRDS(
  list(
    coxite_theta = coxite_theta,
    coxite_covariance_eigenvalues = coxite_eigen,
    coxite_leave_one_out = coxite_loo,
    wind_theta = wind_theta,
    wind_observed_scores = wind_scores,
    wind_observed_score_opg = wind_opg,
    wind_leave_one_out = wind_loo,
    wind_fast_multiplier_diagnostics = wind_fast_diagnostics
  ),
  file.path(output_dir, "fitted_parameter_stability_audit.rds")
)

message("Wrote fitted-parameter diagnostics to: ", output_dir)
