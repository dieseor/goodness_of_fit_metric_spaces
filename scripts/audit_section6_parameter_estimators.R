#!/usr/bin/env Rscript

# Parameter-estimation audit for the Section 6 Normal, logistic-Gaussian and
# vMF models.  This script deliberately does not call any bootstrap routine.
# It checks the production MLE implementations against their defining sample
# equations and compares their repeated-sampling behaviour with exact (Normal
# and LG) or first-order (vMF) theory.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_audit_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[pieces[[1L]]]] <- if (length(pieces) == 1L) "TRUE" else {
      paste(pieces[-1L], collapse = "=")
    }
  }
  out
}

parse_csv_integer <- function(value, default) {
  if (is.null(value) || !nzchar(value)) return(as.integer(default))
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
}

args <- parse_audit_args(commandArgs(trailingOnly = TRUE))
M <- as.integer(args$M %||% 500L)
cores <- as.integer(args$cores %||% 2L)
dimensions <- parse_csv_integer(args$dimensions, c(2L, 5L, 10L))
n_values <- parse_csv_integer(args$n_values, c(50L, 100L, 200L, 400L))
base_seed <- as.integer(args$seed %||% 20260803L)
output_dir <- args$output_dir %||% file.path(
  "simulation_results", "section6_new_scenarios",
  sprintf("parameter_estimator_audit_M%d", M)
)

if (length(M) != 1L || !is.finite(M) || M < 2L) {
  stop("`--M` must be an integer of at least two.")
}
if (length(cores) != 1L || !is.finite(cores) || cores < 1L) {
  stop("`--cores` must be a positive integer.")
}
if (any(!is.finite(dimensions)) || any(dimensions < 2L) ||
    any(!is.finite(n_values)) || any(n_values < 2L)) {
  stop("Dimensions must be at least two and sample sizes must be at least two.")
}

source(file.path("bootstrap", "multiplier_bootstrap.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

e1 <- function(d) c(1, rep.int(0, d - 1L))

section6_sigma <- function(d, sign = c("plus", "minus")) {
  sign <- match.arg(sign)
  sigma <- diag(d)
  rho <- if (identical(sign, "plus")) 0.75 else -0.75
  sigma[1L, 2L] <- rho
  sigma[2L, 1L] <- rho
  sigma
}

normal_vech_index <- function(d) {
  which(lower.tri(matrix(0, d, d), diag = TRUE), arr.ind = TRUE)
}

normal_parameter_names <- function(d, prefix_mu = "mu", prefix_sigma = "Sigma") {
  idx <- normal_vech_index(d)
  c(
    sprintf("%s[%d]", prefix_mu, seq_len(d)),
    sprintf("%s[%d,%d]", prefix_sigma, idx[, 1L], idx[, 2L])
  )
}

normal_theory_covariance <- function(Sigma, n) {
  d <- nrow(Sigma)
  idx <- normal_vech_index(d)
  p_sigma <- nrow(idx)
  out <- matrix(0, d + p_sigma, d + p_sigma)
  out[seq_len(d), seq_len(d)] <- Sigma / n
  for (a in seq_len(p_sigma)) {
    i <- idx[a, 1L]
    j <- idx[a, 2L]
    for (b in seq_len(p_sigma)) {
      k <- idx[b, 1L]
      ell <- idx[b, 2L]
      out[d + a, d + b] <- (n - 1) / n^2 *
        (Sigma[i, k] * Sigma[j, ell] + Sigma[i, ell] * Sigma[j, k])
    }
  }
  out
}

normal_finite_sample_target <- function(mu, Sigma, n) {
  c(mu, fast_multiplier_vech((n - 1) / n * Sigma))
}

normal_exact_sd <- function(Sigma, n) {
  sqrt(diag(normal_theory_covariance(Sigma, n)))
}

safe_condition_number <- function(x) {
  values <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(values)) || min(values) <= 0) return(Inf)
  max(values) / min(values)
}

relative_frobenius_error <- function(x, y) {
  denominator <- sqrt(sum(y^2))
  if (denominator == 0) return(NA_real_)
  sqrt(sum((x - y)^2)) / denominator
}

bind_data_frames_with_missing <- function(rows) {
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  completed <- lapply(rows, function(x) {
    missing <- setdiff(columns, names(x))
    for (name in missing) x[[name]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, completed)
}

summarize_parameter_matrix <- function(estimates, target, theory_sd,
                                       names, family, scenario, d, n,
                                       target_type) {
  empirical_mean <- colMeans(estimates)
  empirical_sd <- apply(estimates, 2L, stats::sd)
  mcse_mean <- empirical_sd / sqrt(nrow(estimates))
  finite_bias <- empirical_mean - target
  z_bias <- finite_bias / mcse_mean
  data.frame(
    family = family,
    scenario = scenario,
    d = as.integer(d),
    n = as.integer(n),
    M = nrow(estimates),
    parameter = names,
    target_type = target_type,
    target = as.numeric(target),
    empirical_mean = as.numeric(empirical_mean),
    empirical_sd = as.numeric(empirical_sd),
    theory_sd = as.numeric(theory_sd),
    variance_ratio = as.numeric((empirical_sd / theory_sd)^2),
    bias_from_target = as.numeric(finite_bias),
    mcse_of_mean = as.numeric(mcse_mean),
    z_bias = as.numeric(z_bias),
    stringsAsFactors = FALSE
  )
}

normal_score_max_abs_mean <- function(x, theta) {
  Sigma_inv <- solve(theta$Sigma)
  centered <- sweep(x, 2L, theta$mu, FUN = "-")
  score_mu <- centered %*% t(Sigma_inv)
  d <- ncol(x)
  score_sigma <- t(vapply(seq_len(nrow(x)), function(i) {
    rr <- centered[i, , drop = FALSE]
    matrix_score <- 0.5 * (Sigma_inv %*% crossprod(rr) %*% Sigma_inv - Sigma_inv)
    fast_multiplier_sym_score_to_vech(matrix_score)
  }, numeric(d * (d + 1L) / 2L)))
  max(abs(colMeans(cbind(score_mu, score_sigma))))
}

fit_normal_once <- function(n, mu, Sigma, seed) {
  set.seed(seed)
  x <- mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
  theta <- fit_mvnormal_theta(
    data = x, null = list(type = "composite"), unknown_param = "both"
  )
  mu_manual <- colMeans(x)
  centered <- sweep(x, 2L, mu_manual, FUN = "-")
  Sigma_manual <- crossprod(centered) / n
  list(
    estimate = c(theta$mu, fast_multiplier_vech(theta$Sigma)),
    max_abs_manual_fit_error = max(abs(c(theta$mu - mu_manual,
                                         theta$Sigma - Sigma_manual))),
    max_abs_mean_score = normal_score_max_abs_mean(x, theta),
    min_eigenvalue = min(eigen(theta$Sigma, symmetric = TRUE, only.values = TRUE)$values),
    covariance_condition_number = safe_condition_number(theta$Sigma)
  )
}

fit_lg_once <- function(n, mu, Sigma, seed) {
  set.seed(seed)
  simplex <- rlogistic_gaussian_simplex(n, mu_ilr = mu, Sigma_ilr = Sigma)
  normalized <- normalize_logistic_gaussian_data(simplex)
  z <- normalized$ilr
  theta <- fit_logistic_gaussian_theta(
    data = simplex, null = list(type = "composite"), unknown_param = "both"
  )
  mu_manual <- colMeans(z)
  centered <- sweep(z, 2L, mu_manual, FUN = "-")
  Sigma_manual <- crossprod(centered) / n
  simplex_back <- logistic_gaussian_ilr_to_simplex(
    z, ambient_dim = ncol(simplex)
  )
  list(
    estimate = c(theta$mu_ilr, fast_multiplier_vech(theta$Sigma_ilr)),
    max_abs_manual_fit_error = max(abs(c(theta$mu_ilr - mu_manual,
                                         theta$Sigma_ilr - Sigma_manual))),
    max_abs_ilr_roundtrip_error = max(abs(simplex - simplex_back)),
    max_abs_mean_score = normal_score_max_abs_mean(
      z, list(mu = theta$mu_ilr, Sigma = theta$Sigma_ilr)
    ),
    min_eigenvalue = min(eigen(theta$Sigma_ilr, symmetric = TRUE, only.values = TRUE)$values),
    covariance_condition_number = safe_condition_number(theta$Sigma_ilr)
  )
}

vmf_information_xi <- function(mu, kappa, q) {
  A <- A_q(kappa, q)
  Aprime <- A_q_prime(kappa, q)
  Aprime * tcrossprod(mu) + (A / kappa) * (diag(q + 1L) - tcrossprod(mu))
}

vmf_aprime_identity <- function(kappa, q) {
  A <- A_q(kappa, q)
  1 - A^2 - q * A / kappa
}

fit_vmf_once <- function(n, mu, kappa, seed) {
  set.seed(seed)
  x <- normalize_vmf_data(rotasym::r_vMF(n, mu = mu, kappa = kappa))
  theta <- fit_vmf_theta(
    data = x, null = list(type = "composite"), unknown_param = "xi"
  )
  resultant <- colMeans(x)
  r_bar <- sqrt(sum(resultant^2))
  theta_direct <- if (r_bar <= 1e-12) {
    normalize_vmf_theta(list(mu = e1(ncol(x)), kappa = 0), ambient_dim = ncol(x))
  } else {
    normalize_vmf_theta(
      list(
        mu = resultant / r_bar,
        kappa = solve_vmf_kappa_from_rbar(r_bar, q = ncol(x) - 1L)
      ),
      ambient_dim = ncol(x)
    )
  }
  q <- ncol(x) - 1L
  scores <- t(vapply(seq_len(nrow(x)), function(i) {
    psi_xi(x[i, ], theta$xi, q)
  }, numeric(ncol(x))))
  A_residual <- abs(A_q(theta$kappa, q) - r_bar)
  list(
    estimate = c(theta$kappa, theta$mu, theta$xi),
    max_abs_direct_fit_error = max(abs(theta$xi - theta_direct$xi)),
    max_abs_mean_score = max(abs(colMeans(scores))),
    A_equation_residual = A_residual,
    kappa_hit_upper_boundary = theta$kappa >= 0.999 * 1e6,
    kappa_zero = theta$kappa == 0
  )
}

run_replicates <- function(M, cores, worker) {
  seeds <- as.integer((base_seed + seq_len(M) * 1009L) %% 2147483646L + 1L)
  if (.Platform$OS.type == "unix" && cores > 1L) {
    return(parallel::mclapply(
      seq_len(M), function(i) worker(seeds[[i]]),
      mc.cores = min(cores, M), mc.set.seed = FALSE
    ))
  }
  lapply(seq_len(M), function(i) worker(seeds[[i]]))
}

normal_cases <- function(d) {
  list(
    normal_1_mixture_null = list(mu = 0.5 * e1(d), Sigma = section6_sigma(d, "plus")),
    normal_2_t3_null = list(mu = rep.int(0, d), Sigma = diag(d))
  )
}

vmf_cases <- function(d) {
  list(
    vmf_1_antipodal_null = list(mu = e1(d + 1L), kappa = as.numeric(d)),
    vmf_2_projected_normal_null = list(mu = e1(d + 1L), kappa = 1.5 * d)
  )
}

parameter_rows <- list()
matrix_rows <- list()
identity_rows <- list()
row_index <- 1L

append_gaussian_audit <- function(family, scenario, d, n, mu, Sigma) {
  replicate_worker <- if (identical(family, "normal")) {
    function(seed) fit_normal_once(n, mu, Sigma, seed)
  } else {
    function(seed) fit_lg_once(n, mu, Sigma, seed)
  }
  fits <- run_replicates(M, cores, replicate_worker)
  estimates <- do.call(rbind, lapply(fits, `[[`, "estimate"))
  target <- normal_finite_sample_target(mu, Sigma, n)
  theory_cov <- normal_theory_covariance(Sigma, n)
  theory_sd <- sqrt(diag(theory_cov))
  prefix_mu <- if (identical(family, "normal")) "mu" else "mu_ilr"
  prefix_sigma <- if (identical(family, "normal")) "Sigma" else "Sigma_ilr"
  parameter_rows[[length(parameter_rows) + 1L]] <<- summarize_parameter_matrix(
    estimates = estimates, target = target, theory_sd = theory_sd,
    names = normal_parameter_names(d, prefix_mu, prefix_sigma), family = family,
    scenario = scenario, d = d, n = n,
    target_type = "exact_finite_sample"
  )
  empirical_cov <- stats::cov(estimates)
  matrix_rows[[length(matrix_rows) + 1L]] <<- data.frame(
    family = family, scenario = scenario, d = d, n = n, M = M,
    parameter_dimension = ncol(estimates),
    covariance_relative_frobenius_error = relative_frobenius_error(empirical_cov, theory_cov),
    covariance_diagonal_median_ratio = stats::median(diag(empirical_cov) / diag(theory_cov)),
    stringsAsFactors = FALSE
  )
  identity_rows[[length(identity_rows) + 1L]] <<- data.frame(
    family = family, scenario = scenario, d = d, n = n, M = M,
    max_abs_manual_fit_error = max(vapply(fits, `[[`, numeric(1), "max_abs_manual_fit_error")),
    max_abs_mean_score = max(vapply(fits, `[[`, numeric(1), "max_abs_mean_score")),
    max_abs_ilr_roundtrip_error = if (identical(family, "lg")) {
      max(vapply(fits, `[[`, numeric(1), "max_abs_ilr_roundtrip_error"))
    } else NA_real_,
    min_fitted_covariance_eigenvalue = min(vapply(fits, `[[`, numeric(1), "min_eigenvalue")),
    max_fitted_covariance_condition_number = max(vapply(fits, `[[`, numeric(1), "covariance_condition_number")),
    stringsAsFactors = FALSE
  )
}

append_vmf_audit <- function(scenario, d, n, mu, kappa) {
  fits <- run_replicates(M, cores, function(seed) fit_vmf_once(n, mu, kappa, seed))
  estimates <- do.call(rbind, lapply(fits, `[[`, "estimate"))
  q <- d
  p <- d + 1L
  information <- vmf_information_xi(mu, kappa, q)
  xi_cov <- solve(information) / n
  A <- A_q(kappa, q)
  Aprime <- A_q_prime(kappa, q)
  theory_sd <- c(
    sqrt(1 / (n * Aprime)),
    c(NA_real_, rep.int(sqrt(1 / (n * kappa * A)), d)),
    sqrt(diag(xi_cov))
  )
  target <- c(kappa, mu, kappa * mu)
  parameter_rows[[length(parameter_rows) + 1L]] <<- summarize_parameter_matrix(
    estimates = estimates, target = target, theory_sd = theory_sd,
    names = c("kappa", sprintf("mu[%d]", seq_len(p)), sprintf("xi[%d]", seq_len(p))),
    family = "vmf", scenario = scenario, d = d, n = n,
    target_type = "asymptotic_for_kappa_and_xi"
  )
  empirical_xi_cov <- stats::cov(estimates[, (p + 2L):(2L * p + 1L), drop = FALSE])
  matrix_rows[[length(matrix_rows) + 1L]] <<- data.frame(
    family = "vmf", scenario = scenario, d = d, n = n, M = M,
    parameter_dimension = p,
    covariance_relative_frobenius_error = relative_frobenius_error(empirical_xi_cov, xi_cov),
    covariance_diagonal_median_ratio = stats::median(diag(empirical_xi_cov) / diag(xi_cov)),
    stringsAsFactors = FALSE
  )
  finite_difference_step <- 1e-5 * max(1, kappa)
  numerical_aprime <- (
    A_q(kappa + finite_difference_step, q) -
      A_q(kappa - finite_difference_step, q)
  ) / (2 * finite_difference_step)
  identity_rows[[length(identity_rows) + 1L]] <<- data.frame(
    family = "vmf", scenario = scenario, d = d, n = n, M = M,
    max_abs_manual_fit_error = max(vapply(fits, `[[`, numeric(1), "max_abs_direct_fit_error")),
    max_abs_mean_score = max(vapply(fits, `[[`, numeric(1), "max_abs_mean_score")),
    max_abs_A_equation_residual = max(vapply(fits, `[[`, numeric(1), "A_equation_residual")),
    A_prime_minus_identity = Aprime - vmf_aprime_identity(kappa, q),
    A_prime_minus_finite_difference = Aprime - numerical_aprime,
    kappa_upper_boundary_hits = sum(vapply(fits, `[[`, logical(1), "kappa_hit_upper_boundary")),
    kappa_zero_hits = sum(vapply(fits, `[[`, logical(1), "kappa_zero")),
    stringsAsFactors = FALSE
  )
}

total_cells <- length(dimensions) * length(n_values) * 6L
completed_cells <- 0L
for (d in dimensions) {
  for (n in n_values) {
    for (scenario in names(normal_cases(d))) {
      case <- normal_cases(d)[[scenario]]
      append_gaussian_audit("normal", scenario, d, n, case$mu, case$Sigma)
      completed_cells <- completed_cells + 1L
      message(sprintf("completed %d/%d cells", completed_cells, total_cells))
    }
    for (scenario in names(normal_cases(d))) {
      case <- normal_cases(d)[[scenario]]
      append_gaussian_audit("lg", sub("^normal", "lg", scenario), d, n, case$mu, case$Sigma)
      completed_cells <- completed_cells + 1L
      message(sprintf("completed %d/%d cells", completed_cells, total_cells))
    }
    for (scenario in names(vmf_cases(d))) {
      case <- vmf_cases(d)[[scenario]]
      append_vmf_audit(scenario, d, n, case$mu, case$kappa)
      completed_cells <- completed_cells + 1L
      message(sprintf("completed %d/%d cells", completed_cells, total_cells))
    }
  }
}

parameter_summary <- bind_data_frames_with_missing(parameter_rows)
matrix_summary <- bind_data_frames_with_missing(matrix_rows)
identity_summary <- bind_data_frames_with_missing(identity_rows)

utils::write.csv(parameter_summary, file.path(output_dir, "parameter_summary.csv"), row.names = FALSE)
utils::write.csv(matrix_summary, file.path(output_dir, "matrix_summary.csv"), row.names = FALSE)
utils::write.csv(identity_summary, file.path(output_dir, "identity_summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(
  M = M, cores = cores, dimensions = paste(dimensions, collapse = ","),
  n_values = paste(n_values, collapse = ","), seed = base_seed,
  stringsAsFactors = FALSE
), file.path(output_dir, "config.csv"), row.names = FALSE)

message("Wrote parameter-estimation audit to: ", normalizePath(output_dir))
