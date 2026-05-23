#!/usr/bin/env Rscript
# Test helpers moved from `utils.R` to live alongside tests so they are clearly
# scoped for testing and diagnostics. These functions were previously used only
# for debugging or entrywise comparisons and are not part of the core library.

source(file.path("utils.R"))

## compute_joint_probability_normal
compute_joint_probability_normal <- function(omega1, t1, omega2, t2, mu, sigma) {
  if (t1 <= 0 || t2 <= 0) return(0)
  lower1 <- omega1 - t1
  upper1 <- omega1 + t1
  lower2 <- omega2 - t2
  upper2 <- omega2 + t2
  intersection_lower <- max(lower1, lower2)
  intersection_upper <- min(upper1, upper2)
  if (intersection_lower >= intersection_upper) return(0)
  return(pnorm((intersection_upper - mu) / sigma) - pnorm((intersection_lower - mu) / sigma))
}


## compute_distance_matrix_vmf
compute_distance_matrix_vmf <- function(sample_matrix, omega_grid, distance_type = "chordal") {
  dot_products <- sample_matrix %*% t(omega_grid)
  dot_products <- check_dot_products(dot_products)
  if (distance_type == "chordal") {
    return(sqrt(2 * (1 - dot_products)))
  } else if (distance_type == "geodesic") {
    dot_products <- check_dot_products(dot_products)
    return(acos(dot_products))
  } else {
    stop("distance_type must be 'chordal' or 'geodesic'")
  }
}


## compute_joint_probability_vmf
compute_joint_probability_vmf <- function(omega1, t1, omega2, t2, mu, kappa,
                                         distance_type, mc_samples) {
  if (t1 <= 0 || t2 <= 0) return(0)
  if (sum((omega1 - omega2)^2) > (t1 + t2)^2) {
    # simple check for chordal distances; use sphere_distance for geodesic
    if (distance_type == 'chordal') {
      if (sqrt(sum((omega1 - omega2)^2)) > t1 + t2) return(0)
    } else {
      if (sphere_distance(omega1, omega2, distance_type) > t1 + t2) return(0)
    }
  }
  dot1 <- mc_samples %*% omega1
  dot2 <- mc_samples %*% omega2
  if (distance_type == 'chordal') {
    dist1 <- sqrt(2 * (1 - dot1)); dist2 <- sqrt(2 * (1 - dot2))
  } else {
    dot1 <- check_dot_products(dot1); dot2 <- check_dot_products(dot2)
    dist1 <- acos(dot1); dist2 <- acos(dot2)
  }
  joint_indicator <- (dist1 <= t1) & (dist2 <= t2)
  return(mean(joint_indicator))
}


## compute_covariance_entry_vmf
compute_covariance_entry_vmf <- function(omega1, t1, omega2, t2, mu, kappa,
                                        distance_type, mc_samples, A_q_kappa, var_X, h0 = c('simple','composite'), unknown_param = NULL) {
  h0 <- match.arg(h0)
  q <- length(mu) - 1
  F1 <- theoretical_distance_profile_vmf(omega1, mu, kappa, t1, distance_type)
  F2 <- theoretical_distance_profile_vmf(omega2, mu, kappa, t2, distance_type)
  P_joint <- compute_joint_probability_vmf(omega1, t1, omega2, t2, mu, kappa, distance_type, mc_samples)
  E1 <- compute_conditional_expectation_vmf(omega1, t1, mu, kappa, distance_type, mc_samples)
  E2 <- compute_conditional_expectation_vmf(omega2, t2, mu, kappa, distance_type, mc_samples)
  m1 <- E1 - A_q_kappa * mu
  m2 <- E2 - A_q_kappa * mu
  inv_var_X <- tryCatch(solve(var_X), error = function(e) MASS::ginv(var_X))
  quadratic_form <- as.numeric(t(m1) %*% inv_var_X %*% m2)
  C_base <- P_joint - F1 * F2
  if (h0 == 'simple') return(C_base)
  composite_corr <- -F1 * F2 * quadratic_form
  cov_value <- C_base + composite_corr
  return(cov_value)
}


## compare_processes_normal (test helper)
compare_processes_normal <- function(omega_grid, t_grid, mu, sigma, M = 1000, n_cores = 2) {
  limit_values <- simulate_limit_gaussian_normal(omega_grid, t_grid, mu, sigma, M = M, n_cores = n_cores)
  empirical_values <- simulate_empirical_process_normal(omega_grid, t_grid, n = 100, mu, sigma, M = M, n_cores = n_cores)
  return(list(limit_values = limit_values, empirical_values = empirical_values))
}


## compare_processes_vmf (test helper)
compare_processes_vmf <- function(omega_grid, t_grid, mu, kappa, M = 1000, n_mc_samples = 1000, n_cores = 2) {
  limit_values <- simulate_limit_gaussian_vmf(omega_grid, t_grid, mu, kappa, M = M, n_mc_samples = n_mc_samples, n_cores = n_cores)
  empirical_values <- simulate_empirical_process_vmf(omega_grid, t_grid, n = 100, mu, kappa, M = M, n_cores = n_cores)
  return(list(limit_values = limit_values, empirical_values = empirical_values))
}


## Simulate empirical process using loops (baseline slow version)
simulate_empirical_process_loops <- function(omega_grid, t_grid, n, mu, sigma, M = 1000) {
  supremum_values <- numeric(M)
  for (i in 1:M) {
    # Generate sample
    sample_x <- rnorm(n, mean = mu, sd = sigma)
    max_difference <- 0
    for (omega in omega_grid) {
      for (t in t_grid) {
        empirical_distances <- abs(sample_x - omega)
        F_hat_omega_t <- mean(empirical_distances <= t)
        F_omega_t <- theoretical_distance_profile_normal(omega, mu, sigma, t)
        scaled_difference <- sqrt(n) * abs(F_hat_omega_t - F_omega_t)
        max_difference <- max(max_difference, scaled_difference)
      }
    }
    supremum_values[i] <- max_difference
  }
  return(supremum_values)
}


## Covariance function used by baseline nested loops (one element at a time)
covariance_gaussian_process <- function(omega1, t1, omega2, t2, mu, sigma) {
  if (t1 <= 0 || t2 <= 0) return(0)
  joint_prob <- compute_joint_probability_normal(omega1, t1, omega2, t2, mu, sigma)
  f_omega1_t1 <- theoretical_distance_profile_normal(omega1, mu, sigma, t1)
  f_omega2_t2 <- theoretical_distance_profile_normal(omega2, mu, sigma, t2)
  return(joint_prob - f_omega1_t1 * f_omega2_t2)
}


## Create covariance matrix using nested loops (baseline slow method)
create_covariance_matrix_loops <- function(omega_grid, t_grid, mu, sigma) {
  n_omega <- length(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  grid_combinations <- expand.grid(omega = omega_grid, t = t_grid)
  cov_matrix <- matrix(0, n_total, n_total)
  total_iterations <- n_total * n_total
  progress_interval <- max(1, floor(total_iterations / 5))
  iteration_count <- 0
  for (i in 1:n_total) {
    for (j in 1:n_total) {
      iteration_count <- iteration_count + 1
      if (iteration_count %% progress_interval == 0) {
        progress_pct <- round(100 * iteration_count / total_iterations, 1)
        cat("  Progress:", progress_pct, "%\n")
      }
      cov_matrix[i, j] <- covariance_gaussian_process(
        grid_combinations$omega[i],
        grid_combinations$t[i],
        grid_combinations$omega[j],
        grid_combinations$t[j], mu, sigma
      )
    }
  }
  cat("Covariance matrix created successfully\n")
  return(cov_matrix)
}


## compute_joint_probability_vmf_no_opt
compute_joint_probability_vmf_no_opt <- function(omega1, t1, omega2, t2, mu, kappa, distance_type, mc_samples) {
  if (t1 <= 0 || t2 <= 0) return(0)
  dot1 <- mc_samples %*% omega1
  dot2 <- mc_samples %*% omega2
  if (distance_type == "chordal") {
    dist1 <- sqrt(2 * (1 - dot1)); dist2 <- sqrt(2 * (1 - dot2))
  } else {
    dist1 <- acos(pmax(-1, pmin(1, dot1)))
    dist2 <- acos(pmax(-1, pmin(1, dot2)))
  }
  joint_indicator <- (dist1 <= t1) & (dist2 <= t2)
  return(mean(joint_indicator))
}


## QQ plot helpers
## Compare empirical process supremums vs limit (Gaussian) quantiles for Normal
qqplot_empirical_vs_limit_normal <- function(omega_grid, t_grid, mu, sigma, n_values = c(10,50,100), M = 500, n_cores = 1, seed = NULL, save_plot = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2 to generate QQ plots")
  if (!is.null(seed)) set.seed(seed)
  # Simulate limit process
  limit_supremums <- simulate_limit_gaussian_normal(omega_grid, t_grid, mu, sigma, M = M, n_cores = n_cores)
  # Prepare plotting dataframe
  df_all <- data.frame()
  probs <- ppoints(M)
  limit_qs <- as.numeric(quantile(limit_supremums, probs = probs, type = 8))
  for (n in n_values) {
    empirical_supremums <- simulate_empirical_process_normal(omega_grid, t_grid, n = n, mu = mu, sigma = sigma, M = M, n_cores = n_cores)
    empirical_qs <- as.numeric(quantile(empirical_supremums, probs = probs, type = 8))
    df <- data.frame(
      sample_size = as.factor(n),
      p = probs,
      theoretical = limit_qs,
      empirical = empirical_qs
    )
    df_all <- rbind(df_all, df)
  }
  # Plot
  library(ggplot2)
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = theoretical, y = empirical, color = sample_size)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(x = "Limit quantiles", y = "Empirical quantiles", color = "n") +
    ggplot2::theme_minimal()
  if (!is.null(save_plot)) ggplot2::ggsave(save_plot, plot = p)
  return(p)
}


## QQ plot helper for vMF
qqplot_empirical_vs_limit_vmf <- function(omega_grid, t_grid, mu, kappa, n_values = c(10,50,100), M = 500, n_mc_samples = 1000, n_cores = 1, seed = NULL, save_plot = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2 to generate QQ plots")
  if (!is.null(seed)) set.seed(seed)
  # Simulate limit process for vMF
  limit_supremums <- simulate_limit_gaussian_vmf(omega_grid, t_grid, mu, kappa, M = M, n_mc_samples = n_mc_samples, n_cores = n_cores)
  probs <- ppoints(M)
  limit_qs <- as.numeric(quantile(limit_supremums, probs = probs, type = 8))
  df_all <- data.frame()
  for (n in n_values) {
    empirical_supremums <- simulate_empirical_process_vmf(omega_grid, t_grid, n = n, mu = mu, kappa = kappa, M = M, n_cores = n_cores)
    empirical_qs <- as.numeric(quantile(empirical_supremums, probs = probs, type = 8))
    df <- data.frame(
      sample_size = as.factor(n),
      p = probs,
      theoretical = limit_qs,
      empirical = empirical_qs
    )
    df_all <- rbind(df_all, df)
  }
  library(ggplot2)
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = theoretical, y = empirical, color = sample_size)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(x = "Limit quantiles", y = "Empirical quantiles", color = "n") +
    ggplot2::theme_minimal()
  if (!is.null(save_plot)) ggplot2::ggsave(save_plot, plot = p)
  return(p)
}


# Basic smoke tests for moved helpers (run when executing this test file directly)
if (identical(environment(), globalenv())) {
  cat("Running basic smoke checks for test helpers...\n")
  # 1) compute_joint_probability_normal basic check
  p <- compute_joint_probability_normal(0, 0.5, 0.5, 0.5, mu = 0, sigma = 1)
  cat("compute_joint_probability_normal ->", p, "\n")

  # 2) compute_distance_matrix_vmf basic check
  sample_matrix <- matrix(c(0, 0, 1, 1, 0, 0), nrow = 2, byrow = TRUE)
  grid <- matrix(c(0, 0, 1), nrow = 1)
  dmat <- compute_distance_matrix_vmf(sample_matrix, grid, distance_type = "chordal")
  cat("compute_distance_matrix_vmf dims ->", paste(dim(dmat), collapse = ","), "\n")

  # 3) compute_covariance_entry_vmf smoke check (very small MC sample)
  mc_samples <- rotasym::r_vMF(20, mu = c(0,0,1), kappa = 1)
  mu <- c(0, 0, 1)
  kappa <- 1
  q <- length(mu) - 1
  A_q_kappa <- besselI(kappa, nu = (q + 1) / 2, expon.scaled = TRUE) / besselI(kappa, nu = (q - 1) / 2, expon.scaled = TRUE)
  scalar_coef <- 1 - A_q_kappa^2 - ((q + 1) * A_q_kappa / kappa)
  var_X <- (A_q_kappa / kappa) * diag(q + 1) + scalar_coef * outer(mu, mu)
  val <- compute_covariance_entry_vmf(c(0,0,1), 0.5, c(0,0,1), 0.5, mu, kappa, "chordal", mc_samples, A_q_kappa, var_X)
  cat("compute_covariance_entry_vmf value ->", val, "\n")
  cat("Test helpers OK.\n")
}

# Unit test to check compute_covariance_entry_vmf against vectorized cov_vmf for a small grid
if (identical(environment(), globalenv())) {
  cat("Running sanity comparison of compute_covariance_entry_vmf vs cov_vmf...\n")
  omega_grid <- generate_canonical_lattice(2, dim = 3)
  t_grid <- c(0.3, 0.6)
  mu <- c(0, 0, 1)
  kappa <- 2
  mc_samples <- rotasym::r_vMF(500, mu = mu, kappa = kappa)
  q <- length(mu) - 1
  A_q_kappa <- besselI(kappa, nu = (q + 1) / 2, expon.scaled = TRUE) / besselI(kappa, nu = (q - 1) / 2, expon.scaled = TRUE)
  scalar_coef <- 1 - A_q_kappa^2 - ((q + 1) * A_q_kappa / kappa)
  var_X <- (A_q_kappa / kappa) * diag(q + 1) + scalar_coef * outer(mu, mu)
  source(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))
  cov_mat <- cov_vmf(omega_grid, t_grid, mu, kappa, distance_type = "chordal", n_mc_samples = 500, n_cores = 1, mc_samples = mc_samples, seed = 123, upper_triangle = FALSE)
  # pick entry (1,2)
  entry_idx <- 2
  omega1_idx <- ((entry_idx - 1) %% nrow(omega_grid)) + 1
  t1_idx <- ((entry_idx - 1) %/% nrow(omega_grid)) + 1
  omega2_idx <- ((1 - 1) %% nrow(omega_grid)) + 1
  t2_idx <- ((1 - 1) %/% nrow(omega_grid)) + 1
  omega1 <- omega_grid[omega1_idx, ]; t1 <- t_grid[t1_idx]
  omega2 <- omega_grid[omega2_idx, ]; t2 <- t_grid[t2_idx]
  manual_val <- compute_covariance_entry_vmf(omega1, t1, omega2, t2, mu, kappa, "chordal", mc_samples, A_q_kappa, var_X, h0 = 'simple')
  cat("cov_vmf[1,2] ->", cov_mat[1,2], "manual ->", manual_val, "\n")
}
# Tests for utility functions in R/utils.R

cat("=== Running utils.R helper tests ===\n")

# Source utilities
source("utils.R")

# ------------------------------------------------------------------
# Test 1: generate_canonical_lattice
# ------------------------------------------------------------------
cat("Test 1: generate_canonical_lattice\n")
mat <- generate_canonical_lattice(20, dim = 3)
stopifnot(is.matrix(mat))
stopifnot(nrow(mat) == 20 && ncol(mat) == 3)
# Check unit length of rows
norms <- sqrt(rowSums(mat^2))
if (all(abs(norms - 1) < 1e-8)) {
  cat("  ✓ Rows are unit vectors\n")
} else {
  cat("  ✗ FAIL: Some rows are not unit vectors (\n")
  print(head(norms))
  stop("generate_canonical_lattice produced non-unit vectors")
}

# ------------------------------------------------------------------
# Test 2: check_dot_products
# ------------------------------------------------------------------
cat("Test 2: check_dot_products\n")
# Valid inputs
v <- c(-1, -0.5, 0, 0.5, 1)
r <- check_dot_products(v)
stopifnot(all.equal(v, r))
# Slightly out-of-range: should error when exceeding tolerance
bad_v <- c(1.1)
err <- tryCatch({ check_dot_products(bad_v); FALSE }, error = function(e) TRUE)
if (isTRUE(err)) {
  cat("  ✓ check_dot_products correctly errors on invalid dot products\n")
} else {
  cat("  ✗ FAIL: check_dot_products did not error on invalid input\n")
  stop("check_dot_products should error on dot products outside [-1 - tol, 1 + tol]")
}

# ------------------------------------------------------------------
# Test 3: compute_precomp_vmf (structural checks only)
# ------------------------------------------------------------------
cat("Test 3: compute_precomp_vmf (structural checks)\n")
# Skip if rotasym not installed — function uses theoretical_distance_profile_vmf which calls rotasym::c_vMF
if (!requireNamespace("rotasym", quietly = TRUE)) {
  cat("  ⚠ Skipping compute_precomp_vmf tests because 'rotasym' not installed.\n")
} else {
  set.seed(123)
  mu <- c(0, 0, 1)
  kappa <- 1
  n_mc <- 30
  # Generate random unit vectors for mc_samples
  mc_samples <- matrix(rnorm(n_mc * 3), nrow = n_mc, ncol = 3)
  mc_samples <- t(apply(mc_samples, 1, function(v) v / sqrt(sum(v^2))))

  omega_grid <- generate_canonical_lattice(5, dim = 3) # 5 omegas
  t_grid <- c(0.1, 0.5)

  # A_q_kappa simple computation
  q <- length(mu) - 1
  A_q_kappa <- besselI(kappa, nu = (q + 1) / 2, expon.scaled = TRUE) / 
               besselI(kappa, nu = (q - 1) / 2, expon.scaled = TRUE)

  res <- compute_precomp_vmf(mc_samples, omega_grid, t_grid, mu, kappa, "chordal", A_q_kappa)
  stopifnot(is.list(res))
  stopifnot("dists_all" %in% names(res))
  stopifnot("F2_matrix" %in% names(res))
  stopifnot("E2_array" %in% names(res))
  stopifnot("E2_mat" %in% names(res))
  stopifnot("F2_vec" %in% names(res))
  stopifnot("m2_mat" %in% names(res))
  stopifnot("in_ball_list" %in% names(res))
  # dists_all should be n_mc x n_omega
  stopifnot(nrow(res$dists_all) == n_mc)
  stopifnot(ncol(res$dists_all) == nrow(omega_grid))
  # F2_matrix should be n_omega x n_t
  stopifnot(nrow(res$F2_matrix) == nrow(omega_grid))
  stopifnot(ncol(res$F2_matrix) == length(t_grid))
  # E2_array should be n_omega x n_t x (q+1)
  stopifnot(dim(res$E2_array)[1] == nrow(omega_grid))
  stopifnot(dim(res$E2_array)[2] == length(t_grid))
  stopifnot(dim(res$E2_array)[3] == (length(mu)))
  cat("  ✓ compute_precomp_vmf returned expected structure\n")
}

cat("=== utils.R helper tests finished ===\n")
