# Utility Functions for Goodness of Fit in Metric Spaces
# This file contains helper functions for the project

#' Calculate Euclidean distance between two points
#' @param x1 First point (vector)
#' @param x2 Second point (vector)
#' @return Euclidean distance
euclidean_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}

## ---------------------------------------------------------------------------
#' Compute MLE estimator xi = kappa * mu for vMF using movMF
#' @param data matrix of sample rows (n x q) with each row a point on sphere
#' @return numeric vector xi_hat (length q)
compute_mle_xi <- function(data) {
  # Require movMF package for MLE
  if (!requireNamespace('movMF', quietly = TRUE)) {
    stop('movMF package is required to compute vMF MLE. Please install movMF.')
  }
  # movMF expects rows as observations; it returns theta matrix with estimated parameters
  fit <- movMF::movMF(data, k = 1)
  xi_hat <- as.vector(fit$theta[1, ])
  return(xi_hat)
}


#' Calculate pairwise distances for a dataset
#' @param data Matrix or data frame with observations in rows
#' @return Distance matrix
pairwise_distances <- function(data) {
  n <- nrow(data)
  dist_matrix <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        dist <- euclidean_distance(data[i, ], data[j, ])
        dist_matrix[i, j] <- dist
        dist_matrix[j, i] <- dist
      }
    }
  }
  
  return(dist_matrix)
}

#' Generate random points in a unit ball
#' @param n Number of points to generate
#' @param dim Dimension of the space
#' @return Matrix with points as rows
generate_ball_points <- function(n, dim = 2) {
  # Generate points using the standard method for uniform distribution in a ball
  points <- matrix(rnorm(n * dim), nrow = n, ncol = dim)
  
  # Normalize to unit sphere
  norms <- sqrt(rowSums(points^2))
  points <- points / norms
  
  # Scale by random radii to get uniform distribution in ball
  radii <- runif(n)^(1/dim)
  points <- points * radii
  
  return(points)
}

#' Kolmogorov-Smirnov test for distances
#' @param distances Vector of observed distances
#' @param reference_cdf Function representing the reference CDF
#' @return KS test result
ks_test_distances <- function(distances, reference_cdf) {
  ks.test(distances, reference_cdf)
}

#' Safe KS test wrapper that returns a result with NA p.value for invalid inputs
#' @param x Numeric vector
#' @param y Numeric vector
#' @return List(list(statistic, p.value, success, message)) where success=TRUE if test ran
safe_ks_test <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) {
    warning('safe_ks_test: inputs are not numeric; returning NA p.value')
    return(list(statistic = NA_real_, p.value = NA_real_, success = FALSE, message = 'non-numeric inputs'))
  }
  if (any(is.na(x)) || any(is.na(y))) {
    warning('safe_ks_test: inputs contain NA; returning NA p.value')
    return(list(statistic = NA_real_, p.value = NA_real_, success = FALSE, message = 'missing values'))
  }
  res <- tryCatch({
    ks <- ks.test(x, y)
    list(statistic = ks$statistic, p.value = ks$p.value, success = TRUE)
  }, error = function(e) {
    warning(sprintf('safe_ks_test: ks.test failed: %s', e$message))
    list(statistic = NA_real_, p.value = NA_real_, success = FALSE, message = e$message)
  })
  return(res)
}

## ---------------------------------------------------------------------------
## ENTRYWISE (vMF) helpers moved from gaussian_process_vmf.R for testing only
## ---------------------------------------------------------------------------
## Note: These functions are provided for testing and diagnostics only. They are
## not used in the vectorized workflow by default.
##
## compute_conditional_expectation_vmf
compute_conditional_expectation_vmf <- function(omega, t, mu, kappa,
                                               distance_type, mc_samples) {
  if (t <= 0) return(rep(0, length(mu)))
  dots <- mc_samples %*% omega
  if (distance_type == "chordal") {
    dists <- sqrt(2 * (1 - dots))
  } else {
    dots <- check_dot_products(dots)
    dists <- acos(dots)
  }
  in_ball <- dists <= t
  if (sum(in_ball) == 0) return(rep(0, length(mu)))
  filtered_samples <- mc_samples[in_ball, , drop = FALSE]
  return(colMeans(filtered_samples))
}

## compute_joint_probability_vmf
compute_joint_probability_vmf <- function(omega1, t1, omega2, t2, mu, kappa,
                                         distance_type, mc_samples) {
  if (t1 <= 0 || t2 <= 0) return(0)
  if (sum((omega1 - omega2)^2) > (t1 + t2)^2) {
    # simple check for chordal distances; use sphere_distance for geodesic
    # Note: we convert to exact computation for accurate results when needed
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

## Note: safe_ks_test removed — use ks.test directly. The codebase now prefers to error on invalid inputs

# ============================================================================
# GAUSSIAN PROCESS GENERIC FUNCTIONS
# Functions that work for ANY distribution (Normal, vMF, etc.)
# ============================================================================

#' Simulate the limiting Gaussian process
#' Works for any distribution once covariance matrix is computed. 
#' If the covariance matrix is invalid (e.g., not numerically positive semi-definite) sampling may fail. There are no internal automatic covariance corrections performed.
#' @param cov_matrix Pre-computed covariance matrix
#' @param M Number of Monte Carlo simulations
#' @return Vector of M supremum values from the Gaussian process
simulate_limit_gaussian <- function(cov_matrix, M = 10000, seed = NULL, tol = 1e-10) {
  n_total <- nrow(cov_matrix)

  # No automatic PSD correction. We will attempt to sample and warn on failure.
  cat("Generating", M, "multivariate normal samples from", n_total, "dimensional process...\n")
  if (!is.null(seed)) set.seed(seed)
  # Sample directly without any automatic PSD correction.
  # If mvtnorm::rmvnorm errors due to a non-PD sigma, we will stop() and raise an error.
  gaussian_samples <- tryCatch({
    mvtnorm::rmvnorm(M, mean = rep(0, n_total), sigma = cov_matrix)
  }, error = function(e) {
    warning(sprintf("simulate_limit_gaussian: Sampling failed due to non-PD covariance or other error: %s. No PSD correction applied; returning NA vector.", e$message))
    return(rep(NA_real_, M))
  })
  
  # Compute supremum for each sample
  supremum_values <- apply(gaussian_samples, 1, function(row) max(abs(row)))
  # No PSD correction applied; return supremum values as-is. No correction attributes are set.
  
  # cat("Supremum statistics, limiting Gaussian process:\n")
  # cat("  Mean:", round(mean(supremum_values), 4), "\n")
  # cat("  Median:", round(median(supremum_values), 4), "\n")
  # cat("  Max:", round(max(supremum_values), 4), "\n")
  # cat("  Min:", round(min(supremum_values), 4), "\n\n")
  
  return(supremum_values)
}

##' Simulate the limiting Gaussian process using eigen-decomposition (alternate sampler)
##' This avoids potential differences due to mvtnorm internal choices and uses a clear eigen-based sampling
## simulate_limit_gaussian_eig removed — use simulate_limit_gaussian which is the canonical sampler.
## The eigen-based alternate sampler used to be here, but it was unused in the repo
## and caused maintenance burden. If you need an alternate eigen-based sampler,
## we can add it back in a scoped utility function and wire it into the test scripts.


## NOTE: Automatic covariance correction tools were intentionally removed.
## If you need to apply a covariance correction, do it externally and pass the corrected matrix to simulate_limit_gaussian().

# ---------------------------------------------------------------------------
# Normal distribution helper functions
# Move these here so they can be reused across scripts and exported to
# parallel workers rather than redefining in many files.
# ---------------------------------------------------------------------------

#' Compute theoretical distance profile for univariate normal distribution
#' Vectorized helper for N(mu, sigma^2)
theoretical_distance_profile_normal <- function(omega, mu, sigma, t_values) {
  n <- max(length(omega), length(t_values))
  omega <- rep(omega, length.out = n)
  t_values <- rep(t_values, length.out = n)
  result <- numeric(n)
  valid_t <- t_values > 0
  if (any(valid_t)) {
    upper_bound <- (omega[valid_t] + t_values[valid_t] - mu) / sigma
    lower_bound <- (omega[valid_t] - t_values[valid_t] - mu) / sigma
    result[valid_t] <- pnorm(upper_bound) - pnorm(lower_bound)
  }
  return(result)
}

#' Compute joint probability for normal distribution
#' P(|X-omega1|<=t1, |X-omega2|<=t2) for X ~ N(mu, sigma^2)
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

## ---------------------------------------------------------------------------
#' Derivatives of distance-profile for normal distribution
#' These are used in constructing the covariance under composite nulls.
#' d/dmu F(omega,t; mu, sigma)
dotF_mu_normal <- function(omega, mu, sigma, t) {
  -(1/sigma) * dnorm((omega - mu + t)/sigma) +
    (1/sigma) * dnorm((omega - mu - t)/sigma)
}

#' d/dsigma F(omega,t; mu, sigma)
dotF_sigma_normal <- function(omega, mu, sigma, t) {
  -(omega - mu + t)/sigma^2 * dnorm((omega - mu + t)/sigma) +
    (omega - mu - t)/sigma^2 * dnorm((omega - mu - t)/sigma)
}


## Diagnostic helpers removed by user request. If you need to inspect eigenvalues, use eigen() directly.

#' Setup parallel cluster for computation (GENERIC)
#' @param n_cores Number of cores to use
#' @param export_vars Character vector of variable names to export
#' @param export_funs List of functions to export to workers
#' @param envir Environment from which to export variables
#' @return Cluster object
setup_parallel_cluster <- function(n_cores, export_vars = NULL, export_funs = NULL, envir = parent.frame()) {
  library(parallel)
  cl <- makeCluster(n_cores)
  
  # Set seed for reproducibility
  clusterSetRNGStream(cl, iseed = 42)
  
  # Export functions if provided
  if (!is.null(export_funs)) {
    clusterExport(cl, names(export_funs), envir = list2env(export_funs))
  }
  
  # Export variables if provided
  if (!is.null(export_vars)) {
    clusterExport(cl, export_vars, envir = envir)
  }
  
  return(cl)
}

#' Distribute work in round-robin fashion (GENERIC)
#' @param n_total Total number of tasks
#' @param n_cores Number of cores
#' @return List of vectors, each containing task indices for one core
distribute_work_round_robin <- function(n_total, n_cores) {
  chunks <- vector("list", n_cores)
  for (i in 1:n_cores) {
    chunks[[i]] <- seq(from = i, to = n_total, by = n_cores)
  }
  return(chunks)
}

# ============================================================================
# DISTANCE PROFILE FOR vMF DISTRIBUTION
# ============================================================================

#' Theoretical distance profile for von Mises-Fisher distribution
#' 
#' This is the ORIGINAL implementation from vmf_distance_profile_analysis.R
#' that was tested and verified to work correctly.
#' 
#' @param omega Reference point on the sphere S^q
#' @param mu Mean direction of vMF distribution (unit vector in R^{q+1})
#' @param kappa Concentration parameter (κ > 0)
#' @param t_values Vector of distance thresholds
#' @param distance_type Either "chordal" or "geodesic"
#' @return P(d(X, omega) <= t) for each t in t_values, where X ~ vMF(μ, κ)
theoretical_distance_profile_vmf <- function(omega, mu, kappa, t_values, distance_type = "chordal") {
  # If omega is a matrix, vectorize over rows
  if (is.matrix(omega)) {
    n <- nrow(omega)
    # If t_values is scalar, repeat
    if (length(t_values) == 1) t_values <- rep(t_values, n)
    stopifnot(length(t_values) == n)
    return(sapply(1:n, function(i) {
      theoretical_distance_profile_vmf(omega[i, ], mu, kappa, t_values[i], distance_type)
    }))
  }
  # Original code for single omega
  rho <- sum(mu * omega)  # μ'ω (cosine of angle between μ and ω)
  q <- length(mu)  # Ambient dimension (for S^{q-1} embedded in R^q)
  sapply(t_values, function(t) {
    threshold <- ifelse(distance_type == "chordal", 1 - (t^2)/2, cos(t))
    density_T <- function(s) {
      log_c_q <- rotasym::c_vMF(p = q, kappa = kappa, log = TRUE)
      log_exp_term <- kappa * rho * s
      log_power_term <- ((q - 3)/2) * log(1 - s^2)
      log_numerator <- log_c_q + log_exp_term + log_power_term
      kappa_term <- kappa * sqrt((1 - s^2) * (1 - rho^2))
      log_denominator <- rotasym::c_vMF(p = q - 1, kappa = kappa_term, log = TRUE)
      result <- exp(log_numerator - log_denominator)
      return(result)
    }
    cdf_at_threshold <- integrate(density_T, lower = -1 + 1e-8, upper = threshold, 
                                   rel.tol = 1e-8, abs.tol = 1e-10)$value
    return(1 - cdf_at_threshold)
  })
}

cat("Utility functions loaded successfully!\n")