# Utility Functions for Goodness of Fit in Metric Spaces
# This file contains helper functions for the project

# ============================================================================
# DISTANCE CALCULATION
# ============================================================================

#' Calculate Euclidean distance between two points
#' @param x1 First point (vector)
#' @param x2 Second point (vector)
#' @return Euclidean distance
euclidean_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}

#' Compute distances between two points on sphere (helper)
#' @param omega1 point on sphere
#' @param omega2 point on sphere
#' @param distance_type either 'chordal' or 'geodesic'
#' @return scalar distance
sphere_distance <- function(omega1, omega2, distance_type = "chordal") {
  omega1 <- omega1 / sqrt(sum(omega1^2))
  omega2 <- omega2 / sqrt(sum(omega2^2))
  dot_product <- sum(omega1 * omega2)
  dot_product <- pmax(pmin(dot_product, 1), -1)
  if (distance_type == "chordal") {
    return(sqrt(2 * (1 - dot_product)))
  } else if (distance_type == "geodesic") {
    return(acos(dot_product))
  } else {
    stop("distance_type must be 'chordal' or 'geodesic'")
  }
}

#' Validate dot products are within [-1, 1]
#' Ensures numerical stability before calling acos; useful for spherical geometry computations
#' @param dot_products Numeric scalar or vector/matrix of dot products
#' @param tolerance Numeric tolerance for out-of-range checks (default 1e-10)
#' @return The input dot_products unchanged if valid; throws error if invalid values are present
check_dot_products <- function(dot_products, tolerance = 1e-10) {
  if (any(dot_products < -1 - tolerance) || any(dot_products > 1 + tolerance)) {
    bad_values <- dot_products[dot_products < -1 - tolerance | dot_products > 1 + tolerance]
    stop(paste0("Invalid dot products detected: range [", 
                round(min(bad_values), 10), ", ", round(max(bad_values), 10), 
                "]. Check that inputs are unit vectors."))
  }
  return(dot_products)
}

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
  
  # cat("Supremum statistics, limiting Gaussian process:\n")
  # cat("  Mean:", round(mean(supremum_values), 4), "\n")
  # cat("  Median:", round(median(supremum_values), 4), "\n")
  # cat("  Max:", round(max(supremum_values), 4), "\n")
  # cat("  Min:", round(min(supremum_values), 4), "\n\n")
  
  return(supremum_values)
}

# ============================================================================
# SPECIFIC FUNCTIONS FOR NORMAL DISTRIBUTION
# ============================================================================

# Compute theoretical distance profile for univariate normal distribution - vectorized
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


# ============================================================================
# SPECIFIC FUNCTIONS FOR vMF DISTRIBUTION
# ============================================================================


# ----------------------------------------------------------------------------
#' Generate canonical lattice on the unit sphere (Fibonacci sphere / golden spiral)
#' @param n Number of points to generate
#' @param dim Dimension of the ambient space (default 3)
#' @return Matrix of size (n, dim) with unit vectors on the sphere
generate_canonical_lattice <- function(n, dim = 3) {
  if (dim == 3) {
    golden_ratio <- (1 + sqrt(5)) / 2
    i <- seq(0, n - 1)
    theta <- 2 * pi * i / golden_ratio
    phi <- acos(1 - 2 * (i + 0.5) / n)
    x <- cos(theta) * sin(phi)
    y <- sin(theta) * sin(phi)
    z <- cos(phi)
    mat <- cbind(x, y, z)
    return(mat)
  } else if (dim > 1) {
    mat <- matrix(rnorm(n * dim), nrow = n, ncol = dim)
    mat <- t(apply(mat, 1, function(v) v / sqrt(sum(v^2))))
    return(mat)
  } else {
    stop("Dimension must be >= 2")
  }
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

## ---------------------------------------------------------------------------
#' Theoretical distance profile for von Mises-Fisher distribution

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

## ---------------------------------------------------------------------------
## Conditional Expectation for vMF distribution
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


#' Precompute common intermediates for vMF vectorized covariance
#' @param mc_samples Monte Carlo sample (n_mc × q)
#' @param omega_grid Matrix (n_omega × q)
#' @param t_grid Vector (n_t)
#' @param mu Mean direction (q+1)
#' @param kappa Concentration parameter
#' @param distance_type Either 'chordal' or 'geodesic'
#' @param A_q_kappa Precomputed A_q(kappa)
#' @return A named list with dists_all, F2_matrix, E2_array, E2_mat, F2_vec, m2_mat, in_ball_list
compute_precomp_vmf <- function(mc_samples, omega_grid, t_grid, mu, kappa, distance_type, A_q_kappa) {
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  q <- length(mu) - 1
  n_mc <- nrow(mc_samples)

  # Dot products and distances from all MC samples to all omegas
  dots_all <- mc_samples %*% t(omega_grid)
  if (distance_type != "chordal") dots_all <- check_dot_products(dots_all)
  dists_all <- if (distance_type == "chordal") sqrt(2 * (1 - dots_all)) else acos(dots_all)

  # F2_matrix: marginal probabilities for each omega2 and t
  F2_matrix <- t(apply(omega_grid, 1, function(omega2) {
    theoretical_distance_profile_vmf(omega2, mu, kappa, t_grid, distance_type)
  }))

  # Build E2_mat (n_total x q+1) using a simpler and explicit rbind-by-t approach
  in_ball_list <- vector("list", n_t)
  E2_rows_by_t <- lapply(seq_len(n_t), function(k) {
    in_ball_k <- dists_all <= t_grid[k]
    in_ball_list[[k]] <<- in_ball_k
    counts_k <- colSums(in_ball_k)
    if (all(counts_k == 0)) return(matrix(0, nrow = n_omega, ncol = q + 1))
    sums_mat <- t(mc_samples) %*% (in_ball_k * 1) # (q+1) x n_omega
    counts_k_safe <- counts_k
    counts_k_safe[counts_k_safe == 0] <- 1
    means_mat <- sweep(sums_mat, 2, counts_k_safe, "/")
    if (any(counts_k == 0)) means_mat[, counts_k == 0] <- 0
    return(t(means_mat))
  })
  E2_mat <- do.call(rbind, E2_rows_by_t)
  E2_array <- array(NA, dim = c(n_omega, n_t, q + 1))
  for (k in seq_len(n_t)) {
    E2_array[, k, ] <- E2_rows_by_t[[k]]
  }

  F2_vec <- as.vector(F2_matrix)
  m2_mat <- sweep(E2_mat, 2, A_q_kappa * mu, "-")

  return(list(
    dists_all = dists_all,
    F2_matrix = F2_matrix,
    E2_array = E2_array,
    E2_mat = E2_mat,
    F2_vec = F2_vec,
    m2_mat = m2_mat,
    in_ball_list = in_ball_list
  ))
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

# ============================================================================
# FUNCTIONS FOR BAHADUR REPRESENTATION
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(movMF)
  library(sphunif)
  library(pracma)
  library(rotasym)
})

#' Calculate A_q function for vMF distribution
A_q <- function(kappa, q) {
  if (kappa == 0) return(0)
  nu1 <- (q + 1) / 2
  nu2 <- (q - 1) / 2
  I_nu1 <- besselI(kappa, nu = nu1, expon.scaled = TRUE)
  I_nu2 <- besselI(kappa, nu = nu2, expon.scaled = TRUE)
  return(I_nu1 / I_nu2)
}

#' Calculate derivative of A_q function
A_q_prime <- function(kappa, q) {
  if (kappa == 0) return(0)
  nu1 <- (q + 1) / 2
  nu2 <- (q - 1) / 2
  nu3 <- (q - 3) / 2
  I_nu1 <- besselI(kappa, nu = nu1, expon.scaled = TRUE)
  I_nu2 <- besselI(kappa, nu = nu2, expon.scaled = TRUE)
  I_nu3 <- besselI(kappa, nu = nu3, expon.scaled = TRUE)
  I_nu1_prime <- I_nu2 - (nu1 / kappa) * I_nu1
  I_nu2_prime <- I_nu3 - (nu2 / kappa) * I_nu2
  numerator <- I_nu1_prime * I_nu2 - I_nu1 * I_nu2_prime
  denominator <- I_nu2^2
  return(numerator / denominator)
}

#' Score function psi_xi
psi_xi <- function(x, xi, q) {
  kappa <- norm(xi, type = "2")
  if (kappa == 0) return(x)
  mu <- xi / kappa
  A_q_val <- A_q(kappa, q)
  return(x - A_q_val * mu)
}

#' Derivative of score function
dot_psi_xi <- function(xi, q) {
  kappa <- norm(xi, type = "2")
  if (kappa == 0) {
    return(-diag(q + 1))
  }
  mu <- xi / kappa
  A_q_val <- A_q(kappa, q)
  A_q_prime_val <- A_q_prime(kappa, q)
  I_q <- diag(q + 1)
  mu_outer <- outer(mu, mu)
  return(-(A_q_prime_val * mu_outer + (A_q_val / kappa) * (I_q - mu_outer)))
}

## NOTE: `compute_mle_xi` was moved to `R/utils.R` to centralize utilities.
## Do not re-define it here — scripts should source `R/utils.R` to import it.

#' Single trajectory analysis
analyze_single_trajectory <- function(mu_true, kappa_true, sample_sizes, trajectory_id) {
  q <- length(mu_true) - 1
  xi_true <- kappa_true * mu_true
  max_n <- max(sample_sizes)
  full_sample <- r_vMF(max_n, mu_true, kappa_true)
  
  results <- data.frame()
  
  for (n in sample_sizes) {
    current_sample <- full_sample[1:n, , drop = FALSE]
    xi_hat <- compute_mle_xi(current_sample)
    quantity_1 <- sqrt(n) * (xi_hat - xi_true)
    
    score_sum <- rep(0, q + 1)
    for (i in 1:n) {
      score_sum <- score_sum + psi_xi(current_sample[i, ], xi_true, q)
    }
    
    dot_psi_matrix <- dot_psi_xi(xi_true, q)
    
    if (det(dot_psi_matrix) != 0) {
      dot_psi_inv <- solve(dot_psi_matrix)
      quantity_2 <- -dot_psi_inv %*% (score_sum / sqrt(n))
    } else {
      stop(sprintf("Singular matrix for trajectory %d at n = %d", trajectory_id, n))
    }
    
    difference_norm <- norm(quantity_1 - as.vector(quantity_2), type = "2")
    
    # Print progress messages at key sample sizes
    if (n %in% c(10000, 25000, 50000, 75000)) {
      cat("  >> Trajectory", trajectory_id, "reached n =", n, "- Difference norm:", 
          round(difference_norm, 6), "\n")
    }
    
    results <- rbind(results, data.frame(
      trajectory = trajectory_id,
      n = n,
      difference_norm = difference_norm
    ))
  }
  
  return(results)
}





cat("Utility functions loaded successfully!\n")