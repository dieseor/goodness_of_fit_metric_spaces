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

#' Validate a covariance matrix before simulation
#' @param cov_matrix Numeric square matrix
#' @param symmetry_tol Maximum accepted symmetry gap
#' @param psd_tol Maximum accepted negativity in the minimum eigenvalue
#' @param stop_on_failure Whether to stop when validation fails
#' @return Named list with validation diagnostics
validate_covariance_matrix <- function(cov_matrix,
                                       symmetry_tol = 1e-8,
                                       psd_tol = 1e-8,
                                       stop_on_failure = TRUE) {
  cov_matrix <- as.matrix(cov_matrix)
  if (nrow(cov_matrix) != ncol(cov_matrix)) {
    stop("Covariance matrix must be square.")
  }
  if (!all(is.finite(cov_matrix))) {
    stop("Covariance matrix contains non-finite entries.")
  }

  symmetry_gap <- max(abs(cov_matrix - t(cov_matrix)))
  eigenvalues <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  min_eigenvalue <- min(eigenvalues)
  max_eigenvalue <- max(eigenvalues)
  negative_eigenvalues <- sum(eigenvalues < -psd_tol)
  is_valid <- symmetry_gap <= symmetry_tol && min_eigenvalue >= -psd_tol

  diagnostics <- list(
    valid = is_valid,
    symmetry_gap = symmetry_gap,
    min_eigenvalue = min_eigenvalue,
    max_eigenvalue = max_eigenvalue,
    negative_eigenvalues = negative_eigenvalues,
    symmetry_tol = symmetry_tol,
    psd_tol = psd_tol
  )

  if (!is_valid && isTRUE(stop_on_failure)) {
    reasons <- character(0)
    if (symmetry_gap > symmetry_tol) {
      reasons <- c(
        reasons,
        sprintf("symmetry gap %.3e exceeds tolerance %.3e", symmetry_gap, symmetry_tol)
      )
    }
    if (min_eigenvalue < -psd_tol) {
      reasons <- c(
        reasons,
        sprintf("minimum eigenvalue %.3e is below tolerance %.3e", min_eigenvalue, -psd_tol)
      )
    }
    stop(
      sprintf(
        "Covariance matrix rejected for simulation: %s.",
        paste(reasons, collapse = "; ")
      )
    )
  }

  diagnostics
}

#' Simulate the limiting Gaussian process
#' Works for any distribution once covariance matrix is computed.
#' The covariance matrix is validated before sampling, and execution stops if it
#' is rejected. No automatic corrections are applied.
#' @param cov_matrix Pre-computed covariance matrix
#' @param M Number of Monte Carlo simulations
#' @return Vector of M supremum values from the Gaussian process
simulate_limit_gaussian <- function(cov_matrix, M = 10000, seed = NULL, tol = 1e-10) {
  cov_matrix <- as.matrix(cov_matrix)
  n_total <- nrow(cov_matrix)

  validate_covariance_matrix(
    cov_matrix,
    symmetry_tol = tol,
    psd_tol = tol,
    stop_on_failure = TRUE
  )

  cat("Generating", M, "multivariate normal samples from", n_total, "dimensional process...\n")
  if (!is.null(seed)) set.seed(seed)
  gaussian_samples <- tryCatch(
    mvtnorm::rmvnorm(M, mean = rep(0, n_total), sigma = cov_matrix),
    error = function(e) {
      stop(sprintf("simulate_limit_gaussian failed: %s", e$message))
    }
  )

  supremum_values <- apply(gaussian_samples, 1, function(row) max(abs(row)))
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
# Hyperboloid helpers and HvMF maximum-likelihood estimation on H^2
# ----------------------------------------------------------------------------

#' Compute the Minkowski pseudo-inner product in R^3
#' @param x Numeric vector of length 3
#' @param y Numeric vector of length 3
#' @return Scalar \eqn{-x_1 y_1 + x_2 y_2 + x_3 y_3}
minkowski_inner_product <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (length(x) != 3L || length(y) != 3L) {
    stop("`x` and `y` must both have length 3.")
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("`x` and `y` must be finite.")
  }

  -x[[1]] * y[[1]] + sum(x[-1L] * y[-1L])
}

#' Validate data on the upper sheet of the two-dimensional hyperboloid
#' @param data Numeric matrix with rows in H^2
#' @param tol Tolerance for checking \eqn{\langle x, x \rangle_M = -1}
#' @return Numeric matrix with 3 columns
normalize_hvmf_h2_data <- function(data, tol = 1e-10) {
  if (is.vector(data)) {
    data <- matrix(as.numeric(data), nrow = 1L)
  } else {
    data <- as.matrix(data)
  }

  if (nrow(data) == 0L || ncol(data) != 3L) {
    stop("`data` must be a non-empty n x 3 matrix with rows in H^2.")
  }
  if (any(!is.finite(data))) {
    stop("HvMF data must be finite.")
  }
  if (any(data[, 1L] <= 0)) {
    stop("HvMF data rows must satisfy x[1] > 0.")
  }

  minkowski_norms <- -data[, 1L]^2 + rowSums(data[, -1L, drop = FALSE]^2)
  if (any(abs(minkowski_norms + 1) > tol)) {
    stop("HvMF data rows must satisfy <x, x>_M = -1 up to tolerance.")
  }

  data
}

#' Generate HvMF samples on H^2 using the Nielsen-Okamura GIG mixture
#'
#' This sampler targets the hyperboloid model under the Minkowski convention
#' \eqn{\langle x, y \rangle_M = -x_1 y_1 + x_2 y_2 + x_3 y_3}.
#'
#' @param n Number of samples to draw
#' @param mu Mean direction in H^2
#' @param kappa Strictly positive concentration parameter
#' @param ... Unused; reserved for future extensions
#' @param check Logical flag to validate inputs and outputs
#' @return An n x 3 matrix with rows in H^2
rhvmf_h2_gig <- function(n, mu, kappa, ..., check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  if (isTRUE(check)) {
    mu <- as.numeric(normalize_hvmf_h2_data(mu, tol = 1e-10)[1L, , drop = TRUE])
  } else {
    mu <- as.numeric(mu)
  }

  kappa <- as.numeric(kappa)
  if (length(mu) != 3L || any(!is.finite(mu))) {
    stop("`mu` must be a finite numeric vector of length 3.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be a strictly positive finite scalar.")
  }

  if (!requireNamespace("GIGrvg", quietly = TRUE)) {
    stop("Package 'GIGrvg' is required for rhvmf_h2_gig(). Install it with install.packages(\"GIGrvg\").")
  }

  s <- GIGrvg::rgig(n, lambda = 0.5, chi = 1, psi = kappa^2)
  if (any(!is.finite(s)) || any(s <= 0)) {
    stop("The GIG sampler returned nonpositive or non-finite values.")
  }

  x1 <- stats::rnorm(n, mean = s * kappa * mu[[2L]], sd = sqrt(s))
  x2 <- stats::rnorm(n, mean = s * kappa * mu[[3L]], sd = sqrt(s))
  x0 <- sqrt(1 + x1^2 + x2^2)
  x <- cbind(x0, x1, x2)

  if (isTRUE(check)) {
    minkowski_norms <- -x[, 1L]^2 + rowSums(x[, -1L, drop = FALSE]^2)
    if (any(!is.finite(x)) || any(abs(minkowski_norms + 1) > 1e-8)) {
      stop("Sampler output does not satisfy <x, x>_M = -1 up to tolerance.")
    }
    if (any(x[, 1L] <= 0)) {
      stop("Sampler output must satisfy x[1] > 0.")
    }
  }

  x
}

#' Closed-form MLE for HvMF on H^2, with optional nonnegative weights
#'
#' This implements the Jensen d = 3 closed form on the two-dimensional
#' hyperboloid embedded in R^3 under the Minkowski convention
#' \eqn{\langle x, y \rangle_M = -x_1 y_1 + x_2 y_2 + x_3 y_3}. The returned
#' `xi` is the location point on H^2 (not the canonical parameter `kappa * xi`).
#'
#' @param data Numeric n x 3 matrix with rows in H^2
#' @param weights Optional nonnegative weight vector of length nrow(data)
#' @param tol Numerical tolerance for degeneracy and hyperboloid checks
#' @return List containing the MLE and auxiliary hyperbolic coordinates
hvmf_mle_h2 <- function(data, weights = NULL, tol = 1e-10) {
  data <- normalize_hvmf_h2_data(data, tol = tol)
  n <- nrow(data)

  if (is.null(weights)) {
    weights <- rep.int(1, n)
  }

  weights <- as.numeric(weights)
  if (length(weights) != n) {
    stop("`weights` must have length nrow(data).")
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and nonnegative.")
  }

  W <- sum(weights)
  if (!is.finite(W) || W <= 0) {
    stop("sum(weights) must be strictly positive.")
  }

  S <- colSums(data * weights)
  S_inner <- minkowski_inner_product(S, S)
  R_sq <- -S_inner
  if (!is.finite(R_sq) || R_sq <= 0) {
    stop("The weighted resultant has non-positive or non-finite hyperbolic norm.")
  }

  R <- sqrt(R_sq)
  if (!is.finite(R)) {
    stop("The weighted resultant has non-finite hyperbolic length.")
  }
  if (R <= W + tol) {
    stop("Degenerate or near-degenerate HvMF MLE: R/W must be greater than 1.")
  }

  xi_hat <- S / R
  xi_inner <- minkowski_inner_product(xi_hat, xi_hat)
  if (!is.finite(xi_inner) || abs(xi_inner + 1) > 100 * tol) {
    stop("Internal HvMF MLE check failed: the estimated `xi` is not on H^2.")
  }
  if (xi_hat[[1]] <= 0) {
    stop("Internal HvMF MLE check failed: the estimated `xi` does not satisfy x[1] > 0.")
  }

  kappa_hat <- W / (R - W)
  if (!is.finite(kappa_hat) || kappa_hat <= 0) {
    stop("Internal HvMF MLE check failed: the estimated `kappa` is not positive and finite.")
  }

  sinh_chi_hat <- sqrt(sum(xi_hat[-1L]^2))
  chi_hat <- asinh(sinh_chi_hat)
  theta_hat <- atan2(xi_hat[[3L]], xi_hat[[2L]])
  theta_hat_deg <- (theta_hat * 180 / pi) %% 360

  list(
    xi = xi_hat,
    mu = xi_hat,
    kappa = kappa_hat,
    chi = chi_hat,
    sinh_chi = sinh_chi_hat,
    theta = theta_hat,
    theta_deg = theta_hat_deg,
    R = R,
    W = W,
    xi_inner = xi_inner,
    resultant = S
  )
}

#' Geodesic distance on H^2 under the Minkowski convention (-,+,+)
#' @param x Numeric vector of length 3 in H^2
#' @param y Numeric vector of length 3 in H^2
#' @param tol Numerical tolerance used in the acosh clipping
#' @return Nonnegative geodesic distance \eqn{\operatorname{acosh}(-\langle x, y \rangle_M)}
hyperbolic_geodesic_distance_h2 <- function(x, y, tol = 1e-12) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (length(x) != 3L || length(y) != 3L) {
    stop("`x` and `y` must both have length 3.")
  }

  inner_xy <- minkowski_inner_product(x, y)
  acosh(pmax(-inner_xy, 1))
}

#' Pairwise geodesic distances between H^2 grid points and data
#' @param omega_grid Numeric matrix whose rows are reference points in H^2
#' @param data Numeric matrix whose rows are observations in H^2
#' @param tol Numerical tolerance used in the acosh clipping
#' @return Matrix with `nrow(omega_grid)` rows and `nrow(data)` columns
hvmf_distance_matrix <- function(omega_grid, data, tol = 1e-12) {
  omega_matrix <- normalize_hvmf_h2_data(omega_grid, tol = tol)
  data_matrix <- normalize_hvmf_h2_data(data, tol = tol)

  minkowski_weighted_data <- data_matrix
  minkowski_weighted_data[, 1L] <- -minkowski_weighted_data[, 1L]

  inner_products <- omega_matrix %*% t(minkowski_weighted_data)
  acosh(pmax(-inner_products, 1))
}

#' Theoretical distance profile for HvMF on H^2
#' @param omega Reference point in H^2
#' @param mu Mean direction in H^2
#' @param kappa Positive concentration parameter
#' @param t_values Nonnegative evaluation points
#' @return Numeric vector of probabilities \eqn{P\{d_H(X,\omega)\le t\}}
theoretical_distance_profile_hvmf <- function(omega, mu, kappa, t_values) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  t_values <- as.numeric(t_values)

  if (length(omega) != 3L || length(mu) != 3L) {
    stop("`omega` and `mu` must both have length 3.")
  }

  mu_omega <- minkowski_inner_product(mu, omega)
  q <- length(mu) - 1L

  sapply(t_values, function(t) {
    density_R <- function(r) {
      nu <- (q - 1) / 2
      log_bessel_k <- log(besselK(kappa, nu = nu, expon.scaled = TRUE)) - kappa

      log_sinh_power_term <- (q - 1) * log(sinh(r))
      log_exp_term <- kappa * mu_omega * cosh(r)
      log_numerator <- nu * log(kappa) - nu * log(2 * pi) - log(2) - log_bessel_k +
        log_sinh_power_term + log_exp_term

      mu_omega_sq_minus_1 <- mu_omega^2 - 1
      if (mu_omega_sq_minus_1 < 0) {
        if (abs(mu_omega_sq_minus_1) < 1e-12) {
          mu_omega_sq_term <- 0
        } else {
          stop(paste(
            "Mathematical error: mu_omega^2 - 1 =",
            mu_omega_sq_minus_1,
            "unexpected negative value inside square root."
          ))
        }
      } else {
        mu_omega_sq_term <- mu_omega_sq_minus_1
      }

      kappa_vmf <- kappa * sinh(r) * sqrt(mu_omega_sq_term)
      log_denominator <- rotasym::c_vMF(p = q, kappa = kappa_vmf, log = TRUE)

      exp(log_numerator - log_denominator)
    }

    cdf_result <- integrate(
      density_R,
      lower = 0,
      upper = t,
      rel.tol = 1e-6,
      abs.tol = 1e-8,
      stop.on.error = FALSE
    )

    cdf_result$value
  })
}

#' Projected density on H^2 for Y = -<omega, X>_M under HvMF(mu, kappa)
#' @param y Scalar or vector with values in [1, \infty)
#' @param alpha Scalar \eqn{\alpha = -\langle \mu, \omega \rangle_M \ge 1}
#' @param kappa Positive concentration parameter
#' @return Density values \eqn{f_{\alpha,\kappa}(y)}
hvmf_projection_density_h2 <- function(y, alpha, kappa) {
  y <- as.numeric(y)
  alpha <- as.numeric(alpha)
  kappa <- as.numeric(kappa)

  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 1) {
    stop("`alpha` must be a finite scalar with alpha >= 1.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be a strictly positive finite scalar.")
  }

  output <- numeric(length(y))
  valid <- is.finite(y) & (y >= 1)
  if (!any(valid)) {
    return(output)
  }

  beta_sq <- max(alpha^2 - 1, 0)
  beta <- sqrt(beta_sq)
  y_valid <- y[valid]
  radial_term <- sqrt(pmax(y_valid^2 - 1, 0))

  if (beta == 0) {
    log_density <- log(kappa) + kappa - kappa * y_valid
  } else {
    z <- kappa * beta * radial_term
    log_i0 <- log(besselI(z, nu = 0, expon.scaled = TRUE)) + z
    log_density <- log(kappa) + kappa - kappa * alpha * y_valid + log_i0
  }

  output[valid] <- exp(log_density)
  output
}

#' Tabulate the projected CDF on H^2 for Y = -<omega, X>_M under HvMF(mu, kappa)
#' @param alpha Scalar \eqn{\alpha = -\langle \mu, \omega \rangle_M \ge 1}
#' @param kappa Positive concentration parameter
#' @param y_max Maximum y-value to cover in the table
#' @param grid_size Number of tabulation points when `grid` is not supplied
#' @param grid Optional user-supplied grid in [1, y_max]
#' @return Data frame with columns `y` and `cdf`
hvmf_tabulate_projection_cdf_h2 <- function(alpha,
                                            kappa,
                                            y_max,
                                            grid_size = 4097L,
                                            grid = NULL) {
  alpha <- as.numeric(alpha)
  kappa <- as.numeric(kappa)
  y_max <- as.numeric(y_max)

  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 1) {
    stop("`alpha` must be a finite scalar with alpha >= 1.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be a strictly positive finite scalar.")
  }
  if (length(y_max) != 1L || !is.finite(y_max) || y_max < 1) {
    stop("`y_max` must be a finite scalar with y_max >= 1.")
  }

  if (is.null(grid)) {
    grid_size <- as.integer(grid_size)
    if (!is.finite(grid_size) || grid_size < 2L) {
      stop("`grid_size` must be an integer >= 2.")
    }
    y_grid <- seq(1, y_max, length.out = grid_size)
  } else {
    y_grid <- sort(as.numeric(grid))
    if (length(y_grid) < 2L) {
      stop("`grid` must contain at least two points.")
    }
    if (any(!is.finite(y_grid))) {
      stop("`grid` must be finite.")
    }
    if (min(y_grid) < 1 || max(y_grid) > y_max) {
      stop("`grid` must be contained in [1, y_max].")
    }
    if (y_grid[[1L]] > 1) {
      y_grid <- c(1, y_grid)
    }
    if (y_grid[[length(y_grid)]] < y_max) {
      y_grid <- c(y_grid, y_max)
    }
    y_grid <- unique(y_grid)
  }

  density_values <- hvmf_projection_density_h2(
    y = y_grid,
    alpha = alpha,
    kappa = kappa
  )

  cdf_values <- numeric(length(y_grid))
  if (length(y_grid) >= 2L) {
    increments <- diff(y_grid) * (density_values[-1L] + density_values[-length(y_grid)]) / 2
    cdf_values[-1L] <- cumsum(increments)
  }

  cdf_values <- pmin(pmax(cdf_values, 0), 1)
  cdf_values[[1L]] <- 0

  data.frame(
    y = y_grid,
    cdf = cdf_values
  )
}

#' Evaluate a tabulated projected CDF on H^2 by linear interpolation
#' @param y Scalar or vector with values in [1, \infty)
#' @param cdf_table Output of `hvmf_tabulate_projection_cdf_h2()`
#' @return Interpolated CDF values clipped to [0, 1]
hvmf_eval_projection_cdf_tabulated <- function(y, cdf_table) {
  y <- as.numeric(y)

  if (!is.data.frame(cdf_table) || !all(c("y", "cdf") %in% names(cdf_table))) {
    stop("`cdf_table` must be a data frame with columns `y` and `cdf`.")
  }

  y_grid <- as.numeric(cdf_table$y)
  cdf_grid <- as.numeric(cdf_table$cdf)
  if (length(y_grid) == 0L) {
    stop("`cdf_table` cannot be empty.")
  }
  if (length(y_grid) == 1L) {
    output <- rep(cdf_grid[[1L]], length(y))
    output[y <= y_grid[[1L]]] <- 0
    return(pmin(pmax(output, 0), 1))
  }

  output <- approx(
    x = y_grid,
    y = cdf_grid,
    xout = y,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y

  output[y <= 1] <- 0
  pmin(pmax(output, 0), 1)
}

#' Tabulated CvM theoretical profile matrix for HvMF on H^2
#' @param data Numeric n x 3 matrix with rows in H^2
#' @param theta List with entries `mu` and `kappa`
#' @param grid_size Number of tabulation points per row
#' @param distance_matrix Optional matrix of pairwise geodesic distances
#' @param tol Numerical tolerance for hyperboloid checks
#' @return Numeric n x n matrix with entries F_{X_i}^{theta}(d_H(X_i, X_j))
hvmf_cvm_profile_matrix_tabulated <- function(data,
                                              theta,
                                              grid_size = 4097L,
                                              distance_matrix = NULL,
                                              tol = 1e-10) {
  x <- normalize_hvmf_h2_data(data, tol = tol)
  if (!is.list(theta)) {
    stop("`theta` must be a list containing `mu` and `kappa`.")
  }
  mu_matrix <- normalize_hvmf_h2_data(theta$mu, tol = tol)
  mu <- as.numeric(mu_matrix[1L, , drop = TRUE])
  kappa <- as.numeric(theta$kappa)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("HvMF theta requires a strictly positive finite scalar `kappa`.")
  }
  n <- nrow(x)

  if (is.null(distance_matrix)) {
    y_matrix <- pmax(cosh(t(hvmf_distance_matrix(x, x))), 1)
  } else {
    distance_matrix <- as.matrix(distance_matrix)
    if (!all(dim(distance_matrix) == c(n, n))) {
      stop("`distance_matrix` must be an n x n matrix compatible with `data`.")
    }
    y_matrix <- pmax(cosh(distance_matrix), 1)
  }

  alpha_values <- as.numeric(x %*% c(mu[[1L]], -mu[[2L]], -mu[[3L]]))
  alpha_values <- pmax(alpha_values, 1)

  output <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    row_y <- pmax(y_matrix[i, ], 1)
    cdf_table <- hvmf_tabulate_projection_cdf_h2(
      alpha = alpha_values[[i]],
      kappa = kappa,
      y_max = max(row_y),
      grid_size = grid_size
    )
    output[i, ] <- hvmf_eval_projection_cdf_tabulated(row_y, cdf_table)
  }

  output
}


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

# ----------------------------------------------------------------------------
# S^1-specific helpers for deterministic angular computations
# ----------------------------------------------------------------------------

#' Wrap angles to [0, 2*pi)
#' @param theta Numeric scalar or vector of angles
#' @return Angles reduced modulo 2*pi
wrap_angle_2pi <- function(theta) {
  wrapped <- theta %% (2 * pi)
  wrapped[wrapped < 0] <- wrapped[wrapped < 0] + 2 * pi
  wrapped
}

#' Generate a deterministic angular grid on the unit circle
#' @param n_angles Number of angles on S^1
#' @param theta0 Starting angle offset
#' @return Data frame with columns theta, x, y, label
generate_circle_grid <- function(n_angles, theta0 = 0) {
  if (!is.numeric(n_angles) || length(n_angles) != 1 || n_angles < 2) {
    stop("`n_angles` must be a single integer >= 2.")
  }
  theta <- wrap_angle_2pi(theta0 + 2 * pi * (0:(n_angles - 1)) / n_angles)
  theta <- sort(theta)
  data.frame(
    theta = theta,
    x = cos(theta),
    y = sin(theta),
    label = sprintf("%.2f", theta),
    stringsAsFactors = FALSE
  )
}

#' Convert a point on S^1 to its angular coordinate
#' @param omega Numeric vector of length 2
#' @param tol Tolerance for the unit-norm check
#' @return Angle in [0, 2*pi)
circle_angle_from_point <- function(omega, tol = 1e-8) {
  omega <- as.numeric(omega)
  if (length(omega) != 2) {
    stop("`omega` must have length 2 for S^1 computations.")
  }
  omega_norm <- sqrt(sum(omega^2))
  if (abs(omega_norm - 1) > tol) {
    stop("`omega` must have unit norm for S^1 computations.")
  }
  wrap_angle_2pi(atan2(omega[2], omega[1]))
}

#' Convert a matrix of points on S^1 to their angular coordinates
#' @param omega_grid Matrix with 2 columns
#' @return Numeric vector of angles in [0, 2*pi)
circle_angles_from_matrix <- function(omega_grid) {
  omega_grid <- as.matrix(omega_grid)
  if (ncol(omega_grid) != 2) {
    stop("`omega_grid` must have exactly 2 columns for S^1 computations.")
  }
  apply(omega_grid, 1, circle_angle_from_point)
}

#' Angular half-width for chordal balls on S^1
#' @param t Chordal distance threshold
#' @return Half-width delta(t) = 2 * asin(t / 2)
s1_chordal_half_width <- function(t) {
  t <- pmin(pmax(t, 0), 2)
  2 * asin(t / 2)
}

#' Convert a wrapped angular interval to non-wrapping segments
#' @param start Interval start
#' @param end Interval end
#' @param tol Tolerance for detecting empty/full intervals
#' @return Matrix with 0, 1, or 2 rows and columns start/end
s1_interval_to_segments <- function(start, end, tol = 1e-12) {
  width <- end - start
  if (width <= tol) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }
  if (width >= 2 * pi - tol) {
    return(matrix(c(0, 2 * pi), ncol = 2, byrow = TRUE,
                  dimnames = list(NULL, c("start", "end"))))
  }

  start_wrapped <- wrap_angle_2pi(start)
  end_wrapped <- wrap_angle_2pi(end)
  if (start_wrapped < end_wrapped) {
    return(matrix(c(start_wrapped, end_wrapped), ncol = 2, byrow = TRUE,
                  dimnames = list(NULL, c("start", "end"))))
  }

  rbind(
    c(0, end_wrapped),
    c(start_wrapped, 2 * pi)
  )
}

#' Wrapped arc associated with a chordal ball on S^1
#' @param omega Point on S^1
#' @param t Chordal distance threshold
#' @param tol Numerical tolerance
#' @return Matrix of non-wrapping segments covering the event arc
s1_event_segments_chordal <- function(omega, t, tol = 1e-12) {
  if (t <= tol) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }
  if (t >= 2 - tol) {
    return(matrix(c(0, 2 * pi), ncol = 2, byrow = TRUE,
                  dimnames = list(NULL, c("start", "end"))))
  }
  theta <- circle_angle_from_point(omega)
  delta <- s1_chordal_half_width(t)
  s1_interval_to_segments(theta - delta, theta + delta, tol = tol)
}

#' Intersect two segment unions on S^1
#' @param seg1 Matrix of interval segments
#' @param seg2 Matrix of interval segments
#' @param tol Numerical tolerance
#' @return Matrix of intersected segments
s1_intersect_segments <- function(seg1, seg2, tol = 1e-12) {
  if (length(seg1) == 0 || length(seg2) == 0 || nrow(seg1) == 0 || nrow(seg2) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }

  intersections <- list()
  idx <- 0
  for (i in seq_len(nrow(seg1))) {
    for (j in seq_len(nrow(seg2))) {
      lower <- max(seg1[i, 1], seg2[j, 1])
      upper <- min(seg1[i, 2], seg2[j, 2])
      if (upper - lower > tol) {
        idx <- idx + 1
        intersections[[idx]] <- c(lower, upper)
      }
    }
  }

  if (length(intersections) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }

  do.call(rbind, intersections)
}

#' Angular density of the vMF distribution on S^1
#' @param phi Angle(s) in radians
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param log Whether to return the log-density
#' @return Density values
vmf_s1_angle_density <- function(phi, mu, kappa, log = FALSE) {
  mu <- as.numeric(mu)
  if (length(mu) != 2) {
    stop("`mu` must have length 2 for S^1 computations.")
  }
  mu <- mu / sqrt(sum(mu^2))
  mu_angle <- circle_angle_from_point(mu)
  log_i0 <- log(besselI(kappa, nu = 0, expon.scaled = TRUE)) + kappa
  log_density <- kappa * cos(phi - mu_angle) - log(2 * pi) - log_i0
  if (log) {
    return(log_density)
  }
  exp(log_density)
}

#' Build a deterministic angular CDF approximation for vMF on S^1
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param n_grid Number of points in the angular grid
#' @return List with phi, density, and cdf
build_vmf_s1_cdf <- function(mu, kappa, n_grid = 16385) {
  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid < 3) {
    stop("`n_grid` must be an integer >= 3.")
  }
  phi <- seq(0, 2 * pi, length.out = n_grid)
  density <- vmf_s1_angle_density(phi, mu, kappa)
  dx <- diff(phi)
  cdf <- numeric(n_grid)
  cdf[-1] <- cumsum((density[-n_grid] + density[-1]) * dx / 2)
  total_mass <- cdf[n_grid]
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("Failed to build a deterministic angular CDF for vMF on S^1.")
  }
  cdf <- cdf / total_mass
  cdf[n_grid] <- 1
  list(
    mu = as.numeric(mu) / sqrt(sum(mu^2)),
    kappa = kappa,
    phi = phi,
    density = density,
    cdf = cdf
  )
}

#' Evaluate a deterministic angular CDF approximation on S^1
#' @param phi Angle(s)
#' @param cdf_object Output of build_vmf_s1_cdf()
#' @param tol Tolerance for identifying 2*pi
#' @return CDF values
evaluate_vmf_s1_cdf <- function(phi, cdf_object, tol = 1e-12) {
  phi <- as.numeric(phi)
  phi_wrapped <- phi %% (2 * pi)
  is_full_turn <- abs(phi_wrapped) < tol & phi > 0
  phi_wrapped[is_full_turn] <- 2 * pi
  approx(
    x = cdf_object$phi,
    y = cdf_object$cdf,
    xout = phi_wrapped,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}

#' Probability of a union of angular segments under vMF on S^1
#' @param segments Matrix of segments
#' @param cdf_object Output of build_vmf_s1_cdf()
#' @return Probability mass of the segment union
vmf_s1_segments_probability <- function(segments, cdf_object) {
  if (length(segments) == 0 || nrow(segments) == 0) {
    return(0)
  }
  sum(vapply(seq_len(nrow(segments)), function(i) {
    evaluate_vmf_s1_cdf(segments[i, 2], cdf_object) -
      evaluate_vmf_s1_cdf(segments[i, 1], cdf_object)
  }, numeric(1)))
}

#' Exact deterministic distance profile for vMF on S^1 with chordal distance
#' @param omega Reference point on S^1
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param t_values Chordal thresholds
#' @param cdf_object Optional deterministic angular CDF
#' @param cdf_grid_size Grid size used when cdf_object is not supplied
#' @return Vector of probabilities
theoretical_distance_profile_vmf_s1_chordal <- function(omega,
                                                        mu,
                                                        kappa,
                                                        t_values,
                                                        cdf_object = NULL,
                                                        cdf_grid_size = 16385) {
  if (is.matrix(omega)) {
    n <- nrow(omega)
    if (length(t_values) == 1) t_values <- rep(t_values, n)
    stopifnot(length(t_values) == n)
    return(vapply(seq_len(n), function(i) {
      theoretical_distance_profile_vmf_s1_chordal(
        omega = omega[i, ],
        mu = mu,
        kappa = kappa,
        t_values = t_values[i],
        cdf_object = cdf_object,
        cdf_grid_size = cdf_grid_size
      )
    }, numeric(1)))
  }

  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }

  t_values <- as.numeric(t_values)
  vapply(t_values, function(t) {
    if (t <= 0) return(0)
    if (t >= 2) return(1)
    segments <- s1_event_segments_chordal(omega, t)
    vmf_s1_segments_probability(segments, cdf_object)
  }, numeric(1))
}

#' Exact deterministic joint probability for chordal balls on S^1
#' @param omega1 First point on S^1
#' @param t1 First threshold
#' @param omega2 Second point on S^1
#' @param t2 Second threshold
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param cdf_object Optional deterministic angular CDF
#' @param cdf_grid_size Grid size used when cdf_object is not supplied
#' @return Joint probability
joint_probability_vmf_s1_chordal_exact <- function(omega1,
                                                   t1,
                                                   omega2,
                                                   t2,
                                                   mu,
                                                   kappa,
                                                   cdf_object = NULL,
                                                   cdf_grid_size = 16385) {
  if (t1 <= 0 || t2 <= 0) return(0)
  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }
  seg1 <- s1_event_segments_chordal(omega1, t1)
  seg2 <- s1_event_segments_chordal(omega2, t2)
  overlap <- s1_intersect_segments(seg1, seg2)
  vmf_s1_segments_probability(overlap, cdf_object)
}

#' Deterministic inversion of the S^1 chordal distance profile
#' @param omega Reference point on S^1
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param u_values Probabilities in [0, 1]
#' @param cdf_object Optional deterministic angular CDF
#' @param cdf_grid_size Grid size used when cdf_object is not supplied
#' @param tol Root-finding tolerance
#' @return Vector of chordal thresholds t in [0, 2]
invert_distance_profile_vmf_s1_chordal <- function(omega,
                                                   mu,
                                                   kappa,
                                                   u_values,
                                                   cdf_object = NULL,
                                                   cdf_grid_size = 16385,
                                                   tol = 1e-8) {
  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }

  u_values <- as.numeric(u_values)
  vapply(u_values, function(u) {
    if (u <= 0) return(0)
    if (u >= 1) return(2)
    root_fun <- function(t) {
      theoretical_distance_profile_vmf_s1_chordal(
        omega = omega,
        mu = mu,
        kappa = kappa,
        t_values = t,
        cdf_object = cdf_object,
        cdf_grid_size = cdf_grid_size
      ) - u
    }
    uniroot(root_fun, interval = c(0, 2), tol = tol)$root
  }, numeric(1))
}

#' Deterministic covariance matrix for the simple vMF process on S^1
#' @param omega_grid Matrix of points on S^1
#' @param t_grid Vector of chordal thresholds
#' @param mu Mean direction on S^1
#' @param kappa Concentration parameter
#' @param cdf_object Optional deterministic angular CDF
#' @param cdf_grid_size Grid size used when cdf_object is not supplied
#' @return Covariance matrix on the flattened grid
cov_vmf_s1_simple_exact <- function(omega_grid,
                                    t_grid,
                                    mu,
                                    kappa,
                                    cdf_object = NULL,
                                    cdf_grid_size = 16385) {
  omega_grid <- as.matrix(omega_grid)
  t_grid <- as.numeric(t_grid)
  if (ncol(omega_grid) != 2) {
    stop("`omega_grid` must have exactly 2 columns for exact S^1 covariance.")
  }
  if (any(t_grid < 0) || any(t_grid > 2)) {
    stop("`t_grid` must lie in [0, 2] for chordal distance on S^1.")
  }
  if (is.null(cdf_object)) {
    cdf_object <- build_vmf_s1_cdf(mu, kappa, n_grid = cdf_grid_size)
  }

  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  F_matrix <- t(vapply(seq_len(n_omega), function(i) {
    theoretical_distance_profile_vmf_s1_chordal(
      omega = omega_grid[i, ],
      mu = mu,
      kappa = kappa,
      t_values = t_grid,
      cdf_object = cdf_object,
      cdf_grid_size = cdf_grid_size
    )
  }, numeric(n_t)))
  F_vec <- as.vector(F_matrix)

  segment_list <- vector("list", n_total)
  for (t_idx in seq_len(n_t)) {
    for (omega_idx in seq_len(n_omega)) {
      flat_idx <- omega_idx + (t_idx - 1) * n_omega
      segment_list[[flat_idx]] <- s1_event_segments_chordal(
        omega = omega_grid[omega_idx, ],
        t = t_grid[t_idx]
      )
    }
  }

  cov_matrix <- matrix(0, nrow = n_total, ncol = n_total)
  for (i in seq_len(n_total)) {
    for (j in i:n_total) {
      joint_prob <- vmf_s1_segments_probability(
        s1_intersect_segments(segment_list[[i]], segment_list[[j]]),
        cdf_object = cdf_object
      )
      cov_value <- joint_prob - F_vec[i] * F_vec[j]
      cov_matrix[i, j] <- cov_value
      cov_matrix[j, i] <- cov_value
    }
  }

  cov_matrix
}

#' Convert a distance threshold on the sphere into a dot-product threshold
#' @param t Distance threshold
#' @param distance_type Either "chordal" or "geodesic"
#' @return Threshold a such that d(x, omega) <= t iff omega^T x >= a
sphere_distance_to_dot_threshold <- function(t, distance_type = "chordal") {
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  t <- as.numeric(t)
  if (distance_type == "chordal") {
    threshold <- 1 - (t^2) / 2
  } else {
    threshold <- cos(t)
  }
  pmin(pmax(threshold, -1), 1)
}

#' Exact projected density on S^2 for U = omega^T X under vMF(mu, kappa)
#' @param u Scalar or vector in [-1, 1]
#' @param omega Unit vector in R^3
#' @param mu Mean direction in R^3
#' @param kappa Concentration parameter
#' @param log Whether to return log-density
#' @return Density values of U
vmf_s2_projected_density <- function(u, omega, mu, kappa, log = FALSE) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3 || length(mu) != 3) {
    stop("`omega` and `mu` must both have length 3 for S^2 projected densities.")
  }
  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))

  u <- as.numeric(u)
  out <- rep(if (log) -Inf else 0, length(u))
  valid <- u >= -1 & u <= 1
  if (!any(valid)) {
    return(out)
  }

  m1 <- sum(mu * omega)
  m1 <- pmin(pmax(m1, -1), 1)
  radial_sq <- pmax(0, 1 - u[valid]^2)
  tangential_sq <- max(0, 1 - m1^2)
  lambda <- kappa * sqrt(radial_sq * tangential_sq)

  log_density <- rotasym::c_vMF(p = 3, kappa = kappa, log = TRUE) +
    kappa * m1 * u[valid] -
    rotasym::c_vMF(p = 2, kappa = lambda, log = TRUE)

  out[valid] <- if (log) log_density else exp(log_density)
  out
}

#' Exact symmetric-arc probability for a von Mises law on S^1
#' @param lambda Concentration parameter of the S^1 conditional law
#' @param alpha Angle between the mean direction and the arc center
#' @param b Dot-product threshold defining the arc
#' @param rel.tol Relative tolerance for integrate()
#' @param abs.tol Absolute tolerance for integrate()
#' @param subdivisions Maximum subdivisions for integrate()
#' @param tol Numerical tolerance
#' @return P(cos(Phi) >= b) for Phi distributed as vM(alpha, lambda)
vmf_s1_cap_probability_integral <- function(lambda,
                                            alpha,
                                            b,
                                            rel.tol = 1e-8,
                                            abs.tol = 1e-10,
                                            subdivisions = 200L,
                                            tol = 1e-10) {
  lambda <- as.numeric(lambda)
  alpha <- as.numeric(alpha)
  b <- as.numeric(b)

  if (b <= -1 + tol) return(1)
  if (b >= 1 - tol) return(0)

  delta <- acos(pmin(pmax(b, -1), 1))
  if (delta >= pi - tol) return(1)
  if (lambda <= tol) {
    return(delta / pi)
  }

  i0_scaled <- besselI(lambda, nu = 0, expon.scaled = TRUE)
  integrand <- function(phi) {
    exp(lambda * (cos(phi - alpha) - 1)) / (2 * pi * i0_scaled)
  }
  integrate(
    f = integrand,
    lower = -delta,
    upper = delta,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions
  )$value
}

#' Exact marginal distance profile on S^2 under vMF by 1D integration
#' @param omega Reference point on S^2
#' @param mu Mean direction in R^3
#' @param kappa Concentration parameter
#' @param t_values Distance thresholds
#' @param distance_type Either "chordal" or "geodesic"
#' @param rel.tol Relative tolerance for integrate()
#' @param abs.tol Absolute tolerance for integrate()
#' @param subdivisions Maximum subdivisions for integrate()
#' @return Vector of probabilities P(d(X, omega) <= t)
distance_profile_vmf_s2_integral <- function(omega,
                                             mu,
                                             kappa,
                                             t_values,
                                             distance_type = "chordal",
                                             rel.tol = 1e-8,
                                             abs.tol = 1e-10,
                                             subdivisions = 200L) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3 || length(mu) != 3) {
    stop("`omega` and `mu` must both have length 3 for S^2 distance profiles.")
  }

  if (is.matrix(omega)) {
    n <- nrow(omega)
    if (length(t_values) == 1) t_values <- rep(t_values, n)
    stopifnot(length(t_values) == n)
    return(vapply(seq_len(n), function(i) {
      distance_profile_vmf_s2_integral(
        omega = omega[i, ],
        mu = mu,
        kappa = kappa,
        t_values = t_values[i],
        distance_type = distance_type,
        rel.tol = rel.tol,
        abs.tol = abs.tol,
        subdivisions = subdivisions
      )
    }, numeric(1)))
  }

  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  t_values <- as.numeric(t_values)

  vapply(t_values, function(t) {
    a <- sphere_distance_to_dot_threshold(t, distance_type)
    if (a <= -1) return(1)
    if (a >= 1) return(0)
    integrate(
      f = function(u) vmf_s2_projected_density(u, omega = omega, mu = mu, kappa = kappa),
      lower = a,
      upper = 1,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      subdivisions = subdivisions
    )$value
  }, numeric(1))
}

#' Exact joint probability for two distance balls on S^2 under vMF by conditioning
#' @param omega1 First point on S^2
#' @param t1 First distance threshold
#' @param omega2 Second point on S^2
#' @param t2 Second distance threshold
#' @param mu Mean direction in R^3
#' @param kappa Concentration parameter
#' @param distance_type Either "chordal" or "geodesic"
#' @param rel.tol_outer Relative tolerance for the outer integral
#' @param abs.tol_outer Absolute tolerance for the outer integral
#' @param rel.tol_inner Relative tolerance for the inner S^1 integral
#' @param abs.tol_inner Absolute tolerance for the inner S^1 integral
#' @param subdivisions_outer Maximum subdivisions for the outer integral
#' @param subdivisions_inner Maximum subdivisions for the inner integral
#' @param tol Numerical tolerance
#' @return Joint probability P(d(X, omega1) <= t1, d(X, omega2) <= t2)
joint_probability_vmf_s2_simple_integral <- function(omega1,
                                                     t1,
                                                     omega2,
                                                     t2,
                                                     mu,
                                                     kappa,
                                                     distance_type = "chordal",
                                                     rel.tol_outer = 1e-7,
                                                     abs.tol_outer = 1e-9,
                                                     rel.tol_inner = 1e-8,
                                                     abs.tol_inner = 1e-10,
                                                     subdivisions_outer = 200L,
                                                     subdivisions_inner = 200L,
                                                     tol = 1e-10) {
  omega1 <- as.numeric(omega1)
  omega2 <- as.numeric(omega2)
  mu <- as.numeric(mu)
  if (length(omega1) != 3 || length(omega2) != 3 || length(mu) != 3) {
    stop("`omega1`, `omega2`, and `mu` must all have length 3 for S^2 exact probabilities.")
  }
  omega1 <- omega1 / sqrt(sum(omega1^2))
  omega2 <- omega2 / sqrt(sum(omega2^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))

  a1 <- sphere_distance_to_dot_threshold(t1, distance_type)
  a2 <- sphere_distance_to_dot_threshold(t2, distance_type)
  if (a1 >= 1 || a2 >= 1) return(0)
  if (a1 <= -1 && a2 <= -1) return(1)

  rho <- sum(omega1 * omega2)
  rho <- pmin(pmax(rho, -1), 1)
  m1 <- sum(mu * omega1)
  m1 <- pmin(pmax(m1, -1), 1)

  density_u <- function(u) {
    vmf_s2_projected_density(u, omega = omega1, mu = mu, kappa = kappa)
  }

  if (rho >= 1 - tol) {
    lower <- max(a1, a2, -1)
    if (lower >= 1) return(0)
    return(integrate(
      f = density_u,
      lower = lower,
      upper = 1,
      rel.tol = rel.tol_outer,
      abs.tol = abs.tol_outer,
      subdivisions = subdivisions_outer
    )$value)
  }

  if (rho <= -1 + tol) {
    lower <- max(a1, -1)
    upper <- min(1, -a2)
    if (upper <= lower + tol) return(0)
    return(integrate(
      f = density_u,
      lower = lower,
      upper = upper,
      rel.tol = rel.tol_outer,
      abs.tol = abs.tol_outer,
      subdivisions = subdivisions_outer
    )$value)
  }

  denom_rho <- sqrt(max(0, 1 - rho^2))
  denom_mu <- sqrt(max(0, 1 - m1^2))
  cos_alpha <- if (denom_mu <= tol) {
    1
  } else {
    (sum(mu * omega2) - m1 * rho) / (denom_mu * denom_rho)
  }
  alpha <- acos(pmin(pmax(cos_alpha, -1), 1))

  lower <- max(a1, -1)
  if (lower >= 1) return(0)

  integrand <- function(u_vec) {
    vapply(u_vec, function(u) {
      radial <- sqrt(max(0, 1 - u^2))
      if (radial <= tol) {
        inner_prob <- as.numeric(rho * u >= a2 - tol)
      } else {
        b <- (a2 - rho * u) / (denom_rho * radial)
        lambda <- kappa * radial * denom_mu
        inner_prob <- vmf_s1_cap_probability_integral(
          lambda = lambda,
          alpha = alpha,
          b = b,
          rel.tol = rel.tol_inner,
          abs.tol = abs.tol_inner,
          subdivisions = subdivisions_inner,
          tol = tol
        )
      }
      density_u(u) * inner_prob
    }, numeric(1))
  }

  integrate(
    f = integrand,
    lower = lower,
    upper = 1,
    rel.tol = rel.tol_outer,
    abs.tol = abs.tol_outer,
    subdivisions = subdivisions_outer
  )$value
}

#' Exact covariance matrix on S^2 for the simple vMF process by nested integration
#' @param omega_grid Matrix of points on S^2
#' @param t_grid Vector of distance thresholds
#' @param mu Mean direction in R^3
#' @param kappa Concentration parameter
#' @param distance_type Either "chordal" or "geodesic"
#' @param rel.tol_outer Relative tolerance for outer integrals
#' @param abs.tol_outer Absolute tolerance for outer integrals
#' @param rel.tol_inner Relative tolerance for inner S^1 integrals
#' @param abs.tol_inner Absolute tolerance for inner S^1 integrals
#' @param subdivisions_outer Maximum subdivisions for outer integrals
#' @param subdivisions_inner Maximum subdivisions for inner integrals
#' @param tol Numerical tolerance
#' @return Covariance matrix under the simple null
cov_vmf_s2_simple_integral <- function(omega_grid,
                                       t_grid,
                                       mu,
                                       kappa,
                                       distance_type = "chordal",
                                       rel.tol_outer = 1e-7,
                                       abs.tol_outer = 1e-9,
                                       rel.tol_inner = 1e-8,
                                       abs.tol_inner = 1e-10,
                                       subdivisions_outer = 200L,
                                       subdivisions_inner = 200L,
                                       tol = 1e-10) {
  omega_grid <- as.matrix(omega_grid)
  t_grid <- as.numeric(t_grid)
  mu <- as.numeric(mu)
  if (ncol(omega_grid) != 3 || length(mu) != 3) {
    stop("`omega_grid` must have 3 columns and `mu` must have length 3 for S^2 exact covariance.")
  }

  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t

  F_matrix <- t(vapply(seq_len(n_omega), function(i) {
    distance_profile_vmf_s2_integral(
      omega = omega_grid[i, ],
      mu = mu,
      kappa = kappa,
      t_values = t_grid,
      distance_type = distance_type,
      rel.tol = rel.tol_outer,
      abs.tol = abs.tol_outer,
      subdivisions = subdivisions_outer
    )
  }, numeric(n_t)))
  F_vec <- as.vector(F_matrix)

  omega_idx_vec <- ((0:(n_total - 1)) %% n_omega) + 1
  t_idx_vec <- ((0:(n_total - 1)) %/% n_omega) + 1
  cov_matrix <- matrix(0, nrow = n_total, ncol = n_total)

  for (i in seq_len(n_total)) {
    omega1 <- omega_grid[omega_idx_vec[i], ]
    t1 <- t_grid[t_idx_vec[i]]
    for (j in i:n_total) {
      joint_prob <- joint_probability_vmf_s2_simple_integral(
        omega1 = omega1,
        t1 = t1,
        omega2 = omega_grid[omega_idx_vec[j], ],
        t2 = t_grid[t_idx_vec[j]],
        mu = mu,
        kappa = kappa,
        distance_type = distance_type,
        rel.tol_outer = rel.tol_outer,
        abs.tol_outer = abs.tol_outer,
        rel.tol_inner = rel.tol_inner,
        abs.tol_inner = abs.tol_inner,
        subdivisions_outer = subdivisions_outer,
        subdivisions_inner = subdivisions_inner,
        tol = tol
      )
      cov_value <- joint_prob - F_vec[i] * F_vec[j]
      cov_matrix[i, j] <- cov_value
      cov_matrix[j, i] <- cov_value
    }
  }

  cov_matrix
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


# -------------------------------------------------------------------------
# BAHADUR helpers (MLE, psi, dot_psi and trajectory analysis)
# -------------------------------------------------------------------------

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
  if (length(mu) == 2 && distance_type == "chordal") {
    return(theoretical_distance_profile_vmf_s1_chordal(omega, mu, kappa, t_values))
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

#' Fast exact distance profile on S^2 under vMF by direct projected-density integration
#'
#' This is mathematically equivalent to `theoretical_distance_profile_vmf()` in the
#' S^2 case, but it hoists all threshold-independent terms out of the integrand and
#' replaces repeated evaluations of the S^1 normalizing constant by `besselI`.
#' The goal is purely computational: preserve the same exact one-dimensional
#' integration while reducing the per-call overhead inside the multiplier bootstrap.
#'
#' @param omega Reference point on S^2
#' @param mu Mean direction in R^3
#' @param kappa Concentration parameter
#' @param t_values Vector of distance thresholds
#' @param distance_type Either "chordal" or "geodesic"
#' @param rel.tol Relative tolerance for integrate()
#' @param abs.tol Absolute tolerance for integrate()
#' @return Vector of probabilities P(d(X, omega) <= t)
theoretical_distance_profile_vmf_s2_fast <- function(omega,
                                                     mu,
                                                     kappa,
                                                     t_values,
                                                     distance_type = "chordal",
                                                     rel.tol = 1e-8,
                                                     abs.tol = 1e-10) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3L || length(mu) != 3L) {
    stop("`omega` and `mu` must both have length 3 for S^2 vMF distance profiles.")
  }

  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  t_values <- as.numeric(t_values)

  rho <- sum(mu * omega)
  rho <- pmin(pmax(rho, -1), 1)
  tangential_norm <- sqrt(pmax(0, 1 - rho^2))
  log_c3 <- rotasym::c_vMF(p = 3, kappa = kappa, log = TRUE)

  vapply(t_values, function(t) {
    threshold <- if (identical(distance_type, "chordal")) {
      1 - (t^2) / 2
    } else {
      cos(t)
    }

    if (!is.finite(threshold) || threshold <= -1) {
      return(1)
    }
    if (threshold >= 1) {
      return(0)
    }

    density_u <- function(u) {
      lambda <- kappa * tangential_norm * sqrt(pmax(0, 1 - u^2))
      log_i0 <- ifelse(
        lambda <= 0,
        0,
        log(besselI(lambda, nu = 0, expon.scaled = TRUE)) + lambda
      )
      exp(log_c3 + kappa * rho * u + log(2 * pi) + log_i0)
    }

    cdf_at_threshold <- integrate(
      f = density_u,
      lower = -1 + 1e-8,
      upper = threshold,
      rel.tol = rel.tol,
      abs.tol = abs.tol
    )$value

    1 - cdf_at_threshold
  }, numeric(1))
}

build_vmf_s2_u_grid <- function(n_u = 4097L) {
  n_u <- as.integer(n_u)
  if (!is.finite(n_u) || n_u < 17L || n_u %% 2L != 1L) {
    stop("`n_u` must be an odd integer >= 17.")
  }

  u_grid <- seq(-1, 1, length.out = n_u)
  list(
    u = u_grid,
    du = u_grid[[2]] - u_grid[[1]],
    sqrt_one_minus_u2 = sqrt(pmax(0, 1 - u_grid^2))
  )
}

vmf_s2_projected_density_matrix <- function(rho,
                                            kappa,
                                            u_grid_object) {
  rho <- as.numeric(rho)
  rho <- pmin(pmax(rho, -1), 1)
  n_u <- length(u_grid_object$u)
  n_rho <- length(rho)

  if (n_rho == 0L) {
    return(matrix(0, nrow = n_u, ncol = 0L))
  }

  if (kappa <= 1e-12) {
    return(matrix(0.5, nrow = n_u, ncol = n_rho))
  }

  sqrt_one_minus_rho2 <- sqrt(pmax(0, 1 - rho^2))
  lambda <- kappa * outer(u_grid_object$sqrt_one_minus_u2, sqrt_one_minus_rho2)
  log_i0 <- matrix(0, nrow = n_u, ncol = n_rho)
  positive_lambda <- lambda > 1e-14
  log_i0[positive_lambda] <- log(besselI(
    lambda[positive_lambda],
    nu = 0,
    expon.scaled = TRUE
  )) + lambda[positive_lambda]

  log_sinh_kappa <- kappa + log1p(-exp(-2 * kappa)) - log(2)
  log_const <- log(kappa) - log(2) - log_sinh_kappa
  log_density <- log_const +
    kappa * outer(u_grid_object$u, rho) +
    log_i0

  exp(log_density)
}

build_vmf_s2_projected_cdf_matrix <- function(rho,
                                              kappa,
                                              n_u = 4097L) {
  u_grid_object <- build_vmf_s2_u_grid(n_u)
  density_matrix <- vmf_s2_projected_density_matrix(
    rho = rho,
    kappa = kappa,
    u_grid_object = u_grid_object
  )
  n_u <- nrow(density_matrix)
  n_rho <- ncol(density_matrix)

  if (n_rho == 0L) {
    return(list(
      u = u_grid_object$u,
      cdf = matrix(0, nrow = n_u, ncol = 0L)
    ))
  }

  increments <- 0.5 * (
    density_matrix[-1, , drop = FALSE] +
      density_matrix[-n_u, , drop = FALSE]
  ) * u_grid_object$du
  cdf_matrix <- rbind(rep(0, n_rho), apply(increments, 2, cumsum))
  cdf_matrix <- sweep(cdf_matrix, 2, cdf_matrix[n_u, ], "/", check.margin = FALSE)
  cdf_matrix[n_u, ] <- 1

  list(
    u = u_grid_object$u,
    cdf = cdf_matrix
  )
}

interpolate_vmf_s2_upper_tail <- function(cdf_values,
                                          u_grid,
                                          a_values) {
  a_values <- pmin(pmax(as.numeric(a_values), -1), 1)
  n_u <- length(u_grid)
  if (n_u < 2L) {
    stop("`u_grid` must have length at least 2.")
  }

  du <- u_grid[[2]] - u_grid[[1]]
  scaled_position <- (a_values - u_grid[[1]]) / du
  left_index <- floor(scaled_position) + 1L
  left_index <- pmin(pmax(left_index, 1L), n_u - 1L)
  lambda <- scaled_position - (left_index - 1L)
  lambda <- pmin(pmax(lambda, 0), 1)

  cdf_left <- cdf_values[left_index]
  cdf_right <- cdf_values[left_index + 1L]
  cdf_at_a <- (1 - lambda) * cdf_left + lambda * cdf_right

  cdf_at_a[a_values <= -1] <- 0
  cdf_at_a[a_values >= 1] <- 1
  1 - cdf_at_a
}

distance_profile_vmf_s2_grid <- function(omega_grid,
                                         mu,
                                         kappa,
                                         t_grid,
                                         distance_type = "geodesic",
                                         n_u = 4097L) {
  omega_grid <- as.matrix(omega_grid)
  mu <- as.numeric(mu)

  if (ncol(omega_grid) != 3L || length(mu) != 3L) {
    stop("`omega_grid` must have 3 columns and `mu` must have length 3.")
  }

  omega_norms <- sqrt(rowSums(omega_grid^2))
  if (any(omega_norms <= 0)) {
    stop("`omega_grid` must contain nonzero rows.")
  }
  omega_grid <- omega_grid / omega_norms
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  t_grid <- as.numeric(t_grid)

  rho <- as.numeric(omega_grid %*% mu)
  cdf_object <- build_vmf_s2_projected_cdf_matrix(
    rho = rho,
    kappa = kappa,
    n_u = n_u
  )
  a_values <- if (identical(distance_type, "geodesic")) {
    cos(t_grid)
  } else {
    1 - (t_grid^2) / 2
  }
  a_values <- pmin(pmax(a_values, -1), 1)

  output <- t(vapply(seq_along(rho), function(i) {
    interpolate_vmf_s2_upper_tail(
      cdf_values = cdf_object$cdf[, i],
      u_grid = cdf_object$u,
      a_values = a_values
    )
  }, numeric(length(t_grid))))

  output[, t_grid <= 0] <- 0
  if (identical(distance_type, "geodesic")) {
    output[, t_grid >= pi] <- 1
  } else {
    output[, t_grid >= 2] <- 1
  }

  output
}

distance_profile_vmf_s2_cvm_grid <- function(X,
                                             mu,
                                             kappa,
                                             n_u = 4097L) {
  X <- as.matrix(X)
  mu <- as.numeric(mu)

  if (ncol(X) != 3L || length(mu) != 3L) {
    stop("`X` must have 3 columns and `mu` must have length 3.")
  }

  X_norms <- sqrt(rowSums(X^2))
  if (any(X_norms <= 0)) {
    stop("`X` must contain nonzero rows.")
  }
  X <- X / X_norms
  mu <- mu / sqrt(sum(mu^2))

  rho <- as.numeric(X %*% mu)
  a_matrix <- X %*% t(X)
  a_matrix <- pmin(pmax(a_matrix, -1), 1)

  cdf_object <- build_vmf_s2_projected_cdf_matrix(
    rho = rho,
    kappa = kappa,
    n_u = n_u
  )

  t(vapply(seq_along(rho), function(i) {
    interpolate_vmf_s2_upper_tail(
      cdf_values = cdf_object$cdf[, i],
      u_grid = cdf_object$u,
      a_values = a_matrix[i, ]
    )
  }, numeric(nrow(X))))
}

validate_vmf_s2_grid_profile <- function(n_checks = 200,
                                         kappa_values = c(0.5, 2, 5),
                                         n_u = 4097L,
                                         distance_type = "geodesic",
                                         seed = 1) {
  set.seed(seed)
  errors <- numeric(0)

  for (kappa in kappa_values) {
    mu <- stats::rnorm(3)
    mu <- mu / sqrt(sum(mu^2))
    omega_grid <- matrix(stats::rnorm(n_checks * 3), ncol = 3)
    omega_grid <- omega_grid / sqrt(rowSums(omega_grid^2))
    t_max <- if (identical(distance_type, "geodesic")) pi else 2
    t_grid <- stats::runif(n_checks, 0, t_max)

    fast_matrix <- distance_profile_vmf_s2_grid(
      omega_grid = omega_grid,
      mu = mu,
      kappa = kappa,
      t_grid = t_grid,
      distance_type = distance_type,
      n_u = n_u
    )
    fast_values <- diag(fast_matrix)
    exact_values <- vapply(seq_len(n_checks), function(i) {
      distance_profile_vmf_s2_integral(
        omega = omega_grid[i, ],
        mu = mu,
        kappa = kappa,
        t_values = t_grid[i],
        distance_type = distance_type
      )
    }, numeric(1))

    errors <- c(errors, abs(fast_values - exact_values))
  }

  c(
    max_error = max(errors),
    mean_error = mean(errors),
    q95_error = unname(stats::quantile(errors, probs = 0.95, names = FALSE, type = 8))
  )
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

# -------------------------------------------------------------------------
# BAHADUR: Single trajectory analysis and helpers
# -------------------------------------------------------------------------
#' Single trajectory analysis
analyze_single_trajectory <- function(mu_true, kappa_true, sample_sizes, trajectory_id) {
  q <- length(mu_true) - 1
  xi_true <- kappa_true * mu_true
  
  results <- data.frame()
  
  for (n in sample_sizes) {
    # Draw an independent sample for each n (no cumulative sample growth)
    current_sample <- r_vMF(n, mu_true, kappa_true)
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
