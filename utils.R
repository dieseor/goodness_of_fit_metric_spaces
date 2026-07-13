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
  hvmf_eval_projection_cdf_tabulated_spline(y = y, cdf_table = cdf_table)
}

hvmf_prepare_projection_cdf_interpolator <- function(cdf_table) {
  if (!is.data.frame(cdf_table) || !all(c("y", "cdf") %in% names(cdf_table))) {
    stop("`cdf_table` must be a data frame with columns `y` and `cdf`.")
  }

  y_grid <- as.numeric(cdf_table$y)
  cdf_grid <- as.numeric(cdf_table$cdf)
  if (length(y_grid) == 0L) {
    stop("`cdf_table` cannot be empty.")
  }

  if (length(y_grid) == 1L) {
    return(list(
      type = "constant",
      x_min = y_grid[[1L]],
      x_max = y_grid[[1L]],
      y_min = cdf_grid[[1L]],
      y_max = cdf_grid[[1L]]
    ))
  }

  list(
    type = "spline",
    x_min = y_grid[[1L]],
    x_max = y_grid[[length(y_grid)]],
    y_min = min(cdf_grid),
    y_max = max(cdf_grid),
    fn = stats::splinefun(
      x = y_grid,
      y = cdf_grid,
      method = "monoH.FC"
    )
  )
}

hvmf_eval_projection_cdf_tabulated_linear <- function(y, cdf_table) {
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

hvmf_eval_projection_cdf_tabulated_spline <- function(y,
                                                      cdf_table,
                                                      interpolator = NULL) {
  y <- as.numeric(y)
  interpolator <- interpolator %||% hvmf_prepare_projection_cdf_interpolator(cdf_table)

  if (identical(interpolator$type, "constant")) {
    output <- rep(interpolator$y_min, length(y))
    output[y <= interpolator$x_min] <- 0
    return(pmin(pmax(output, 0), 1))
  }

  y_eval <- pmin(pmax(y, interpolator$x_min), interpolator$x_max)
  output <- interpolator$fn(y_eval)
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
    cdf_interpolator <- hvmf_prepare_projection_cdf_interpolator(cdf_table)
    output[i, ] <- hvmf_eval_projection_cdf_tabulated_spline(
      y = row_y,
      cdf_table = cdf_table,
      interpolator = cdf_interpolator
    )
  }

  output
}

hvmf_cvm_profile_matrix_tabulated_linear <- function(data,
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
    output[i, ] <- hvmf_eval_projection_cdf_tabulated_linear(row_y, cdf_table)
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
  input_is_matrix <- is.matrix(t)
  t_dims <- if (input_is_matrix) dim(t) else NULL
  t <- as.numeric(t)
  if (distance_type == "chordal") {
    threshold <- 1 - (t^2) / 2
  } else {
    threshold <- cos(t)
  }
  threshold <- pmin(pmax(threshold, -1), 1)
  if (isTRUE(input_is_matrix)) {
    matrix(threshold, nrow = t_dims[[1L]], ncol = t_dims[[2L]])
  } else {
    threshold
  }
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

#' Validate spherical Cauchy parameters on S^2
#' @param mu Mean direction on S^2
#' @param rho Concentration parameter in [0, 1)
#' @return Normalized parameter list
spherical_cauchy_validate_parameters <- function(mu, rho) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Spherical Cauchy utilities currently support only S^2 (ambient dimension 3).")
  }

  rho <- as.numeric(rho)
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) {
    stop("`rho` must be a finite scalar in [0, 1).")
  }

  list(
    mu = mu,
    rho = rho,
    phi = rho * mu,
    q = 2L,
    ambient_dim = 3L
  )
}

spherical_cauchy_projected_cdf_from_inner_products <- function(r,
                                                               x,
                                                               rho,
                                                               tail_tol = 1e-10,
                                                               l_max = NULL,
                                                               max_l_max = NULL,
                                                               warn = TRUE) {
  r <- as.numeric(r)
  if (length(r) == 0L || any(!is.finite(r))) {
    stop("`r` must be a non-empty finite vector.")
  }
  r <- pmin(pmax(r, -1), 1)

  x_is_matrix <- is.matrix(x)
  if (x_is_matrix) {
    x_mat <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
    if (nrow(x_mat) != length(r)) {
      stop("When `x` is a matrix, `nrow(x)` must equal `length(r)`.")
    }
  } else {
    x_vec <- as.numeric(x)
    x_mat <- matrix(
      rep(x_vec, each = length(r)),
      nrow = length(r),
      ncol = length(x_vec)
    )
  }

  out <- matrix(0, nrow = nrow(x_mat), ncol = ncol(x_mat))
  x_clamped <- pmin(pmax(x_mat, -1), 1)
  out[x_mat <= -1] <- 0
  out[x_mat >= 1] <- 1

  interior_mask <- x_mat > -1 & x_mat < 1
  if (!any(interior_mask)) {
    return(out)
  }

  if (rho <= 1e-14) {
    out[interior_mask] <- ((x_clamped + 1) / 2)[interior_mask]
    return(pmin(pmax(out, 0), 1))
  }

  exact_axis_rows <- abs(r - 1) <= 1e-12
  if (any(exact_axis_rows)) {
    out[exact_axis_rows, ] <- matrix(
      spherical_cauchy_axis_projected_cdf(
        x_clamped[exact_axis_rows, , drop = FALSE],
        rho = rho
      ),
      nrow = sum(exact_axis_rows),
      ncol = ncol(x_clamped)
    )
  }

  general_rows <- which(!exact_axis_rows)
  if (length(general_rows) == 0L) {
    return(pmin(pmax(out, 0), 1))
  }

  if (is.null(l_max)) {
    l_max <- if (rho > 0.95) 8192L else 2048L
  }
  if (is.null(max_l_max)) {
    max_l_max <- max(l_max, if (rho > 0.95) 32768L else 8192L)
  }
  l_max <- as.integer(l_max)
  max_l_max <- as.integer(max_l_max)
  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 1L) {
    stop("`l_max` must be a positive integer.")
  }
  if (length(max_l_max) != 1L || !is.finite(max_l_max) || max_l_max < l_max) {
    stop("`max_l_max` must be an integer >= `l_max`.")
  }

  x_general <- x_clamped[general_rows, , drop = FALSE]
  cdf_general <- (x_general + 1) / 2
  r_general <- r[general_rows]
  p_r_prev <- rep(1, length(general_rows))
  p_r_curr <- r_general
  p_x_prev <- matrix(1, nrow = nrow(x_general), ncol = ncol(x_general))
  p_x_curr <- x_general
  rho_power <- rho
  ell <- 1L
  current_cap <- l_max

  repeat {
    while (ell <= current_cap) {
      p_r_next <- ((2 * ell + 1) * r_general * p_r_curr - ell * p_r_prev) / (ell + 1)
      p_x_next <- ((2 * ell + 1) * x_general * p_x_curr - ell * p_x_prev) / (ell + 1)

      cdf_general <- cdf_general + 0.5 * rho_power * sweep(
        p_x_next - p_x_prev,
        1L,
        p_r_curr,
        "*"
      )
      tail_bound <- rho_power * rho / (1 - rho)
      if (tail_bound <= tail_tol) {
        break
      }

      p_r_prev <- p_r_curr
      p_r_curr <- p_r_next
      p_x_prev <- p_x_curr
      p_x_curr <- p_x_next
      rho_power <- rho_power * rho
      ell <- ell + 1L
    }

    if (tail_bound <= tail_tol) {
      break
    }

    if (current_cap >= max_l_max) {
      if (isTRUE(warn)) {
        warning(sprintf(
          "Spherical Cauchy Legendre truncation reached L_max = %d at rho = %.4f with tail bound %.3e.",
          current_cap,
          rho,
          tail_bound
        ))
      }
      break
    }

    next_cap <- min(max_l_max, max(current_cap + 1L, 2L * current_cap))
    if (isTRUE(warn) && rho > 0.95 && next_cap > current_cap) {
      warning(sprintf(
        "rho = %.4f is close to 1; increasing spherical Cauchy Legendre L_max from %d to %d.",
        rho,
        current_cap,
        next_cap
      ))
    }
    current_cap <- next_cap
  }

  out[general_rows, ] <- cdf_general
  out[interior_mask] <- pmin(pmax(out[interior_mask], 0), 1)
  out
}

#' Normalize spherical Cauchy theta supplied as {mu, rho} or {phi}
#' @param theta Parameter list
#' @param ambient_dim Expected ambient dimension
#' @return Normalized parameter list with mu, rho, phi, and q
spherical_cauchy_normalize_theta <- function(theta, ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Spherical Cauchy theta must be a list containing either `phi` or (`mu`, `rho`).")
  }

  if (!is.null(theta$phi)) {
    phi <- as.numeric(theta$phi)
    if (length(phi) != ambient_dim || any(!is.finite(phi))) {
      stop("`theta$phi` must be a finite vector with the expected ambient dimension.")
    }
    rho <- sqrt(sum(phi^2))
    if (!is.finite(rho) || rho >= 1) {
      stop("`theta$phi` must satisfy ||phi|| < 1.")
    }
    mu <- if (rho <= 1e-14) c(0, 0, 1) else phi / rho
    return(list(mu = mu, rho = rho, phi = phi, q = ambient_dim - 1L, ambient_dim = ambient_dim))
  }

  spherical_cauchy_validate_parameters(mu = theta$mu, rho = theta$rho)
}

#' Lower projected CDF on the axis omega = mu for spherical Cauchy on S^2
#' @param x Scalar or vector in [-1, 1]
#' @param rho Concentration parameter in [0, 1)
#' @return P(mu^T X <= x)
spherical_cauchy_axis_projected_cdf <- function(x, rho) {
  rho <- as.numeric(rho)
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) {
    stop("`rho` must be a finite scalar in [0, 1).")
  }

  x <- pmin(pmax(as.numeric(x), -1), 1)
  if (rho <= 1e-14) {
    return((x + 1) / 2)
  }

  coeff <- (1 - rho^2) / (2 * rho)
  coeff * (
    1 / sqrt(1 + rho^2 - 2 * rho * x) -
      1 / (1 + rho)
  )
}

#' Axis projected density for spherical Cauchy on S^2
#' @param x Scalar or vector in [-1, 1]
#' @param rho Concentration parameter in [0, 1)
#' @param log Whether to return the log-density
#' @return Density values
spherical_cauchy_axis_projected_density <- function(x, rho, log = FALSE) {
  rho <- as.numeric(rho)
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) {
    stop("`rho` must be a finite scalar in [0, 1).")
  }

  x <- as.numeric(x)
  out <- rep(if (log) -Inf else 0, length(x))
  valid <- x >= -1 & x <= 1
  if (!any(valid)) {
    return(out)
  }

  log_density <- log1p(-rho^2) - log(2) - 1.5 * log(1 + rho^2 - 2 * rho * x[valid])
  out[valid] <- if (log) log_density else exp(log_density)
  out
}

#' Log-density on S^2 under the spherical Cauchy / Poisson-kernel law
#' @param x Matrix of observations on S^2
#' @param mu Mean direction on S^2
#' @param rho Concentration parameter in [0, 1)
#' @param phi Optional Euclidean parameter with ||phi|| < 1
#' @param log Whether to return log-density
#' @return Density or log-density values
d_spherical_cauchy_s2 <- function(x, mu = NULL, rho = NULL, phi = NULL, log = FALSE) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Spherical Cauchy density currently supports only S^2 data.")
  }

  theta <- if (!is.null(phi)) {
    spherical_cauchy_normalize_theta(list(phi = phi), ambient_dim = 3L)
  } else {
    spherical_cauchy_validate_parameters(mu = mu, rho = rho)
  }

  phi_vec <- theta$phi
  rho_sq <- sum(phi_vec^2)
  denom <- as.numeric(1 - 2 * (x %*% phi_vec) + rho_sq)
  if (any(denom <= 0)) {
    stop("Encountered a non-positive denominator in the spherical Cauchy density.")
  }

  log_density <- log1p(-rho_sq) - log(4 * pi) - 1.5 * log(denom)
  if (log) log_density else exp(log_density)
}

#' Lower projected CDF for spherical Cauchy on S^2 via Legendre expansion
#' @param x Scalar or vector in [-1, 1]
#' @param omega Reference point on S^2
#' @param mu Mean direction on S^2
#' @param rho Concentration parameter in [0, 1)
#' @param tail_tol Geometric tail tolerance for truncation
#' @param l_max Initial truncation cap
#' @param max_l_max Maximum truncation cap after adaptive growth
#' @param warn Whether to emit warnings on aggressive truncation growth
#' @return Lower CDF values P(omega^T X <= x)
spherical_cauchy_projected_cdf <- function(x,
                                           omega,
                                           mu,
                                           rho,
                                           tail_tol = 1e-10,
                                           l_max = NULL,
                                           max_l_max = NULL,
                                           warn = TRUE) {
  params <- spherical_cauchy_validate_parameters(mu = mu, rho = rho)
  x <- as.numeric(x)

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (ncol(omega) != 3L) {
      stop("`omega` must have three columns for spherical Cauchy profiles on S^2.")
    }
    if (length(x) == 1L) {
      x <- rep(x, nrow(omega))
    }
    if (length(x) != nrow(omega)) {
      stop("When `omega` is a matrix, `x` must have length 1 or nrow(omega).")
    }

    return(vapply(seq_len(nrow(omega)), function(i) {
      spherical_cauchy_projected_cdf(
        x = x[i],
        omega = omega[i, ],
        mu = params$mu,
        rho = params$rho,
        tail_tol = tail_tol,
        l_max = l_max,
        max_l_max = max_l_max,
        warn = warn
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  if (length(omega) != 3L) {
    stop("`omega` must have length 3 for spherical Cauchy profiles on S^2.")
  }

  rho <- params$rho
  x_clamped <- pmin(pmax(x, -1), 1)
  out <- numeric(length(x_clamped))
  out[x <= -1] <- 0
  out[x >= 1] <- 1
  active <- which(x > -1 & x < 1)
  if (length(active) == 0L) {
    return(out)
  }

  if (rho <= 1e-14) {
    out[active] <- (x_clamped[active] + 1) / 2
    return(pmin(pmax(out, 0), 1))
  }

  r <- pmin(pmax(sum(params$mu * omega), -1), 1)
  out[active] <- drop(spherical_cauchy_projected_cdf_from_inner_products(
    r = r,
    x = x_clamped[active],
    rho = rho,
    tail_tol = tail_tol,
    l_max = l_max,
    max_l_max = max_l_max,
    warn = warn
  ))
  pmin(pmax(out, 0), 1)
}

#' Theoretical distance profile for spherical Cauchy on S^2
#' @param omega Reference point on S^2
#' @param t_values Distance thresholds
#' @param mu Mean direction on S^2
#' @param rho Concentration parameter in [0, 1)
#' @param distance_type Either "geodesic" or "chordal"
#' @param tail_tol Geometric tail tolerance for Legendre truncation
#' @param l_max Initial truncation cap
#' @param max_l_max Maximum truncation cap after adaptive growth
#' @param warn Whether to emit warnings when the truncation budget is stressed
#' @return Vector of probabilities P(d(X, omega) <= t)
distance_profile_spherical_cauchy <- function(omega,
                                              t_values,
                                              mu,
                                              rho,
                                              distance_type = c("geodesic", "chordal"),
                                              tail_tol = 1e-10,
                                              l_max = NULL,
                                              max_l_max = NULL,
                                              warn = TRUE) {
  distance_type <- match.arg(distance_type)
  params <- spherical_cauchy_validate_parameters(mu = mu, rho = rho)
  t_values <- as.numeric(t_values)

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (ncol(omega) != 3L) {
      stop("`omega` must have three columns for spherical Cauchy profiles on S^2.")
    }
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }

    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_spherical_cauchy(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = params$mu,
        rho = params$rho,
        distance_type = distance_type,
        tail_tol = tail_tol,
        l_max = l_max,
        max_l_max = max_l_max,
        warn = warn
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  out <- numeric(length(t_values))
  if (identical(distance_type, "geodesic")) {
    out[t_values <= 0] <- 0
    out[t_values >= pi] <- 1
    active <- which(t_values > 0 & t_values < pi)
  } else {
    out[t_values <= 0] <- 0
    out[t_values >= 2] <- 1
    active <- which(t_values > 0 & t_values < 2)
  }

  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  out[active] <- 1 - spherical_cauchy_projected_cdf(
    x = thresholds,
    omega = omega,
    mu = params$mu,
    rho = params$rho,
    tail_tol = tail_tol,
    l_max = l_max,
    max_l_max = max_l_max,
    warn = warn
  )

  pmin(pmax(out, 0), 1)
}

#' Profile matrix on S^2 for spherical Cauchy
#' @param omega_grid Matrix of reference points on S^2
#' @param mu Mean direction on S^2
#' @param rho Concentration parameter in [0, 1)
#' @param t_grid Distance thresholds
#' @param distance_type Either "geodesic" or "chordal"
#' @param tail_tol Geometric tail tolerance for Legendre truncation
#' @param l_max Initial truncation cap
#' @param max_l_max Maximum truncation cap after adaptive growth
#' @param warn Whether to emit warnings when the truncation budget is stressed
#' @return Matrix with one row per omega and one column per t
distance_profile_spherical_cauchy_grid <- function(omega_grid,
                                                   mu,
                                                   rho,
                                                   t_grid,
                                                   distance_type = c("geodesic", "chordal"),
                                                   tail_tol = 1e-10,
                                                   l_max = NULL,
                                                   max_l_max = NULL,
                                                   warn = TRUE) {
  distance_type <- match.arg(distance_type)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  t_grid <- as.numeric(t_grid)
  params <- spherical_cauchy_validate_parameters(mu = mu, rho = rho)

  out <- matrix(0, nrow = nrow(omega_grid), ncol = length(t_grid))
  if (identical(distance_type, "geodesic")) {
    out[, t_grid <= 0] <- 0
    out[, t_grid >= pi] <- 1
    active <- which(t_grid > 0 & t_grid < pi)
  } else {
    out[, t_grid <= 0] <- 0
    out[, t_grid >= 2] <- 1
    active <- which(t_grid > 0 & t_grid < 2)
  }

  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_grid[active], distance_type = distance_type)
  r_vec <- as.numeric(omega_grid %*% params$mu)
  out[, active] <- 1 - spherical_cauchy_projected_cdf_from_inner_products(
    r = r_vec,
    x = thresholds,
    rho = params$rho,
    tail_tol = tail_tol,
    l_max = l_max,
    max_l_max = max_l_max,
    warn = warn
  )

  pmin(pmax(out, 0), 1)
}

distance_profile_spherical_cauchy_cvm_grid <- function(data,
                                                       mu,
                                                       rho,
                                                       distance_matrix,
                                                       distance_type = c("geodesic", "chordal"),
                                                       tail_tol = 1e-10,
                                                       l_max = NULL,
                                                       max_l_max = NULL,
                                                       warn = TRUE) {
  distance_type <- match.arg(distance_type)
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Spherical Cauchy CvM profile matrix currently supports only S^2 data.")
  }

  distance_matrix <- as.matrix(distance_matrix)
  if (nrow(distance_matrix) != nrow(x) || ncol(distance_matrix) != nrow(x)) {
    stop("`distance_matrix` must be a square matrix with `nrow(data)` rows.")
  }

  params <- spherical_cauchy_validate_parameters(mu = mu, rho = rho)
  thresholds <- sphere_distance_to_dot_threshold(distance_matrix, distance_type = distance_type)
  r_vec <- as.numeric(x %*% params$mu)

  out <- 1 - spherical_cauchy_projected_cdf_from_inner_products(
    r = r_vec,
    x = thresholds,
    rho = params$rho,
    tail_tol = tail_tol,
    l_max = l_max,
    max_l_max = max_l_max,
    warn = warn
  )

  pmin(pmax(out, 0), 1)
}

#' Convert unconstrained Euclidean parameter u into phi with ||phi|| < 1
#' @param u Unconstrained vector in R^3
#' @return Phi inside the open unit ball
spherical_cauchy_phi_from_u <- function(u) {
  u <- as.numeric(u)
  if (length(u) != 3L || any(!is.finite(u))) {
    stop("`u` must be a finite vector of length 3.")
  }

  u / sqrt(1 + sum(u^2))
}

#' Convert spherical Cauchy phi into unconstrained Euclidean parameter u
#' @param phi Euclidean parameter with ||phi|| < 1
#' @return Unconstrained vector u
spherical_cauchy_u_from_phi <- function(phi) {
  phi <- as.numeric(phi)
  if (length(phi) != 3L || any(!is.finite(phi))) {
    stop("`phi` must be a finite vector of length 3.")
  }
  phi_norm_sq <- sum(phi^2)
  if (!is.finite(phi_norm_sq) || phi_norm_sq >= 1) {
    stop("`phi` must satisfy ||phi|| < 1.")
  }

  phi / sqrt(1 - phi_norm_sq)
}

#' Weighted log-likelihood for spherical Cauchy on S^2 in phi parametrization
#' @param phi Euclidean parameter with ||phi|| < 1
#' @param x Matrix of observations on S^2
#' @param prob_weights Probability weights summing to one
#' @param include_constant Whether to include the additive -log(4*pi) constant
#' @return Weighted average log-likelihood
spherical_cauchy_weighted_loglik_phi <- function(phi,
                                                 x,
                                                 prob_weights,
                                                 include_constant = TRUE) {
  phi <- as.numeric(phi)
  if (length(phi) != 3L || any(!is.finite(phi))) {
    stop("`phi` must be a finite vector of length 3.")
  }
  phi_norm_sq <- sum(phi^2)
  if (!is.finite(phi_norm_sq) || phi_norm_sq >= 1) {
    return(-Inf)
  }

  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- jp_normalize_probability_weights(prob_weights, nrow(x))
  denom <- as.numeric(1 - 2 * (x %*% phi) + phi_norm_sq)
  if (any(!is.finite(denom)) || any(denom <= 0)) {
    return(-Inf)
  }

  constant_term <- if (isTRUE(include_constant)) -log(4 * pi) else 0
  log1p(-phi_norm_sq) + constant_term - 1.5 * sum(prob_weights * log(denom))
}

#' Gradient of the weighted log-likelihood for spherical Cauchy on S^2 in phi
#' @param phi Euclidean parameter with ||phi|| < 1
#' @param x Matrix of observations on S^2
#' @param prob_weights Probability weights summing to one
#' @return Gradient vector with respect to phi
spherical_cauchy_weighted_loglik_phi_grad <- function(phi, x, prob_weights) {
  phi <- as.numeric(phi)
  phi_norm_sq <- sum(phi^2)
  if (length(phi) != 3L || any(!is.finite(phi)) || !is.finite(phi_norm_sq) || phi_norm_sq >= 1) {
    return(rep(NaN, 3L))
  }

  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- jp_normalize_probability_weights(prob_weights, nrow(x))
  denom <- as.numeric(1 - 2 * (x %*% phi) + phi_norm_sq)
  if (any(!is.finite(denom)) || any(denom <= 0)) {
    return(rep(NaN, 3L))
  }

  -2 * phi / (1 - phi_norm_sq) +
    3 * colSums(x * (prob_weights / denom)) -
    3 * phi * sum(prob_weights / denom)
}

spherical_cauchy_negloglik_u <- function(u, x, prob_weights) {
  phi <- spherical_cauchy_phi_from_u(u)
  -spherical_cauchy_weighted_loglik_phi(phi = phi, x = x, prob_weights = prob_weights)
}

spherical_cauchy_negloglik_u_grad <- function(u, x, prob_weights) {
  phi <- spherical_cauchy_phi_from_u(u)
  grad_phi <- spherical_cauchy_weighted_loglik_phi_grad(phi = phi, x = x, prob_weights = prob_weights)
  scale <- sqrt(1 + sum(u^2))
  jacobian <- diag(1 / scale, nrow = 3L, ncol = 3L) - outer(u, u) / (scale^3)
  drop(-crossprod(jacobian, grad_phi))
}

#' Weighted spherical Cauchy MLE on S^2
#' @param data Matrix of observations on S^2
#' @param weights Optional nonnegative weights
#' @param control Optimizer and warm-start controls
#' @return Estimated theta with mu, rho, phi, q, ll, opt, weighted_mle
spherical_cauchy_mle_s2_weighted <- function(data, weights = NULL, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("`spherical_cauchy_mle_s2_weighted()` currently supports only S^2 data.")
  }

  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  theta_start <- control$spherical_cauchy_start_theta %||%
    control$cauchy_start_theta %||%
    control$theta_start %||%
    NULL
  maxit <- as.integer(control$spherical_cauchy_maxit %||% 500L)
  reltol <- as.numeric(control$spherical_cauchy_reltol %||% 1e-10)
  method <- as.character(control$spherical_cauchy_optim_method %||% "BFGS")
  use_gradient <- isTRUE(control$spherical_cauchy_use_gradient %||% TRUE)

  if (length(maxit) != 1L || !is.finite(maxit) || maxit < 1L) {
    stop("`control$spherical_cauchy_maxit` must be a positive integer.")
  }
  if (length(reltol) != 1L || !is.finite(reltol) || reltol <= 0) {
    stop("`control$spherical_cauchy_reltol` must be a positive finite scalar.")
  }

  resultant <- colSums(x * prob_weights)
  resultant_norm <- sqrt(sum(resultant^2))
  mu0 <- if (resultant_norm <= 1e-12) c(0, 0, 1) else resultant / resultant_norm
  rho0 <- min(max(resultant_norm, 1e-4), 0.95)
  phi0 <- rho0 * mu0

  if (!is.null(theta_start)) {
    theta_start <- spherical_cauchy_normalize_theta(theta_start, ambient_dim = 3L)
    phi0 <- theta_start$phi
  }

  u0 <- spherical_cauchy_u_from_phi(phi0)
  optim_control <- list(maxit = maxit, reltol = reltol)

  fit <- stats::optim(
    par = u0,
    fn = spherical_cauchy_negloglik_u,
    gr = if (isTRUE(use_gradient)) spherical_cauchy_negloglik_u_grad else NULL,
    x = x,
    prob_weights = prob_weights,
    method = method,
    control = optim_control
  )

  phi_hat <- spherical_cauchy_phi_from_u(fit$par)
  rho_hat <- sqrt(sum(phi_hat^2))
  mu_hat <- if (rho_hat <= 1e-12) c(0, 0, 1) else phi_hat / rho_hat
  ll_hat <- spherical_cauchy_weighted_loglik_phi(
    phi = phi_hat,
    x = x,
    prob_weights = prob_weights,
    include_constant = TRUE
  )

  list(
    mu = mu_hat,
    rho = rho_hat,
    phi = phi_hat,
    q = 2L,
    ll = ll_hat,
    opt = fit,
    weighted_mle = !is.null(weights)
  )
}

#' Sample from the spherical Cauchy / Poisson-kernel law on S^2
#' @param n Sample size
#' @param mu Mean direction on S^2
#' @param rho Concentration parameter in [0, 1)
#' @param check Whether to validate output norms
#' @return n x 3 matrix with rows on S^2
r_sph_spherical_cauchy <- function(n, mu, rho, check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  params <- spherical_cauchy_validate_parameters(mu = mu, rho = rho)
  if (params$rho <= 1e-14) {
    return(jp_uniform_sphere(n = n, ambient_dim = 3L))
  }

  u_samples <- stats::runif(n)
  y <- 1 / (1 + params$rho) + 2 * params$rho * u_samples / (1 - params$rho^2)
  t_samples <- (1 + params$rho^2 - y^(-2)) / (2 * params$rho)
  t_samples <- pmin(pmax(t_samples, -1), 1)

  angles <- stats::runif(n, min = 0, max = 2 * pi)
  xi <- cbind(cos(angles), sin(angles))
  b_mu <- jp_orthonormal_complement(params$mu)
  tangential <- xi %*% t(b_mu)
  sample <- t_samples * matrix(params$mu, nrow = n, ncol = 3L, byrow = TRUE) +
    sqrt(pmax(0, 1 - t_samples^2)) * tangential

  if (isTRUE(check)) {
    sample <- sample / sqrt(rowSums(sample^2))
  }

  sample
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

vmf_s2_legendre_cache <- new.env(parent = emptyenv())

vmf_s2_legendre_coefficients <- function(kappa, l_max) {
  kappa <- as.numeric(kappa)
  l_max <- as.integer(l_max)

  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite nonnegative scalar.")
  }
  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  key <- paste0(formatC(kappa, digits = 17, format = "fg"), "::", l_max)
  if (exists(key, envir = vmf_s2_legendre_cache, inherits = FALSE)) {
    return(get(key, envir = vmf_s2_legendre_cache, inherits = FALSE))
  }

  coeffs <- numeric(l_max + 1L)
  coeffs[[1L]] <- 1
  if (l_max > 0L && kappa > 0) {
    ell <- 0:l_max
    ive <- besselI(kappa, nu = ell + 0.5, expon.scaled = TRUE)
    i_ell <- sqrt(pi / (2 * kappa)) * exp(kappa) * ive
    coeffs <- (kappa / sinh(kappa)) * (2 * ell + 1) * i_ell
    coeffs[[1L]] <- 1
  }

  assign(key, coeffs, envir = vmf_s2_legendre_cache)
  coeffs
}

select_vmf_s2_legendre_l_max <- function(kappa,
                                         tail_tol = 1e-10,
                                         min_l = 10L,
                                         max_l = 200L,
                                         tail_length = 5L) {
  kappa <- as.numeric(kappa)
  tail_tol <- as.numeric(tail_tol)
  min_l <- as.integer(min_l)
  max_l <- as.integer(max_l)
  tail_length <- as.integer(tail_length)

  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite nonnegative scalar.")
  }
  if (length(tail_tol) != 1L || !is.finite(tail_tol) || tail_tol <= 0) {
    stop("`tail_tol` must be a strictly positive scalar.")
  }
  if (length(min_l) != 1L || !is.finite(min_l) || min_l < 0L) {
    stop("`min_l` must be a nonnegative integer.")
  }
  if (length(max_l) != 1L || !is.finite(max_l) || max_l < min_l) {
    stop("`max_l` must be an integer not smaller than `min_l`.")
  }
  if (length(tail_length) != 1L || !is.finite(tail_length) || tail_length < 1L) {
    stop("`tail_length` must be a strictly positive integer.")
  }

  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = max_l + tail_length)
  abs_coeffs <- abs(coeffs)

  for (l_max_candidate in min_l:max_l) {
    tail_idx <- (l_max_candidate + 2L):(l_max_candidate + tail_length + 1L)
    if (max(abs_coeffs[tail_idx]) <= tail_tol) {
      return(l_max_candidate)
    }
  }

  max_l
}

distance_profile_vmf_s2_legendre <- function(omega,
                                             mu,
                                             kappa,
                                             t_values,
                                             distance_type = "geodesic",
                                             l_max = NULL,
                                             tail_tol = 1e-10) {
  omega <- as.numeric(omega)
  mu <- as.numeric(mu)
  if (length(omega) != 3L || length(mu) != 3L) {
    stop("`omega` and `mu` must both have length 3 for S^2 Legendre vMF distance profiles.")
  }

  omega <- omega / sqrt(sum(omega^2))
  mu <- mu / sqrt(sum(mu^2))
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  if (is.null(l_max)) {
    l_max <- select_vmf_s2_legendre_l_max(kappa = kappa, tail_tol = tail_tol)
  }
  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = l_max)

  rotational_profile_legendre(
    t = t_values,
    omega = omega,
    mu = mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )
}

distance_profile_vmf_s2_legendre_grid <- function(omega_grid,
                                                  mu,
                                                  kappa,
                                                  t_grid,
                                                  distance_type = "geodesic",
                                                  l_max = NULL,
                                                  tail_tol = 1e-10) {
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
  if (is.null(l_max)) {
    l_max <- select_vmf_s2_legendre_l_max(kappa = kappa, tail_tol = tail_tol)
  }
  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = l_max)

  rotational_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )
}

distance_profile_vmf_s2_legendre_cvm_grid <- function(X,
                                                      mu,
                                                      kappa,
                                                      l_max = NULL,
                                                      tail_tol = 1e-10) {
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
  if (is.null(l_max)) {
    l_max <- select_vmf_s2_legendre_l_max(kappa = kappa, tail_tol = tail_tol)
  }
  coeffs <- vmf_s2_legendre_coefficients(kappa = kappa, l_max = l_max)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  r_values <- as.numeric(X %*% mu)

  1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )
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
## Jones-Pewsey distribution on S^q

jp_cache_env <- new.env(parent = emptyenv())

#' Surface area of the q-dimensional unit sphere S^q
#' @param q Sphere dimension
#' @return Surface area of S^q
sphere_surface_area <- function(q) {
  q <- as.integer(q)
  if (length(q) != 1L || !is.finite(q) || q < 0L) {
    stop("`q` must be a nonnegative integer.")
  }

  2 * pi^((q + 1) / 2) / gamma((q + 1) / 2)
}

jp_normalize_unit_vector <- function(x, arg_name = "`x`", min_length = 3L) {
  x <- as.numeric(x)
  if (length(x) < min_length) {
    stop(sprintf("%s must have length at least %d.", arg_name, min_length))
  }
  if (any(!is.finite(x))) {
    stop(sprintf("%s must be finite.", arg_name))
  }

  x_norm <- sqrt(sum(x^2))
  if (!is.finite(x_norm) || x_norm <= 0) {
    stop(sprintf("%s must have strictly positive norm.", arg_name))
  }

  x / x_norm
}

jp_normalize_unit_matrix <- function(x, arg_name = "`x`", min_ncol = 3L) {
  if (is.vector(x)) {
    x <- matrix(as.numeric(x), nrow = 1L)
  } else {
    x <- as.matrix(x)
  }

  if (nrow(x) == 0L || ncol(x) < min_ncol) {
    stop(sprintf("%s must be a non-empty matrix with at least %d columns.", arg_name, min_ncol))
  }
  if (any(!is.finite(x))) {
    stop(sprintf("%s must be finite.", arg_name))
  }

  norms <- sqrt(rowSums(x^2))
  if (any(!is.finite(norms)) || any(norms <= 0)) {
    stop(sprintf("%s must contain rows with strictly positive norm.", arg_name))
  }

  x / norms
}

jp_validate_parameters <- function(mu, kappa, psi) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  q <- length(mu) - 1L
  if (q < 2L) {
    stop("Jones-Pewsey functions currently support only q >= 2. The case q = 1 is not implemented yet.")
  }

  kappa <- as.numeric(kappa)
  psi <- as.numeric(psi)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a nonnegative finite scalar.")
  }
  if (length(psi) != 1L || !is.finite(psi)) {
    stop("`psi` must be a finite scalar.")
  }

  alpha <- if (psi == 0) 0 else tanh(kappa * psi)
  beta <- if (psi == 0) NA_real_ else 1 / psi

  list(
    mu = mu,
    q = q,
    ambient_dim = length(mu),
    kappa = kappa,
    psi = psi,
    alpha = alpha,
    beta = beta
  )
}

jp_vmf_near_zero_abs_kappa_psi_default <- 1e-3

jp_is_near_zero_vmf_s2 <- function(ambient_dim,
                                   kappa,
                                   psi,
                                   abs_kappa_psi_tol = jp_vmf_near_zero_abs_kappa_psi_default) {
  ambient_dim <- as.integer(ambient_dim)
  kappa <- as.numeric(kappa)
  psi <- as.numeric(psi)
  abs_kappa_psi_tol <- as.numeric(abs_kappa_psi_tol)

  if (length(ambient_dim) != 1L || !is.finite(ambient_dim)) {
    return(FALSE)
  }
  if (length(kappa) != 1L || length(psi) != 1L || !is.finite(kappa) || !is.finite(psi)) {
    return(FALSE)
  }
  if (length(abs_kappa_psi_tol) != 1L || !is.finite(abs_kappa_psi_tol) || abs_kappa_psi_tol <= 0) {
    stop("`abs_kappa_psi_tol` must be a strictly positive finite scalar.")
  }

  # Numerical regularization on S^2:
  # the (kappa, psi) parametrization becomes ill-conditioned near the vMF
  # submodel psi = 0. In that regime we deliberately replace the JP model by
  # its exact vMF limit once |kappa * psi| falls below a user-controlled
  # threshold. This is not an exact JP evaluation for nonzero psi.
  ambient_dim == 3L && abs(kappa * psi) <= abs_kappa_psi_tol
}

jp_log_one_minus_u2_term <- function(u, exponent) {
  u <- as.numeric(u)
  if (abs(exponent) <= 1e-15) {
    return(rep(0, length(u)))
  }

  exponent * log(pmax(0, 1 - u^2))
}

jp_uniform_sphere <- function(n, ambient_dim) {
  z <- matrix(stats::rnorm(n * ambient_dim), nrow = n, ncol = ambient_dim)
  z / sqrt(rowSums(z^2))
}

jp_orthonormal_complement <- function(mu) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  q_complete <- qr.Q(qr(cbind(mu, diag(length(mu)))), complete = TRUE)
  q_complete[, -1L, drop = FALSE]
}

jp_build_u_grid <- function(n_u = 4097L) {
  n_u <- as.integer(n_u)
  if (length(n_u) != 1L || !is.finite(n_u) || n_u < 17L || n_u %% 2L != 1L) {
    stop("`n_u` must be an odd integer >= 17.")
  }

  phi <- seq(pi, 0, length.out = n_u)
  u <- cos(phi)
  list(
    phi = phi,
    u = u,
    du = diff(u),
    sqrt_one_minus_u2 = sqrt(pmax(0, 1 - u^2))
  )
}

jp_cdf_from_density_grid <- function(u, density) {
  density <- as.numeric(density)
  if (length(u) != length(density) || length(u) < 2L) {
    stop("`u` and `density` must have the same length >= 2.")
  }

  increments <- 0.5 * (density[-1L] + density[-length(density)]) * diff(u)
  cdf <- c(0, cumsum(increments))
  total_mass <- cdf[[length(cdf)]]
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("Failed to build a valid CDF table for the Jones-Pewsey distribution.")
  }

  cdf <- cdf / total_mass
  cdf[[1L]] <- 0
  cdf[[length(cdf)]] <- 1
  cdf <- cummax(pmin(pmax(cdf, 0), 1))
  cdf[[length(cdf)]] <- 1
  cdf
}

jp_prepare_inverse_cdf_table <- function(u, cdf) {
  cdf <- cummax(pmin(pmax(as.numeric(cdf), 0), 1))
  cdf[[1L]] <- 0
  cdf[[length(cdf)]] <- 1

  keep <- c(TRUE, diff(cdf) > 0)
  keep[[length(keep)]] <- TRUE

  list(
    u = as.numeric(u)[keep],
    cdf = cdf[keep]
  )
}

jp_prepare_cdf_interpolator <- function(x_grid,
                                        cdf_grid) {
  x_grid <- as.numeric(x_grid)
  cdf_grid <- as.numeric(cdf_grid)
  if (length(x_grid) != length(cdf_grid) || length(x_grid) < 2L) {
    stop("`x_grid` and `cdf_grid` must have the same length >= 2.")
  }

  list(
    x_min = x_grid[[1L]],
    x_max = x_grid[[length(x_grid)]],
    fn = stats::splinefun(
      x = x_grid,
      y = cdf_grid,
      method = "monoH.FC"
    )
  )
}

jp_prepare_inverse_cdf_interpolator <- function(inverse_table) {
  if (!is.list(inverse_table) || !all(c("cdf", "u") %in% names(inverse_table))) {
    stop("`inverse_table` must be a list with entries `cdf` and `u`.")
  }

  jp_prepare_cdf_interpolator(
    x_grid = inverse_table$cdf,
    cdf_grid = inverse_table$u
  )
}

jp_interpolate_cdf_linear <- function(x, x_grid, cdf_grid) {
  x <- pmin(pmax(as.numeric(x), min(x_grid)), max(x_grid))
  stats::approx(
    x = x_grid,
    y = cdf_grid,
    xout = x,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}

jp_interpolate_cdf <- function(x, x_grid, cdf_grid) {
  jp_interpolate_cdf_spline(x = x, x_grid = x_grid, cdf_grid = cdf_grid)
}

jp_interpolate_cdf_spline <- function(x,
                                      x_grid,
                                      cdf_grid,
                                      interpolator = NULL) {
  interpolator <- interpolator %||% jp_prepare_cdf_interpolator(
    x_grid = x_grid,
    cdf_grid = cdf_grid
  )
  x_eval <- pmin(pmax(as.numeric(x), interpolator$x_min), interpolator$x_max)
  pmin(pmax(interpolator$fn(x_eval), 0), 1)
}

jp_interpolate_upper_tail <- function(threshold,
                                      x_grid,
                                      cdf_grid,
                                      interpolator = NULL) {
  1 - jp_interpolate_cdf_spline(
    x = threshold,
    x_grid = x_grid,
    cdf_grid = cdf_grid,
    interpolator = interpolator
  )
}

jp_interpolate_upper_tail_linear <- function(threshold,
                                             x_grid,
                                             cdf_grid) {
  1 - jp_interpolate_cdf_linear(
    x = threshold,
    x_grid = x_grid,
    cdf_grid = cdf_grid
  )
}

jp_interpolate_upper_tail_spline <- function(threshold,
                                             x_grid,
                                             cdf_grid,
                                             interpolator = NULL) {
  1 - jp_interpolate_cdf_spline(
    x = threshold,
    x_grid = x_grid,
    cdf_grid = cdf_grid,
    interpolator = interpolator
  )
}

jp_uniform_s2_distance_profile <- function(t_values,
                                           distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  t_values <- as.numeric(t_values)
  output <- numeric(length(t_values))

  if (identical(distance_type, "geodesic")) {
    output[t_values <= 0] <- 0
    output[t_values >= pi] <- 1
    active <- which(t_values > 0 & t_values < pi)
    output[active] <- (1 - cos(t_values[active])) / 2
    return(output)
  }

  output[t_values <= 0] <- 0
  output[t_values >= 2] <- 1
  active <- which(t_values > 0 & t_values < 2)
  output[active] <- (t_values[active]^2) / 4
  output
}

jp_vmf_axis_density <- function(t, q, kappa, log = FALSE) {
  t <- as.numeric(t)
  out <- rep(if (log) -Inf else 0, length(t))
  valid <- is.finite(t) & t >= -1 & t <= 1
  if (!any(valid)) {
    return(out)
  }

  if (kappa <= 1e-12) {
    log_density <- log(sphere_surface_area(q - 1L)) - log(sphere_surface_area(q)) +
      jp_log_one_minus_u2_term(t[valid], q / 2 - 1)
  } else {
    log_density <- log(sphere_surface_area(q - 1L)) +
      rotasym::c_vMF(p = q + 1L, kappa = kappa, log = TRUE) +
      kappa * t[valid] +
      jp_log_one_minus_u2_term(t[valid], q / 2 - 1)
  }

  out[valid] <- if (log) log_density else exp(log_density)
  out
}

jp_vmf_projected_density <- function(s, rho, ambient_dim, kappa, log = FALSE) {
  s <- as.numeric(s)
  out <- rep(if (log) -Inf else 0, length(s))
  valid <- is.finite(s) & s >= -1 & s <= 1
  if (!any(valid)) {
    return(out)
  }

  rho <- pmin(pmax(as.numeric(rho), -1), 1)
  exponent <- (ambient_dim - 3L) / 2
  radial_term <- pmax(0, 1 - s[valid]^2)
  kappa_term <- kappa * sqrt(radial_term * pmax(0, 1 - rho^2))

  log_density <- rotasym::c_vMF(p = ambient_dim, kappa = kappa, log = TRUE) +
    kappa * rho * s[valid] +
    jp_log_one_minus_u2_term(s[valid], exponent) -
    rotasym::c_vMF(p = ambient_dim - 1L, kappa = kappa_term, log = TRUE)

  out[valid] <- if (log) log_density else exp(log_density)
  out
}

jp_normalize_probability_weights <- function(weights, n_expected) {
  weights <- as.numeric(weights)

  if (length(weights) != n_expected) {
    stop("`weights` has incompatible length.")
  }
  if (any(!is.finite(weights))) {
    stop("`weights` must be finite.")
  }
  if (any(weights < 0)) {
    stop("`weights` must be nonnegative.")
  }

  total_weight <- sum(weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("`weights` must have strictly positive sum.")
  }

  weights / total_weight
}

jp_log_abs_diff_exp <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  max_ab <- pmax(a, b)
  min_ab <- pmin(a, b)
  max_ab + log1p(-exp(min_ab - max_ab))
}

jp_log_norm_constant_s2 <- function(alpha, beta) {
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)

  len <- max(length(alpha), length(beta))
  alpha <- rep_len(alpha, len)
  beta <- rep_len(beta, len)

  if (any(!is.finite(alpha)) || any(abs(alpha) >= 1)) {
    stop("`alpha` must be finite and lie in (-1, 1).")
  }
  if (any(!is.finite(beta))) {
    stop("`beta` must be finite.")
  }

  out <- numeric(len)
  alpha_zero <- abs(alpha) <= 1e-15
  out[alpha_zero] <- -log(4 * pi)

  active <- which(!alpha_zero)
  if (length(active) == 0L) {
    return(out)
  }

  beta_active <- beta[active]
  alpha_active <- alpha[active]
  beta_near_minus_one <- abs(beta_active + 1) <= 1e-8

  if (any(beta_near_minus_one)) {
    idx <- active[beta_near_minus_one]
    log_ratio <- log1p(alpha[idx]) - log1p(-alpha[idx])
    out[idx] <- -(log(2 * pi) + log(abs(log_ratio)) - log(abs(alpha[idx])))
  }

  if (any(!beta_near_minus_one)) {
    idx <- active[!beta_near_minus_one]
    g <- beta[idx] + 1
    a_term <- g * log1p(alpha[idx])
    b_term <- g * log1p(-alpha[idx])
    out[idx] <- -(log(2 * pi) +
      jp_log_abs_diff_exp(a_term, b_term) -
      log(abs(alpha[idx] * g)))
  }

  out
}


jp_weighted_loglik_s2 <- function(mu, kappa, psi, x, prob_weights) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("`jp_weighted_loglik_s2()` currently supports only S^2 data.")
  }

  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  if (params$ambient_dim != 3L) {
    stop("`jp_weighted_loglik_s2()` currently supports only S^2.")
  }

  prob_weights <- as.numeric(prob_weights)
  if (length(prob_weights) != nrow(x)) {
    stop("`prob_weights` must have length nrow(x).")
  }

  z <- as.numeric(x %*% params$mu)

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # Regularized evaluation near the numerically ill-conditioned vMF limit.
    # On S^2, once |kappa * psi| is sufficiently small we replace the JP
    # loglikelihood by the exact vMF one. This stabilizes optimization but is
    # not the exact JP likelihood for nonzero psi.
    return(rotasym::c_vMF(p = 3L, kappa = params$kappa, log = TRUE) +
      params$kappa * sum(prob_weights * z))
  }

  term <- 1 + params$alpha * z
  if (any(term <= 0) || any(!is.finite(term))) {
    return(-Inf)
  }

  jp_log_norm_constant_s2(params$alpha, params$beta) +
    params$beta * sum(prob_weights * log(term))
}

#' Weighted log-likelihood for Jones-Pewsey model on S^2 (pre-validated/prepared)
jp_weighted_loglik_s2_prepared <- function(mu, kappa, psi, x, prob_weights) {
  mu <- as.numeric(mu)
  kappa <- as.numeric(kappa)
  psi <- as.numeric(psi)

  if (length(mu) != 3L || length(kappa) != 1L || length(psi) != 1L) {
    return(-Inf)
  }
  if (any(!is.finite(mu)) || !is.finite(kappa) || !is.finite(psi) || kappa < 0) {
    return(-Inf)
  }

  mu_norm <- sqrt(sum(mu^2))
  if (!is.finite(mu_norm) || mu_norm <= 0) {
    return(-Inf)
  }
  if (abs(mu_norm - 1) > 1e-8) {
    mu <- mu / mu_norm
  }

  prob_weights <- as.numeric(prob_weights)
  if (length(prob_weights) != nrow(x) || any(!is.finite(prob_weights))) {
    return(-Inf)
  }

  z <- as.numeric(x %*% mu)

  if (jp_is_near_zero_vmf_s2(ambient_dim = 3L, kappa = kappa, psi = psi)) {
    # Same S^2 regularization as above, but with already prepared inputs:
    # use the exact vMF objective in the nearly-singular JP regime.
    return(rotasym::c_vMF(p = 3L, kappa = kappa, log = TRUE) +
      kappa * sum(prob_weights * z))
  }

  alpha <- tanh(kappa * psi)
  beta <- 1 / psi
  term <- 1 + alpha * z
  if (any(term <= 0) || any(!is.finite(term))) {
    return(-Inf)
  }

  jp_log_norm_constant_s2(alpha, beta) + beta * sum(prob_weights * log(term))
}

jp_solve_vmf_kappa_from_rbar <- function(r_bar,
                                         q,
                                         tol = 1e-10,
                                         max_kappa = 1e6) {
  r_bar <- as.numeric(r_bar)
  q <- as.numeric(q)

  if (length(r_bar) != 1L || !is.finite(r_bar) || r_bar < 0 || r_bar > 1) {
    stop("`r_bar` must be a scalar in [0, 1].")
  }
  if (length(q) != 1L || !is.finite(q) || q < 1) {
    stop("`q` must be a positive scalar.")
  }

  if (r_bar <= tol) {
    return(0)
  }
  if (r_bar >= 1 - tol) {
    return(max_kappa)
  }

  objective <- function(kappa) A_q(kappa, q) - r_bar
  upper <- 1
  while (objective(upper) < 0 && upper < max_kappa) {
    upper <- upper * 2
  }

  if (upper >= max_kappa && objective(upper) < 0) {
    return(max_kappa)
  }

  stats::uniroot(objective, interval = c(0, upper), tol = tol)$root
}

jp_weighted_vmf_mle_s2 <- function(x, prob_weights) {
  resultant <- colSums(x * prob_weights)
  r_bar <- sqrt(sum(resultant^2))

  if (!is.finite(r_bar) || r_bar <= 1e-12) {
    mu_hat <- c(0, 0, 1)
    kappa_hat <- 0
  } else {
    mu_hat <- resultant / r_bar
    kappa_hat <- jp_solve_vmf_kappa_from_rbar(r_bar, q = 2L)
  }

  list(
    mu = mu_hat,
    kappa = kappa_hat,
    psi = 0,
    loglik = jp_weighted_loglik_s2_prepared(mu_hat, kappa_hat, 0, x, prob_weights),
    source = "vmf"
  )
}

jp_mu_s2_to_raw <- function(mu) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  phi <- acos(pmin(pmax(mu[[3L]], -1), 1))
  sin_phi <- sin(phi)
  lambda <- if (abs(sin_phi) <= 1e-10) 0 else atan2(mu[[2L]], mu[[1L]])
  phi_scaled <- pmin(pmax(phi / pi, 1e-8), 1 - 1e-8)

  c(lambda, stats::qlogis(phi_scaled))
}

jp_mu_s2_from_raw <- function(lambda_raw, phi_raw) {
  phi <- pi * stats::plogis(phi_raw)
  sin_phi <- sin(phi)
  c(cos(lambda_raw) * sin_phi, sin(lambda_raw) * sin_phi, cos(phi))
}

jp_params_s2_from_raw <- function(raw, sign_branch, psi_min) {
  raw <- as.numeric(raw)
  if (length(raw) != 4L || any(!is.finite(raw))) {
    return(NULL)
  }
  if (raw[[3L]] > 700 || raw[[4L]] > 700) {
    return(NULL)
  }

  mu <- jp_mu_s2_from_raw(raw[[1L]], raw[[2L]])
  a <- exp(raw[[3L]])
  psi_abs <- psi_min + exp(raw[[4L]])
  psi <- sign_branch * psi_abs
  kappa <- a / psi_abs
  alpha <- sign_branch * tanh(a)
  beta <- 1 / psi

  list(
    mu = mu,
    a = a,
    psi_abs = psi_abs,
    psi = psi,
    kappa = kappa,
    alpha = alpha,
    beta = beta
  )
}

jp_neg_weighted_loglik_s2_raw <- function(raw,
                                          sign_branch,
                                          x,
                                          prob_weights,
                                          psi_min,
                                          max_abs_kappa_psi = Inf) {
  if (!all(is.finite(raw)) || raw[[3L]] > 700) {
    return(Inf)
  }

  a_value <- exp(raw[[3L]])
  if (!is.finite(a_value) || a_value <= 0 || a_value > max_abs_kappa_psi) {
    return(Inf)
  }

  params <- jp_params_s2_from_raw(raw, sign_branch = sign_branch, psi_min = psi_min)
  if (is.null(params)) {
    return(Inf)
  }

  loglik <- jp_weighted_loglik_s2_prepared(
    mu = params$mu,
    kappa = params$kappa,
    psi = params$psi,
    x = x,
    prob_weights = prob_weights
  )
  if (!is.finite(loglik)) {
    return(Inf)
  }

  -loglik
}

#' Weighted numerical MLE for the Jones--Pewsey model on S^2
#'
#' Unlike the composite normal, vMF, and HvMF fits used elsewhere in the bootstrap
#' pipeline, the spherical Jones--Pewsey composite fit does not have a closed-form
#' MLE in the current parametrization. The likelihood involves the normalizing
#' constant of the Jones--Pewsey family and the shape parameter `psi`, so the fit
#' is obtained by numerical maximization.
#'
#' The implementation below uses the vMF MLE as a cheap baseline/start.
#'
#' The control entries `jp_mle_sign_branches`, `jp_mle_psi_abs_starts`,
#' `jp_mle_maxit`, and `jp_mle_reltol` are therefore the main speed/accuracy
#' knobs for calibration studies. The optional entry `jp_mle_start_theta` may
#' contain a list with entries `mu`, `kappa`, and `psi`; when supplied, it is
#' added as a warm start for the numerical search. This is useful for bootstrap
#' refits, where the observed MLE is often a much better starting point than a
#' fresh vMF initialization.
jp_mle_s2_weighted <- function(data, weights = NULL, control = list()) {
  data_already_normalized <- isTRUE(control$jp_data_already_normalized)
  x <- if (data_already_normalized) {
    as.matrix(data)
  } else {
    jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  }
  if (ncol(x) != 3L) {
    stop("`jp_mle_s2_weighted()` currently supports only S^2 data.")
  }

  psi_min <- as.numeric(control$jp_mle_psi_min %||% 1e-3)
  psi_abs_starts <- as.numeric(control$jp_mle_psi_abs_starts %||% c(0.25, 0.5, 1, 2))
  maxit <- as.integer(control$jp_mle_maxit %||% 500L)
  reltol <- as.numeric(control$jp_mle_reltol %||% 1e-10)
  optim_method <- as.character(control$jp_mle_method %||% "BFGS")
  max_abs_kappa_psi <- as.numeric(control$jp_mle_max_abs_kappa_psi %||% Inf)
  branches <- as.integer(control$jp_mle_sign_branches %||% c(-1L, 1L))
  vmf_switch_abs_kappa_psi <- as.numeric(
    control$jp_vmf_switch_abs_kappa_psi %||% jp_vmf_near_zero_abs_kappa_psi_default
  )
  bootstrap_refit <- isTRUE(control$jp_mle_bootstrap_refit)
  bootstrap_allow_global_fallback <- isTRUE(control$jp_mle_bootstrap_allow_global_fallback)
  warm_start_only <- isTRUE(control$jp_mle_warm_start_only) ||
    bootstrap_refit

  if (length(psi_min) != 1L || !is.finite(psi_min) || psi_min <= 0) {
    stop("`control$jp_mle_psi_min` must be a strictly positive scalar.")
  }
  if (length(psi_abs_starts) == 0L || any(!is.finite(psi_abs_starts)) || any(psi_abs_starts <= psi_min)) {
    stop("`control$jp_mle_psi_abs_starts` must contain finite values strictly greater than `psi_min`.")
  }
  if (length(branches) == 0L || any(!branches %in% c(-1L, 1L))) {
    stop("`control$jp_mle_sign_branches` must be a vector with entries in {-1, 1}.")
  }
  if (length(max_abs_kappa_psi) != 1L || is.na(max_abs_kappa_psi) || max_abs_kappa_psi <= 0) {
    stop("`control$jp_mle_max_abs_kappa_psi` must be a strictly positive scalar (finite or Inf).")
  }
  if (length(vmf_switch_abs_kappa_psi) != 1L ||
      !is.finite(vmf_switch_abs_kappa_psi) ||
      vmf_switch_abs_kappa_psi <= 0) {
    stop("`control$jp_vmf_switch_abs_kappa_psi` must be a strictly positive finite scalar.")
  }

  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  vmf_candidate <- NULL
  raw_mu_start <- NULL
  candidates <- list()
  if (!isTRUE(warm_start_only)) {
    vmf_candidate <- jp_weighted_vmf_mle_s2(x, prob_weights)
    raw_mu_start <- jp_mu_s2_to_raw(vmf_candidate$mu)
    candidates <- list(vmf_candidate)
  }

  warm_start_theta <- control$jp_mle_start_theta %||% NULL
  warm_mu <- NULL
  warm_kappa <- NA_real_
  warm_psi <- NA_real_
  warm_is_valid <- FALSE
  warm_is_near_zero_vmf <- FALSE

  if (!is.null(warm_start_theta)) {
    warm_mu_try <- try(
      jp_normalize_unit_vector(
        warm_start_theta$mu,
        arg_name = "`control$jp_mle_start_theta$mu`",
        min_length = 3L
      ),
      silent = TRUE
    )
    warm_kappa_try <- suppressWarnings(as.numeric(warm_start_theta$kappa))
    warm_psi_try <- suppressWarnings(as.numeric(warm_start_theta$psi))

    if (!inherits(warm_mu_try, "try-error") &&
        length(warm_kappa_try) == 1L && is.finite(warm_kappa_try) && warm_kappa_try >= 0 &&
        length(warm_psi_try) == 1L && is.finite(warm_psi_try)) {
      warm_mu <- warm_mu_try
      warm_kappa <- warm_kappa_try
      warm_psi <- warm_psi_try
      warm_is_valid <- TRUE
      warm_is_near_zero_vmf <- jp_is_near_zero_vmf_s2(
        ambient_dim = 3L,
        kappa = warm_kappa,
        psi = warm_psi,
        abs_kappa_psi_tol = vmf_switch_abs_kappa_psi
      )
    }
  }

  if (bootstrap_refit) {
    if (!warm_is_valid) {
      if (bootstrap_allow_global_fallback) {
        fallback_control <- control
        fallback_control$jp_mle_warm_start_only <- FALSE
        fallback_control$jp_mle_bootstrap_refit <- FALSE
        return(jp_mle_s2_weighted(data = x, weights = prob_weights, control = fallback_control))
      }
      stop("JP bootstrap refit requires a valid `jp_mle_start_theta` (theta_hat).")
    }

    if (warm_is_near_zero_vmf) {
      best_vmf <- jp_weighted_vmf_mle_s2(x, prob_weights)
      return(list(
        mu = best_vmf$mu,
        kappa = best_vmf$kappa,
        psi = best_vmf$psi
      ))
    }

    # Stabilized bootstrap refit:
    # once the observed fit selects a nonzero-psi JP branch, the weighted
    # bootstrap refits are kept on that same sign branch unless the caller
    # explicitly overrides this behavior. This avoids unstable cross-branch
    # jumps during local optimization, but it is a computational restriction,
    # not the full unconstrained composite JP re-optimization.
    # TODO: revisit this restriction empirically. In particular, check whether
    # the bootstrap calibration remains acceptable if branch forcing is removed
    # and both JP sign branches are explored again in the refits.
    branches <- as.integer(ifelse(warm_psi >= 0, 1L, -1L))
  }

  raw_start_entries <- list()
  if (warm_is_valid) {
      warm_loglik <- jp_weighted_loglik_s2_prepared(
        mu = warm_mu,
        kappa = warm_kappa,
        psi = warm_psi,
        x = x,
        prob_weights = prob_weights
      )
      warm_abs_kappa_psi <- abs(warm_kappa * warm_psi)
      if (is.finite(warm_loglik) && is.finite(warm_abs_kappa_psi) && warm_abs_kappa_psi <= max_abs_kappa_psi) {
        candidates[[length(candidates) + 1L]] <- list(
          mu = warm_mu,
          kappa = warm_kappa,
          psi = warm_psi,
          loglik = warm_loglik,
          source = "warm_start"
        )
        warm_a <- max(abs(warm_kappa * warm_psi), 1e-4)
        warm_psi_abs <- max(abs(warm_psi), psi_min * (1 + 1e-6))
        raw_start_entries[[length(raw_start_entries) + 1L]] <- list(
          sign_branch = ifelse(warm_psi >= 0, 1L, -1L),
          raw_start = c(
            jp_mu_s2_to_raw(warm_mu),
            log(warm_a),
            log(warm_psi_abs - psi_min)
          )
        )
      }
  }

  if (!isTRUE(warm_start_only)) {
    for (sign_branch in unique(branches)) {
      for (psi_abs0 in psi_abs_starts) {
        a0 <- max(vmf_candidate$kappa * psi_abs0, 1e-4)
        raw_start_entries[[length(raw_start_entries) + 1L]] <- list(
          sign_branch = sign_branch,
          raw_start = c(
            raw_mu_start,
            log(a0),
            log(psi_abs0 - psi_min)
          )
        )
      }
    }
  }

  if (isTRUE(warm_start_only) && length(raw_start_entries) == 0L) {
    if (bootstrap_refit && !bootstrap_allow_global_fallback) {
      stop("JP bootstrap refit failed to build local warm-start entries.")
    }
    fallback_control <- control
    fallback_control$jp_mle_warm_start_only <- FALSE
    fallback_control$jp_mle_bootstrap_refit <- FALSE
    return(jp_mle_s2_weighted(data = x, weights = prob_weights, control = fallback_control))
  }

  for (raw_entry in raw_start_entries) {
    sign_branch <- raw_entry$sign_branch
    raw_start <- raw_entry$raw_start

    fit <- try(
      stats::optim(
        par = raw_start,
        fn = jp_neg_weighted_loglik_s2_raw,
        sign_branch = sign_branch,
        x = x,
        prob_weights = prob_weights,
        psi_min = psi_min,
        max_abs_kappa_psi = max_abs_kappa_psi,
        method = optim_method,
        control = list(maxit = maxit, reltol = reltol)
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error") || !is.list(fit) || fit$convergence != 0L || !is.finite(fit$value)) {
      next
    }

    params <- jp_params_s2_from_raw(fit$par, sign_branch = sign_branch, psi_min = psi_min)
    if (is.null(params)) {
      next
    }

    abs_kappa_psi <- abs(params$kappa * params$psi)
    if (!is.finite(abs_kappa_psi) || abs_kappa_psi > max_abs_kappa_psi) {
      next
    }

    loglik <- jp_weighted_loglik_s2_prepared(
      mu = params$mu,
      kappa = params$kappa,
      psi = params$psi,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(loglik)) {
      next
    }

    candidates[[length(candidates) + 1L]] <- list(
      mu = params$mu,
      kappa = params$kappa,
      psi = params$psi,
      loglik = loglik,
      source = "jp"
    )
  }

  if (isTRUE(warm_start_only)) {
    optimized_candidate <- vapply(candidates, function(candidate) {
      identical(candidate$source, "jp") && is.finite(candidate$loglik)
    }, logical(1))

    if (!any(optimized_candidate)) {
      if (bootstrap_refit && !bootstrap_allow_global_fallback) {
        finite_candidate <- vapply(candidates, function(candidate) {
          is.finite(candidate$loglik)
        }, logical(1))
        if (any(finite_candidate)) {
          candidates <- candidates[finite_candidate]
        } else {
          stop("JP bootstrap refit local optimization did not converge.")
        }
      } else {
        fallback_control <- control
        fallback_control$jp_mle_warm_start_only <- FALSE
        fallback_control$jp_mle_bootstrap_refit <- FALSE
        return(jp_mle_s2_weighted(data = x, weights = prob_weights, control = fallback_control))
      }
    } else {
      candidates <- candidates[optimized_candidate]
    }
  }

  loglik_values <- vapply(candidates, function(candidate) candidate$loglik, numeric(1))
  if (all(!is.finite(loglik_values))) {
    stop("JP composite MLE failed: no finite candidate was found.")
  }

  best <- candidates[[which.max(loglik_values)]]

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = 3L,
    kappa = best$kappa,
    psi = best$psi,
    abs_kappa_psi_tol = vmf_switch_abs_kappa_psi
  )) {
    if (is.null(vmf_candidate)) {
      vmf_candidate <- jp_weighted_vmf_mle_s2(x, prob_weights)
    }
    best <- vmf_candidate
  }

  list(
    mu = best$mu,
    kappa = best$kappa,
    psi = best$psi
  )
}

jp_reciprocal_constant_sphere <- function(q,
                                          alpha,
                                          beta,
                                          rel.tol = 1e-10,
                                          abs.tol = 1e-12,
                                          subdivisions = 2000L) {
  if (q == 0L) {
    return((1 + alpha)^beta + (1 - alpha)^beta)
  }

  integrand <- function(phi) {
    exp(beta * log1p(alpha * cos(phi))) * sin(phi)^(q - 1L)
  }

  result <- stats::integrate(
    f = integrand,
    lower = 0,
    upper = pi,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = as.integer(subdivisions),
    stop.on.error = FALSE
  )

  integration_ok <- is.character(result$message) &&
    length(result$message) == 1L &&
    identical(result$message, "OK")

  if (integration_ok && is.finite(result$value) && result$value > 0) {
    return(sphere_surface_area(q - 1L) * result$value)
  }

  # Fallback for difficult parameter configurations where adaptive quadrature
  # reports roundoff problems or otherwise fails to certify the result.
  n_grid <- 65537L
  phi_grid <- seq(0, pi, length.out = n_grid)
  density_grid <- integrand(phi_grid)
  density_grid[!is.finite(density_grid)] <- 0
  trap_value <- sum(0.5 * (density_grid[-1L] + density_grid[-n_grid]) * diff(phi_grid))

  if (!is.finite(trap_value) || trap_value <= 0) {
    stop("Failed to evaluate the Jones-Pewsey normalizing integral.")
  }

  sphere_surface_area(q - 1L) * trap_value
}

#' Jones-Pewsey normalizing constant on S^q
#' @param q Sphere dimension
#' @param alpha Scalar in (-1, 1)
#' @param beta Scalar exponent
#' @param log Whether to return the log-constant
#' @return Normalizing constant for the Jones-Pewsey density on S^q
c_jp_sphere <- function(q, alpha, beta, log = FALSE) {
  q <- as.integer(q)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)

  if (length(q) != 1L || !is.finite(q) || q < 0L) {
    stop("`q` must be a nonnegative integer.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || abs(alpha) >= 1) {
    stop("`alpha` must be a finite scalar in (-1, 1).")
  }
  if (length(beta) != 1L || !is.finite(beta)) {
    stop("`beta` must be a finite scalar.")
  }

  if (abs(alpha) <= 1e-15) {
    log_const <- -log(sphere_surface_area(q))
    return(if (log) log_const else exp(log_const))
  }

  cache_key <- sprintf("cjp_q%d_a%.17g_b%.17g", q, alpha, beta)
  if (!exists(cache_key, envir = jp_cache_env, inherits = FALSE)) {
    reciprocal_constant <- jp_reciprocal_constant_sphere(q = q, alpha = alpha, beta = beta)
    assign(cache_key, 1 / reciprocal_constant, envir = jp_cache_env)
  }

  constant <- get(cache_key, envir = jp_cache_env, inherits = FALSE)
  if (log) log(constant) else constant
}

get_jp_delta_reciprocal_table <- function(q,
                                          beta,
                                          n_delta = 1025L,
                                          delta_max = 1 - 1e-10) {
  q <- as.integer(q)
  n_delta <- as.integer(n_delta)
  delta_max <- as.numeric(delta_max)

  if (length(q) != 1L || !is.finite(q) || q < 1L) {
    stop("`q` must be an integer >= 1 for the delta table.")
  }
  if (length(beta) != 1L || !is.finite(beta)) {
    stop("`beta` must be finite.")
  }
  if (length(n_delta) != 1L || !is.finite(n_delta) || n_delta < 33L) {
    stop("`n_delta` must be an integer >= 33.")
  }
  if (length(delta_max) != 1L || !is.finite(delta_max) || delta_max <= 0 || delta_max >= 1) {
    stop("`delta_max` must be a scalar in (0, 1).")
  }

  cache_key <- sprintf("jp_delta_recip_q%d_b%.17g_n%d_dmax%.17g", q, beta, n_delta, delta_max)
  if (!exists(cache_key, envir = jp_cache_env, inherits = FALSE)) {
    base_grid <- seq(0, 1, length.out = n_delta)
    delta_grid <- delta_max * (1 - (1 - base_grid)^2)
    reciprocal_grid <- vapply(delta_grid, function(delta) {
      1 / c_jp_sphere(q = q, alpha = delta, beta = beta)
    }, numeric(1))

    assign(
      cache_key,
      data.frame(delta = delta_grid, reciprocal = reciprocal_grid),
      envir = jp_cache_env
    )
  }

  get(cache_key, envir = jp_cache_env, inherits = FALSE)
}

jp_max_delta_from_rho <- function(alpha, rho, n_search = 1025L) {
  alpha_signed <- as.numeric(alpha)
  alpha <- abs(alpha_signed)
  rho <- pmin(pmax(as.numeric(rho), -1), 1)
  if (alpha <= 1e-15 || length(rho) == 0L) {
    return(0)
  }

  u_grid_object <- jp_build_u_grid(n_u = as.integer(n_search))
  max_delta <- 0
  for (rho_value in rho) {
    num0 <- 1 + alpha_signed * rho_value * u_grid_object$u
    delta_values <- alpha * sqrt(pmax(0, 1 - rho_value^2)) * u_grid_object$sqrt_one_minus_u2 / num0
    max_delta <- max(max_delta, delta_values)
  }

  pmin(1 - 1e-10, max_delta + 1e-6)
}

jp_eval_delta_reciprocal <- function(delta, reciprocal_table) {
  delta <- abs(as.numeric(delta))
  max_delta <- reciprocal_table$delta[[nrow(reciprocal_table)]]
  delta <- pmin(pmax(delta, 0), max_delta)

  stats::approx(
    x = reciprocal_table$delta,
    y = reciprocal_table$reciprocal,
    xout = delta,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}

jp_projected_density_matrix <- function(rho,
                                        q,
                                        alpha,
                                        beta,
                                        u_grid_object,
                                        reciprocal_table = NULL) {
  rho <- pmin(pmax(as.numeric(rho), -1), 1)
  n_u <- length(u_grid_object$u)
  n_rho <- length(rho)

  if (n_rho == 0L) {
    return(matrix(0, nrow = n_u, ncol = 0L))
  }

  if (is.null(reciprocal_table)) {
    reciprocal_table <- get_jp_delta_reciprocal_table(
      q = q - 1L,
      beta = beta,
      delta_max = jp_max_delta_from_rho(alpha = alpha, rho = rho)
    )
  }

  log_c_q <- c_jp_sphere(q = q, alpha = alpha, beta = beta, log = TRUE)
  base_log_term <- jp_log_one_minus_u2_term(u_grid_object$u, q / 2 - 1)
  sqrt_one_minus_rho2 <- sqrt(pmax(0, 1 - rho^2))
  density_matrix <- matrix(0, nrow = n_u, ncol = n_rho)

  for (j in seq_along(rho)) {
    num0 <- 1 + alpha * rho[[j]] * u_grid_object$u
    delta <- abs(alpha) * sqrt_one_minus_rho2[[j]] * u_grid_object$sqrt_one_minus_u2 / num0
    reciprocal_values <- jp_eval_delta_reciprocal(delta, reciprocal_table)
    log_density <- log_c_q +
      base_log_term +
      beta * log(num0) +
      log(reciprocal_values)
    density_matrix[, j] <- exp(log_density)
  }

  density_matrix
}

build_jp_projected_cdf_matrix <- function(rho,
                                          q,
                                          alpha,
                                          beta,
                                          n_u = 4097L,
                                          n_delta = 1025L) {
  u_grid_object <- jp_build_u_grid(n_u = n_u)
  reciprocal_table <- get_jp_delta_reciprocal_table(
    q = q - 1L,
    beta = beta,
    n_delta = n_delta,
    delta_max = jp_max_delta_from_rho(alpha = alpha, rho = rho, n_search = n_u)
  )
  density_matrix <- jp_projected_density_matrix(
    rho = rho,
    q = q,
    alpha = alpha,
    beta = beta,
    u_grid_object = u_grid_object,
    reciprocal_table = reciprocal_table
  )

  n_u <- nrow(density_matrix)
  n_rho <- ncol(density_matrix)
  if (n_rho == 0L) {
    return(list(
      u = u_grid_object$u,
      cdf = matrix(0, nrow = n_u, ncol = 0L),
      density = density_matrix
    ))
  }

  increments <- 0.5 * (density_matrix[-1L, , drop = FALSE] + density_matrix[-n_u, , drop = FALSE]) *
    u_grid_object$du
  cdf_matrix <- rbind(rep(0, n_rho), apply(increments, 2, cumsum))
  last_values <- pmax(cdf_matrix[n_u, ], .Machine$double.xmin)
  cdf_matrix <- sweep(cdf_matrix, 2, last_values, "/", check.margin = FALSE)
  cdf_matrix <- apply(cdf_matrix, 2, function(column) {
    column <- cummax(pmin(pmax(column, 0), 1))
    column[[1L]] <- 0
    column[[length(column)]] <- 1
    column
  })
  if (!is.matrix(cdf_matrix)) {
    cdf_matrix <- matrix(cdf_matrix, ncol = n_rho)
  }

  list(
    u = u_grid_object$u,
    cdf = cdf_matrix,
    density = density_matrix
  )
}

build_jp_axis_cdf_table <- function(mu,
                                    kappa,
                                    psi,
                                    grid_size = 8193L) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  u_grid_object <- jp_build_u_grid(n_u = grid_size)
  density <- d_proj_jp(
    t = u_grid_object$u,
    mu = params$mu,
    kappa = params$kappa,
    psi = params$psi,
    log = FALSE
  )
  cdf <- jp_cdf_from_density_grid(u = u_grid_object$u, density = density)

  data.frame(
    u = u_grid_object$u,
    density = density,
    cdf = cdf
  )
}

#' Density of T = mu^T X under the spherical Jones-Pewsey law
#' @param t Scalar or vector in [-1, 1]
#' @param mu Mean direction on S^q
#' @param kappa Nonnegative concentration parameter
#' @param psi Shape parameter
#' @param log Whether to return the log-density
#' @return Density values for T = mu^T X
d_proj_jp <- function(t, mu, kappa, psi, log = FALSE) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  t <- as.numeric(t)
  out <- rep(if (log) -Inf else 0, length(t))
  valid <- is.finite(t) & t >= -1 & t <= 1
  if (!any(valid)) {
    return(out)
  }

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # The axis density inherits the same S^2 regularization: near psi = 0 we
    # evaluate the exact vMF projected density instead of the JP one.
    return(jp_vmf_axis_density(
      t = t,
      q = params$q,
      kappa = params$kappa,
      log = log
    ))
  }

  log_density <- log(sphere_surface_area(params$q - 1L)) +
    c_jp_sphere(q = params$q, alpha = params$alpha, beta = params$beta, log = TRUE) +
    params$beta * log1p(params$alpha * t[valid]) +
    jp_log_one_minus_u2_term(t[valid], params$q / 2 - 1)

  out[valid] <- if (log) log_density else exp(log_density)
  out
}

#' Projected density of S_omega = omega^T X under spherical Jones-Pewsey
#' @param s Scalar or vector in [-1, 1]
#' @param omega Center direction on S^q
#' @param mu Mean direction on S^q
#' @param kappa Nonnegative concentration parameter
#' @param psi Shape parameter
#' @param log Whether to return the log-density
#' @return Density values for S_omega = omega^T X
p_proj_jp <- function(s, omega, mu, kappa, psi, log = FALSE) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = length(params$mu))
  if (length(omega) != length(params$mu)) {
    stop("`omega` and `mu` must have the same length.")
  }

  s <- as.numeric(s)
  out <- rep(if (log) -Inf else 0, length(s))
  valid <- is.finite(s) & s >= -1 & s <= 1
  if (!any(valid)) {
    return(out)
  }

  rho <- pmin(pmax(sum(omega * params$mu), -1), 1)
  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # Same regularized vMF limit for general one-dimensional projections.
    return(jp_vmf_projected_density(
      s = s,
      rho = rho,
      ambient_dim = params$ambient_dim,
      kappa = params$kappa,
      log = log
    ))
  }

  reciprocal_table <- get_jp_delta_reciprocal_table(
    q = params$q - 1L,
    beta = params$beta,
    delta_max = jp_max_delta_from_rho(alpha = params$alpha, rho = rho)
  )
  num0 <- 1 + params$alpha * rho * s[valid]
  delta <- abs(params$alpha) * sqrt(pmax(0, 1 - rho^2)) * sqrt(pmax(0, 1 - s[valid]^2)) / num0
  reciprocal_values <- jp_eval_delta_reciprocal(delta, reciprocal_table)
  log_density <- c_jp_sphere(q = params$q, alpha = params$alpha, beta = params$beta, log = TRUE) +
    jp_log_one_minus_u2_term(s[valid], params$q / 2 - 1) +
    params$beta * log(num0) +
    log(reciprocal_values)

  out[valid] <- if (log) log_density else exp(log_density)
  out
}

#' Theoretical distance profile under the spherical Jones-Pewsey law
#' @param omega Reference point on S^q
#' @param t_values Distance thresholds
#' @param mu Mean direction on S^q
#' @param kappa Nonnegative concentration parameter
#' @param psi Shape parameter
#' @param distance_type Either "geodesic" or "chordal"
#' @return Vector of probabilities P(d(X, omega) <= t)
distance_profile_jp <- function(omega,
                                t_values,
                                mu,
                                kappa,
                                psi,
                                distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  t_values <- as.numeric(t_values)

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = length(params$mu))
    if (ncol(omega) != length(params$mu)) {
      stop("`omega` and `mu` must have compatible ambient dimensions.")
    }
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }

    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_jp(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = params$mu,
        kappa = params$kappa,
        psi = params$psi,
        distance_type = distance_type
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = length(params$mu))
  if (length(omega) != length(params$mu)) {
    stop("`omega` and `mu` must have the same length.")
  }

  if (params$ambient_dim == 3L && params$psi != 0 && abs(params$alpha) <= 1e-15) {
    return(jp_uniform_s2_distance_profile(t_values, distance_type = distance_type))
  }

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # Theoretical profiles near the JP-vMF boundary are evaluated through the
    # exact vMF profile. This preserves numerical stability in S^2 at the cost
    # of replacing JP by its limit model when |kappa * psi| is small.
    if (params$ambient_dim == 3L) {
      return(theoretical_distance_profile_vmf_s2_fast(
        omega = omega,
        mu = params$mu,
        kappa = params$kappa,
        t_values = t_values,
        distance_type = distance_type
      ))
    }

    return(theoretical_distance_profile_vmf(
      omega = omega,
      mu = params$mu,
      kappa = params$kappa,
      t_values = t_values,
      distance_type = distance_type
    ))
  }

  output <- numeric(length(t_values))
  if (identical(distance_type, "geodesic")) {
    output[t_values <= 0] <- 0
    output[t_values >= pi] <- 1
    active <- which(t_values > 0 & t_values < pi)
  } else {
    output[t_values <= 0] <- 0
    output[t_values >= 2] <- 1
    active <- which(t_values > 0 & t_values < 2)
  }

  if (length(active) == 0L) {
    return(output)
  }

  rho <- pmin(pmax(sum(omega * params$mu), -1), 1)
  cdf_object <- build_jp_projected_cdf_matrix(
    rho = rho,
    q = params$q,
    alpha = params$alpha,
    beta = params$beta,
    n_u = 4097L,
    n_delta = 1025L
  )
  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  thresholds <- pmin(pmax(thresholds, -1), 1)
  cdf_interpolator <- jp_prepare_cdf_interpolator(
    x_grid = cdf_object$u,
    cdf_grid = cdf_object$cdf[, 1L]
  )
  output[active] <- jp_interpolate_upper_tail(
    threshold = thresholds,
    x_grid = cdf_object$u,
    cdf_grid = cdf_object$cdf[, 1L],
    interpolator = cdf_interpolator
  )

  output
}

#' Sample from the spherical Jones-Pewsey law by inversion of the projected CDF
#' @param n Sample size
#' @param mu Mean direction on S^q
#' @param kappa Nonnegative concentration parameter
#' @param psi Shape parameter
#' @param grid_size Odd grid size for the tabulated projected CDF
#' @param check Whether to validate the output norms
#' @return An n x (q + 1) matrix with rows on S^q
r_sph_jp <- function(n,
                     mu,
                     kappa,
                     psi,
                     grid_size = 8193L,
                     check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    # The sampler follows the same regularization policy near psi = 0 on S^2,
    # drawing from the exact vMF limit instead of inverting the JP CDF.
    return(rotasym::r_vMF(n = n, mu = params$mu, kappa = params$kappa))
  }
  if (abs(params$alpha) <= 1e-15) {
    return(jp_uniform_sphere(n = n, ambient_dim = params$ambient_dim))
  }

  cdf_table <- build_jp_axis_cdf_table(
    mu = params$mu,
    kappa = params$kappa,
    psi = params$psi,
    grid_size = grid_size
  )
  inverse_table <- jp_prepare_inverse_cdf_table(u = cdf_table$u, cdf = cdf_table$cdf)
  inverse_interpolator <- jp_prepare_inverse_cdf_interpolator(inverse_table)
  u_samples <- stats::runif(n)
  t_samples <- inverse_interpolator$fn(
    pmin(pmax(u_samples, inverse_interpolator$x_min), inverse_interpolator$x_max)
  )
  t_samples <- pmin(pmax(t_samples, -1), 1)

  xi <- jp_uniform_sphere(n = n, ambient_dim = params$q)
  b_mu <- jp_orthonormal_complement(params$mu)
  tangent_part <- t(b_mu %*% t(xi))
  radial_scales <- sqrt(pmax(0, 1 - t_samples^2))
  x <- tcrossprod(t_samples, params$mu) + sweep(tangent_part, 1, radial_scales, "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Jones-Pewsey sampler returned non-unit vectors.")
    }
  }

  x
}

distance_profile_jp_grid <- function(omega_grid,
                                     mu,
                                     kappa,
                                     psi,
                                     t_grid,
                                     distance_type = "geodesic",
                                     n_u = 4097L,
                                     n_delta = 1025L) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  omega_grid <- jp_normalize_unit_matrix(
    omega_grid,
    arg_name = "`omega_grid`",
    min_ncol = length(params$mu)
  )
  if (ncol(omega_grid) != length(params$mu)) {
    stop("`omega_grid` and `mu` must have compatible ambient dimensions.")
  }

  distance_type <- match.arg(distance_type, choices = c("geodesic", "chordal"))
  t_grid <- as.numeric(t_grid)

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  ) && params$ambient_dim == 3L) {
    return(distance_profile_vmf_s2_grid(
      omega_grid = omega_grid,
      mu = params$mu,
      kappa = params$kappa,
      t_grid = t_grid,
      distance_type = distance_type,
      n_u = n_u
    ))
  }

  rho <- as.numeric(omega_grid %*% params$mu)
  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    return(t(vapply(seq_len(nrow(omega_grid)), function(i) {
      theoretical_distance_profile_vmf(
        omega = omega_grid[i, ],
        mu = params$mu,
        kappa = params$kappa,
        t_values = t_grid,
        distance_type = distance_type
      )
    }, numeric(length(t_grid)))))
  }

  cdf_object <- build_jp_projected_cdf_matrix(
    rho = rho,
    q = params$q,
    alpha = params$alpha,
    beta = params$beta,
    n_u = n_u,
    n_delta = n_delta
  )
  thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
  thresholds <- pmin(pmax(thresholds, -1), 1)
  cdf_interpolators <- lapply(seq_along(rho), function(i) {
    jp_prepare_cdf_interpolator(
      x_grid = cdf_object$u,
      cdf_grid = cdf_object$cdf[, i]
    )
  })

  output <- t(vapply(seq_along(rho), function(i) {
    jp_interpolate_upper_tail_spline(
      threshold = thresholds,
      x_grid = cdf_object$u,
      cdf_grid = cdf_object$cdf[, i],
      interpolator = cdf_interpolators[[i]]
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

distance_profile_jp_cvm_grid <- function(X,
                                         mu,
                                         kappa,
                                         psi,
                                         n_u = 4097L,
                                         n_delta = 1025L) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = length(params$mu))
  if (ncol(X) != length(params$mu)) {
    stop("`X` and `mu` must have compatible ambient dimensions.")
  }

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  ) && params$ambient_dim == 3L) {
    return(distance_profile_vmf_s2_cvm_grid(
      X = X,
      mu = params$mu,
      kappa = params$kappa,
      n_u = n_u
    ))
  }

  rho <- as.numeric(X %*% params$mu)
  a_matrix <- X %*% t(X)
  a_matrix <- pmin(pmax(a_matrix, -1), 1)

  if (jp_is_near_zero_vmf_s2(
    ambient_dim = params$ambient_dim,
    kappa = params$kappa,
    psi = params$psi
  )) {
    return(t(vapply(seq_len(nrow(X)), function(i) {
      theoretical_distance_profile_vmf(
        omega = X[i, ],
        mu = params$mu,
        kappa = params$kappa,
        t_values = acos(a_matrix[i, ]),
        distance_type = "geodesic"
      )
    }, numeric(nrow(X)))))
  }

  cdf_object <- build_jp_projected_cdf_matrix(
    rho = rho,
    q = params$q,
    alpha = params$alpha,
    beta = params$beta,
    n_u = n_u,
    n_delta = n_delta
  )
  cdf_interpolators <- lapply(seq_along(rho), function(i) {
    jp_prepare_cdf_interpolator(
      x_grid = cdf_object$u,
      cdf_grid = cdf_object$cdf[, i]
    )
  })

  t(vapply(seq_along(rho), function(i) {
    jp_interpolate_upper_tail_spline(
      threshold = a_matrix[i, ],
      x_grid = cdf_object$u,
      cdf_grid = cdf_object$cdf[, i],
      interpolator = cdf_interpolators[[i]]
    )
  }, numeric(nrow(X))))
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

small_circle_validate_parameters <- function(mu, kappa, nu, allow_negative_nu = FALSE) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("The Small Circle implementation currently supports only S^2, i.e. `mu` of length 3.")
  }

  kappa <- as.numeric(kappa)
  nu <- as.numeric(nu)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite scalar in [0, Inf).")
  }
  if (length(nu) != 1L || !is.finite(nu)) {
    stop("`nu` must be a finite scalar.")
  }

  if (!allow_negative_nu && (nu < 0 || nu >= 1)) {
    stop("`nu` must lie in [0, 1).")
  }
  if (allow_negative_nu && abs(nu) >= 1) {
    stop("`nu` must lie in (-1, 1) before canonicalization.")
  }

  list(mu = mu, kappa = kappa, nu = nu)
}

small_circle_canonicalize_theta <- function(mu, nu, tol = 1e-12) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  nu <- as.numeric(nu)
  if (length(nu) != 1L || !is.finite(nu) || abs(nu) >= 1) {
    stop("`nu` must be a finite scalar in (-1, 1).")
  }

  if (nu < 0) {
    mu <- -mu
    nu <- -nu
  }
  if (abs(nu) <= tol) {
    nu <- 0
  }

  list(mu = mu, nu = nu)
}

small_circle_erf <- function(x) {
  2 * stats::pnorm(sqrt(2) * as.numeric(x)) - 1
}

small_circle_log_norm_constant <- function(kappa, nu) {
  # Keep the scalar-kappa small-circle model while allowing batched axial offsets.
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = 0)
  nu <- as.numeric(nu)
  if (length(nu) == 0L || any(!is.finite(nu)) || any(nu < 0 | nu >= 1)) {
    stop("`nu` must contain finite values in [0, 1).")
  }
  if (params$kappa <= 0) {
    return(rep(0, length(nu)))
  }

  root_kappa <- sqrt(params$kappa)
  erf_sum <- small_circle_erf(root_kappa * (1 - nu)) +
    small_circle_erf(root_kappa * (1 + nu))
  log(sqrt(pi) / (4 * root_kappa)) + log(erf_sum)
}

small_circle_norm_constant <- function(kappa, nu) {
  exp(small_circle_log_norm_constant(kappa = kappa, nu = nu))
}

small_circle_axis_density <- function(z, kappa, nu, log = FALSE) {
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = nu)
  z <- as.numeric(z)
  out <- rep(if (log) -Inf else 0, length(z))
  valid <- is.finite(z) & z >= -1 & z <= 1
  if (!any(valid)) {
    return(out)
  }

  if (params$kappa <= 0) {
    log_density <- rep(log(0.5), sum(valid))
  } else {
    log_density <- -log(2) - small_circle_log_norm_constant(params$kappa, params$nu) -
      params$kappa * (z[valid] - params$nu)^2
  }

  out[valid] <- if (log) log_density else exp(log_density)
  out
}

small_circle_axis_cdf <- function(z, kappa, nu) {
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = nu)
  z <- as.numeric(z)
  out <- numeric(length(z))
  out[z <= -1] <- 0
  out[z >= 1] <- 1
  active <- which(is.finite(z) & z > -1 & z < 1)
  if (length(active) == 0L) {
    return(out)
  }

  if (params$kappa <= 0) {
    out[active] <- (z[active] + 1) / 2
    return(out)
  }

  root_kappa <- sqrt(params$kappa)
  denominator <- small_circle_erf(root_kappa * (1 - params$nu)) +
    small_circle_erf(root_kappa * (1 + params$nu))
  numerator <- small_circle_erf(root_kappa * (z[active] - params$nu)) +
    small_circle_erf(root_kappa * (1 + params$nu))
  out[active] <- numerator / denominator
  pmin(pmax(out, 0), 1)
}

small_circle_gauss_legendre_cache <- new.env(parent = emptyenv())

small_circle_gauss_legendre <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("`n` must be an integer >= 2.")
  }

  key <- as.character(n)
  if (exists(key, envir = small_circle_gauss_legendre_cache, inherits = FALSE)) {
    return(get(key, envir = small_circle_gauss_legendre_cache, inherits = FALSE))
  }

  beta <- seq_len(n - 1L) / sqrt(4 * seq_len(n - 1L)^2 - 1)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- beta
  jacobi[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- beta
  eig <- eigen(jacobi, symmetric = TRUE)
  order_idx <- order(eig$values)
  nodes <- eig$values[order_idx]
  weights <- 2 * (eig$vectors[1L, order_idx]^2)
  out <- list(nodes = nodes, weights = weights)
  assign(key, out, envir = small_circle_gauss_legendre_cache)
  out
}

small_circle_legendre_matrix <- function(x, l_max) {
  x <- pmin(pmax(as.numeric(x), -1), 1)
  l_max <- as.integer(l_max)
  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  out <- matrix(0, nrow = length(x), ncol = l_max + 1L)
  out[, 1L] <- 1
  if (l_max == 0L) {
    return(out)
  }

  out[, 2L] <- x
  if (l_max == 1L) {
    return(out)
  }

  for (ell in seq_len(l_max - 1L)) {
    out[, ell + 2L] <- ((2 * ell + 1) * x * out[, ell + 1L] - ell * out[, ell]) / (ell + 1)
  }
  out
}

small_circle_legendre_coefficients <- function(kappa,
                                               nu,
                                               l_max = 150L,
                                               quad_n = 1000L,
                                               tol = 1e-10) {
  params <- small_circle_validate_parameters(mu = c(0, 0, 1), kappa = kappa, nu = nu)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (params$kappa <= 0) {
    coeffs <- c(1, rep(0, l_max))
    return(list(coefficients = coeffs, a0_error = 0, max_abs_nonzero = 0))
  }

  quad <- small_circle_gauss_legendre(quad_n)
  legendre_matrix <- small_circle_legendre_matrix(quad$nodes, l_max = l_max)
  h_values <- exp(-params$kappa * (quad$nodes - params$nu)^2 -
    small_circle_log_norm_constant(params$kappa, params$nu))
  raw_moments <- as.numeric(crossprod(legendre_matrix, quad$weights * h_values))
  ell <- 0:l_max
  coeffs <- ((2 * ell + 1) / 2) * raw_moments
  coeffs[[1L]] <- 1

  a0_error <- abs(((1 / 2) * raw_moments[[1L]]) - 1)
  if (a0_error > tol) {
    stop(sprintf("Small Circle Legendre coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error,
    max_abs_nonzero = if (l_max >= 1L) max(abs(coeffs[-1L])) else 0
  )
}

small_circle_projection_cdf_legendre <- function(x,
                                                 r,
                                                 coefficients,
                                                 enforce_bounds = TRUE) {
  x <- pmin(pmax(as.numeric(x), -1), 1)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != 1L || !is.finite(r)) {
    stop("`r` must be a finite scalar.")
  }

  out <- (x + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- as.numeric(small_circle_legendre_matrix(r, l_max = l_max))
  p_x <- small_circle_legendre_matrix(x, l_max = l_max + 1L)
  basis_matrix <- matrix(0, nrow = length(x), ncol = l_max)
  for (ell in seq_len(l_max)) {
    basis_matrix[, ell] <- (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1))
  }

  out <- out + as.numeric(basis_matrix %*% (coefficients[-1L] * p_r[-1L]))
  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}

small_circle_projection_cdf_legendre_matrix <- function(x_matrix,
                                                        r,
                                                        coefficients,
                                                        enforce_bounds = TRUE) {
  x_matrix <- as.matrix(x_matrix)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != nrow(x_matrix) || any(!is.finite(r))) {
    stop("`r` must be a finite vector of length nrow(`x_matrix`).")
  }

  x_matrix <- pmin(pmax(x_matrix, -1), 1)
  out <- (x_matrix + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- small_circle_legendre_matrix(r, l_max = l_max)
  p_x <- small_circle_legendre_matrix(as.numeric(x_matrix), l_max = l_max + 1L)

  for (ell in seq_len(l_max)) {
    basis_ell <- matrix(
      (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1)),
      nrow = nrow(x_matrix),
      ncol = ncol(x_matrix)
    )
    out <- out + sweep(
      basis_ell,
      1L,
      coefficients[[ell + 1L]] * p_r[, ell + 1L],
      FUN = "*"
    )
  }

  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}

small_circle_projection_kernel <- function(c_thresholds, z_nodes, r) {
  c_thresholds <- pmin(pmax(as.numeric(c_thresholds), -1), 1)
  z_nodes <- pmin(pmax(as.numeric(z_nodes), -1), 1)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  one_minus_r2 <- pmax(0, 1 - r^2)
  a_values <- r * z_nodes
  b_values <- sqrt(one_minus_r2) * sqrt(pmax(0, 1 - z_nodes^2))
  kernel <- matrix(0, nrow = length(c_thresholds), ncol = length(z_nodes))

  for (j in seq_along(z_nodes)) {
    if (b_values[[j]] <= 1e-15) {
      kernel[, j] <- as.numeric(a_values[[j]] >= c_thresholds)
    } else {
      ratio <- (c_thresholds - a_values[[j]]) / b_values[[j]]
      kernel[, j] <- ifelse(
        ratio <= -1,
        1,
        ifelse(ratio >= 1, 0, acos(pmin(pmax(ratio, -1), 1)) / pi)
      )
    }
  }

  kernel
}

small_circle_distance_profile_integral <- function(omega,
                                                   t_values,
                                                   mu,
                                                   kappa,
                                                   nu,
                                                   distance_type = c("geodesic", "chordal"),
                                                   quad_n = 1000L) {
  distance_type <- match.arg(distance_type)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  if (length(omega) != 3L) {
    stop("The Small Circle implementation currently supports only S^2.")
  }

  t_values <- as.numeric(t_values)
  out <- numeric(length(t_values))
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  if (params$kappa <= 0) {
    out[active] <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_values[active])) / 2
    } else {
      (t_values[active]^2) / 4
    }
    return(out)
  }

  r_value <- sum(omega * params$mu)
  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  if (abs(r_value - 1) <= 1e-12) {
    out[active] <- 1 - small_circle_axis_cdf(thresholds, kappa = params$kappa, nu = params$nu)
    return(pmin(pmax(out, 0), 1))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- small_circle_axis_cdf(-thresholds, kappa = params$kappa, nu = params$nu)
    return(pmin(pmax(out, 0), 1))
  }

  quad <- small_circle_gauss_legendre(as.integer(quad_n))
  density_z <- small_circle_axis_density(quad$nodes, kappa = params$kappa, nu = params$nu)
  kernel <- small_circle_projection_kernel(
    c_thresholds = thresholds,
    z_nodes = quad$nodes,
    r = r_value
  )
  out[active] <- as.numeric(kernel %*% (quad$weights * density_z))
  pmin(pmax(out, 0), 1)
}

small_circle_monotone_clip <- function(t_values, values, upper_bound) {
  if (length(values) <= 1L) {
    return(pmin(pmax(values, 0), 1))
  }

  order_idx <- order(t_values)
  sorted_values <- pmin(pmax(values[order_idx], 0), 1)
  sorted_values <- cummax(sorted_values)
  sorted_values[t_values[order_idx] <= 0] <- 0
  sorted_values[t_values[order_idx] >= upper_bound] <- 1
  out <- numeric(length(sorted_values))
  out[order_idx] <- sorted_values
  out
}

small_circle_monotone_clip_dot <- function(dot_values, values) {
  if (length(values) <= 1L) {
    out <- pmin(pmax(values, 0), 1)
    out[dot_values >= 1] <- 0
    out[dot_values <= -1] <- 1
    return(out)
  }

  order_idx <- order(dot_values, decreasing = TRUE)
  sorted_values <- pmin(pmax(values[order_idx], 0), 1)
  sorted_values <- cummax(sorted_values)
  sorted_dots <- dot_values[order_idx]
  sorted_values[sorted_dots >= 1] <- 0
  sorted_values[sorted_dots <= -1] <- 1
  out <- numeric(length(sorted_values))
  out[order_idx] <- sorted_values
  out
}

small_circle_profile_matrix_legendre <- function(t_grid,
                                                 omega_grid,
                                                 mu,
                                                 coeffs,
                                                 l_max = length(coeffs) - 1L,
                                                 distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(omega_grid))
  coeffs <- as.numeric(coeffs)
  l_max <- as.integer(l_max)
  if (l_max != length(coeffs) - 1L) {
    stop("`l_max` must match `length(coeffs) - 1`.")
  }

  t_grid <- as.numeric(t_grid)
  n_omega <- nrow(omega_grid)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- matrix(0, nrow = n_omega, ncol = length(t_grid))

  if (length(t_grid) == 0L) {
    return(out)
  }

  active <- which(is.finite(t_grid) & t_grid > 0 & t_grid < upper_bound)
  out[, t_grid >= upper_bound] <- 1
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_grid[active], distance_type = distance_type)
  threshold_matrix <- matrix(thresholds, nrow = n_omega, ncol = length(active), byrow = TRUE)
  out[, active] <- 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = threshold_matrix,
    r = as.numeric(omega_grid %*% mu),
    coefficients = coeffs
  )
  out <- pmin(pmax(out, 0), 1)

  for (i in seq_len(n_omega)) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = t_grid,
      values = out[i, ],
      upper_bound = upper_bound
    )
  }
  out
}

small_circle_sample_profile_matrix_legendre <- function(X,
                                                        mu,
                                                        coeffs,
                                                        l_max = length(coeffs) - 1L) {
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(X))
  coeffs <- as.numeric(coeffs)
  l_max <- as.integer(l_max)
  if (l_max != length(coeffs) - 1L) {
    stop("`l_max` must match `length(coeffs) - 1`.")
  }

  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  out <- 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = as.numeric(X %*% mu),
    coefficients = coeffs
  )
  out <- pmin(pmax(out, 0), 1)

  for (i in seq_len(nrow(X))) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = acos(dot_products[i, ]),
      values = out[i, ],
      upper_bound = pi
    )
  }
  out
}

distance_profile_small_circle <- function(omega,
                                          t_values,
                                          mu,
                                          kappa,
                                          nu,
                                          distance_type = c("geodesic", "chordal"),
                                          method = c("legendre", "integral"),
                                          l_max = 150L,
                                          quad_n = 1000L,
                                          tol = 1e-10,
                                          validate_against_integral = FALSE,
                                          validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  t_values <- as.numeric(t_values)

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (ncol(omega) != 3L) {
      stop("`omega` must have three columns for Small Circle profiles on S^2.")
    }
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }

    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_small_circle(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = params$mu,
        kappa = params$kappa,
        nu = params$nu,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol,
        validate_against_integral = validate_against_integral,
        validation_tol = validation_tol
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  if (params$kappa <= 0) {
    out[active] <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_values[active])) / 2
    } else {
      (t_values[active]^2) / 4
    }
    return(out)
  }

  if (identical(method, "integral")) {
    return(small_circle_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = params$mu,
      kappa = params$kappa,
      nu = params$nu,
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  coeffs <- small_circle_legendre_coefficients(
    kappa = params$kappa,
    nu = params$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  cdf_values <- small_circle_projection_cdf_legendre(
    x = thresholds,
    r = sum(omega * params$mu),
    coefficients = coeffs
  )
  out[active] <- 1 - cdf_values
  out <- small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound)

  if (isTRUE(validate_against_integral)) {
    integral_values <- small_circle_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = params$mu,
      kappa = params$kappa,
      nu = params$nu,
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - integral_values))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Small Circle Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}

distance_profile_small_circle_grid <- function(omega_grid,
                                               mu,
                                               kappa,
                                               nu,
                                               t_grid,
                                               distance_type = c("geodesic", "chordal"),
                                               method = c("legendre", "integral"),
                                               l_max = 150L,
                                               quad_n = 1000L,
                                               tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  t_grid <- as.numeric(t_grid)

  if (params$kappa <= 0) {
    base_profile <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_grid)) / 2
    } else {
      (t_grid^2) / 4
    }
    return(matrix(base_profile, nrow = nrow(omega_grid), ncol = length(t_grid), byrow = TRUE))
  }

  if (identical(method, "integral")) {
    quad <- small_circle_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * small_circle_axis_density(
      quad$nodes,
      kappa = params$kappa,
      nu = params$nu
    )
    out <- projection_profile_matrix_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = params$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density,
      distance_type = distance_type
    )
    thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
    r_values <- as.numeric(omega_grid %*% params$mu)
    pos_idx <- which(abs(r_values - 1) <= 1e-12)
    if (length(pos_idx) > 0L) {
      out[pos_idx, ] <- 1 - small_circle_axis_cdf(
        thresholds,
        kappa = params$kappa,
        nu = params$nu
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- small_circle_axis_cdf(
        -thresholds,
        kappa = params$kappa,
        nu = params$nu
      )
    }
    return(out)
  }

  coeffs <- small_circle_legendre_coefficients(
    kappa = params$kappa,
    nu = params$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  small_circle_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = params$mu,
    coeffs = coeffs,
    l_max = l_max,
    distance_type = distance_type
  )
}

distance_profile_small_circle_cvm_grid <- function(X,
                                                   mu,
                                                   kappa,
                                                   nu,
                                                   method = c("legendre", "integral"),
                                                   l_max = 150L,
                                                   quad_n = 1000L,
                                                   tol = 1e-10) {
  method <- match.arg(method)
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)

  if (params$kappa <= 0) {
    dot_products <- pmin(pmax(X %*% t(X), -1), 1)
    return((1 - dot_products) / 2)
  }

  if (identical(method, "integral")) {
    quad <- small_circle_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * small_circle_axis_density(
      quad$nodes,
      kappa = params$kappa,
      nu = params$nu
    )
    return(projection_sample_profile_matrix_integral(
      X = X,
      mu = params$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    ))
  }

  coeffs <- small_circle_legendre_coefficients(
    kappa = params$kappa,
    nu = params$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  small_circle_sample_profile_matrix_legendre(
    X = X,
    mu = params$mu,
    coeffs = coeffs,
    l_max = l_max
  )
}

r_sph_small_circle <- function(n, mu, kappa, nu, check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  if (params$kappa <= 0) {
    x <- jp_uniform_sphere(n = n, ambient_dim = 3L)
    if (isTRUE(check)) {
      expect_norms <- sqrt(rowSums(x^2))
      if (max(abs(expect_norms - 1)) > 1e-8) {
        stop("Small Circle uniform sampler returned non-unit vectors.")
      }
    }
    return(x)
  }

  sigma <- 1 / sqrt(2 * params$kappa)
  lower_prob <- stats::pnorm((-1 - params$nu) / sigma)
  upper_prob <- stats::pnorm((1 - params$nu) / sigma)
  uniforms <- stats::runif(n)
  z <- params$nu + sigma * stats::qnorm(lower_prob + uniforms * (upper_prob - lower_prob))
  z <- pmin(pmax(z, -1), 1)

  phi <- stats::runif(n, min = 0, max = 2 * pi)
  basis <- jp_orthonormal_complement(params$mu)
  tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
  radial <- sqrt(pmax(0, 1 - z^2))
  x <- tcrossprod(z, params$mu) + sweep(tangent, 1, radial, "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Small Circle sampler returned non-unit vectors.")
    }
  }

  x
}

small_circle_weighted_loglik_s2 <- function(mu,
                                            kappa,
                                            nu,
                                            x,
                                            prob_weights = NULL) {
  params <- small_circle_validate_parameters(mu = mu, kappa = kappa, nu = nu)
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  projections <- as.numeric(x %*% params$mu)
  -small_circle_log_norm_constant(params$kappa, params$nu) -
    params$kappa * sum(prob_weights * (projections - params$nu)^2)
}

small_circle_start_theta_s2 <- function(x,
                                        weights = NULL,
                                        nu_min = 1e-6,
                                        nu_eps = 1e-6,
                                        kappa_min = 1e-8,
                                        kappa_max = 1e6) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  xbar <- colSums(x * prob_weights)
  weighted_x <- sweep(x, 1, sqrt(prob_weights), "*")
  s_matrix <- crossprod(weighted_x)
  sigma_matrix <- s_matrix - tcrossprod(xbar)
  eig <- eigen(sigma_matrix, symmetric = TRUE)
  mu0 <- eig$vectors[, which.min(eig$values)]
  nu0 <- sum(mu0 * xbar)
  canonical <- small_circle_canonicalize_theta(mu0, nu0)
  mu0 <- canonical$mu
  nu0 <- canonical$nu
  nu0 <- min(max(nu0, nu_min), 1 - nu_eps)

  t3 <- max(min(eig$values), .Machine$double.eps)
  kappa0 <- min(max(1 / (2 * t3), kappa_min), kappa_max)

  list(mu = mu0, kappa = kappa0, nu = nu0)
}

small_circle_logistic_bounded <- function(x, upper = 1 - 1e-6) {
  upper / (1 + exp(-x))
}

small_circle_inverse_logistic_bounded <- function(y, upper = 1 - 1e-6) {
  y <- min(max(y, 1e-8), upper - 1e-8)
  stats::qlogis(y / upper)
}

small_circle_mle_s2_weighted <- function(x,
                                         weights = NULL,
                                         control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  nu_eps <- as.numeric(control$small_circle_nu_eps %||% 1e-6)
  kappa_min <- as.numeric(control$small_circle_kappa_min %||% 1e-8)
  kappa_max <- as.numeric(control$small_circle_kappa_max %||% 1e6)
  nu_optim <- match.arg(
    as.character(control$small_circle_nu_optim %||% "logistic"),
    choices = c("logistic", "box")
  )
  theta_start <- control$small_circle_mle_start_theta %||% control$theta_start %||% control$jp_mle_start_theta %||% NULL

  if (is.null(theta_start)) {
    theta_start <- small_circle_start_theta_s2(
      x = x,
      weights = prob_weights,
      nu_min = nu_eps,
      nu_eps = nu_eps,
      kappa_min = kappa_min,
      kappa_max = kappa_max
    )
  } else {
    theta_start <- small_circle_validate_parameters(
      mu = theta_start$mu,
      kappa = theta_start$kappa,
      nu = theta_start$nu,
      allow_negative_nu = TRUE
    )
    canonical <- small_circle_canonicalize_theta(theta_start$mu, theta_start$nu)
    theta_start$mu <- canonical$mu
    theta_start$nu <- min(max(canonical$nu, nu_eps), 1 - nu_eps)
    theta_start$kappa <- min(max(theta_start$kappa, kappa_min), kappa_max)
  }

  objective <- function(par) {
    if (identical(nu_optim, "logistic")) {
      kappa_value <- min(max(log1p(exp(par[[1L]])), kappa_min), kappa_max)
      nu_value <- small_circle_logistic_bounded(par[[2L]], upper = 1 - nu_eps)
    } else {
      kappa_value <- min(max(par[[1L]], kappa_min), kappa_max)
      nu_value <- min(max(par[[2L]], 0), 1 - nu_eps)
    }

    mu_raw <- par[3:5]
    mu_norm <- sqrt(sum(mu_raw^2))
    if (!is.finite(mu_norm) || mu_norm <= 0) {
      return(.Machine$double.xmax / 100)
    }

    mu_value <- mu_raw / mu_norm
    value <- -small_circle_weighted_loglik_s2(
      mu = mu_value,
      kappa = kappa_value,
      nu = nu_value,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(value)) {
      return(.Machine$double.xmax / 100)
    }
    value
  }

  optim_control <- control$small_circle_optim_control %||% list(maxit = 500L, reltol = 1e-10)
  if (identical(nu_optim, "logistic")) {
    opt <- stats::optim(
      par = c(
        log(expm1(theta_start$kappa)),
        small_circle_inverse_logistic_bounded(theta_start$nu, upper = 1 - nu_eps),
        theta_start$mu
      ),
      fn = objective,
      method = control$small_circle_optim_method %||% "BFGS",
      control = optim_control
    )
    kappa_hat <- min(max(log1p(exp(opt$par[[1L]])), kappa_min), kappa_max)
    nu_hat <- small_circle_logistic_bounded(opt$par[[2L]], upper = 1 - nu_eps)
  } else {
    optim_method_box <- control$small_circle_optim_method %||% "L-BFGS-B"
    if (identical(optim_method_box, "L-BFGS-B") && !is.null(optim_control$reltol)) {
      reltol_value <- as.numeric(optim_control$reltol)
      if (length(reltol_value) == 1L && is.finite(reltol_value) && reltol_value > 0) {
        optim_control$factr <- reltol_value / .Machine$double.eps
      }
      optim_control$reltol <- NULL
    }
    lower_bounds <- c(kappa_min, 0, rep(-Inf, 3L))
    upper_bounds <- c(kappa_max, 1 - nu_eps, rep(Inf, 3L))
    opt <- stats::optim(
      par = c(theta_start$kappa, theta_start$nu, theta_start$mu),
      fn = objective,
      method = optim_method_box,
      lower = lower_bounds,
      upper = upper_bounds,
      control = optim_control
    )
    kappa_hat <- min(max(opt$par[[1L]], kappa_min), kappa_max)
    nu_hat <- min(max(opt$par[[2L]], 0), 1 - nu_eps)
  }

  mu_hat <- opt$par[3:5] / sqrt(sum(opt$par[3:5]^2))
  canonical <- small_circle_canonicalize_theta(mu_hat, nu_hat)

  list(
    mu = canonical$mu,
    kappa = kappa_hat,
    nu = canonical$nu,
    loglik = -opt$value,
    opt = opt,
    weighted_mle = TRUE,
    start_theta = theta_start,
    nu_optim = nu_optim
  )
}

small_circle_benchmark_nu_optimization <- function(n = 200L,
                                                   n_rep = 40L,
                                                   true_mu = c(0, 0, 1),
                                                   true_kappa = 8,
                                                   true_nu = 0.15,
                                                   seed = 1,
                                                   control_logistic = list(),
                                                   control_box = list()) {
  n <- as.integer(n)
  n_rep <- as.integer(n_rep)
  if (length(n) != 1L || !is.finite(n) || n < 5L) {
    stop("`n` must be an integer >= 5.")
  }
  if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1L) {
    stop("`n_rep` must be a strictly positive integer.")
  }

  true_params <- small_circle_validate_parameters(mu = true_mu, kappa = true_kappa, nu = true_nu)
  set.seed(as.integer(seed))

  evaluate_fit <- function(method_name,
                           x,
                           theta_start,
                           base_control) {
    control_method <- utils::modifyList(
      list(
        small_circle_nu_optim = method_name,
        small_circle_mle_start_theta = theta_start
      ),
      base_control
    )

    start_time <- proc.time()[[3L]]
    fit <- try(small_circle_mle_s2_weighted(x = x, control = control_method), silent = TRUE)
    elapsed <- proc.time()[[3L]] - start_time

    if (inherits(fit, "try-error")) {
      return(data.frame(
        method = method_name,
        elapsed_sec = elapsed,
        converged = FALSE,
        loglik = NA_real_,
        abs_kappa_error = NA_real_,
        abs_nu_error = NA_real_,
        mu_angle_error = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    dot_mu <- pmin(pmax(sum(fit$mu * true_params$mu), -1), 1)
    data.frame(
      method = method_name,
      elapsed_sec = elapsed,
      converged = is.list(fit$opt) && isTRUE(fit$opt$convergence == 0L),
      loglik = fit$loglik,
      abs_kappa_error = abs(fit$kappa - true_params$kappa),
      abs_nu_error = abs(fit$nu - true_params$nu),
      mu_angle_error = acos(dot_mu),
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(seq_len(n_rep), function(rep_id) {
    x <- r_sph_small_circle(
      n = n,
      mu = true_params$mu,
      kappa = true_params$kappa,
      nu = true_params$nu,
      check = FALSE
    )
    theta_start <- small_circle_start_theta_s2(x = x)

    row_logistic <- evaluate_fit(
      method_name = "logistic",
      x = x,
      theta_start = theta_start,
      base_control = control_logistic
    )
    row_box <- evaluate_fit(
      method_name = "box",
      x = x,
      theta_start = theta_start,
      base_control = control_box
    )

    cbind(rep = rep_id, rbind(row_logistic, row_box), stringsAsFactors = FALSE)
  })

  details <- do.call(rbind, rows)
  split_details <- split(details, details$method)

  summary <- do.call(rbind, lapply(names(split_details), function(method_name) {
    d <- split_details[[method_name]]
    converged_idx <- which(d$converged)
    d_conv <- if (length(converged_idx) == 0L) d[FALSE, , drop = FALSE] else d[converged_idx, , drop = FALSE]

    data.frame(
      method = method_name,
      n_rep = nrow(d),
      convergence_rate = mean(d$converged),
      mean_elapsed_sec = mean(d$elapsed_sec, na.rm = TRUE),
      median_elapsed_sec = stats::median(d$elapsed_sec, na.rm = TRUE),
      mean_loglik = if (nrow(d_conv) > 0L) mean(d_conv$loglik, na.rm = TRUE) else NA_real_,
      mean_abs_kappa_error = if (nrow(d_conv) > 0L) mean(d_conv$abs_kappa_error, na.rm = TRUE) else NA_real_,
      mean_abs_nu_error = if (nrow(d_conv) > 0L) mean(d_conv$abs_nu_error, na.rm = TRUE) else NA_real_,
      mean_mu_angle_error = if (nrow(d_conv) > 0L) mean(d_conv$mu_angle_error, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(
    setup = list(
      n = n,
      n_rep = n_rep,
      true_mu = true_params$mu,
      true_kappa = true_params$kappa,
      true_nu = true_params$nu,
      seed = seed
    ),
    summary = summary,
    details = details
  )
}

small_circle_compare_profile_methods <- function(mu,
                                                 kappa,
                                                 nu,
                                                 omega_list,
                                                 t_grid,
                                                 distance_type = c("geodesic", "chordal"),
                                                 l_max = 150L,
                                                 quad_n = 1000L,
                                                 tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  comparison_rows <- lapply(seq_along(omega_list), function(i) {
    omega <- omega_list[[i]]
    legendre <- distance_profile_small_circle(
      omega = omega,
      t_values = t_grid,
      mu = mu,
      kappa = kappa,
      nu = nu,
      distance_type = distance_type,
      method = "legendre",
      l_max = l_max,
      quad_n = quad_n,
      tol = tol
    )
    integral <- distance_profile_small_circle(
      omega = omega,
      t_values = t_grid,
      mu = mu,
      kappa = kappa,
      nu = nu,
      distance_type = distance_type,
      method = "integral",
      quad_n = quad_n,
      tol = tol
    )

    data.frame(
      omega_id = i,
      max_abs_diff = max(abs(legendre - integral)),
      mean_abs_diff = mean(abs(legendre - integral)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, comparison_rows)
}

rotational_logsumexp2 <- function(log_x, log_y) {
  m <- pmax(log_x, log_y)
  m + log(exp(log_x - m) + exp(log_y - m))
}

rotational_clamp_unit_interval <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

rotational_bounded_weight <- function(eta, weight_eps = 1e-6) {
  pmin(pmax(stats::plogis(eta), weight_eps), 1 - weight_eps)
}

rotational_positive_parameter <- function(log_value,
                                          lower = 1e-6,
                                          upper = 1e6) {
  pmin(pmax(exp(log_value), lower), upper)
}

rotational_weighted_mean <- function(x, weights) {
  sum(as.numeric(weights) * as.numeric(x))
}

rotational_weighted_variance <- function(x, weights, center = NULL) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  if (is.null(center)) {
    center <- rotational_weighted_mean(x, weights)
  }
  sum(weights * (x - center)^2)
}

rotational_weighted_quantile <- function(x, weights, prob) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  prob <- as.numeric(prob)

  if (length(prob) != 1L || !is.finite(prob) || prob < 0 || prob > 1) {
    stop("`prob` must be a finite scalar in [0, 1].")
  }

  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- weights[ord]
  cum_weights <- cumsum(w_sorted)
  total <- cum_weights[[length(cum_weights)]]
  x_sorted[[which(cum_weights >= prob * total)[[1L]]]]
}

rotational_weighted_covariance_matrix <- function(x, weights) {
  x <- as.matrix(x)
  weights <- as.numeric(weights)
  center <- colSums(x * weights)
  centered <- sweep(x, 2L, center, FUN = "-")
  crossprod(sweep(centered, 1L, sqrt(weights), FUN = "*"))
}

rotational_unit_vector_fallback <- function(x,
                                            fallback = c(0, 0, 1),
                                            tol = 1e-12) {
  x <- as.numeric(x)
  norm_x <- sqrt(sum(x^2))
  if (!is.finite(norm_x) || norm_x <= tol) {
    fallback <- as.numeric(fallback)
    fallback / sqrt(sum(fallback^2))
  } else {
    x / norm_x
  }
}

rotational_unique_mu_candidates <- function(mu_candidates,
                                            tol = 1e-8) {
  out <- list()
  for (mu in mu_candidates) {
    mu_vec <- rotational_unit_vector_fallback(mu)
    keep <- TRUE
    if (length(out) > 0L) {
      for (existing in out) {
        if (max(abs(mu_vec - existing)) <= tol || max(abs(mu_vec + existing)) <= tol) {
          keep <- FALSE
          break
        }
      }
    }
    if (keep) {
      out[[length(out) + 1L]] <- mu_vec
    }
  }
  out
}

rotational_gauss_legendre <- function(n) {
  small_circle_gauss_legendre(n)
}

rotational_gauss_hermite_cache <- new.env(parent = emptyenv())

rotational_gauss_hermite <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  key <- as.character(n)
  if (exists(key, envir = rotational_gauss_hermite_cache, inherits = FALSE)) {
    return(get(key, envir = rotational_gauss_hermite_cache, inherits = FALSE))
  }

  off_diag <- sqrt(seq_len(n - 1L) / 2)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off_diag
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off_diag
  eig <- eigen(jacobi, symmetric = TRUE)
  order_idx <- order(eig$values)

  out <- list(
    nodes = eig$values[order_idx],
    weights = sqrt(pi) * (eig$vectors[1L, order_idx]^2)
  )
  assign(key, out, envir = rotational_gauss_hermite_cache)
  out
}

rotational_legendre_matrix <- function(x, l_max) {
  small_circle_legendre_matrix(x, l_max)
}

rotational_legendre_coefficients <- function(density_h,
                                             Lmax,
                                             quad_n = 1000L,
                                             tol = 1e-10) {
  if (!is.function(density_h)) {
    stop("`density_h` must be a function.")
  }

  Lmax <- as.integer(Lmax)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(Lmax) != 1L || !is.finite(Lmax) || Lmax < 0L) {
    stop("`Lmax` must be a nonnegative integer.")
  }

  quad <- rotational_gauss_legendre(quad_n)
  legendre_matrix <- rotational_legendre_matrix(quad$nodes, l_max = Lmax)
  h_values <- as.numeric(density_h(quad$nodes))
  if (length(h_values) != length(quad$nodes) || any(!is.finite(h_values)) || any(h_values < 0)) {
    stop("`density_h` must return finite nonnegative values of the same length as its input.")
  }

  raw_moments <- as.numeric(crossprod(legendre_matrix, quad$weights * h_values))
  ell <- 0:Lmax
  coeffs <- ((2 * ell + 1) / 2) * raw_moments
  a0_error <- abs(raw_moments[[1L]] / 2 - 1)
  coeffs[[1L]] <- 1

  if (a0_error > tol) {
    stop(sprintf("Rotational Legendre coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error
  )
}

rotational_projection_cdf_legendre <- function(x,
                                               r,
                                               coefficients,
                                               enforce_bounds = TRUE) {
  x <- pmin(pmax(as.numeric(x), -1), 1)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != 1L || !is.finite(r)) {
    stop("`r` must be a finite scalar.")
  }

  out <- (x + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- as.numeric(rotational_legendre_matrix(r, l_max = l_max))
  p_x <- rotational_legendre_matrix(x, l_max = l_max + 1L)
  basis_matrix <- matrix(0, nrow = length(x), ncol = l_max)
  for (ell in seq_len(l_max)) {
    basis_matrix[, ell] <- (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1))
  }

  out <- out + as.numeric(basis_matrix %*% (coefficients[-1L] * p_r[-1L]))
  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}

rotational_projection_cdf_legendre_matrix <- function(x_matrix,
                                                      r,
                                                      coefficients,
                                                      enforce_bounds = TRUE) {
  x_matrix <- as.matrix(x_matrix)
  r <- pmin(pmax(as.numeric(r), -1), 1)
  coefficients <- as.numeric(coefficients)
  l_max <- length(coefficients) - 1L

  if (length(r) != nrow(x_matrix) || any(!is.finite(r))) {
    stop("`r` must be a finite vector of length nrow(`x_matrix`).")
  }

  x_matrix <- pmin(pmax(x_matrix, -1), 1)
  out <- (x_matrix + 1) / 2
  if (l_max < 1L) {
    return(out)
  }

  p_r <- rotational_legendre_matrix(r, l_max = l_max)
  p_x <- rotational_legendre_matrix(as.numeric(x_matrix), l_max = l_max + 1L)

  for (ell in seq_len(l_max)) {
    basis_ell <- matrix(
      (p_x[, ell + 2L] - p_x[, ell]) / (2 * (2 * ell + 1)),
      nrow = nrow(x_matrix),
      ncol = ncol(x_matrix)
    )
    out <- out + sweep(
      basis_ell,
      1L,
      coefficients[[ell + 1L]] * p_r[, ell + 1L],
      FUN = "*"
    )
  }

  if (isTRUE(enforce_bounds)) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}

rotational_profile_legendre <- function(t,
                                        omega,
                                        mu,
                                        coeffs,
                                        Lmax = length(coeffs) - 1L,
                                        distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = length(omega))
  coeffs <- as.numeric(coeffs)
  Lmax <- as.integer(Lmax)
  if (Lmax != length(coeffs) - 1L) {
    stop("`Lmax` must match `length(coeffs) - 1`.")
  }

  t <- as.numeric(t)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t))
  out[t <= 0] <- 0
  out[t >= upper_bound] <- 1
  active <- which(is.finite(t) & t > 0 & t < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t[active], distance_type = distance_type)
  out[active] <- 1 - rotational_projection_cdf_legendre(
    x = thresholds,
    r = sum(omega * mu),
    coefficients = coeffs
  )
  small_circle_monotone_clip(t_values = t, values = out, upper_bound = upper_bound)
}

rotational_profile_matrix_legendre <- function(t_grid,
                                               omega_grid,
                                               mu,
                                               coeffs,
                                               Lmax = length(coeffs) - 1L,
                                               distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(omega_grid))
  coeffs <- as.numeric(coeffs)
  Lmax <- as.integer(Lmax)
  if (Lmax != length(coeffs) - 1L) {
    stop("`Lmax` must match `length(coeffs) - 1`.")
  }

  t_grid <- as.numeric(t_grid)
  n_omega <- nrow(omega_grid)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- matrix(0, nrow = n_omega, ncol = length(t_grid))

  if (length(t_grid) == 0L) {
    return(out)
  }

  active <- which(is.finite(t_grid) & t_grid > 0 & t_grid < upper_bound)
  out[, t_grid >= upper_bound] <- 1
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_grid[active], distance_type = distance_type)
  p_r <- rotational_legendre_matrix(as.numeric(omega_grid %*% mu), l_max = Lmax)
  if (Lmax == 0L) {
    out[, active] <- matrix((1 - thresholds) / 2, nrow = n_omega, ncol = length(active), byrow = TRUE)
    return(out)
  }

  p_x <- rotational_legendre_matrix(thresholds, l_max = Lmax + 1L)
  basis <- matrix(0, nrow = length(active), ncol = Lmax)
  for (ell in seq_len(Lmax)) {
    basis[, ell] <- (p_x[, ell] - p_x[, ell + 2L]) / (2 * (2 * ell + 1))
  }

  out[, active] <- matrix((1 - thresholds) / 2, nrow = n_omega, ncol = length(active), byrow = TRUE) +
    p_r[, -1L, drop = FALSE] %*% t(sweep(basis, 2L, coeffs[-1L], FUN = "*"))
  out <- pmin(pmax(out, 0), 1)

  for (i in seq_len(n_omega)) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = t_grid,
      values = out[i, ],
      upper_bound = upper_bound
    )
  }
  out
}

projection_profile_matrix_integral <- function(omega_grid,
                                               t_grid,
                                               mu,
                                               z_nodes,
                                               weighted_density,
                                               distance_type = c("geodesic", "chordal")) {
  distance_type <- match.arg(distance_type)
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(omega_grid))
  t_grid <- as.numeric(t_grid)
  z_nodes <- as.numeric(z_nodes)
  weighted_density <- as.numeric(weighted_density)

  if (length(z_nodes) != length(weighted_density)) {
    stop("`z_nodes` and `weighted_density` must have the same length.")
  }

  n_omega <- nrow(omega_grid)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- matrix(0, nrow = n_omega, ncol = length(t_grid))

  if (length(t_grid) == 0L) {
    return(out)
  }

  active <- which(is.finite(t_grid) & t_grid > 0 & t_grid < upper_bound)
  out[, t_grid >= upper_bound] <- 1
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_grid[active], distance_type = distance_type)
  r_values <- as.numeric(omega_grid %*% mu)

  for (i in seq_len(n_omega)) {
    kernel <- small_circle_projection_kernel(
      c_thresholds = thresholds,
      z_nodes = z_nodes,
      r = r_values[[i]]
    )
    out[i, active] <- as.numeric(kernel %*% weighted_density)
  }

  pmin(pmax(out, 0), 1)
}

projection_sample_profile_matrix_integral <- function(X,
                                                      mu,
                                                      z_nodes,
                                                      weighted_density) {
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = ncol(X))
  z_nodes <- as.numeric(z_nodes)
  weighted_density <- as.numeric(weighted_density)

  if (length(z_nodes) != length(weighted_density)) {
    stop("`z_nodes` and `weighted_density` must have the same length.")
  }

  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  r_values <- as.numeric(X %*% mu)
  out <- matrix(0, nrow = nrow(X), ncol = nrow(X))

  for (i in seq_len(nrow(X))) {
    kernel <- small_circle_projection_kernel(
      c_thresholds = dot_products[i, ],
      z_nodes = z_nodes,
      r = r_values[[i]]
    )
    out[i, ] <- as.numeric(kernel %*% weighted_density)
  }

  pmin(pmax(out, 0), 1)
}

rotational_distance_profile_integral <- function(omega,
                                                 t_values,
                                                 mu,
                                                 density_gz,
                                                 distance_type = c("geodesic", "chordal"),
                                                 quad_n = 1000L) {
  distance_type <- match.arg(distance_type)
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = length(omega))
  if (!is.function(density_gz)) {
    stop("`density_gz` must be a function.")
  }

  t_values <- as.numeric(t_values)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  quad <- rotational_gauss_legendre(as.integer(quad_n))
  density_z <- as.numeric(density_gz(quad$nodes))
  if (length(density_z) != length(quad$nodes) || any(!is.finite(density_z)) || any(density_z < 0)) {
    stop("`density_gz` must return finite nonnegative values of the same length as its input.")
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  kernel <- small_circle_projection_kernel(
    c_thresholds = thresholds,
    z_nodes = quad$nodes,
    r = sum(omega * mu)
  )
  out[active] <- as.numeric(kernel %*% (quad$weights * density_z))
  pmin(pmax(out, 0), 1)
}

rotational_sample_profile_matrix <- function(X,
                                             mu,
                                             profile_fun,
                                             ...) {
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)
  t(vapply(seq_len(nrow(X)), function(i) {
    as.numeric(profile_fun(
      omega = X[i, ],
      t_values = acos(dot_products[i, ]),
      mu = mu,
      ...
    ))
  }, numeric(nrow(X))))
}

debug_memory_log <- function(control = list(),
                             label,
                             objects = list(),
                             force = FALSE) {
  enabled <- isTRUE(force) || isTRUE(control$cvm_memory_debug %||% FALSE)
  if (!enabled) {
    return(invisible(NULL))
  }

  gc_info <- tryCatch(gc(), error = function(e) NULL)
  gc_summary <- if (is.null(gc_info)) {
    "gc=unavailable"
  } else {
    used_col <- grep("^used", colnames(gc_info), value = TRUE)[1L]
    trigger_col <- grep("^gc trigger", colnames(gc_info), value = TRUE)[1L]
    paste(
      sprintf(
        "%s used=%sMb gc_trigger=%sMb",
        rownames(gc_info),
        format(round(gc_info[, used_col], 1), trim = TRUE),
        format(round(gc_info[, trigger_col], 1), trim = TRUE)
      ),
      collapse = " | "
    )
  }

  message(sprintf("[CvM debug] %s | %s", label, gc_summary))
  if (length(objects) > 0L) {
    for (object_name in names(objects)) {
      obj <- objects[[object_name]]
      dims <- dim(obj)
      dim_text <- if (is.null(dims)) {
        sprintf("length=%d", length(obj))
      } else {
        sprintf("dim=%s", paste(dims, collapse = "x"))
      }
      size_mb <- as.numeric(utils::object.size(obj)) / 1024^2
      message(sprintf(
        "[CvM debug]   %s: %s, size=%.2f Mb, class=%s",
        object_name,
        dim_text,
        size_mb,
        paste(class(obj), collapse = "/")
      ))
    }
  }

  invisible(NULL)
}

small_circle_symmetric_mixture2_validate_parameters <- function(mu, kappa, nu) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Symmetric small-circle mixture utilities currently support only S^2.")
  }

  kappa <- as.numeric(kappa)
  nu <- as.numeric(nu)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be a finite scalar in [0, Inf).")
  }
  if (length(nu) != 1L || !is.finite(nu) || nu < 0 || nu >= 1) {
    stop("`nu` must be a finite scalar in [0, 1).")
  }

  list(mu = mu, kappa = kappa, nu = nu, ambient_dim = 3L)
}

small_circle_symmetric_mixture2_canonicalize_theta <- function(theta, tol = 1e-12) {
  if (!is.list(theta)) {
    stop("Symmetric small-circle-mixture theta must be a list.")
  }

  mu <- jp_normalize_unit_vector(theta$mu, arg_name = "`theta$mu`", min_length = 3L)
  kappa <- as.numeric(theta$kappa)
  nu <- abs(as.numeric(theta$nu))
  tol <- as.numeric(tol)

  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("`theta$kappa` must be a finite scalar in [0, Inf).")
  }
  if (length(nu) != 1L || !is.finite(nu) || nu >= 1) {
    stop("`theta$nu` must be a finite scalar with absolute value in [0, 1).")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a finite nonnegative scalar.")
  }

  first_nonzero <- which(abs(mu) > tol)[1L]
  if (!is.na(first_nonzero) && mu[[first_nonzero]] < 0) {
    mu <- -mu
  }
  if (nu <= tol) {
    nu <- 0
  }

  list(mu = mu, kappa = kappa, nu = nu, ambient_dim = length(mu))
}

small_circle_symmetric_mixture2_normalize_theta <- function(theta,
                                                            ambient_dim = 3L,
                                                            tol = 1e-12) {
  params <- small_circle_symmetric_mixture2_canonicalize_theta(theta, tol = tol)
  if (params$ambient_dim != ambient_dim) {
    stop("Symmetric small-circle-mixture theta has incompatible ambient dimension.")
  }
  params
}

small_circle_weighted_mixture2_validate_parameters <- function(mu,
                                                               pi,
                                                               kappa1,
                                                               nu1,
                                                               kappa2,
                                                               nu2) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Weighted small-circle mixture utilities currently support only S^2.")
  }

  pi <- as.numeric(pi)
  kappa1 <- as.numeric(kappa1)
  nu1 <- as.numeric(nu1)
  kappa2 <- as.numeric(kappa2)
  nu2 <- as.numeric(nu2)

  if (length(pi) != 1L || !is.finite(pi) || pi <= 0 || pi >= 1) {
    stop("`pi` must be a finite scalar in (0, 1).")
  }
  if (length(kappa1) != 1L || !is.finite(kappa1) || kappa1 < 0) {
    stop("`kappa1` must be a finite scalar in [0, Inf).")
  }
  if (length(nu1) != 1L || !is.finite(nu1) || nu1 < 0 || nu1 >= 1) {
    stop("`nu1` must be a finite scalar in [0, 1).")
  }
  if (length(kappa2) != 1L || !is.finite(kappa2) || kappa2 < 0) {
    stop("`kappa2` must be a finite scalar in [0, Inf).")
  }
  if (length(nu2) != 1L || !is.finite(nu2) || nu2 < 0 || nu2 >= 1) {
    stop("`nu2` must be a finite scalar in [0, 1).")
  }

  list(
    mu = mu,
    pi = pi,
    kappa1 = kappa1,
    nu1 = nu1,
    kappa2 = kappa2,
    nu2 = nu2,
    ambient_dim = 3L
  )
}

small_circle_weighted_mixture2_canonicalize_theta <- function(theta, tol = 1e-12) {
  if (!is.list(theta)) {
    stop("Weighted small-circle-mixture theta must be a list.")
  }

  params <- small_circle_weighted_mixture2_validate_parameters(
    mu = theta$mu,
    pi = theta$pi,
    kappa1 = theta$kappa1,
    nu1 = theta$nu1,
    kappa2 = theta$kappa2,
    nu2 = theta$nu2
  )
  tol <- as.numeric(tol)
  if (length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a finite nonnegative scalar.")
  }

  mu <- params$mu
  do_flip <- FALSE
  if (mu[[3L]] < -tol) {
    do_flip <- TRUE
  } else if (abs(mu[[3L]]) <= tol) {
    first_nonzero <- which(abs(mu) > tol)[1L]
    if (!is.na(first_nonzero) && mu[[first_nonzero]] < 0) {
      do_flip <- TRUE
    }
  }

  if (!do_flip) {
    return(params)
  }

  list(
    mu = -mu,
    pi = 1 - params$pi,
    kappa1 = params$kappa2,
    nu1 = params$nu2,
    kappa2 = params$kappa1,
    nu2 = params$nu1,
    ambient_dim = params$ambient_dim
  )
}

small_circle_weighted_mixture2_normalize_theta <- function(theta,
                                                           ambient_dim = 3L,
                                                           tol = 1e-12) {
  params <- small_circle_weighted_mixture2_canonicalize_theta(theta, tol = tol)
  if (params$ambient_dim != ambient_dim) {
    stop("Weighted small-circle-mixture theta has incompatible ambient dimension.")
  }
  params
}

d_sph_small_circle_weighted_mixture2_s2 <- function(x,
                                                    mu,
                                                    pi,
                                                    kappa1,
                                                    nu1,
                                                    kappa2,
                                                    nu2,
                                                    log = FALSE) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- small_circle_weighted_mixture2_validate_parameters(
    mu = mu,
    pi = pi,
    kappa1 = kappa1,
    nu1 = nu1,
    kappa2 = kappa2,
    nu2 = nu2
  )

  z <- pmin(pmax(as.numeric(x %*% params$mu), -1), 1)
  log_density <- -log(2 * base::pi) + rotational_logsumexp2(
    log(params$pi) + small_circle_axis_density(z, kappa = params$kappa1, nu = params$nu1, log = TRUE),
    log1p(-params$pi) + small_circle_axis_density(-z, kappa = params$kappa2, nu = params$nu2, log = TRUE)
  )
  if (log) log_density else exp(log_density)
}

small_circle_weighted_mixture2_weighted_loglik_s2 <- function(mu,
                                                              pi,
                                                              kappa1,
                                                              nu1,
                                                              kappa2,
                                                              nu2,
                                                              x,
                                                              prob_weights = NULL) {
  params <- small_circle_weighted_mixture2_validate_parameters(
    mu = mu,
    pi = pi,
    kappa1 = kappa1,
    nu1 = nu1,
    kappa2 = kappa2,
    nu2 = nu2
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  sum(prob_weights * d_sph_small_circle_weighted_mixture2_s2(
    x = x,
    mu = params$mu,
    pi = params$pi,
    kappa1 = params$kappa1,
    nu1 = params$nu1,
    kappa2 = params$kappa2,
    nu2 = params$nu2,
    log = TRUE
  ))
}

small_circle_weighted_mixture2_start_thetas_s2 <- function(x,
                                                           weights = NULL,
                                                           control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  nu_eps <- as.numeric(control$small_circle_weighted_mixture2_nu_eps %||% 1e-6)
  kappa_min <- as.numeric(control$small_circle_weighted_mixture2_kappa_min %||% 1e-8)
  kappa_max <- as.numeric(control$small_circle_weighted_mixture2_kappa_max %||% 1e6)

  warm_start <- control$small_circle_weighted_mixture2_start_theta %||%
    control$theta_start %||%
    control$small_circle_symmetric_mixture2_start_theta %||%
    NULL
  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- small_circle_weighted_mixture2_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$small_circle_weighted_mixture2_warm_start_only %||% FALSE)) {
      return(list(warm_start))
    }
  }

  symmetric_start <- control$small_circle_weighted_mixture2_symmetric_start_theta %||% NULL
  if (!is.null(symmetric_start)) {
    symmetric_start <- small_circle_symmetric_mixture2_normalize_theta(symmetric_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- small_circle_weighted_mixture2_canonicalize_theta(list(
      mu = symmetric_start$mu,
      pi = 0.5,
      kappa1 = symmetric_start$kappa,
      nu1 = symmetric_start$nu,
      kappa2 = symmetric_start$kappa,
      nu2 = symmetric_start$nu
    ))
  }

  resultant <- colSums(x * prob_weights)
  weighted_x <- sweep(x, 1L, sqrt(prob_weights), FUN = "*")
  second_moment <- crossprod(weighted_x)
  eig <- eigen(second_moment, symmetric = TRUE)
  mu_candidates <- rotational_unique_mu_candidates(list(
    resultant,
    -resultant,
    eig$vectors[, which.max(eig$values)],
    -eig$vectors[, which.max(eig$values)],
    eig$vectors[, which.min(eig$values)],
    -eig$vectors[, which.min(eig$values)],
    c(0, 0, 1),
    c(0, 0, -1)
  ))

  derive_component_start <- function(values, value_weights) {
    if (length(values) == 0L) {
      return(list(kappa = 10, nu = 0.5))
    }
    normalized_component_weights <- value_weights / sum(value_weights)
    nu0 <- min(max(rotational_weighted_mean(values, normalized_component_weights), nu_eps), 1 - nu_eps)
    var0 <- rotational_weighted_variance(values, normalized_component_weights, center = nu0)
    kappa0 <- min(max(1 / (2 * max(var0, 1 / (2 * kappa_max))), kappa_min), kappa_max)
    list(kappa = kappa0, nu = nu0)
  }

  for (mu0 in mu_candidates) {
    projection <- as.numeric(x %*% mu0)
    north_idx <- projection >= 0
    south_idx <- !north_idx
    if (!any(north_idx) || !any(south_idx)) {
      next
    }

    pi0 <- sum(prob_weights[north_idx])
    if (!is.finite(pi0) || pi0 <= 1e-6 || pi0 >= 1 - 1e-6) {
      next
    }

    north_start <- derive_component_start(projection[north_idx], prob_weights[north_idx])
    south_start <- derive_component_start(-projection[south_idx], prob_weights[south_idx])
    out[[length(out) + 1L]] <- small_circle_weighted_mixture2_canonicalize_theta(list(
      mu = mu0,
      pi = pi0,
      kappa1 = north_start$kappa,
      nu1 = north_start$nu,
      kappa2 = south_start$kappa,
      nu2 = south_start$nu
    ))
  }

  if (length(out) == 0L) {
    out[[1L]] <- small_circle_weighted_mixture2_canonicalize_theta(list(
      mu = c(0, 0, 1),
      pi = 0.5,
      kappa1 = 10,
      nu1 = 0.5,
      kappa2 = 10,
      nu2 = 0.5
    ))
  }

  out
}

small_circle_weighted_mixture2_pack_par <- function(theta,
                                                   control = list()) {
  theta <- small_circle_weighted_mixture2_normalize_theta(theta, ambient_dim = 3L)
  nu_upper <- 1 - as.numeric(control$small_circle_weighted_mixture2_nu_eps %||% 1e-6)
  list(
    par = c(
      stats::qlogis(theta$pi),
      log(pmax(expm1(theta$kappa1), .Machine$double.eps)),
      small_circle_inverse_logistic_bounded(theta$nu1, upper = nu_upper),
      log(pmax(expm1(theta$kappa2), .Machine$double.eps)),
      small_circle_inverse_logistic_bounded(theta$nu2, upper = nu_upper),
      theta$mu
    ),
    theta = theta
  )
}

small_circle_weighted_mixture2_unpack_par <- function(par,
                                                     control = list()) {
  kappa_min <- as.numeric(control$small_circle_weighted_mixture2_kappa_min %||% 1e-8)
  kappa_max <- as.numeric(control$small_circle_weighted_mixture2_kappa_max %||% 1e6)
  nu_eps <- as.numeric(control$small_circle_weighted_mixture2_nu_eps %||% 1e-6)
  weight_eps <- as.numeric(control$small_circle_weighted_mixture2_weight_eps %||% 1e-6)

  mu_hat <- rotational_unit_vector_fallback(par[6:8])
  small_circle_weighted_mixture2_canonicalize_theta(list(
    mu = mu_hat,
    pi = rotational_bounded_weight(par[[1L]], weight_eps = weight_eps),
    kappa1 = min(max(log1p(exp(par[[2L]])), kappa_min), kappa_max),
    nu1 = small_circle_logistic_bounded(par[[3L]], upper = 1 - nu_eps),
    kappa2 = min(max(log1p(exp(par[[4L]])), kappa_min), kappa_max),
    nu2 = small_circle_logistic_bounded(par[[5L]], upper = 1 - nu_eps)
  ))
}

small_circle_weighted_mixture2_mle_s2_weighted <- function(x,
                                                           weights = NULL,
                                                           control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- small_circle_weighted_mixture2_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(
    length(candidate_thetas),
    as.integer(control$small_circle_weighted_mixture2_n_starts %||% 12L)
  ))]

  objective <- function(par) {
    theta <- small_circle_weighted_mixture2_unpack_par(par, control = control)
    value <- -small_circle_weighted_mixture2_weighted_loglik_s2(
      mu = theta$mu,
      pi = theta$pi,
      kappa1 = theta$kappa1,
      nu1 = theta$nu1,
      kappa2 = theta$kappa2,
      nu2 = theta$nu2,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$small_circle_weighted_mixture2_optim_method %||% "BFGS"
  optim_control <- control$small_circle_weighted_mixture2_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- small_circle_weighted_mixture2_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- small_circle_weighted_mixture2_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$small_circle_weighted_mixture2_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$small_circle_weighted_mixture2_warm_start_only <- FALSE
    return(small_circle_weighted_mixture2_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Weighted small-circle-mixture weighted MLE failed for all starting values.")
  }

  c(
    best$theta,
    list(
      loglik = -best$opt$value,
      opt = best$opt,
      weighted_mle = TRUE,
      start_theta = best$start_theta
    )
  )
}

distance_profile_small_circle_weighted_mixture2 <- function(omega,
                                                           t_values,
                                                           mu,
                                                           pi,
                                                           kappa1,
                                                           nu1,
                                                           kappa2,
                                                           nu2,
                                                           distance_type = c("geodesic", "chordal"),
                                                           method = c("legendre", "integral"),
                                                           l_max = 150L,
                                                           quad_n = 1000L,
                                                           tol = 1e-10,
                                                           validate_against_integral = FALSE,
                                                           validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- small_circle_weighted_mixture2_validate_parameters(
    mu = mu,
    pi = pi,
    kappa1 = kappa1,
    nu1 = nu1,
    kappa2 = kappa2,
    nu2 = nu2
  )

  theta$pi * distance_profile_small_circle(
    omega = omega,
    t_values = t_values,
    mu = theta$mu,
    kappa = theta$kappa1,
    nu = theta$nu1,
    distance_type = distance_type,
    method = method,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    validate_against_integral = validate_against_integral,
    validation_tol = validation_tol
  ) + (1 - theta$pi) * distance_profile_small_circle(
    omega = omega,
    t_values = t_values,
    mu = -theta$mu,
    kappa = theta$kappa2,
    nu = theta$nu2,
    distance_type = distance_type,
    method = method,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    validate_against_integral = validate_against_integral,
    validation_tol = validation_tol
  )
}

distance_profile_small_circle_weighted_mixture2_grid <- function(omega_grid,
                                                                mu,
                                                                pi,
                                                                kappa1,
                                                                nu1,
                                                                kappa2,
                                                                nu2,
                                                                t_grid,
                                                                distance_type = c("geodesic", "chordal"),
                                                                method = c("legendre", "integral"),
                                                                l_max = 150L,
                                                                quad_n = 1000L,
                                                                tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- small_circle_weighted_mixture2_validate_parameters(
    mu = mu,
    pi = pi,
    kappa1 = kappa1,
    nu1 = nu1,
    kappa2 = kappa2,
    nu2 = nu2
  )

  theta$pi * distance_profile_small_circle_grid(
    omega_grid = omega_grid,
    mu = theta$mu,
    kappa = theta$kappa1,
    nu = theta$nu1,
    t_grid = t_grid,
    distance_type = distance_type,
    method = method,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  ) + (1 - theta$pi) * distance_profile_small_circle_grid(
    omega_grid = omega_grid,
    mu = -theta$mu,
    kappa = theta$kappa2,
    nu = theta$nu2,
    t_grid = t_grid,
    distance_type = distance_type,
    method = method,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )
}

r_sph_small_circle_weighted_mixture2 <- function(n,
                                                 mu,
                                                 pi,
                                                 kappa1,
                                                 nu1,
                                                 kappa2,
                                                 nu2,
                                                 check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  params <- small_circle_weighted_mixture2_validate_parameters(
    mu = mu,
    pi = pi,
    kappa1 = kappa1,
    nu1 = nu1,
    kappa2 = kappa2,
    nu2 = nu2
  )
  component_one <- stats::runif(n) <= params$pi
  x <- matrix(0, nrow = n, ncol = 3L)
  if (any(component_one)) {
    x[component_one, ] <- r_sph_small_circle(
      n = sum(component_one),
      mu = params$mu,
      kappa = params$kappa1,
      nu = params$nu1,
      check = check
    )
  }
  if (any(!component_one)) {
    x[!component_one, ] <- r_sph_small_circle(
      n = sum(!component_one),
      mu = -params$mu,
      kappa = params$kappa2,
      nu = params$nu2,
      check = check
    )
  }

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Weighted small-circle-mixture sampler returned non-unit vectors.")
    }
  }

  x
}

small_circle_symmetric_mixture2_axis_density <- function(z, kappa, nu) {
  params <- small_circle_symmetric_mixture2_validate_parameters(
    mu = c(0, 0, 1),
    kappa = kappa,
    nu = nu
  )
  z <- as.numeric(z)
  out <- numeric(length(z))
  valid <- is.finite(z) & z >= -1 & z <= 1
  if (!any(valid)) {
    return(out)
  }

  out[valid] <- 0.5 * small_circle_axis_density(
    z = z[valid],
    kappa = params$kappa,
    nu = params$nu
  ) + 0.5 * small_circle_axis_density(
    z = -z[valid],
    kappa = params$kappa,
    nu = params$nu
  )
  out
}

small_circle_symmetric_mixture2_axis_cdf <- function(z, kappa, nu) {
  params <- small_circle_symmetric_mixture2_validate_parameters(
    mu = c(0, 0, 1),
    kappa = kappa,
    nu = nu
  )
  z <- as.numeric(z)
  out <- numeric(length(z))
  out[z <= -1] <- 0
  out[z >= 1] <- 1
  active <- which(is.finite(z) & z > -1 & z < 1)
  if (length(active) == 0L) {
    return(out)
  }

  out[active] <- 0.5 * small_circle_axis_cdf(
    z = z[active],
    kappa = params$kappa,
    nu = params$nu
  ) + 0.5 * (
    1 - small_circle_axis_cdf(
      z = -z[active],
      kappa = params$kappa,
      nu = params$nu
    )
  )
  pmin(pmax(out, 0), 1)
}

d_sph_small_circle_symmetric_mixture2_s2 <- function(x,
                                                     mu,
                                                     kappa,
                                                     nu,
                                                     log = FALSE) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- small_circle_symmetric_mixture2_validate_parameters(
    mu = mu,
    kappa = kappa,
    nu = nu
  )

  z <- pmin(pmax(as.numeric(x %*% params$mu), -1), 1)
  log_density <- -log(2 * pi) + rotational_logsumexp2(
    log(0.5) + small_circle_axis_density(z, kappa = params$kappa, nu = params$nu, log = TRUE),
    log(0.5) + small_circle_axis_density(-z, kappa = params$kappa, nu = params$nu, log = TRUE)
  )
  if (log) log_density else exp(log_density)
}

small_circle_symmetric_mixture2_weighted_loglik_s2 <- function(mu,
                                                               kappa,
                                                               nu,
                                                               x,
                                                               prob_weights = NULL) {
  params <- small_circle_symmetric_mixture2_validate_parameters(
    mu = mu,
    kappa = kappa,
    nu = nu
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  sum(prob_weights * d_sph_small_circle_symmetric_mixture2_s2(
    x = x,
    mu = params$mu,
    kappa = params$kappa,
    nu = params$nu,
    log = TRUE
  ))
}

small_circle_symmetric_mixture2_start_thetas_s2 <- function(x,
                                                            weights = NULL,
                                                            control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  nu_eps <- as.numeric(control$small_circle_symmetric_mixture2_nu_eps %||% 1e-6)
  kappa_min <- as.numeric(control$small_circle_symmetric_mixture2_kappa_min %||% 1e-8)
  kappa_max <- as.numeric(control$small_circle_symmetric_mixture2_kappa_max %||% 1e6)

  warm_start <- control$small_circle_symmetric_mixture2_start_theta %||%
    control$theta_start %||%
    control$jp_mle_start_theta %||%
    NULL

  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- small_circle_symmetric_mixture2_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$small_circle_symmetric_mixture2_warm_start_only %||% FALSE)) {
      return(list(warm_start))
    }
  }

  resultant <- colSums(x * prob_weights)
  weighted_x <- sweep(x, 1L, sqrt(prob_weights), FUN = "*")
  second_moment <- crossprod(weighted_x)
  eig <- eigen(second_moment, symmetric = TRUE)
  mu_candidates <- rotational_unique_mu_candidates(list(
    resultant,
    -resultant,
    eig$vectors[, which.max(eig$values)],
    -eig$vectors[, which.max(eig$values)],
    eig$vectors[, which.min(eig$values)],
    -eig$vectors[, which.min(eig$values)],
    c(0, 0, 1),
    c(0, 0, -1)
  ))

  for (mu0 in mu_candidates) {
    folded_u <- abs(pmin(pmax(as.numeric(x %*% mu0), -1), 1))
    nu_candidates <- unique(pmin(pmax(c(
      rotational_weighted_mean(folded_u, prob_weights),
      rotational_weighted_quantile(folded_u, prob_weights, prob = 0.50),
      rotational_weighted_quantile(folded_u, prob_weights, prob = 0.70),
      rotational_weighted_quantile(folded_u, prob_weights, prob = 0.85)
    ), nu_eps), 1 - nu_eps))

    for (nu0 in nu_candidates) {
      var0 <- rotational_weighted_variance(folded_u, prob_weights, center = nu0)
      kappa0 <- min(max(1 / (2 * max(var0, 1 / (2 * kappa_max))), kappa_min), kappa_max)
      out[[length(out) + 1L]] <- small_circle_symmetric_mixture2_canonicalize_theta(list(
        mu = mu0,
        kappa = kappa0,
        nu = nu0
      ))
    }
  }

  if (length(out) == 0L) {
    out[[1L]] <- small_circle_symmetric_mixture2_canonicalize_theta(list(
      mu = c(0, 0, 1),
      kappa = 10,
      nu = 0.5
    ))
  }

  out
}

small_circle_symmetric_mixture2_pack_par <- function(theta,
                                                     control = list()) {
  theta <- small_circle_symmetric_mixture2_normalize_theta(theta, ambient_dim = 3L)
  nu_upper <- 1 - as.numeric(control$small_circle_symmetric_mixture2_nu_eps %||% 1e-6)
  list(
    par = c(
      log(pmax(expm1(theta$kappa), .Machine$double.eps)),
      small_circle_inverse_logistic_bounded(theta$nu, upper = nu_upper),
      theta$mu
    ),
    theta = theta
  )
}

small_circle_symmetric_mixture2_unpack_par <- function(par,
                                                       control = list()) {
  kappa_min <- as.numeric(control$small_circle_symmetric_mixture2_kappa_min %||% 1e-8)
  kappa_max <- as.numeric(control$small_circle_symmetric_mixture2_kappa_max %||% 1e6)
  nu_eps <- as.numeric(control$small_circle_symmetric_mixture2_nu_eps %||% 1e-6)

  mu_hat <- rotational_unit_vector_fallback(par[3:5])
  small_circle_symmetric_mixture2_canonicalize_theta(list(
    mu = mu_hat,
    kappa = min(max(log1p(exp(par[[1L]])), kappa_min), kappa_max),
    nu = small_circle_logistic_bounded(par[[2L]], upper = 1 - nu_eps)
  ))
}

small_circle_symmetric_mixture2_mle_s2_weighted <- function(x,
                                                            weights = NULL,
                                                            control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- small_circle_symmetric_mixture2_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(
    length(candidate_thetas),
    as.integer(control$small_circle_symmetric_mixture2_n_starts %||% 16L)
  ))]

  objective <- function(par) {
    theta <- small_circle_symmetric_mixture2_unpack_par(par, control = control)
    value <- -small_circle_symmetric_mixture2_weighted_loglik_s2(
      mu = theta$mu,
      kappa = theta$kappa,
      nu = theta$nu,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$small_circle_symmetric_mixture2_optim_method %||% "BFGS"
  optim_control <- control$small_circle_symmetric_mixture2_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- small_circle_symmetric_mixture2_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- small_circle_symmetric_mixture2_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$small_circle_symmetric_mixture2_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$small_circle_symmetric_mixture2_warm_start_only <- FALSE
    return(small_circle_symmetric_mixture2_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Symmetric small-circle-mixture weighted MLE failed for all starting values.")
  }

  c(
    best$theta,
    list(
      loglik = -best$opt$value,
      opt = best$opt,
      weighted_mle = TRUE,
      start_theta = best$start_theta
    )
  )
}

distance_profile_small_circle_symmetric_mixture2 <- function(omega,
                                                            t_values,
                                                            mu,
                                                            kappa,
                                                            nu,
                                                            distance_type = c("geodesic", "chordal"),
                                                            method = c("legendre", "integral"),
                                                            l_max = 150L,
                                                            quad_n = 1000L,
                                                            tol = 1e-10,
                                                            validate_against_integral = FALSE,
                                                            validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- small_circle_symmetric_mixture2_validate_parameters(
    mu = mu,
    kappa = kappa,
    nu = nu
  )

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }
    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_small_circle_symmetric_mixture2(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = theta$mu,
        kappa = theta$kappa,
        nu = theta$nu,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol,
        validate_against_integral = validate_against_integral,
        validation_tol = validation_tol
      )
    }, numeric(1)))
  }

  t_values <- as.numeric(t_values)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  if (theta$kappa <= 0) {
    out[active] <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_values[active])) / 2
    } else {
      (t_values[active]^2) / 4
    }
    return(out)
  }

  profile_plus <- distance_profile_small_circle(
    omega = omega,
    t_values = t_values,
    mu = theta$mu,
    kappa = theta$kappa,
    nu = theta$nu,
    distance_type = distance_type,
    method = method,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    validate_against_integral = FALSE,
    validation_tol = validation_tol
  )
  profile_minus <- distance_profile_small_circle(
    omega = omega,
    t_values = t_values,
    mu = -theta$mu,
    kappa = theta$kappa,
    nu = theta$nu,
    distance_type = distance_type,
    method = method,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    validate_against_integral = FALSE,
    validation_tol = validation_tol
  )
  out <- 0.5 * (profile_plus + profile_minus)
  out <- small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound)

  if (isTRUE(validate_against_integral) && identical(method, "legendre")) {
    out_integral <- distance_profile_small_circle_symmetric_mixture2(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      kappa = theta$kappa,
      nu = theta$nu,
      distance_type = distance_type,
      method = "integral",
      l_max = l_max,
      quad_n = quad_n,
      tol = tol,
      validate_against_integral = FALSE,
      validation_tol = validation_tol
    )
    discrepancy <- max(abs(out - out_integral))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Symmetric small-circle-mixture Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}

distance_profile_small_circle_symmetric_mixture2_grid <- function(omega_grid,
                                                                 mu,
                                                                 kappa,
                                                                 nu,
                                                                 t_grid,
                                                                 distance_type = c("geodesic", "chordal"),
                                                                 method = c("legendre", "integral"),
                                                                 l_max = 150L,
                                                                 quad_n = 1000L,
                                                                 tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- small_circle_symmetric_mixture2_validate_parameters(
    mu = mu,
    kappa = kappa,
    nu = nu
  )
  omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
  t_grid <- as.numeric(t_grid)

  if (theta$kappa <= 0) {
    base_profile <- if (identical(distance_type, "geodesic")) {
      (1 - cos(t_grid)) / 2
    } else {
      (t_grid^2) / 4
    }
    return(matrix(base_profile, nrow = nrow(omega_grid), ncol = length(t_grid), byrow = TRUE))
  }

  0.5 * (
    distance_profile_small_circle_grid(
      omega_grid = omega_grid,
      mu = theta$mu,
      kappa = theta$kappa,
      nu = theta$nu,
      t_grid = t_grid,
      distance_type = distance_type,
      method = method,
      l_max = l_max,
      quad_n = quad_n,
      tol = tol
    ) +
      distance_profile_small_circle_grid(
        omega_grid = omega_grid,
        mu = -theta$mu,
        kappa = theta$kappa,
        nu = theta$nu,
        t_grid = t_grid,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol
      )
  )
}

distance_profile_small_circle_symmetric_mixture2_cvm_grid <- function(X,
                                                                     mu,
                                                                     kappa,
                                                                     nu,
                                                                     distance_matrix = NULL,
                                                                     distance_type = c("geodesic", "chordal"),
                                                                     control = list(),
                                                                     method = c("legendre", "integral"),
                                                                     l_max = 150L,
                                                                     quad_n = 1000L,
                                                                     tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- small_circle_symmetric_mixture2_validate_parameters(
    mu = mu,
    kappa = kappa,
    nu = nu
  )
  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  if (!is.null(distance_matrix)) {
    distance_matrix <- as.matrix(distance_matrix)
    if (!identical(dim(distance_matrix), c(nrow(X), nrow(X)))) {
      stop("`distance_matrix` must be an n x n matrix matching `X`.")
    }
  }
  debug_memory_log(
    control = control,
    label = "enter distance_profile_small_circle_symmetric_mixture2_cvm_grid",
    objects = list(
      X = X,
      distance_matrix = distance_matrix
    )
  )

  if (theta$kappa <= 0) {
    dot_products <- if (is.null(distance_matrix)) {
      pmin(pmax(X %*% t(X), -1), 1)
    } else {
      if (identical(distance_type, "geodesic")) {
        pmin(pmax(cos(distance_matrix), -1), 1)
      } else {
        pmin(pmax(1 - (distance_matrix^2) / 2, -1), 1)
      }
    }
    return((1 - dot_products) / 2)
  }

  if (identical(method, "integral")) {
    debug_memory_log(control, "cvm_grid integral branch before first component")
    out <- distance_profile_small_circle_cvm_grid(
      X = X,
      mu = theta$mu,
      kappa = theta$kappa,
      nu = theta$nu,
      method = method,
      l_max = l_max,
      quad_n = quad_n,
      tol = tol
    )
    debug_memory_log(control, "cvm_grid integral branch after first component", list(out = out))
    out <- 0.5 * out
    out <- out + 0.5 * distance_profile_small_circle_cvm_grid(
      X = X,
      mu = -theta$mu,
      kappa = theta$kappa,
      nu = theta$nu,
      method = method,
      l_max = l_max,
      quad_n = quad_n,
      tol = tol
    )
    debug_memory_log(control, "cvm_grid integral branch after second component", list(out = out))
    return(out)
  }

  dot_products <- if (is.null(distance_matrix)) {
    pmin(pmax(X %*% t(X), -1), 1)
  } else {
    if (identical(distance_type, "geodesic")) {
      pmin(pmax(cos(distance_matrix), -1), 1)
    } else {
      pmin(pmax(1 - (distance_matrix^2) / 2, -1), 1)
    }
  }
  debug_memory_log(control, "cvm_grid legendre branch after dot_products", list(dot_products = dot_products))
  coeffs <- small_circle_legendre_coefficients(
    kappa = theta$kappa,
    nu = theta$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  r_values <- as.numeric(X %*% theta$mu)
  debug_memory_log(control, "cvm_grid legendre branch after coeffs/r_values", list(coeffs = coeffs, r_values = r_values))
  out <- 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )
  debug_memory_log(control, "cvm_grid legendre branch after first component", list(out = out))
  out <- out + 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = -r_values,
    coefficients = coeffs
  )
  debug_memory_log(control, "cvm_grid legendre branch after second component", list(out = out))
  out <- 0.5 * out
  out <- pmin(pmax(out, 0), 1)

  for (i in seq_len(nrow(X))) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = acos(dot_products[i, ]),
      values = out[i, ],
      upper_bound = pi
    )
  }

  out
}

small_circle_symmetric_mixture2_cvm_profile_block <- function(X_block,
                                                             dot_threshold_block,
                                                             mu,
                                                             kappa,
                                                             nu,
                                                             distance_type = c("geodesic", "chordal"),
                                                             method = c("legendre", "integral"),
                                                             l_max = 150L,
                                                             quad_n = 1000L,
                                                             tol = 1e-10,
                                                             control = list()) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  X_block <- jp_normalize_unit_matrix(X_block, arg_name = "`X_block`", min_ncol = 3L)
  dot_threshold_block <- pmin(pmax(as.matrix(dot_threshold_block), -1), 1)
  theta <- small_circle_symmetric_mixture2_validate_parameters(mu = mu, kappa = kappa, nu = nu)

  if (!identical(dim(dot_threshold_block), c(nrow(X_block), ncol(dot_threshold_block)))) {
    stop("`dot_threshold_block` must have nrow(`X_block`) rows.")
  }

  if (theta$kappa <= 0) {
    return((1 - dot_threshold_block) / 2)
  }

  if (identical(method, "integral")) {
    distance_block <- if (identical(distance_type, "geodesic")) {
      acos(dot_threshold_block)
    } else {
      sqrt(pmax(0, 2 * (1 - dot_threshold_block)))
    }
    out_plus <- t(vapply(seq_len(nrow(X_block)), function(i) {
      distance_profile_small_circle(
        omega = X_block[i, ],
        t_values = distance_block[i, ],
        mu = theta$mu,
        kappa = theta$kappa,
        nu = theta$nu,
        distance_type = distance_type,
        method = "integral",
        quad_n = quad_n,
        tol = tol
      )
    }, numeric(ncol(dot_threshold_block))))
    out_minus <- t(vapply(seq_len(nrow(X_block)), function(i) {
      distance_profile_small_circle(
        omega = X_block[i, ],
        t_values = distance_block[i, ],
        mu = -theta$mu,
        kappa = theta$kappa,
        nu = theta$nu,
        distance_type = distance_type,
        method = "integral",
        quad_n = quad_n,
        tol = tol
      )
    }, numeric(ncol(dot_threshold_block))))
    return(0.5 * (out_plus + out_minus))
  }

  coeffs <- small_circle_legendre_coefficients(
    kappa = theta$kappa,
    nu = theta$nu,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  r_values <- as.numeric(X_block %*% theta$mu)

  out_plus <- 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = dot_threshold_block,
    r = r_values,
    coefficients = coeffs
  )
  out_minus <- 1 - small_circle_projection_cdf_legendre_matrix(
    x_matrix = dot_threshold_block,
    r = -r_values,
    coefficients = coeffs
  )
  out <- 0.5 * (out_plus + out_minus)
  out <- pmin(pmax(out, 0), 1)

  for (i in seq_len(nrow(X_block))) {
    out[i, ] <- small_circle_monotone_clip_dot(dot_threshold_block[i, ], out[i, ])
  }

  debug_memory_log(
    control = control,
    label = "small_circle_symmetric_mixture2_cvm_profile_block",
    objects = list(
      X_block = X_block,
      dot_threshold_block = dot_threshold_block,
      out = out
    )
  )

  out
}

compute_weighted_sample_profile_block <- function(order_matrix_block,
                                                  rank_linear_index_block,
                                                  normalized_weights,
                                                  n_total) {
  order_matrix_block <- as.matrix(order_matrix_block)
  rank_linear_index_block <- as.matrix(rank_linear_index_block)
  block_rows <- nrow(order_matrix_block)
  n <- ncol(order_matrix_block)

  if (length(normalized_weights) != n || n_total != n) {
    stop("Incompatible dimensions in `compute_weighted_sample_profile_block()`.")
  }
  if (any(!is.finite(normalized_weights)) || any(normalized_weights < 0)) {
    stop("`normalized_weights` must be finite and nonnegative.")
  }

  total_weight <- sum(normalized_weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("`normalized_weights` must have strictly positive finite sum.")
  }

  ordered_weights_matrix <- matrix(
    normalized_weights[order_matrix_block],
    nrow = block_rows,
    ncol = n
  )
  cumulative_weights_matrix <- ordered_weights_matrix

  if (n >= 2L) {
    for (j in 2:n) {
      cumulative_weights_matrix[, j] <- cumulative_weights_matrix[, j] +
        cumulative_weights_matrix[, j - 1L]
    }
  }

  global_linear_index <- as.integer(rank_linear_index_block)
  global_row_index <- ((global_linear_index - 1L) %% n_total) + 1L
  global_rank_index <- ((global_linear_index - global_row_index) %/% n_total) + 1L
  local_row_index <- matrix(rep.int(seq_len(block_rows), n), nrow = block_rows, ncol = n)
  local_linear_index <- local_row_index + (matrix(global_rank_index, nrow = block_rows, ncol = n) - 1L) * block_rows

  out <- matrix(cumulative_weights_matrix[local_linear_index] / total_weight, nrow = block_rows, ncol = n)
  if (any(!is.finite(out)) || any(out < -1e-12) || any(out > 1 + 1e-12)) {
    stop("`compute_weighted_sample_profile_block()` produced values outside [0, 1].")
  }
  out
}

r_sph_small_circle_symmetric_mixture2 <- function(n,
                                                  mu,
                                                  kappa,
                                                  nu,
                                                  check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- small_circle_symmetric_mixture2_validate_parameters(
    mu = mu,
    kappa = kappa,
    nu = nu
  )
  component_plus <- stats::runif(n) <= 0.5
  x <- matrix(0, nrow = n, ncol = 3L)

  n_plus <- sum(component_plus)
  n_minus <- n - n_plus
  if (n_plus > 0L) {
    x[component_plus, ] <- r_sph_small_circle(
      n = n_plus,
      mu = theta$mu,
      kappa = theta$kappa,
      nu = theta$nu,
      check = FALSE
    )
  }
  if (n_minus > 0L) {
    x[!component_plus, ] <- r_sph_small_circle(
      n = n_minus,
      mu = -theta$mu,
      kappa = theta$kappa,
      nu = theta$nu,
      check = FALSE
    )
  }

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Symmetric small-circle-mixture sampler returned non-unit vectors.")
    }
  }

  x
}

beta_mixture2_validate_parameters <- function(mu,
                                                         weight1,
                                                         alpha1,
                                                         beta1,
                                                         alpha2,
                                                         beta2) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Rotational beta-mixture utilities currently support only S^2.")
  }

  weight1 <- as.numeric(weight1)
  alpha1 <- as.numeric(alpha1)
  beta1 <- as.numeric(beta1)
  alpha2 <- as.numeric(alpha2)
  beta2 <- as.numeric(beta2)

  if (length(weight1) != 1L || !is.finite(weight1) || weight1 <= 0 || weight1 >= 1) {
    stop("`weight1` must be a finite scalar in (0, 1).")
  }
  positive_shapes <- c(alpha1, beta1, alpha2, beta2)
  if (any(!is.finite(positive_shapes)) || any(positive_shapes <= 0)) {
    stop("Beta-mixture shape parameters must be finite and strictly positive.")
  }

  list(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2,
    ambient_dim = 3L
  )
}

beta_mixture2_canonicalize_theta <- function(theta) {
  params <- beta_mixture2_validate_parameters(
    mu = theta$mu,
    weight1 = theta$weight1,
    alpha1 = theta$alpha1,
    beta1 = theta$beta1,
    alpha2 = theta$alpha2,
    beta2 = theta$beta2
  )

  mean1 <- params$alpha1 / (params$alpha1 + params$beta1)
  mean2 <- params$alpha2 / (params$alpha2 + params$beta2)
  if (mean1 <= mean2) {
    return(params)
  }

  list(
    mu = params$mu,
    weight1 = 1 - params$weight1,
    alpha1 = params$alpha2,
    beta1 = params$beta2,
    alpha2 = params$alpha1,
    beta2 = params$beta1,
    ambient_dim = params$ambient_dim
  )
}

beta_mixture2_normalize_theta <- function(theta,
                                                     ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Beta-mixture theta must be a list.")
  }

  params <- beta_mixture2_canonicalize_theta(theta)
  if (params$ambient_dim != ambient_dim) {
    stop("Beta-mixture theta has incompatible ambient dimension.")
  }
  params
}

beta_mixture2_density_y <- function(y,
                                               weight1,
                                               alpha1,
                                               beta1,
                                               alpha2,
                                               beta2,
                                               log = FALSE,
                                               eps = 1e-12) {
  y <- rotational_clamp_unit_interval(as.numeric(y), eps = eps)
  density1 <- stats::dbeta(y, shape1 = alpha1, shape2 = beta1, log = TRUE)
  density2 <- stats::dbeta(y, shape1 = alpha2, shape2 = beta2, log = TRUE)
  log_density <- rotational_logsumexp2(
    log(weight1) + density1,
    log1p(-weight1) + density2
  )
  if (log) log_density else exp(log_density)
}

beta_mixture2_cdf_y <- function(y,
                                           weight1,
                                           alpha1,
                                           beta1,
                                           alpha2,
                                           beta2) {
  y <- as.numeric(y)
  out <- numeric(length(y))
  out[y <= 0] <- 0
  out[y >= 1] <- 1
  active <- which(y > 0 & y < 1)
  if (length(active) == 0L) {
    return(out)
  }

  out[active] <- weight1 * stats::pbeta(y[active], alpha1, beta1) +
    (1 - weight1) * stats::pbeta(y[active], alpha2, beta2)
  pmin(pmax(out, 0), 1)
}

beta_mixture2_density_h <- function(z,
                                               weight1,
                                               alpha1,
                                               beta1,
                                               alpha2,
                                               beta2,
                                               eps = 1e-12) {
  z <- as.numeric(z)
  y <- rotational_clamp_unit_interval((z + 1) / 2, eps = eps)
  out <- numeric(length(z))
  valid <- z >= -1 & z <= 1
  if (!any(valid)) {
    return(out)
  }

  out[valid] <- beta_mixture2_density_y(
    y = y[valid],
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2,
    log = FALSE,
    eps = eps
  )
  out
}

beta_mixture2_density_gz <- function(z,
                                                weight1,
                                                alpha1,
                                                beta1,
                                                alpha2,
                                                beta2) {
  0.5 * beta_mixture2_density_h(
    z = z,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )
}

d_sph_beta_mixture2_s2 <- function(x,
                                              mu,
                                              weight1,
                                              alpha1,
                                              beta1,
                                              alpha2,
                                              beta2,
                                              log = FALSE,
                                              eps = 1e-12) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )
  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  log_density <- -log(4 * pi) + beta_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    alpha1 = params$alpha1,
    beta1 = params$beta1,
    alpha2 = params$alpha2,
    beta2 = params$beta2,
    log = TRUE,
    eps = eps
  )
  if (log) log_density else exp(log_density)
}

beta_mixture2_weighted_loglik_s2 <- function(mu,
                                                        weight1,
                                                        alpha1,
                                                        beta1,
                                                        alpha2,
                                                        beta2,
                                                        x,
                                                        prob_weights = NULL,
                                                        eps = 1e-12) {
  params <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  sum(prob_weights * beta_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    alpha1 = params$alpha1,
    beta1 = params$beta1,
    alpha2 = params$alpha2,
    beta2 = params$beta2,
    log = TRUE,
    eps = eps
  ))
}

beta_mixture2_moment_match <- function(y,
                                                  weights,
                                                  shape_floor = 1e-3,
                                                  shape_ceiling = 1e4) {
  y <- rotational_clamp_unit_interval(y, eps = 1e-8)
  weights <- jp_normalize_probability_weights(weights, length(y))
  m <- min(max(rotational_weighted_mean(y, weights), 1e-4), 1 - 1e-4)
  v <- rotational_weighted_variance(y, weights, center = m)
  max_var <- m * (1 - m) * (1 - 1e-6)
  v <- min(max(v, 1e-6), max_var)
  precision <- min(max(m * (1 - m) / v - 1, 2 * shape_floor), shape_ceiling)

  list(
    alpha = min(max(m * precision, shape_floor), shape_ceiling),
    beta = min(max((1 - m) * precision, shape_floor), shape_ceiling)
  )
}

beta_mixture2_start_thetas_s2 <- function(x,
                                                      weights = NULL,
                                                      control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  warm_start <- control$beta_mixture2_start_theta %||%
    control$theta_start %||%
    control$jp_mle_start_theta %||%
    NULL
  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- beta_mixture2_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$beta_mixture2_warm_start_only %||% FALSE)) {
      return(list(warm_start))
    }
  }

  resultant <- colSums(x * prob_weights)
  second_moment <- rotational_weighted_covariance_matrix(x, prob_weights)
  eig <- eigen(second_moment, symmetric = TRUE)
  mu_candidates <- rotational_unique_mu_candidates(list(
    resultant,
    -resultant,
    eig$vectors[, which.max(eig$values)],
    -eig$vectors[, which.max(eig$values)],
    eig$vectors[, which.min(eig$values)],
    -eig$vectors[, which.min(eig$values)],
    c(0, 0, 1),
    c(0, 0, -1)
  ))

  split_probs <- c(0.35, 0.5, 0.65)
  for (mu0 in mu_candidates) {
    y <- rotational_clamp_unit_interval((pmin(pmax(as.numeric(x %*% mu0), -1), 1) + 1) / 2, eps = 1e-8)
    for (split_prob in split_probs) {
      threshold <- rotational_weighted_quantile(y, prob_weights, prob = split_prob)
      group1 <- y <= threshold
      if (all(group1) || !any(group1)) {
        next
      }

      w1 <- prob_weights[group1]
      w2 <- prob_weights[!group1]
      p1 <- sum(w1)
      if (p1 <= 1e-6 || p1 >= 1 - 1e-6) {
        next
      }

      comp1 <- beta_mixture2_moment_match(y[group1], w1 / sum(w1))
      comp2 <- beta_mixture2_moment_match(y[!group1], w2 / sum(w2))
      out[[length(out) + 1L]] <- beta_mixture2_canonicalize_theta(list(
        mu = mu0,
        weight1 = p1,
        alpha1 = comp1$alpha,
        beta1 = comp1$beta,
        alpha2 = comp2$alpha,
        beta2 = comp2$beta
      ))
    }
  }

  if (length(out) == 0L) {
    out[[1L]] <- beta_mixture2_canonicalize_theta(list(
      mu = c(0, 0, 1),
      weight1 = 0.5,
      alpha1 = 2,
      beta1 = 5,
      alpha2 = 5,
      beta2 = 2
    ))
  }

  out
}

beta_mixture2_pack_par <- function(theta,
                                              control = list()) {
  theta <- beta_mixture2_normalize_theta(theta, ambient_dim = 3L)
  list(
    par = c(
      stats::qlogis(theta$weight1),
      log(theta$alpha1),
      log(theta$beta1),
      log(theta$alpha2),
      log(theta$beta2),
      theta$mu
    ),
    theta = theta
  )
}

beta_mixture2_unpack_par <- function(par,
                                                control = list()) {
  shape_lower <- as.numeric(control$beta_mixture2_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$beta_mixture2_shape_upper %||% 1e3)
  weight_eps <- as.numeric(control$beta_mixture2_weight_eps %||% 0.01)

  mu_raw <- par[6:8]
  mu_hat <- rotational_unit_vector_fallback(mu_raw)
  beta_mixture2_canonicalize_theta(list(
    mu = mu_hat,
    weight1 = rotational_bounded_weight(par[[1L]], weight_eps = weight_eps),
    alpha1 = rotational_positive_parameter(par[[2L]], lower = shape_lower, upper = shape_upper),
    beta1 = rotational_positive_parameter(par[[3L]], lower = shape_lower, upper = shape_upper),
    alpha2 = rotational_positive_parameter(par[[4L]], lower = shape_lower, upper = shape_upper),
    beta2 = rotational_positive_parameter(par[[5L]], lower = shape_lower, upper = shape_upper)
  ))
}

beta_mixture2_mle_s2_weighted <- function(x,
                                                      weights = NULL,
                                                      control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- beta_mixture2_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(length(candidate_thetas), as.integer(control$beta_mixture2_n_starts %||% 12L)))]

  objective <- function(par) {
    theta <- beta_mixture2_unpack_par(par, control = control)
    value <- -beta_mixture2_weighted_loglik_s2(
      mu = theta$mu,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2,
      x = x,
      prob_weights = prob_weights
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$beta_mixture2_optim_method %||% "BFGS"
  optim_control <- control$beta_mixture2_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- beta_mixture2_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- beta_mixture2_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$beta_mixture2_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$beta_mixture2_warm_start_only <- FALSE
    return(beta_mixture2_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Beta-mixture weighted MLE failed for all starting values.")
  }

  c(
    best$theta,
    list(
      loglik = -best$opt$value,
      opt = best$opt,
      weighted_mle = TRUE,
      start_theta = best$start_theta
    )
  )
}

beta_mixture2_gauss_jacobi <- function(n, alpha, beta) {
  n <- as.integer(n)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= -1 ||
      length(beta) != 1L || !is.finite(beta) || beta <= -1) {
    stop("Gauss-Jacobi parameters must be finite and greater than -1.")
  }

  ab <- alpha + beta
  total_mass <- exp((ab + 1) * log(2) + lgamma(alpha + 1) +
    lgamma(beta + 1) - lgamma(ab + 2))

  diagonal <- numeric(n)
  diagonal[[1L]] <- (beta - alpha) / (ab + 2)
  if (n > 1L) {
    k_diag <- seq_len(n - 1L)
    diagonal[-1L] <- (beta^2 - alpha^2) /
      ((2 * k_diag + ab) * (2 * k_diag + ab + 2))

    k_off <- seq_len(n - 1L)
    offdiag <- sqrt(
      4 * k_off * (k_off + alpha) * (k_off + beta) * (k_off + ab) /
        ((2 * k_off + ab)^2 * ((2 * k_off + ab)^2 - 1))
    )
  } else {
    offdiag <- numeric(0L)
  }

  jacobi_matrix <- matrix(0, nrow = n, ncol = n)
  diag(jacobi_matrix) <- diagonal
  if (n > 1L) {
    jacobi_matrix[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- offdiag
    jacobi_matrix[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- offdiag
  }

  eig <- eigen(jacobi_matrix, symmetric = TRUE)
  order_idx <- order(eig$values)
  nodes <- eig$values[order_idx]
  weights <- total_mass * eig$vectors[1L, order_idx]^2

  list(nodes = nodes, weights = weights, total_mass = total_mass)
}

beta_mixture2_legendre_expectations_one_beta <- function(alpha, beta, l_max, quad_n) {
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 ||
      length(beta) != 1L || !is.finite(beta) || beta <= 0) {
    stop("Beta shape parameters must be finite and strictly positive.")
  }
  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  n_nodes <- max(quad_n, as.integer(ceiling((l_max + 1L) / 2)))
  jacobi <- beta_mixture2_gauss_jacobi(
    n = n_nodes,
    alpha = beta - 1,
    beta = alpha - 1
  )
  normalized_weights <- jacobi$weights / jacobi$total_mass
  legendre_matrix <- rotational_legendre_matrix(jacobi$nodes, l_max = l_max)
  expectations <- as.numeric(crossprod(legendre_matrix, normalized_weights))

  list(
    expectations = expectations,
    mass_error = abs(sum(normalized_weights) - 1),
    n_nodes = n_nodes
  )
}

beta_mixture2_legendre_coefficients <- function(theta,
                                                l_max = 150L,
                                                quad_n = 1000L,
                                                tol = 1e-10) {
  theta <- beta_mixture2_normalize_theta(theta, ambient_dim = 3L)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }
  if (length(quad_n) != 1L || !is.finite(quad_n) || quad_n < 1L) {
    stop("`quad_n` must be a strictly positive integer.")
  }

  comp1 <- beta_mixture2_legendre_expectations_one_beta(
    alpha = theta$alpha1,
    beta = theta$beta1,
    l_max = l_max,
    quad_n = quad_n
  )
  comp2 <- beta_mixture2_legendre_expectations_one_beta(
    alpha = theta$alpha2,
    beta = theta$beta2,
    l_max = l_max,
    quad_n = quad_n
  )

  expectations <- theta$weight1 * comp1$expectations +
    (1 - theta$weight1) * comp2$expectations
  ell <- 0:l_max
  coeffs <- (2 * ell + 1) * expectations
  coeffs[[1L]] <- 1

  if (any(!is.finite(coeffs))) {
    stop("Beta-mixture Gauss-Jacobi Legendre coefficient computation produced nonfinite coefficients.")
  }

  a0_error <- abs(expectations[[1L]] - 1)
  if (a0_error > tol) {
    stop(sprintf("Beta-mixture Gauss-Jacobi coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error,
    mass_error = max(comp1$mass_error, comp2$mass_error),
    method = "gauss_jacobi",
    quad_n = max(comp1$n_nodes, comp2$n_nodes)
  )
}

distance_profile_beta_mixture2 <- function(omega,
                                                      t_values,
                                                      mu,
                                                      weight1,
                                                      alpha1,
                                                      beta1,
                                                      alpha2,
                                                      beta2,
                                                      distance_type = c("geodesic", "chordal"),
                                                      method = c("legendre", "integral"),
                                                      l_max = 150L,
                                                      quad_n = 1000L,
                                                      tol = 1e-10,
                                                      validate_against_integral = FALSE,
                                                      validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }
    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_beta_mixture2(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol,
        validate_against_integral = validate_against_integral,
        validation_tol = validation_tol
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  t_values <- as.numeric(t_values)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  r_value <- sum(omega * theta$mu)
  if (abs(r_value - 1) <= 1e-12) {
    out[active] <- 1 - beta_mixture2_cdf_y(
      y = (thresholds + 1) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- beta_mixture2_cdf_y(
      y = (1 - thresholds) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }

  if (identical(method, "integral")) {
    return(rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        beta_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  coeffs <- beta_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  out[active] <- 1 - rotational_projection_cdf_legendre(
    x = thresholds,
    r = r_value,
    coefficients = coeffs
  )
  out <- small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound)

  if (isTRUE(validate_against_integral)) {
    out_integral <- rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        beta_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - out_integral))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Beta-mixture Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}

distance_profile_beta_mixture2_grid <- function(omega_grid,
                                                           mu,
                                                           weight1,
                                                           alpha1,
                                                           beta1,
                                                           alpha2,
                                                           beta2,
                                                           t_grid,
                                                           distance_type = c("geodesic", "chordal"),
                                                           method = c("legendre", "integral"),
                                                           l_max = 150L,
                                                           quad_n = 1000L,
                                                           tol = 1e-10) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )

  if (identical(method, "integral")) {
    omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * beta_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    out <- projection_profile_matrix_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density,
      distance_type = distance_type
    )
    return(out)
  }

  coeffs <- beta_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients

  out <- rotational_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = theta$mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )

  out
}

distance_profile_beta_mixture2_cvm_grid <- function(X,
                                                               mu,
                                                               weight1,
                                                               alpha1,
                                                               beta1,
                                                               alpha2,
                                                               beta2,
                                                               method = c("legendre", "integral"),
                                                               l_max = 150L,
                                                               quad_n = 1000L,
                                                               tol = 1e-10) {
  method <- match.arg(method)
  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )

  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)

  if (identical(method, "integral")) {
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * beta_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
    return(projection_sample_profile_matrix_integral(
      X = X,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    ))
  }

  coeffs <- beta_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  r_values <- as.numeric(X %*% theta$mu)
  out <- 1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )

  pos_idx <- which(abs(r_values - 1) <= 1e-12)
  if (length(pos_idx) > 0L) {
    out[pos_idx, ] <- 1 - beta_mixture2_cdf_y(
      y = (dot_products[pos_idx, , drop = FALSE] + 1) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
  }

  neg_idx <- which(abs(r_values + 1) <= 1e-12)
  if (length(neg_idx) > 0L) {
    out[neg_idx, ] <- beta_mixture2_cdf_y(
      y = (1 - dot_products[neg_idx, , drop = FALSE]) / 2,
      weight1 = theta$weight1,
      alpha1 = theta$alpha1,
      beta1 = theta$beta1,
      alpha2 = theta$alpha2,
      beta2 = theta$beta2
    )
  }

  out <- pmin(pmax(out, 0), 1)
  for (i in seq_len(nrow(X))) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = acos(dot_products[i, ]),
      values = out[i, ],
      upper_bound = pi
    )
  }
  out
}

r_sph_beta_mixture2 <- function(n,
                                           mu,
                                           weight1,
                                           alpha1,
                                           beta1,
                                           alpha2,
                                           beta2,
                                           check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- beta_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )

  component1 <- stats::runif(n) <= theta$weight1
  y <- numeric(n)
  n1 <- sum(component1)
  n2 <- n - n1
  if (n1 > 0L) {
    y[component1] <- stats::rbeta(n1, theta$alpha1, theta$beta1)
  }
  if (n2 > 0L) {
    y[!component1] <- stats::rbeta(n2, theta$alpha2, theta$beta2)
  }

  z <- pmin(pmax(2 * y - 1, -1), 1)
  phi <- stats::runif(n, min = 0, max = 2 * pi)
  basis <- jp_orthonormal_complement(theta$mu)
  tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
  radial <- sqrt(pmax(0, 1 - z^2))
  x <- tcrossprod(z, theta$mu) + sweep(tangent, 1L, radial, FUN = "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Beta-mixture sampler returned non-unit vectors.")
    }
  }
  x
}

uniform_beta_mixture_validate_parameters <- function(mu,
                                                     weight_uniform,
                                                     alpha,
                                                     beta) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Rotational uniform-beta-mixture utilities currently support only S^2.")
  }

  weight_uniform <- as.numeric(weight_uniform)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)

  if (length(weight_uniform) != 1L || !is.finite(weight_uniform) ||
      weight_uniform <= 0 || weight_uniform >= 1) {
    stop("`weight_uniform` must be a finite scalar in (0, 1).")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 ||
      length(beta) != 1L || !is.finite(beta) || beta <= 0) {
    stop("Uniform-beta-mixture shape parameters must be finite and strictly positive.")
  }

  list(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta,
    ambient_dim = 3L
  )
}

uniform_beta_mixture_normalize_theta <- function(theta,
                                                 ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Uniform-beta-mixture theta must be a list.")
  }

  params <- uniform_beta_mixture_validate_parameters(
    mu = theta$mu,
    weight_uniform = theta$weight_uniform,
    alpha = theta$alpha,
    beta = theta$beta
  )
  if (params$ambient_dim != ambient_dim) {
    stop("Uniform-beta-mixture theta has incompatible ambient dimension.")
  }
  params
}

uniform_beta_mixture_density_y <- function(y,
                                           weight_uniform,
                                           alpha,
                                           beta,
                                           log = FALSE,
                                           eps = 1e-12) {
  y <- rotational_clamp_unit_interval(as.numeric(y), eps = eps)
  log_density <- rotational_logsumexp2(
    log(weight_uniform),
    log1p(-weight_uniform) + stats::dbeta(y, shape1 = alpha, shape2 = beta, log = TRUE)
  )
  if (log) log_density else exp(log_density)
}

uniform_beta_mixture_cdf_y <- function(y,
                                       weight_uniform,
                                       alpha,
                                       beta) {
  y <- as.numeric(y)
  out <- numeric(length(y))
  out[y <= 0] <- 0
  out[y >= 1] <- 1
  active <- which(y > 0 & y < 1)
  if (length(active) == 0L) {
    return(out)
  }

  out[active] <- weight_uniform * y[active] +
    (1 - weight_uniform) * stats::pbeta(y[active], shape1 = alpha, shape2 = beta)
  pmin(pmax(out, 0), 1)
}

uniform_beta_mixture_density_h <- function(z,
                                           weight_uniform,
                                           alpha,
                                           beta,
                                           eps = 1e-12) {
  z <- as.numeric(z)
  y <- rotational_clamp_unit_interval((z + 1) / 2, eps = eps)
  out <- numeric(length(z))
  valid <- z >= -1 & z <= 1
  if (!any(valid)) {
    return(out)
  }

  out[valid] <- uniform_beta_mixture_density_y(
    y = y[valid],
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta,
    log = FALSE,
    eps = eps
  )
  out
}

uniform_beta_mixture_density_gz <- function(z,
                                            weight_uniform,
                                            alpha,
                                            beta,
                                            eps = 1e-12) {
  0.5 * uniform_beta_mixture_density_h(
    z = z,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta,
    eps = eps
  )
}

d_sph_uniform_beta_mixture_s2 <- function(x,
                                          mu,
                                          weight_uniform,
                                          alpha,
                                          beta,
                                          log = FALSE,
                                          eps = 1e-12) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )
  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  log_density <- -log(4 * pi) + uniform_beta_mixture_density_y(
    y = y,
    weight_uniform = params$weight_uniform,
    alpha = params$alpha,
    beta = params$beta,
    log = TRUE,
    eps = eps
  )
  if (log) log_density else exp(log_density)
}

uniform_beta_mixture_weighted_loglik_s2 <- function(mu,
                                                    weight_uniform,
                                                    alpha,
                                                    beta,
                                                    x,
                                                    prob_weights = NULL,
                                                    eps = 1e-12) {
  params <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  y <- rotational_clamp_unit_interval(
    (pmin(pmax(as.numeric(x %*% params$mu), -1), 1) + 1) / 2,
    eps = eps
  )
  sum(prob_weights * uniform_beta_mixture_density_y(
    y = y,
    weight_uniform = params$weight_uniform,
    alpha = params$alpha,
    beta = params$beta,
    log = TRUE,
    eps = eps
  ))
}

uniform_beta_mixture_start_thetas_s2 <- function(x,
                                                 weights = NULL,
                                                 control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  warm_start <- control$uniform_beta_mixture_start_theta %||%
    control$theta_start %||%
    control$jp_mle_start_theta %||%
    NULL
  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- uniform_beta_mixture_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$uniform_beta_mixture_warm_start_only %||% FALSE)) {
      return(list(warm_start))
    }
  }

  resultant <- colSums(x * prob_weights)
  second_moment <- rotational_weighted_covariance_matrix(x, prob_weights)
  eig <- eigen(second_moment, symmetric = TRUE)
  mu_candidates <- rotational_unique_mu_candidates(list(
    resultant,
    -resultant,
    eig$vectors[, which.max(eig$values)],
    -eig$vectors[, which.max(eig$values)],
    eig$vectors[, which.min(eig$values)],
    -eig$vectors[, which.min(eig$values)],
    c(0, 0, 1),
    c(0, 0, -1)
  ))

  split_probs <- c(0.15, 0.25, 0.35, 0.5)
  for (mu0 in mu_candidates) {
    y <- rotational_clamp_unit_interval((pmin(pmax(as.numeric(x %*% mu0), -1), 1) + 1) / 2, eps = 1e-8)
    for (split_prob in split_probs) {
      threshold <- rotational_weighted_quantile(y, prob_weights, prob = split_prob)
      beta_group <- y > threshold
      if (all(beta_group) || !any(beta_group)) {
        next
      }

      w_beta <- prob_weights[beta_group]
      weight_uniform0 <- sum(prob_weights[!beta_group])
      if (weight_uniform0 <= 1e-6 || weight_uniform0 >= 1 - 1e-6) {
        next
      }

      beta_start <- beta_mixture2_moment_match(y[beta_group], w_beta / sum(w_beta))
      out[[length(out) + 1L]] <- uniform_beta_mixture_normalize_theta(list(
        mu = mu0,
        weight_uniform = weight_uniform0,
        alpha = beta_start$alpha,
        beta = beta_start$beta
      ))
    }
  }

  if (length(out) == 0L) {
    out[[1L]] <- uniform_beta_mixture_normalize_theta(list(
      mu = c(0, 0, 1),
      weight_uniform = 0.2,
      alpha = 8,
      beta = 2
    ))
  }

  out
}

uniform_beta_mixture_pack_par <- function(theta,
                                          control = list()) {
  theta <- uniform_beta_mixture_normalize_theta(theta, ambient_dim = 3L)
  list(
    par = c(
      stats::qlogis(theta$weight_uniform),
      log(theta$alpha),
      log(theta$beta),
      theta$mu
    ),
    theta = theta
  )
}

uniform_beta_mixture_unpack_par <- function(par,
                                            control = list()) {
  shape_lower <- as.numeric(control$uniform_beta_mixture_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$uniform_beta_mixture_shape_upper %||% 1e3)
  weight_eps <- as.numeric(control$uniform_beta_mixture_weight_eps %||% 0.01)

  mu_raw <- par[4:6]
  mu_hat <- rotational_unit_vector_fallback(mu_raw)
  uniform_beta_mixture_normalize_theta(list(
    mu = mu_hat,
    weight_uniform = rotational_bounded_weight(par[[1L]], weight_eps = weight_eps),
    alpha = rotational_positive_parameter(par[[2L]], lower = shape_lower, upper = shape_upper),
    beta = rotational_positive_parameter(par[[3L]], lower = shape_lower, upper = shape_upper)
  ))
}

uniform_beta_mixture_mle_s2_weighted <- function(x,
                                                 weights = NULL,
                                                 control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- uniform_beta_mixture_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(length(candidate_thetas), as.integer(control$uniform_beta_mixture_n_starts %||% 12L)))]

  objective <- function(par) {
    theta <- uniform_beta_mixture_unpack_par(par, control = control)
    value <- -uniform_beta_mixture_weighted_loglik_s2(
      mu = theta$mu,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      x = x,
      prob_weights = prob_weights,
      eps = as.numeric(control$uniform_beta_mixture_eps %||% 1e-12)
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$uniform_beta_mixture_optim_method %||% "BFGS"
  optim_control <- control$uniform_beta_mixture_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- uniform_beta_mixture_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- uniform_beta_mixture_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$uniform_beta_mixture_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$uniform_beta_mixture_warm_start_only <- FALSE
    return(uniform_beta_mixture_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Uniform-beta-mixture weighted MLE failed for all starting values.")
  }

  c(
    best$theta,
    list(
      loglik = -best$opt$value,
      opt = best$opt,
      weighted_mle = TRUE,
      start_theta = best$start_theta
    )
  )
}

uniform_beta_mixture_legendre_coefficients <- function(theta,
                                                       l_max = 150L,
                                                       quad_n = 1000L,
                                                       tol = 1e-10) {
  theta <- uniform_beta_mixture_normalize_theta(theta, ambient_dim = 3L)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }
  if (length(quad_n) != 1L || !is.finite(quad_n) || quad_n < 1L) {
    stop("`quad_n` must be a strictly positive integer.")
  }

  beta_part <- beta_mixture2_legendre_expectations_one_beta(
    alpha = theta$alpha,
    beta = theta$beta,
    l_max = l_max,
    quad_n = quad_n
  )
  ell <- 0:l_max
  coeffs <- numeric(l_max + 1L)
  coeffs[[1L]] <- 1
  if (l_max >= 1L) {
    coeffs[-1L] <- (1 - theta$weight_uniform) *
      (2 * ell[-1L] + 1) * beta_part$expectations[-1L]
  }

  if (any(!is.finite(coeffs))) {
    stop("Uniform-beta-mixture Legendre coefficient computation produced nonfinite coefficients.")
  }

  a0_error <- abs(beta_part$expectations[[1L]] - 1)
  if (a0_error > tol) {
    stop(sprintf("Uniform-beta-mixture coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error,
    mass_error = beta_part$mass_error,
    method = "gauss_jacobi",
    quad_n = beta_part$n_nodes
  )
}

distance_profile_uniform_beta_mixture <- function(omega,
                                                  t_values,
                                                  mu,
                                                  weight_uniform,
                                                  alpha,
                                                  beta,
                                                  distance_type = c("geodesic", "chordal"),
                                                  method = c("legendre", "integral"),
                                                  l_max = 150L,
                                                  quad_n = 1000L,
                                                  tol = 1e-10,
                                                  eps = 1e-12,
                                                  validate_against_integral = FALSE,
                                                  validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }
    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_uniform_beta_mixture(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = theta$mu,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol,
        eps = eps,
        validate_against_integral = validate_against_integral,
        validation_tol = validation_tol
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  t_values <- as.numeric(t_values)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  r_value <- sum(omega * theta$mu)
  if (abs(r_value - 1) <= 1e-12) {
    out[active] <- 1 - uniform_beta_mixture_cdf_y(
      y = (thresholds + 1) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- uniform_beta_mixture_cdf_y(
      y = (1 - thresholds) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }

  if (identical(method, "integral")) {
    return(rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        uniform_beta_mixture_density_gz(
          z = z,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  coeffs <- uniform_beta_mixture_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  out[active] <- 1 - rotational_projection_cdf_legendre(
    x = thresholds,
    r = r_value,
    coefficients = coeffs
  )
  out <- small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound)

  if (isTRUE(validate_against_integral)) {
    out_integral <- rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        uniform_beta_mixture_density_gz(
          z = z,
          weight_uniform = theta$weight_uniform,
          alpha = theta$alpha,
          beta = theta$beta,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - out_integral))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Uniform-beta-mixture Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}

distance_profile_uniform_beta_mixture_grid <- function(omega_grid,
                                                       mu,
                                                       weight_uniform,
                                                       alpha,
                                                       beta,
                                                       t_grid,
                                                       distance_type = c("geodesic", "chordal"),
                                                       method = c("legendre", "integral"),
                                                       l_max = 150L,
                                                       quad_n = 1000L,
                                                       tol = 1e-10,
                                                       eps = 1e-12) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )

  if (identical(method, "integral")) {
    omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * uniform_beta_mixture_density_gz(
      z = quad$nodes,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      eps = eps
    )
    out <- projection_profile_matrix_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density,
      distance_type = distance_type
    )
    thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
    r_values <- as.numeric(omega_grid %*% theta$mu)
    pos_idx <- which(abs(r_values - 1) <= 1e-12)
    if (length(pos_idx) > 0L) {
      out[pos_idx, ] <- 1 - uniform_beta_mixture_cdf_y(
        y = (thresholds + 1) / 2,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- uniform_beta_mixture_cdf_y(
        y = (1 - thresholds) / 2,
        weight_uniform = theta$weight_uniform,
        alpha = theta$alpha,
        beta = theta$beta
      )
    }
    return(out)
  }

  coeffs <- uniform_beta_mixture_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  rotational_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = theta$mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )
}

distance_profile_uniform_beta_mixture_cvm_grid <- function(X,
                                                           mu,
                                                           weight_uniform,
                                                           alpha,
                                                           beta,
                                                           distance_matrix = NULL,
                                                           method = c("legendre", "integral"),
                                                           l_max = 150L,
                                                           quad_n = 1000L,
                                                           tol = 1e-10,
                                                           eps = 1e-12) {
  method <- match.arg(method)
  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )

  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  if (is.null(distance_matrix)) {
    dot_products <- pmin(pmax(X %*% t(X), -1), 1)
    distance_matrix <- acos(dot_products)
  } else {
    distance_matrix <- as.matrix(distance_matrix)
    if (!all(dim(distance_matrix) == c(nrow(X), nrow(X)))) {
      stop("`distance_matrix` must be an n x n matrix compatible with `X`.")
    }
    dot_products <- cos(pmin(pmax(distance_matrix, 0), pi))
  }

  if (identical(method, "integral")) {
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * uniform_beta_mixture_density_gz(
      z = quad$nodes,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta,
      eps = eps
    )
    return(projection_sample_profile_matrix_integral(
      X = X,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    ))
  }

  coeffs <- uniform_beta_mixture_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )$coefficients
  r_values <- as.numeric(X %*% theta$mu)
  out <- 1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )

  pos_idx <- which(abs(r_values - 1) <= 1e-12)
  if (length(pos_idx) > 0L) {
    out[pos_idx, ] <- 1 - uniform_beta_mixture_cdf_y(
      y = (dot_products[pos_idx, , drop = FALSE] + 1) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
  }

  neg_idx <- which(abs(r_values + 1) <= 1e-12)
  if (length(neg_idx) > 0L) {
    out[neg_idx, ] <- uniform_beta_mixture_cdf_y(
      y = (1 - dot_products[neg_idx, , drop = FALSE]) / 2,
      weight_uniform = theta$weight_uniform,
      alpha = theta$alpha,
      beta = theta$beta
    )
  }

  out <- pmin(pmax(out, 0), 1)
  for (i in seq_len(nrow(X))) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = distance_matrix[i, ],
      values = out[i, ],
      upper_bound = pi
    )
  }
  out
}

r_sph_uniform_beta_mixture <- function(n,
                                       mu,
                                       weight_uniform,
                                       alpha,
                                       beta,
                                       check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- uniform_beta_mixture_validate_parameters(
    mu = mu,
    weight_uniform = weight_uniform,
    alpha = alpha,
    beta = beta
  )

  x <- matrix(0, nrow = n, ncol = 3L)
  component_uniform <- stats::runif(n) <= theta$weight_uniform
  n_uniform <- sum(component_uniform)
  n_beta <- n - n_uniform

  if (n_uniform > 0L) {
    x[component_uniform, ] <- jp_uniform_sphere(n = n_uniform, ambient_dim = 3L)
  }
  if (n_beta > 0L) {
    y <- stats::rbeta(n_beta, shape1 = theta$alpha, shape2 = theta$beta)
    z <- pmin(pmax(2 * y - 1, -1), 1)
    phi <- stats::runif(n_beta, min = 0, max = 2 * pi)
    basis <- jp_orthonormal_complement(theta$mu)
    tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
    radial <- sqrt(pmax(0, 1 - z^2))
    x[!component_uniform, ] <- tcrossprod(z, theta$mu) + sweep(tangent, 1L, radial, FUN = "*")
  }

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Uniform-beta-mixture sampler returned non-unit vectors.")
    }
  }

  x
}

logitnormal_mixture2_validate_parameters <- function(mu,
                                                                weight1,
                                                                mean1,
                                                                sd1,
                                                                mean2,
                                                                sd2) {
  mu <- jp_normalize_unit_vector(mu, arg_name = "`mu`", min_length = 3L)
  if (length(mu) != 3L) {
    stop("Rotational logit-normal-mixture utilities currently support only S^2.")
  }

  weight1 <- as.numeric(weight1)
  mean1 <- as.numeric(mean1)
  mean2 <- as.numeric(mean2)
  sd1 <- as.numeric(sd1)
  sd2 <- as.numeric(sd2)

  if (length(weight1) != 1L || !is.finite(weight1) || weight1 <= 0 || weight1 >= 1) {
    stop("`weight1` must be a finite scalar in (0, 1).")
  }
  if (length(mean1) != 1L || !is.finite(mean1) || length(mean2) != 1L || !is.finite(mean2)) {
    stop("Logit-normal mixture means must be finite scalars.")
  }
  if (length(sd1) != 1L || !is.finite(sd1) || sd1 <= 0 || length(sd2) != 1L || !is.finite(sd2) || sd2 <= 0) {
    stop("Logit-normal mixture standard deviations must be strictly positive finite scalars.")
  }

  list(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2,
    ambient_dim = 3L
  )
}

logitnormal_mixture2_canonicalize_theta <- function(theta) {
  params <- logitnormal_mixture2_validate_parameters(
    mu = theta$mu,
    weight1 = theta$weight1,
    mean1 = theta$mean1,
    sd1 = theta$sd1,
    mean2 = theta$mean2,
    sd2 = theta$sd2
  )

  if (params$mean1 <= params$mean2) {
    return(params)
  }

  list(
    mu = params$mu,
    weight1 = 1 - params$weight1,
    mean1 = params$mean2,
    sd1 = params$sd2,
    mean2 = params$mean1,
    sd2 = params$sd1,
    ambient_dim = params$ambient_dim
  )
}

logitnormal_mixture2_normalize_theta <- function(theta,
                                                            ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Logit-normal-mixture theta must be a list.")
  }
  params <- logitnormal_mixture2_canonicalize_theta(theta)
  if (params$ambient_dim != ambient_dim) {
    stop("Logit-normal-mixture theta has incompatible ambient dimension.")
  }
  params
}

logitnormal_mixture2_density_y <- function(y,
                                                      weight1,
                                                      mean1,
                                                      sd1,
                                                      mean2,
                                                      sd2,
                                                      log = FALSE,
                                                      eps = 1e-12) {
  y <- as.numeric(y)
  out <- rep(if (log) -Inf else 0, length(y))
  active <- y > 0 & y < 1
  if (!any(active)) {
    return(out)
  }

  y_active <- rotational_clamp_unit_interval(y[active], eps = eps)
  x_active <- stats::qlogis(y_active)
  log_jacobian <- -log(y_active) - log1p(-y_active)
  log_density <- rotational_logsumexp2(
    log(weight1) + stats::dnorm(x_active, mean = mean1, sd = sd1, log = TRUE) + log_jacobian,
    log1p(-weight1) + stats::dnorm(x_active, mean = mean2, sd = sd2, log = TRUE) + log_jacobian
  )
  out[active] <- if (log) log_density else exp(log_density)
  out
}

logitnormal_mixture2_cdf_y <- function(y,
                                                  weight1,
                                                  mean1,
                                                  sd1,
                                                  mean2,
                                                  sd2,
                                                  eps = 1e-12) {
  y <- as.numeric(y)
  out <- numeric(length(y))
  out[y <= 0] <- 0
  out[y >= 1] <- 1
  active <- which(y > 0 & y < 1)
  if (length(active) == 0L) {
    return(out)
  }

  y_active <- rotational_clamp_unit_interval(y[active], eps = eps)
  x_active <- stats::qlogis(y_active)
  out[active] <- weight1 * stats::pnorm((x_active - mean1) / sd1) +
    (1 - weight1) * stats::pnorm((x_active - mean2) / sd2)
  pmin(pmax(out, 0), 1)
}

logitnormal_mixture2_density_h <- function(z,
                                                      weight1,
                                                      mean1,
                                                      sd1,
                                                      mean2,
                                                      sd2,
                                                      eps = 1e-12) {
  z <- as.numeric(z)
  y <- (z + 1) / 2
  out <- numeric(length(z))
  valid <- z > -1 & z < 1
  if (!any(valid)) {
    return(out)
  }

  out[valid] <- logitnormal_mixture2_density_y(
    y = y[valid],
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2,
    log = FALSE,
    eps = eps
  )
  out
}

logitnormal_mixture2_density_gz <- function(z,
                                                       weight1,
                                                       mean1,
                                                       sd1,
                                                       mean2,
                                                       sd2,
                                                       eps = 1e-12) {
  0.5 * logitnormal_mixture2_density_h(
    z = z,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2,
    eps = eps
  )
}

d_sph_logitnormal_mixture2_s2 <- function(x,
                                                     mu,
                                                     weight1,
                                                     mean1,
                                                     sd1,
                                                     mean2,
                                                     sd2,
                                                     log = FALSE,
                                                     eps = 1e-12) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  params <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )
  z <- pmin(pmax(as.numeric(x %*% params$mu), -1 + eps, 1 - eps), 1 - eps)
  y <- (z + 1) / 2
  log_density <- -log(4 * pi) + logitnormal_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    mean1 = params$mean1,
    sd1 = params$sd1,
    mean2 = params$mean2,
    sd2 = params$sd2,
    log = TRUE,
    eps = eps
  )
  if (log) log_density else exp(log_density)
}

logitnormal_mixture2_weighted_loglik_s2 <- function(mu,
                                                               weight1,
                                                               mean1,
                                                               sd1,
                                                               mean2,
                                                               sd2,
                                                               x,
                                                               prob_weights = NULL,
                                                               eps = 1e-12) {
  params <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(prob_weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(prob_weights, nrow(x))
  }

  z <- pmin(pmax(as.numeric(x %*% params$mu), -1 + eps, 1 - eps), 1 - eps)
  y <- (z + 1) / 2
  sum(prob_weights * logitnormal_mixture2_density_y(
    y = y,
    weight1 = params$weight1,
    mean1 = params$mean1,
    sd1 = params$sd1,
    mean2 = params$mean2,
    sd2 = params$sd2,
    log = TRUE,
    eps = eps
  ))
}

logitnormal_mixture2_start_thetas_s2 <- function(x,
                                                            weights = NULL,
                                                            control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  warm_start <- control$logitnormal_mixture2_start_theta %||%
    control$theta_start %||%
    control$jp_mle_start_theta %||%
    NULL
  out <- list()
  if (!is.null(warm_start)) {
    warm_start <- logitnormal_mixture2_normalize_theta(warm_start, ambient_dim = 3L)
    out[[length(out) + 1L]] <- warm_start
    if (isTRUE(control$logitnormal_mixture2_warm_start_only %||% FALSE)) {
      return(list(warm_start))
    }
  }

  resultant <- colSums(x * prob_weights)
  second_moment <- rotational_weighted_covariance_matrix(x, prob_weights)
  eig <- eigen(second_moment, symmetric = TRUE)
  mu_candidates <- rotational_unique_mu_candidates(list(
    resultant,
    -resultant,
    eig$vectors[, which.max(eig$values)],
    -eig$vectors[, which.max(eig$values)],
    eig$vectors[, which.min(eig$values)],
    -eig$vectors[, which.min(eig$values)],
    c(0, 0, 1),
    c(0, 0, -1)
  ))

  split_probs <- c(0.35, 0.5, 0.65)
  for (mu0 in mu_candidates) {
    y <- rotational_clamp_unit_interval((pmin(pmax(as.numeric(x %*% mu0), -1 + 1e-10), 1 - 1e-10) + 1) / 2, eps = 1e-8)
    x_logit <- stats::qlogis(y)
    for (split_prob in split_probs) {
      threshold <- rotational_weighted_quantile(x_logit, prob_weights, prob = split_prob)
      group1 <- x_logit <= threshold
      if (all(group1) || !any(group1)) {
        next
      }

      w1 <- prob_weights[group1]
      w2 <- prob_weights[!group1]
      p1 <- sum(w1)
      if (p1 <= 1e-6 || p1 >= 1 - 1e-6) {
        next
      }

      mean1 <- rotational_weighted_mean(x_logit[group1], w1 / sum(w1))
      mean2 <- rotational_weighted_mean(x_logit[!group1], w2 / sum(w2))
      sd1 <- sqrt(max(rotational_weighted_variance(x_logit[group1], w1 / sum(w1), center = mean1), 0.05^2))
      sd2 <- sqrt(max(rotational_weighted_variance(x_logit[!group1], w2 / sum(w2), center = mean2), 0.05^2))
      out[[length(out) + 1L]] <- logitnormal_mixture2_canonicalize_theta(list(
        mu = mu0,
        weight1 = p1,
        mean1 = mean1,
        sd1 = sd1,
        mean2 = mean2,
        sd2 = sd2
      ))
    }
  }

  if (length(out) == 0L) {
    out[[1L]] <- logitnormal_mixture2_canonicalize_theta(list(
      mu = c(0, 0, 1),
      weight1 = 0.5,
      mean1 = -1,
      sd1 = 0.7,
      mean2 = 1,
      sd2 = 0.7
    ))
  }

  out
}

logitnormal_mixture2_pack_par <- function(theta,
                                                     control = list()) {
  theta <- logitnormal_mixture2_normalize_theta(theta, ambient_dim = 3L)
  list(
    par = c(
      stats::qlogis(theta$weight1),
      theta$mean1,
      log(theta$sd1),
      theta$mean2,
      log(theta$sd2),
      theta$mu
    ),
    theta = theta
  )
}

logitnormal_mixture2_unpack_par <- function(par,
                                                       control = list()) {
  clip_means <- isTRUE(control$logitnormal_mixture2_clip_means %||% FALSE)
  mean_lower <- as.numeric(control$logitnormal_mixture2_mean_lower %||% -8)
  mean_upper <- as.numeric(control$logitnormal_mixture2_mean_upper %||% 8)
  sd_lower <- as.numeric(control$logitnormal_mixture2_sd_lower %||% 0.05)
  sd_upper <- as.numeric(control$logitnormal_mixture2_sd_upper %||% 5)
  weight_eps <- as.numeric(control$logitnormal_mixture2_weight_eps %||% 0.01)

  mu_raw <- par[6:8]
  mu_hat <- rotational_unit_vector_fallback(mu_raw)
  mean1_hat <- par[[2L]]
  mean2_hat <- par[[4L]]
  if (clip_means) {
    mean1_hat <- pmin(pmax(mean1_hat, mean_lower), mean_upper)
    mean2_hat <- pmin(pmax(mean2_hat, mean_lower), mean_upper)
  }
  logitnormal_mixture2_canonicalize_theta(list(
    mu = mu_hat,
    weight1 = rotational_bounded_weight(par[[1L]], weight_eps = weight_eps),
    mean1 = mean1_hat,
    sd1 = rotational_positive_parameter(par[[3L]], lower = sd_lower, upper = sd_upper),
    mean2 = mean2_hat,
    sd2 = rotational_positive_parameter(par[[5L]], lower = sd_lower, upper = sd_upper)
  ))
}

logitnormal_mixture2_mle_s2_weighted <- function(x,
                                                            weights = NULL,
                                                            control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  prob_weights <- if (is.null(weights)) {
    rep(1 / nrow(x), nrow(x))
  } else {
    jp_normalize_probability_weights(weights, nrow(x))
  }

  candidate_thetas <- logitnormal_mixture2_start_thetas_s2(
    x = x,
    weights = prob_weights,
    control = control
  )
  candidate_thetas <- candidate_thetas[seq_len(min(length(candidate_thetas), as.integer(control$logitnormal_mixture2_n_starts %||% 12L)))]

  objective <- function(par) {
    theta <- logitnormal_mixture2_unpack_par(par, control = control)
    value <- -logitnormal_mixture2_weighted_loglik_s2(
      mu = theta$mu,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      x = x,
      prob_weights = prob_weights,
      eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12)
    )
    if (!is.finite(value)) {
      .Machine$double.xmax / 100
    } else {
      value
    }
  }

  optim_method <- control$logitnormal_mixture2_optim_method %||% "BFGS"
  optim_control <- control$logitnormal_mixture2_optim_control %||% list(maxit = 400L, reltol = 1e-9)

  best <- NULL
  for (theta0 in candidate_thetas) {
    par0 <- logitnormal_mixture2_pack_par(theta0, control = control)$par
    opt <- try(stats::optim(
      par = par0,
      fn = objective,
      method = optim_method,
      control = optim_control
    ), silent = TRUE)
    if (inherits(opt, "try-error")) {
      next
    }

    theta_hat <- logitnormal_mixture2_unpack_par(opt$par, control = control)
    if (is.null(best) || opt$value < best$opt$value) {
      best <- list(theta = theta_hat, opt = opt, start_theta = theta0)
    }
  }

  if ((is.null(best) || isTRUE(best$opt$convergence != 0L)) &&
      isTRUE(control$logitnormal_mixture2_warm_start_only %||% FALSE)) {
    fallback_control <- control
    fallback_control$logitnormal_mixture2_warm_start_only <- FALSE
    return(logitnormal_mixture2_mle_s2_weighted(
      x = x,
      weights = prob_weights,
      control = fallback_control
    ))
  }

  if (is.null(best)) {
    stop("Logit-normal-mixture weighted MLE failed for all starting values.")
  }

  c(
    best$theta,
    list(
      loglik = -best$opt$value,
      opt = best$opt,
      weighted_mle = TRUE,
      start_theta = best$start_theta
    )
  )
}

logitnormal_mixture2_legendre_coefficients <- function(theta,
                                                                  l_max = 150L,
                                                                  quad_n = 1000L,
                                                                  tol = 1e-10,
                                                                  eps = 1e-12) {
  theta <- logitnormal_mixture2_normalize_theta(theta, ambient_dim = 3L)
  l_max <- as.integer(l_max)
  quad_n <- as.integer(quad_n)
  tol <- as.numeric(tol)

  if (length(l_max) != 1L || !is.finite(l_max) || l_max < 0L) {
    stop("`l_max` must be a nonnegative integer.")
  }

  quad <- rotational_gauss_hermite(quad_n)
  x1 <- theta$mean1 + theta$sd1 * sqrt(2) * quad$nodes
  x2 <- theta$mean2 + theta$sd2 * sqrt(2) * quad$nodes
  z1 <- tanh(x1 / 2)
  z2 <- tanh(x2 / 2)

  legendre_1 <- rotational_legendre_matrix(z1, l_max = l_max)
  legendre_2 <- rotational_legendre_matrix(z2, l_max = l_max)
  expectation_1 <- as.numeric(crossprod(legendre_1, quad$weights / sqrt(pi)))
  expectation_2 <- as.numeric(crossprod(legendre_2, quad$weights / sqrt(pi)))
  expectations <- theta$weight1 * expectation_1 + (1 - theta$weight1) * expectation_2

  ell <- 0:l_max
  coeffs <- (2 * ell + 1) * expectations
  coeffs[[1L]] <- 1
  a0_error <- abs(expectations[[1L]] - 1)
  if (a0_error > tol) {
    stop(sprintf("Rotational logit-normal Legendre coefficient check failed: |a0 - 1| = %.3e.", a0_error))
  }

  list(
    coefficients = coeffs,
    a0_error = a0_error
  )
}

distance_profile_logitnormal_mixture2 <- function(omega,
                                                             t_values,
                                                             mu,
                                                             weight1,
                                                             mean1,
                                                             sd1,
                                                             mean2,
                                                             sd2,
                                                             distance_type = c("geodesic", "chordal"),
                                                             method = c("legendre", "integral"),
                                                             l_max = 150L,
                                                             quad_n = 1000L,
                                                             tol = 1e-10,
                                                             eps = 1e-12,
                                                             validate_against_integral = FALSE,
                                                             validation_tol = 5e-6) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  if (is.matrix(omega)) {
    omega <- jp_normalize_unit_matrix(omega, arg_name = "`omega`", min_ncol = 3L)
    if (length(t_values) == 1L) {
      t_values <- rep(t_values, nrow(omega))
    }
    if (length(t_values) != nrow(omega)) {
      stop("When `omega` is a matrix, `t_values` must have length 1 or nrow(omega).")
    }
    return(vapply(seq_len(nrow(omega)), function(i) {
      distance_profile_logitnormal_mixture2(
        omega = omega[i, ],
        t_values = t_values[i],
        mu = theta$mu,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        distance_type = distance_type,
        method = method,
        l_max = l_max,
        quad_n = quad_n,
        tol = tol,
        eps = eps,
        validate_against_integral = validate_against_integral,
        validation_tol = validation_tol
      )
    }, numeric(1)))
  }

  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  t_values <- as.numeric(t_values)
  upper_bound <- if (identical(distance_type, "geodesic")) pi else 2
  out <- numeric(length(t_values))
  out[t_values <= 0] <- 0
  out[t_values >= upper_bound] <- 1
  active <- which(is.finite(t_values) & t_values > 0 & t_values < upper_bound)
  if (length(active) == 0L) {
    return(out)
  }

  thresholds <- sphere_distance_to_dot_threshold(t_values[active], distance_type = distance_type)
  r_value <- sum(omega * theta$mu)
  if (abs(r_value - 1) <= 1e-12) {
    out[active] <- 1 - logitnormal_mixture2_cdf_y(
      y = (thresholds + 1) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }
  if (abs(r_value + 1) <= 1e-12) {
    out[active] <- logitnormal_mixture2_cdf_y(
      y = (1 - thresholds) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    return(small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound))
  }

  if (identical(method, "integral")) {
    return(rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        logitnormal_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    ))
  }

  coeffs <- logitnormal_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    eps = eps
  )$coefficients
  out[active] <- 1 - rotational_projection_cdf_legendre(
    x = thresholds,
    r = r_value,
    coefficients = coeffs
  )
  out <- small_circle_monotone_clip(t_values = t_values, values = out, upper_bound = upper_bound)

  if (isTRUE(validate_against_integral)) {
    out_integral <- rotational_distance_profile_integral(
      omega = omega,
      t_values = t_values,
      mu = theta$mu,
      density_gz = function(z) {
        logitnormal_mixture2_density_gz(
          z = z,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          eps = eps
        )
      },
      distance_type = distance_type,
      quad_n = quad_n
    )
    discrepancy <- max(abs(out - out_integral))
    if (discrepancy > validation_tol) {
      stop(sprintf(
        "Logit-normal-mixture Legendre profile validation failed: max discrepancy %.3e exceeds %.3e.",
        discrepancy,
        validation_tol
      ))
    }
  }

  out
}

distance_profile_logitnormal_mixture2_grid <- function(omega_grid,
                                                                  mu,
                                                                  weight1,
                                                                  mean1,
                                                                  sd1,
                                                                  mean2,
                                                                  sd2,
                                                                  t_grid,
                                                                  distance_type = c("geodesic", "chordal"),
                                                                  method = c("legendre", "integral"),
                                                                  l_max = 150L,
                                                                  quad_n = 1000L,
                                                                  tol = 1e-10,
                                                                  eps = 1e-12) {
  distance_type <- match.arg(distance_type)
  method <- match.arg(method)
  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  if (identical(method, "integral")) {
    omega_grid <- jp_normalize_unit_matrix(omega_grid, arg_name = "`omega_grid`", min_ncol = 3L)
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * logitnormal_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    out <- projection_profile_matrix_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density,
      distance_type = distance_type
    )
    thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
    r_values <- as.numeric(omega_grid %*% theta$mu)
    pos_idx <- which(abs(r_values - 1) <= 1e-12)
    if (length(pos_idx) > 0L) {
      out[pos_idx, ] <- 1 - logitnormal_mixture2_cdf_y(
        y = (thresholds + 1) / 2,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        eps = eps
      )
    }
    neg_idx <- which(abs(r_values + 1) <= 1e-12)
    if (length(neg_idx) > 0L) {
      out[neg_idx, ] <- logitnormal_mixture2_cdf_y(
        y = (1 - thresholds) / 2,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        eps = eps
      )
    }
    return(out)
  }

  coeffs <- logitnormal_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    eps = eps
  )$coefficients
  rotational_profile_matrix_legendre(
    t_grid = t_grid,
    omega_grid = omega_grid,
    mu = theta$mu,
    coeffs = coeffs,
    Lmax = l_max,
    distance_type = distance_type
  )
}

distance_profile_logitnormal_mixture2_cvm_grid <- function(X,
                                                                      mu,
                                                                      weight1,
                                                                      mean1,
                                                                      sd1,
                                                                      mean2,
                                                                      sd2,
                                                                      method = c("legendre", "integral"),
                                                                      l_max = 150L,
                                                                      quad_n = 1000L,
                                                                      tol = 1e-10,
                                                                      eps = 1e-12) {
  method <- match.arg(method)
  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  X <- jp_normalize_unit_matrix(X, arg_name = "`X`", min_ncol = 3L)
  dot_products <- pmin(pmax(X %*% t(X), -1), 1)

  if (identical(method, "integral")) {
    quad <- rotational_gauss_legendre(as.integer(quad_n))
    weighted_density <- quad$weights * logitnormal_mixture2_density_gz(
      z = quad$nodes,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
    return(projection_sample_profile_matrix_integral(
      X = X,
      mu = theta$mu,
      z_nodes = quad$nodes,
      weighted_density = weighted_density
    ))
  }

  coeffs <- logitnormal_mixture2_legendre_coefficients(
    theta = theta,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol,
    eps = eps
  )$coefficients
  r_values <- as.numeric(X %*% theta$mu)
  out <- 1 - rotational_projection_cdf_legendre_matrix(
    x_matrix = dot_products,
    r = r_values,
    coefficients = coeffs
  )

  pos_idx <- which(abs(r_values - 1) <= 1e-12)
  if (length(pos_idx) > 0L) {
    out[pos_idx, ] <- 1 - logitnormal_mixture2_cdf_y(
      y = (dot_products[pos_idx, , drop = FALSE] + 1) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
  }

  neg_idx <- which(abs(r_values + 1) <= 1e-12)
  if (length(neg_idx) > 0L) {
    out[neg_idx, ] <- logitnormal_mixture2_cdf_y(
      y = (1 - dot_products[neg_idx, , drop = FALSE]) / 2,
      weight1 = theta$weight1,
      mean1 = theta$mean1,
      sd1 = theta$sd1,
      mean2 = theta$mean2,
      sd2 = theta$sd2,
      eps = eps
    )
  }

  out <- pmin(pmax(out, 0), 1)
  for (i in seq_len(nrow(X))) {
    out[i, ] <- small_circle_monotone_clip(
      t_values = acos(dot_products[i, ]),
      values = out[i, ],
      upper_bound = pi
    )
  }
  out
}

r_sph_logitnormal_mixture2 <- function(n,
                                                  mu,
                                                  weight1,
                                                  mean1,
                                                  sd1,
                                                  mean2,
                                                  sd2,
                                                  check = TRUE) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a strictly positive integer.")
  }

  theta <- logitnormal_mixture2_validate_parameters(
    mu = mu,
    weight1 = weight1,
    mean1 = mean1,
    sd1 = sd1,
    mean2 = mean2,
    sd2 = sd2
  )

  component1 <- stats::runif(n) <= theta$weight1
  t_values <- numeric(n)
  n1 <- sum(component1)
  n2 <- n - n1
  if (n1 > 0L) {
    t_values[component1] <- stats::rnorm(n1, mean = theta$mean1, sd = theta$sd1)
  }
  if (n2 > 0L) {
    t_values[!component1] <- stats::rnorm(n2, mean = theta$mean2, sd = theta$sd2)
  }

  y <- stats::plogis(t_values)
  z <- pmin(pmax(2 * y - 1, -1), 1)
  phi <- stats::runif(n, min = 0, max = 2 * pi)
  basis <- jp_orthonormal_complement(theta$mu)
  tangent <- tcrossprod(cos(phi), basis[, 1L]) + tcrossprod(sin(phi), basis[, 2L])
  radial <- sqrt(pmax(0, 1 - z^2))
  x <- tcrossprod(z, theta$mu) + sweep(tangent, 1L, radial, FUN = "*")

  if (isTRUE(check)) {
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || max(abs(norms - 1)) > 1e-8) {
      stop("Logit-normal-mixture sampler returned non-unit vectors.")
    }
  }
  x
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
