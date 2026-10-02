# Mathematical and numerical helpers: normal.

theoretical_distance_profile_normal <- function(omega, mu, sigma, t_values) {
  n <- max(length(omega), length(t_values))
  omega <- rep(omega, length.out = n)
  t_values <- rep(t_values, length.out = n)
  result <- numeric(n)
  valid_t <- t_values > 0
  if (any(valid_t)) {
    if (exists("distance_profile_backend_current", mode = "function") &&
        identical(distance_profile_backend_current(), "cpp") &&
        length(mu) == 1L && length(sigma) == 1L) {
      return(distance_profile_cpp_call(
        "cpp_dp_normal_profile",
        as.numeric(omega),
        as.numeric(t_values),
        as.numeric(mu),
        as.numeric(sigma)
      ))
    }
    upper_bound <- (omega[valid_t] + t_values[valid_t] - mu) / sigma
    lower_bound <- (omega[valid_t] - t_values[valid_t] - mu) / sigma
    result[valid_t] <- pnorm(upper_bound) - pnorm(lower_bound)
  }
  return(result)
}
dotF_mu_normal <- function(omega, mu, sigma, t) {
  -(1/sigma) * dnorm((omega - mu + t)/sigma) +
    (1/sigma) * dnorm((omega - mu - t)/sigma)
}
dotF_sigma_normal <- function(omega, mu, sigma, t) {
  -(omega - mu + t)/sigma^2 * dnorm((omega - mu + t)/sigma) +
    (omega - mu - t)/sigma^2 * dnorm((omega - mu - t)/sigma)
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
