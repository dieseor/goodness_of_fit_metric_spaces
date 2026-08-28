# Restricted Logistic-Gaussian AR(1) model used by the t-logistic scenario:
#
#   ilr(X) ~ N_d(mu, R_d(rho)),
#   [R_d(rho)]_{jk} = rho^{|j-k|},
#   mu in R^d, rho in (-1, 1).
#
# Both mu and rho are estimated.  This adapter is deliberately separate
# from the unrestricted Logistic-Gaussian and sigma^2 I_d adapters.

if (!exists("make_logistic_gaussian_spec", mode = "function")) {
  lg_ar1_candidates <- c(
    file.path("bootstrap", "logistic_gaussian_model_spec.R"),
    "logistic_gaussian_model_spec.R",
    file.path("..", "bootstrap", "logistic_gaussian_model_spec.R"),
    file.path("..", "..", "bootstrap", "logistic_gaussian_model_spec.R")
  )
  lg_ar1_path <- lg_ar1_candidates[file.exists(lg_ar1_candidates)][1L]
  if (is.na(lg_ar1_path)) {
    stop("Could not locate bootstrap/logistic_gaussian_model_spec.R.")
  }
  source(lg_ar1_path)
}

logistic_gaussian_ar1_covariance <- function(ilr_dim, rho) {
  d <- as.integer(ilr_dim)
  rho <- as.numeric(rho)

  if (length(d) != 1L || !is.finite(d) || d < 1L) {
    stop("`ilr_dim` must be a strictly positive integer.")
  }
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("AR(1) parameter `rho` must be a finite scalar in (-1, 1).")
  }

  idx <- seq_len(d)
  rho ^ abs(outer(idx, idx, "-"))
}

logistic_gaussian_ar1_covariance_derivative <- function(ilr_dim, rho) {
  d <- as.integer(ilr_dim)
  rho <- as.numeric(rho)

  if (length(d) != 1L || !is.finite(d) || d < 1L) {
    stop("`ilr_dim` must be a strictly positive integer.")
  }
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("AR(1) parameter `rho` must be a finite scalar in (-1, 1).")
  }

  idx <- seq_len(d)
  lag <- abs(outer(idx, idx, "-"))
  out <- matrix(0, nrow = d, ncol = d)
  nz <- lag > 0
  out[nz] <- lag[nz] * rho ^ (lag[nz] - 1)
  out
}

normalize_logistic_gaussian_ar1_theta <- function(theta,
                                                   ambient_dim = NULL,
                                                   control = list()) {
  if (!is.list(theta)) {
    stop("Logistic-Gaussian AR(1) theta must contain `mu_ilr` and `rho`.")
  }

  mu_ilr <- as.numeric(theta$mu_ilr %||% theta$mu)
  rho <- suppressWarnings(as.numeric(theta$rho))

  if (!length(mu_ilr) || any(!is.finite(mu_ilr))) {
    stop("Logistic-Gaussian AR(1) theta requires a finite `mu_ilr` vector.")
  }
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("Logistic-Gaussian AR(1) `rho` must lie in (-1, 1).")
  }

  ambient_dim <- as.integer(
    ambient_dim %||% theta$ambient_dim %||% (length(mu_ilr) + 1L)
  )
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) ||
      ambient_dim < 2L || length(mu_ilr) != ambient_dim - 1L) {
    stop("Logistic-Gaussian AR(1) theta has incompatible dimension.")
  }

  d <- ambient_dim - 1L
  Sigma_ilr <- logistic_gaussian_ar1_covariance(d, rho)

  base <- normalize_logistic_gaussian_theta(
    list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr),
    ambient_dim = ambient_dim,
    control = control
  )

  c(base, list(
    rho = rho,
    dSigma_drho = logistic_gaussian_ar1_covariance_derivative(d, rho)
  ))
}

rlogistic_gaussian_ar1 <- function(n, mu_ilr, rho, control = list()) {
  theta <- normalize_logistic_gaussian_ar1_theta(
    list(mu_ilr = mu_ilr, rho = rho),
    control = control
  )

  rlogistic_gaussian_simplex(
    n = n,
    mu_ilr = theta$mu_ilr,
    Sigma_ilr = theta$Sigma_ilr,
    control = control
  )
}

logistic_gaussian_ar1_loglik <- function(data, theta, weights = NULL,
                                          control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr

  theta <- normalize_logistic_gaussian_ar1_theta(
    theta,
    ambient_dim = normalized$ambient_dim,
    control = control
  )

  probability_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(z), nrow(z))
  } else {
    normalize_probability_weights(weights, nrow(z))
  }

  centered <- sweep(z, 2L, theta$mu_ilr, FUN = "-")
  Sigma_inv <- solve(theta$Sigma_ilr)
  quadratic <- rowSums((centered %*% Sigma_inv) * centered)
  logdet <- as.numeric(
    determinant(theta$Sigma_ilr, logarithm = TRUE)$modulus
  )

  sum(probability_weights * (
    -0.5 * theta$ilr_dim * log(2 * pi) -
      0.5 * logdet -
      0.5 * quadratic
  ))
}

logistic_gaussian_ar1_score_matrix <- function(data, theta, control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr

  theta <- normalize_logistic_gaussian_ar1_theta(
    theta,
    ambient_dim = normalized$ambient_dim,
    control = control
  )

  centered <- sweep(z, 2L, theta$mu_ilr, FUN = "-")
  Sigma_inv <- solve(theta$Sigma_ilr)
  D <- theta$dSigma_drho

  score_mu <- centered %*% Sigma_inv

  A <- Sigma_inv %*% D %*% Sigma_inv
  trace_term <- sum(diag(Sigma_inv %*% D))
  quadratic_term <- rowSums((centered %*% A) * centered)

  score_rho <- -0.5 * trace_term + 0.5 * quadratic_term

  out <- cbind(score_mu, rho = score_rho)
  colnames(out)[seq_len(theta$ilr_dim)] <-
    paste0("mu", seq_len(theta$ilr_dim))
  out
}

logistic_gaussian_ar1_fisher_information <- function(theta, control = list()) {
  theta <- normalize_logistic_gaussian_ar1_theta(theta, control = control)

  d <- theta$ilr_dim
  Sigma_inv <- solve(theta$Sigma_ilr)
  D <- theta$dSigma_drho

  information <- matrix(0, nrow = d + 1L, ncol = d + 1L)
  information[seq_len(d), seq_len(d)] <- Sigma_inv

  K <- Sigma_inv %*% D
  information[d + 1L, d + 1L] <-
    0.5 * sum(diag(K %*% K))

  information
}

# Paper convention: V = E[dot psi] = - Fisher information for the score.
logistic_gaussian_ar1_score_jacobian_V <- function(theta, control = list()) {
  -logistic_gaussian_ar1_fisher_information(theta, control = control)
}

logistic_gaussian_ar1_profile_moments <- function(centered,
                                                     probability_weights) {
  centered <- as.matrix(centered)
  probability_weights <- as.numeric(probability_weights)

  n <- nrow(centered)
  d <- ncol(centered)

  if (d < 2L) {
    stop("The AR(1) correlation parameter requires ilr dimension d >= 2.")
  }
  if (length(probability_weights) != n ||
      any(!is.finite(probability_weights)) ||
      any(probability_weights < 0)) {
    stop("Invalid probability weights in AR(1) profile likelihood.")
  }

  weight_sum <- sum(probability_weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    stop("AR(1) profile-likelihood weights must have positive total mass.")
  }
  probability_weights <- probability_weights / weight_sum

  A <- sum(
    probability_weights * rowSums(centered^2)
  )

  B <- sum(
    probability_weights *
      rowSums(
        centered[, seq_len(d - 1L), drop = FALSE] *
          centered[, 2:d, drop = FALSE]
      )
  )

  C <- if (d > 2L) {
    sum(
      probability_weights *
        rowSums(centered[, 2:(d - 1L), drop = FALSE]^2)
    )
  } else {
    0
  }

  list(
    d = d,
    A = A,
    B = B,
    C = C
  )
}

logistic_gaussian_ar1_profile_loglik_from_moments <- function(rho, moments) {
  rho <- as.numeric(rho)

  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    return(-Inf)
  }

  d <- moments$d
  A <- moments$A
  B <- moments$B
  C <- moments$C

  one_minus_rho2 <- 1 - rho^2

  -0.5 * (
    (d - 1L) * log1p(-rho^2) +
      (A + rho^2 * C - 2 * rho * B) / one_minus_rho2
  )
}

# Legacy implementation retained only for validation/benchmarking.
# This is the method previously used by fit_logistic_gaussian_ar1_theta():
# one-dimensional stats::optimize() with explicit covariance inversions.
fit_logistic_gaussian_ar1_rho_optimize <- function(
    centered,
    probability_weights,
    rho_bound = 0.995,
    rho_tol = 1e-9) {

  centered <- as.matrix(centered)
  d <- ncol(centered)

  profile_loglik <- function(rho) {
    Sigma <- logistic_gaussian_ar1_covariance(d, rho)
    Sigma_inv <- solve(Sigma)
    logdet <- as.numeric(
      determinant(Sigma, logarithm = TRUE)$modulus
    )
    quadratic <- rowSums((centered %*% Sigma_inv) * centered)

    sum(
      probability_weights *
        (-0.5 * logdet - 0.5 * quadratic)
    )
  }

  opt <- stats::optimize(
    f = profile_loglik,
    interval = c(-rho_bound, rho_bound),
    maximum = TRUE,
    tol = rho_tol
  )

  list(
    rho = opt$maximum,
    objective = opt$objective,
    method = "profile_optimize_1d"
  )
}

# Global maximization of the AR(1) profile likelihood.
#
# The stationary equation is
#
# -(d-1) rho^3
# + B rho^2
# + {(d-1) - (A+C)} rho
# + B = 0.
#
# We evaluate every real stationary point in the admissible interval,
# together with both interval endpoints. Therefore the returned solution
# is the global maximizer on [-rho_bound, rho_bound].
fit_logistic_gaussian_ar1_rho_global <- function(
    centered,
    probability_weights,
    rho_bound = 0.995,
    root_imag_tol = 1e-9) {

  centered <- as.matrix(centered)
  d <- ncol(centered)

  if (d < 2L) {
    stop("The AR(1) correlation parameter requires ilr dimension d >= 2.")
  }

  rho_bound <- as.numeric(rho_bound)
  if (length(rho_bound) != 1L || !is.finite(rho_bound) ||
      rho_bound <= 0 || rho_bound >= 1) {
    stop("`rho_bound` must lie in (0,1).")
  }

  root_imag_tol <- as.numeric(root_imag_tol)
  if (length(root_imag_tol) != 1L || !is.finite(root_imag_tol) ||
      root_imag_tol <= 0) {
    stop("`root_imag_tol` must be strictly positive.")
  }

  moments <- logistic_gaussian_ar1_profile_moments(
    centered = centered,
    probability_weights = probability_weights
  )

  A <- moments$A
  B <- moments$B
  C <- moments$C

  # polyroot() expects coefficients in increasing powers:
  # a0 + a1*rho + a2*rho^2 + a3*rho^3.
  coefficients <- c(
    B,
    (d - 1L) - (A + C),
    B,
    -(d - 1L)
  )

  roots <- polyroot(coefficients)

  real_roots <- Re(roots[
    abs(Im(roots)) <= root_imag_tol *
      pmax(1, abs(Re(roots)))
  ])

  interior_roots <- real_roots[
    real_roots > -rho_bound &
      real_roots < rho_bound
  ]

  candidates <- unique(c(
    -rho_bound,
    interior_roots,
    rho_bound
  ))

  objectives <- vapply(
    candidates,
    logistic_gaussian_ar1_profile_loglik_from_moments,
    numeric(1),
    moments = moments
  )

  best <- which.max(objectives)

  list(
    rho = candidates[[best]],
    objective = objectives[[best]],
    method = "global_profile_cubic_roots",
    stationary_roots = sort(interior_roots),
    candidates = candidates,
    candidate_objectives = objectives,
    at_boundary = abs(candidates[[best]]) >= rho_bound - 1e-12,
    moments = moments
  )
}

fit_logistic_gaussian_ar1_theta <- function(data, weights = NULL, null,
                                             control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_logistic_gaussian_ar1_theta(
      null$theta,
      ambient_dim = normalized$ambient_dim,
      control = control
    ))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (normalized$ilr_dim < 2L) {
    stop("The Logistic-Gaussian AR(1) model with unknown rho requires d >= 2.")
  }

  probability_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(z), nrow(z))
  } else {
    normalize_probability_weights(weights, nrow(z))
  }

  # Exact profiled MLE:
  # for every fixed rho, the Gaussian MLE of mu is the weighted mean.
  mu_hat <- colSums(z * probability_weights)
  centered <- sweep(z, 2L, mu_hat, FUN = "-")

  rho_bound <- as.numeric(
    control$logistic_gaussian_ar1_fit_bound %||% 0.995
  )
  if (length(rho_bound) != 1L || !is.finite(rho_bound) ||
      rho_bound <= 0 || rho_bound >= 1) {
    stop("`control$logistic_gaussian_ar1_fit_bound` must lie in (0,1).")
  }

  root_imag_tol <- as.numeric(
    control$logistic_gaussian_ar1_root_imag_tol %||% 1e-9
  )
  if (length(root_imag_tol) != 1L || !is.finite(root_imag_tol) ||
      root_imag_tol <= 0) {
    stop("`control$logistic_gaussian_ar1_root_imag_tol` must be positive.")
  }

  rho_fit <- fit_logistic_gaussian_ar1_rho_global(
    centered = centered,
    probability_weights = probability_weights,
    rho_bound = rho_bound,
    root_imag_tol = root_imag_tol
  )

  rho_hat <- rho_fit$rho

  fit <- normalize_logistic_gaussian_ar1_theta(
    list(mu_ilr = mu_hat, rho = rho_hat),
    ambient_dim = normalized$ambient_dim,
    control = control
  )

  fit$loglik <- logistic_gaussian_ar1_loglik(
    normalized,
    fit,
    probability_weights,
    control
  )

  fit$fit_diagnostics <- list(
    weighted = !is.null(weights),
    optimizer = rho_fit$method,
    covariance_shape = "AR1",
    rho_bound = rho_bound,
    rho_at_boundary = rho_fit$at_boundary,
    rho_stationary_roots = rho_fit$stationary_roots,
    rho_candidate_objectives = rho_fit$candidate_objectives
  )

  fit
}

prepare_logistic_gaussian_ar1_fast_multiplier <- function(spec, data, theta_hat,
                                                           ks_prep = NULL,
                                                           cvm_prep = NULL,
                                                           control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)

  theta_hat <- normalize_logistic_gaussian_ar1_theta(
    theta_hat,
    ambient_dim = normalized$ambient_dim,
    control = control
  )

  d <- theta_hat$ilr_dim
  par0 <- c(theta_hat$mu_ilr, theta_hat$rho)

  state_from_par <- function(par) {
    par <- as.numeric(par)
    if (length(par) != d + 1L) {
      stop("Logistic-Gaussian AR(1) parameter vector has incompatible dimension.")
    }

    normalize_logistic_gaussian_ar1_theta(
      list(
        mu_ilr = par[seq_len(d)],
        rho = par[[d + 1L]]
      ),
      ambient_dim = normalized$ambient_dim,
      control = control
    )
  }

  score_matrix_fn <- function(sample, par) {
    logistic_gaussian_ar1_score_matrix(
      sample,
      state_from_par(par),
      control = control
    )
  }

  sample_fn <- function(n_aux, par) {
    state <- state_from_par(par)

    rlogistic_gaussian_ar1(
      n = n_aux,
      mu_ilr = state$mu_ilr,
      rho = state$rho,
      control = control
    )
  }

  fisher_hat <- logistic_gaussian_ar1_fisher_information(
    theta_hat,
    control
  )

  prepared <- prepare_fast_multiplier_score_model(
    spec = spec,
    data = normalized$simplex,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn,
    vhat_fn = function(data, par0, S_obs, aux_sample, Psi_aux) fisher_hat
  )

  prepared$vhat_method <- "logistic_gaussian_ar1_analytic_fisher"
  prepared$paper_Vhat <- -fisher_hat
  prepared$paper_Vhat_inverse <- solve(-fisher_hat)
  prepared$paper_Vhat_method <- "analytic_ar1_score_jacobian"
  prepared$paper_Vhat_diagnostics <-
    fast_multiplier_matrix_condition_diagnostics(-fisher_hat)
  prepared$score_parameterization <- "beta=(mu_ilr,rho)"
  prepared$auxiliary_sampler <- "logistic_gaussian_ar1"
  prepared$covariance_shape <- theta_hat$Sigma_ilr

  prepared
}

make_logistic_gaussian_ar1_spec <- function() {
  # Reuse the existing Aitchison-distance and LG profile machinery.
  spec <- make_logistic_gaussian_spec(unknown_param = "both")

  base_profile_eval <- spec$profile_eval
  base_profile_matrix_eval <- spec$profile_matrix_eval
  base_sample_profile_matrix_eval <- spec$sample_profile_matrix_eval

  spec$name <- "logistic_gaussian_ar1_aitchison"

  spec$fit_theta <- function(data, weights = NULL, null, control = list()) {
    fit_logistic_gaussian_ar1_theta(data, weights, null, control)
  }

  spec$normalize_data <- normalize_logistic_gaussian_data

  spec$profile_eval <- function(omega, t, theta, control = list()) {
    fitted <- normalize_logistic_gaussian_ar1_theta(
      theta,
      control = control
    )

    base_profile_eval(
      omega,
      t,
      list(
        mu_ilr = fitted$mu_ilr,
        Sigma_ilr = fitted$Sigma_ilr
      ),
      control
    )
  }

  if (is.function(base_profile_matrix_eval)) {
    spec$profile_matrix_eval <- function(omega_grid, t_grid, theta,
                                         control = list()) {
      fitted <- normalize_logistic_gaussian_ar1_theta(
        theta,
        control = control
      )

      base_profile_matrix_eval(
        omega_grid,
        t_grid,
        list(
          mu_ilr = fitted$mu_ilr,
          Sigma_ilr = fitted$Sigma_ilr
        ),
        control
      )
    }
  }

  if (is.function(base_sample_profile_matrix_eval)) {
    spec$sample_profile_matrix_eval <- function(data, distance_matrix, theta,
                                                control = list()) {
      fitted <- normalize_logistic_gaussian_ar1_theta(
        theta,
        control = control
      )

      base_sample_profile_matrix_eval(
        data,
        distance_matrix,
        list(
          mu_ilr = fitted$mu_ilr,
          Sigma_ilr = fitted$Sigma_ilr
        ),
        control
      )
    }
  }

  spec$score_matrix <- logistic_gaussian_ar1_score_matrix
  spec$fisher_information <- logistic_gaussian_ar1_fisher_information
  spec$score_jacobian_V <- logistic_gaussian_ar1_score_jacobian_V
  spec$weighted_mle <- TRUE

  spec$fast_multiplier_prepare <- function(data, theta_hat,
                                            ks_prep = NULL,
                                            cvm_prep = NULL,
                                            control = list()) {
    prepare_logistic_gaussian_ar1_fast_multiplier(
      spec = make_logistic_gaussian_ar1_spec(),
      data = data,
      theta_hat = theta_hat,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      control = control
    )
  }

  validate_model_spec(spec)
  spec
}
