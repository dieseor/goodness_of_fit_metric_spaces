# Restricted Gaussian model with a mean-aligned rank-one covariance spike.
#
#   X ~ N_d(theta, I_d + lambda u(theta) u(theta)^T),
#   u(theta) = theta / ||theta||, lambda > 0.
#
# This is deliberately distinct from both the unrestricted multivariate
# Gaussian adapter and PPCA: theta controls the mean and the spike direction.

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_restricted_spiked <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_restricted_spiked <- model_specs_candidates_restricted_spiked[
    file.exists(model_specs_candidates_restricted_spiked)
  ][1L]
  if (is.na(model_specs_path_restricted_spiked)) {
    stop("Could not locate bootstrap/model_specs.R.")
  }
  source(model_specs_path_restricted_spiked)
}

restricted_spiked_normal_radius_tolerance <- function(x) {
  x <- as.matrix(x)
  if (!length(x) || any(!is.finite(x))) {
    stop("Cannot determine the numerical radius tolerance from non-finite data.")
  }
  # This is solely a floating-point diagnostic threshold.  It is not a model
  # constraint or a regularisation: a fit below it is stopped explicitly.
  sqrt(.Machine$double.eps) * max(1, sqrt(mean(rowSums(x^2))))
}

normalize_restricted_spiked_normal_data <- function(data, control = list()) {
  normalize_mvnormal_data(data, control)
}

normalize_restricted_spiked_normal_theta <- function(theta,
                                                      ambient_dim = NULL,
                                                      control = list()) {
  if (!is.list(theta)) {
    stop("Restricted-spiked normal theta must be a list containing `theta` and `lambda`.")
  }
  theta_vector <- as.numeric(theta$theta %||% theta$mu)
  lambda <- suppressWarnings(as.numeric(theta$lambda))
  if (!length(theta_vector) || any(!is.finite(theta_vector))) {
    stop("Restricted-spiked normal theta requires a finite non-empty `theta` vector.")
  }
  ambient_dim <- as.integer(ambient_dim %||% theta$ambient_dim %||% length(theta_vector))
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 1L ||
      length(theta_vector) != ambient_dim) {
    stop("Restricted-spiked normal theta has incompatible ambient dimension.")
  }
  if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("Restricted-spiked normal `lambda` must be a strictly positive finite scalar.")
  }
  radius <- sqrt(sum(theta_vector^2))
  if (!is.finite(radius) || radius == 0) {
    stop("Restricted-spiked normal theta must be non-zero.")
  }
  u <- theta_vector / radius
  Sigma <- diag(ambient_dim) + lambda * tcrossprod(u)
  list(
    theta = theta_vector,
    mu = theta_vector,
    lambda = lambda,
    u = u,
    radius = radius,
    Sigma = Sigma,
    ambient_dim = ambient_dim
  )
}

restricted_spiked_normal_covariance <- function(theta, lambda = NULL,
                                                 ambient_dim = NULL,
                                                 control = list()) {
  theta_input <- if (is.list(theta)) theta else list(theta = theta, lambda = lambda)
  normalize_restricted_spiked_normal_theta(
    theta_input,
    ambient_dim = ambient_dim,
    control = control
  )$Sigma
}

rrestricted_spiked_normal <- function(n, theta, lambda, control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer.")
  }
  normalized <- normalize_restricted_spiked_normal_theta(
    list(theta = theta, lambda = lambda), control = control
  )
  mvtnorm::rmvnorm(n, mean = normalized$theta, sigma = normalized$Sigma)
}

restricted_spiked_normal_loglik <- function(data, theta, weights = NULL,
                                             control = list()) {
  x <- normalize_restricted_spiked_normal_data(data, control)
  normalized <- normalize_restricted_spiked_normal_theta(
    theta, ambient_dim = ncol(x), control = control
  )
  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    normalize_probability_weights(weights, nrow(x))
  }
  centered <- sweep(x, 2L, normalized$theta, FUN = "-")
  projections <- drop(centered %*% normalized$u)
  quadratic <- rowSums(centered^2) -
    normalized$lambda / (1 + normalized$lambda) * projections^2
  sum(prob_weights * (
    -0.5 * (ncol(x) * log(2 * pi) + log1p(normalized$lambda) + quadratic)
  ))
}

restricted_spiked_normal_score_matrix <- function(data, theta,
                                                   control = list()) {
  x <- normalize_restricted_spiked_normal_data(data, control)
  normalized <- normalize_restricted_spiked_normal_theta(
    theta, ambient_dim = ncol(x), control = control
  )
  y <- sweep(x, 2L, normalized$theta, FUN = "-")
  s <- drop(y %*% normalized$u)
  perpendicular_y <- y - tcrossprod(s, normalized$u)
  lambda <- normalized$lambda
  score_theta <- y + lambda / (1 + lambda) * (
    perpendicular_y * (s / normalized$radius) - tcrossprod(s, normalized$u)
  )
  score_lambda <- -0.5 / (1 + lambda) + 0.5 * s^2 / (1 + lambda)^2
  cbind(score_theta, lambda = score_lambda)
}

restricted_spiked_normal_fisher_information <- function(theta,
                                                        control = list()) {
  fitted <- normalize_restricted_spiked_normal_theta(theta, control = control)
  d <- fitted$ambient_dim
  lambda <- fitted$lambda
  uu <- tcrossprod(fitted$u)
  P <- diag(d) - uu
  information_theta <-
    uu / (1 + lambda) +
    (1 + lambda^2 / ((1 + lambda) * fitted$radius^2)) * P
  information <- matrix(0, nrow = d + 1L, ncol = d + 1L)
  information[seq_len(d), seq_len(d)] <- information_theta
  information[d + 1L, d + 1L] <- 1 / (2 * (1 + lambda)^2)
  0.5 * (information + t(information))
}

# V is the derivative of the expected likelihood score, using the paper's
# convention.  The generic fast-multiplier kernel internally stores -V, i.e.
# the positive Fisher information; see prepare_restricted_spiked_normal_fast_multiplier().
restricted_spiked_normal_score_jacobian_V <- function(theta,
                                                       control = list()) {
  -restricted_spiked_normal_fisher_information(theta, control = control)
}

restricted_spiked_normal_fast_state_from_par <- function(par,
                                                          ambient_dim,
                                                          control = list()) {
  par <- as.numeric(par)
  if (length(par) != ambient_dim + 1L) {
    stop("Restricted-spiked fast-multiplier parameters have incompatible dimension.")
  }
  normalize_restricted_spiked_normal_theta(
    list(theta = par[seq_len(ambient_dim)], lambda = par[[ambient_dim + 1L]]),
    ambient_dim = ambient_dim,
    control = control
  )
}

prepare_restricted_spiked_normal_fast_multiplier <- function(spec,
                                                              data,
                                                              theta_hat,
                                                              ks_prep = NULL,
                                                              cvm_prep = NULL,
                                                              control = list()) {
  x <- normalize_restricted_spiked_normal_data(data, control)
  theta_hat <- normalize_restricted_spiked_normal_theta(
    theta_hat, ambient_dim = ncol(x), control = control
  )
  d <- theta_hat$ambient_dim
  par0 <- c(theta_hat$theta, theta_hat$lambda)

  state_from_par <- function(par) {
    restricted_spiked_normal_fast_state_from_par(
      par, ambient_dim = d, control = control
    )
  }
  score_matrix_fn <- function(sample, par) {
    state <- state_from_par(par)
    restricted_spiked_normal_score_matrix(sample, state, control = control)
  }
  sample_fn <- function(n_aux, par) {
    state <- state_from_par(par)
    rrestricted_spiked_normal(
      n = n_aux,
      theta = state$theta,
      lambda = state$lambda,
      control = control
    )
  }
  fisher_hat <- restricted_spiked_normal_fisher_information(theta_hat, control)
  paper_Vhat <- -fisher_hat

  # `prepare_fast_multiplier_score_model()` represents the same correction as
  # the paper with a positive information matrix: its internal expression is
  # y - psi I^{-1} dotF.  Since V_paper = -I, this equals
  # y + dotF^T V_paper^{-1} psi after transposition.
  prepared <- prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn,
    vhat_fn = function(data, par0, S_obs, aux_sample, Psi_aux) fisher_hat
  )
  prepared$vhat_method <- "restricted_spiked_analytic_fisher"
  prepared$paper_Vhat <- paper_Vhat
  prepared$paper_Vhat_inverse <- solve(paper_Vhat)
  prepared$paper_Vhat_method <- "analytic_restricted_score_jacobian"
  prepared$paper_Vhat_diagnostics <- fast_multiplier_matrix_condition_diagnostics(
    paper_Vhat
  )
  prepared$score_parameterization <- "beta=(theta,lambda)"
  prepared$auxiliary_sampler <- "restricted_spiked_normal"
  prepared
}

restricted_spiked_normal_profile_at_tau <- function(tau, xbar, S) {
  tau <- as.numeric(tau)
  if (length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("The restricted-spiked profile coordinate `tau` must be finite and non-negative.")
  }
  c_tau <- -expm1(-tau)
  A <- tcrossprod(xbar) + c_tau * S
  eig <- eigen(0.5 * (A + t(A)), symmetric = TRUE)
  u <- eig$vectors[, 1L]
  projection <- drop(crossprod(u, xbar))
  if (projection < 0) u <- -u
  list(
    tau = tau,
    lambda = expm1(tau),
    value = eig$values[[1L]] - tau,
    u = u,
    r = drop(crossprod(u, xbar)),
    leading_eigenvalue = eig$values[[1L]]
  )
}

restricted_spiked_normal_profile_upper_tau <- function(xbar, S) {
  q0 <- sum(xbar^2)
  C <- eigen(0.5 * (tcrossprod(xbar) + S + t(tcrossprod(xbar) + S)),
             symmetric = TRUE, only.values = TRUE)$values[[1L]]
  # Q(tau) <= C - tau while Q(0) = q0.  The displayed upper endpoint
  # therefore puts the entire remaining tail at least five log-likelihood
  # units below the boundary candidate.
  tau_upper <- max(1, C - q0 + 5)
  if (!is.finite(tau_upper) || tau_upper >= log(.Machine$double.xmax)) {
    stop("The restricted-spiked profile requires a lambda beyond floating-point range.")
  }
  list(tau_upper = tau_upper, q0 = q0, tail_bound_constant = C)
}

restricted_spiked_normal_profile_maximize <- function(xbar, S,
                                                       grid_size = 401L,
                                                       tol = sqrt(.Machine$double.eps)) {
  grid_size <- as.integer(grid_size)
  if (length(grid_size) != 1L || !is.finite(grid_size) || grid_size < 5L) {
    stop("`grid_size` must be an integer of at least five profile points.")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a strictly positive finite scalar.")
  }
  bound <- restricted_spiked_normal_profile_upper_tau(xbar, S)
  # A linear scan resolves interior maxima, while the quadratic scan also
  # diagnoses a profile maximum at the excluded boundary lambda = 0.  This is
  # necessary because a narrow positive-lambda maximum can lie before the
  # first point of a purely linear grid.
  grid_fraction <- seq(0, 1, length.out = grid_size)
  tau_grid <- sort(c(
    bound$tau_upper * grid_fraction,
    bound$tau_upper * grid_fraction^2
  ))
  tau_grid <- tau_grid[c(
    TRUE,
    diff(tau_grid) > sqrt(.Machine$double.eps) * max(1, bound$tau_upper)
  )]
  n_grid <- length(tau_grid)
  profile_grid <- lapply(tau_grid, restricted_spiked_normal_profile_at_tau,
                         xbar = xbar, S = S)
  values <- vapply(profile_grid, `[[`, numeric(1L), "value")
  candidate_idx <- unique(c(
    1L,
    n_grid,
    which(values[2:(n_grid - 1L)] >= values[1:(n_grid - 2L)] &
            values[2:(n_grid - 1L)] >= values[3:n_grid]) + 1L
  ))
  refined <- lapply(candidate_idx, function(index) {
    if (index %in% c(1L, n_grid)) return(profile_grid[[index]])
    optimum <- stats::optimize(
      f = function(tau) -restricted_spiked_normal_profile_at_tau(tau, xbar, S)$value,
      interval = c(tau_grid[[index - 1L]], tau_grid[[index + 1L]]),
      tol = tol
    )
    restricted_spiked_normal_profile_at_tau(optimum$minimum, xbar, S)
  })
  all_candidates <- c(profile_grid[candidate_idx], refined)
  candidate_values <- vapply(all_candidates, `[[`, numeric(1L), "value")
  best <- all_candidates[[which.max(candidate_values)]]
  list(
    best = best,
    tau_grid = tau_grid,
    profile_values = values,
    local_candidate_indices = candidate_idx,
    candidate_values = candidate_values,
    tau_upper = bound$tau_upper,
    q0 = bound$q0,
    tail_bound_constant = bound$tail_bound_constant,
    grid_size = grid_size,
    profile_points_evaluated = n_grid,
    tol = tol
  )
}

fit_restricted_spiked_normal_theta <- function(data, weights = NULL, null,
                                               control = list()) {
  x <- normalize_restricted_spiked_normal_data(data, control)
  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }
  if (identical(null$type, "simple")) {
    return(normalize_restricted_spiked_normal_theta(
      null$theta, ambient_dim = ncol(x), control = control
    ))
  }
  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }
  prob_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    normalize_probability_weights(weights, nrow(x))
  }
  xbar <- colSums(x * prob_weights)
  S <- mvnormal_weighted_covariance(x, prob_weights, xbar)
  profile <- restricted_spiked_normal_profile_maximize(
    xbar = xbar,
    S = S,
    grid_size = control$restricted_spiked_profile_grid_size %||% 401L,
    tol = control$restricted_spiked_profile_tol %||% sqrt(.Machine$double.eps)
  )
  radius_tolerance <- restricted_spiked_normal_radius_tolerance(x)
  if (!is.finite(profile$best$r) || profile$best$r <= radius_tolerance) {
    stop(sprintf(
      paste0(
        "Restricted-spiked normal MLE stopped: the profiled radius r = u^T xbar ",
        "is numerically too close to zero (r = %.6e; tolerance = %.6e)."
      ),
      profile$best$r, radius_tolerance
    ))
  }
  if (!is.finite(profile$best$lambda) || profile$best$lambda <= 0) {
    stop(sprintf(
      paste0(
        "Restricted-spiked normal MLE stopped: the profiled lambda is numerically ",
        "at the excluded boundary lambda = 0 (lambda = %.6e)."
      ),
      profile$best$lambda
    ))
  }
  theta_hat <- profile$best$r * profile$best$u
  fitted <- normalize_restricted_spiked_normal_theta(
    list(theta = theta_hat, lambda = profile$best$lambda),
    ambient_dim = ncol(x), control = control
  )
  fitted$loglik <- restricted_spiked_normal_loglik(x, fitted, prob_weights, control)
  fitted$fit_diagnostics <- list(
    optimizer = "profile_grid_local_refinement",
    profile_tau = profile$best$tau,
    profile_value = profile$best$value,
    profile_grid_size = profile$grid_size,
    profile_tolerance = profile$tol,
    profile_tau_upper = profile$tau_upper,
    profile_tail_bound_constant = profile$tail_bound_constant,
    profile_boundary_value = profile$q0,
    profile_local_candidate_indices = profile$local_candidate_indices,
    profile_candidate_values = profile$candidate_values,
    profiled_radius = profile$best$r,
    radius_tolerance = radius_tolerance,
    weighted = !is.null(weights)
  )
  fitted
}

make_restricted_spiked_normal_spec <- function() {
  new_model_spec(
    name = "restricted_spiked_normal_euclidean",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_restricted_spiked_normal_theta(data, weights, null, control)
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_restricted_spiked_normal_data(data, control)
      omega_matrix <- normalize_restricted_spiked_normal_data(omega, control)
      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }
      sqrt(pmax(
        outer(rowSums(x^2), rowSums(omega_matrix^2), FUN = "+") -
          2 * (x %*% t(omega_matrix)),
        0
      ))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      fitted <- normalize_restricted_spiked_normal_theta(theta, control = control)
      omega <- as.numeric(normalize_restricted_spiked_normal_data(omega, control)[1L, ])
      spectral <- eigen(fitted$Sigma, symmetric = TRUE)
      evaluate_mvnorm_distance_profile(
        shift = fitted$theta - omega,
        t_values = as.numeric(t),
        eigenvalues_full = spectral$values,
        eigenvectors_full = spectral$vectors,
        positive_idx = rep(TRUE, fitted$ambient_dim),
        control = mvnormal_quadform_with_label(control, "restricted-spiked normal distance profile")
      )
    },
    normalize_data = normalize_restricted_spiked_normal_data,
    n_obs = function(data, control = list()) nrow(normalize_restricted_spiked_normal_data(data, control)),
    observation_at = function(data, idx, control = list()) {
      normalize_restricted_spiked_normal_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        fitted <- normalize_restricted_spiked_normal_theta(theta, control = control)
        omega_matrix <- normalize_restricted_spiked_normal_data(omega_grid, control)
        if (ncol(omega_matrix) != fitted$ambient_dim) {
          stop("`omega_grid` has incompatible ambient dimension.")
        }
        spectral <- eigen(fitted$Sigma, symmetric = TRUE)
        shift_matrix <- matrix(
          rep(fitted$theta, each = nrow(omega_matrix)),
          nrow = nrow(omega_matrix), ncol = fitted$ambient_dim
        ) - omega_matrix
        t_matrix <- matrix(
          rep(as.numeric(t_grid), each = nrow(omega_matrix)),
          nrow = nrow(omega_matrix), ncol = length(t_grid)
        )
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = t_matrix,
          eigenvalues_full = spectral$values,
          eigenvectors_full = spectral$vectors,
          positive_idx = rep(TRUE, fitted$ambient_dim),
          control = mvnormal_quadform_with_label(
            control, "restricted-spiked normal distance-profile grid"
          )
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta,
                                             control = list()) {
        fitted <- normalize_restricted_spiked_normal_theta(theta, control = control)
        x <- normalize_restricted_spiked_normal_data(data, control)
        if (ncol(x) != fitted$ambient_dim) {
          stop("`data` has incompatible ambient dimension.")
        }
        spectral <- eigen(fitted$Sigma, symmetric = TRUE)
        shift_matrix <- matrix(
          rep(fitted$theta, each = nrow(x)),
          nrow = nrow(x), ncol = fitted$ambient_dim
        ) - x
        evaluate_mvnorm_distance_profile_matrix(
          shift_matrix = shift_matrix,
          t_matrix = as.matrix(distance_matrix),
          eigenvalues_full = spectral$values,
          eigenvectors_full = spectral$vectors,
          positive_idx = rep(TRUE, fitted$ambient_dim),
          control = mvnormal_quadform_with_label(
            control, "restricted-spiked normal sample distance-profile grid"
          )
        )
      },
      distance_type = "euclidean",
      weighted_mle = TRUE,
      score_matrix = restricted_spiked_normal_score_matrix,
      fisher_information = restricted_spiked_normal_fisher_information,
      score_jacobian_V = restricted_spiked_normal_score_jacobian_V,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_restricted_spiked_normal_fast_multiplier(
          spec = make_restricted_spiked_normal_spec(),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control
        )
      }
    )
  )
}
