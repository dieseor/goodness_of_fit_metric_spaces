# Spherical multivariate Gaussian model.
#
#   X ~ N_d(mu, sigma^2 I_d),  mu in R^d, sigma > 0.
#
# This adapter deliberately reuses the unrestricted Gaussian distance-profile
# evaluator.  Only the fitted parameter space, score and information matrix
# differ from the unrestricted multivariate-normal model.

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_normal_sigma_Id <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_normal_sigma_Id <- model_specs_candidates_normal_sigma_Id[
    file.exists(model_specs_candidates_normal_sigma_Id)
  ][1L]
  if (is.na(model_specs_path_normal_sigma_Id)) {
    stop("Could not locate bootstrap/model_specs.R.")
  }
  source(model_specs_path_normal_sigma_Id)
}

normalize_normal_sigma_Id_data <- function(data, control = list()) {
  normalize_mvnormal_data(data, control)
}

normalize_normal_sigma_Id_theta <- function(theta,
                                              ambient_dim = NULL,
                                              control = list()) {
  if (!is.list(theta)) {
    stop("Normal-sigma-Id theta must be a list containing `mu` and `sigma`.")
  }
  mu <- as.numeric(theta$mu %||% theta$theta)
  sigma <- suppressWarnings(as.numeric(theta$sigma))
  if (!length(mu) || any(!is.finite(mu))) {
    stop("Normal-sigma-Id theta requires a finite, non-empty `mu` vector.")
  }
  ambient_dim <- as.integer(ambient_dim %||% theta$ambient_dim %||% length(mu))
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 1L ||
      length(mu) != ambient_dim) {
    stop("Normal-sigma-Id theta has incompatible ambient dimension.")
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("Normal-sigma-Id `sigma` must be a strictly positive finite scalar.")
  }
  base <- normalize_mvnormal_theta(
    list(mu = mu, Sigma = diag(sigma^2, ambient_dim)),
    ambient_dim = ambient_dim,
    control = control
  )
  c(base, list(sigma = sigma, sigma2 = sigma^2))
}

rnormal_sigma_Id <- function(n, mu, sigma, control = list()) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n <= 0L) {
    stop("`n` must be a strictly positive integer.")
  }
  theta <- normalize_normal_sigma_Id_theta(list(mu = mu, sigma = sigma), control = control)
  mvtnorm::rmvnorm(n, mean = theta$mu, sigma = theta$Sigma)
}

normal_sigma_Id_loglik <- function(data, theta, weights = NULL,
                                    control = list()) {
  x <- normalize_normal_sigma_Id_data(data, control)
  theta <- normalize_normal_sigma_Id_theta(theta, ambient_dim = ncol(x), control = control)
  probability_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    normalize_probability_weights(weights, nrow(x))
  }
  z <- sweep(x, 2L, theta$mu, FUN = "-")
  sum(probability_weights * (
    -0.5 * ncol(x) * log(2 * pi * theta$sigma2) -
      rowSums(z^2) / (2 * theta$sigma2)
  ))
}

normal_sigma_Id_score_matrix <- function(data, theta, control = list()) {
  x <- normalize_normal_sigma_Id_data(data, control)
  theta <- normalize_normal_sigma_Id_theta(theta, ambient_dim = ncol(x), control = control)
  z <- sweep(x, 2L, theta$mu, FUN = "-")
  cbind(
    z / theta$sigma2,
    sigma = -ncol(x) / theta$sigma + rowSums(z^2) / theta$sigma^3
  )
}

normal_sigma_Id_fisher_information <- function(theta, control = list()) {
  theta <- normalize_normal_sigma_Id_theta(theta, control = control)
  d <- theta$ambient_dim
  information <- diag(c(rep.int(1 / theta$sigma2, d), 2 * d / theta$sigma2))
  information
}

# Paper convention: V = E[dot psi] = - Fisher information.
normal_sigma_Id_score_jacobian_V <- function(theta, control = list()) {
  -normal_sigma_Id_fisher_information(theta, control = control)
}

fit_normal_sigma_Id_theta <- function(data, weights = NULL, null,
                                       control = list()) {
  x <- normalize_normal_sigma_Id_data(data, control)
  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }
  if (identical(null$type, "simple")) {
    return(normalize_normal_sigma_Id_theta(
      null$theta, ambient_dim = ncol(x), control = control
    ))
  }
  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }
  probability_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(x), nrow(x))
  } else {
    normalize_probability_weights(weights, nrow(x))
  }
  mu_hat <- colSums(x * probability_weights)
  z <- sweep(x, 2L, mu_hat, FUN = "-")
  sigma2_hat <- sum(probability_weights * rowSums(z^2)) / ncol(x)
  if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
    stop("Normal-sigma-Id MLE stopped: the fitted variance is not strictly positive.")
  }
  fit <- normalize_normal_sigma_Id_theta(
    list(mu = mu_hat, sigma = sqrt(sigma2_hat)),
    ambient_dim = ncol(x), control = control
  )
  fit$loglik <- normal_sigma_Id_loglik(x, fit, probability_weights, control)
  fit$fit_diagnostics <- list(weighted = !is.null(weights), optimizer = "closed_form_mle")
  fit
}

prepare_normal_sigma_Id_fast_multiplier <- function(spec, data, theta_hat,
                                                      ks_prep = NULL,
                                                      cvm_prep = NULL,
                                                      control = list()) {
  x <- normalize_normal_sigma_Id_data(data, control)
  theta_hat <- normalize_normal_sigma_Id_theta(
    theta_hat, ambient_dim = ncol(x), control = control
  )
  d <- theta_hat$ambient_dim
  par0 <- c(theta_hat$mu, theta_hat$sigma)
  state_from_par <- function(par) {
    par <- as.numeric(par)
    if (length(par) != d + 1L) {
      stop("Normal-sigma-Id fast-multiplier parameters have incompatible dimension.")
    }
    normalize_normal_sigma_Id_theta(
      list(mu = par[seq_len(d)], sigma = par[[d + 1L]]),
      ambient_dim = d, control = control
    )
  }
  score_matrix_fn <- function(sample, par) {
    normal_sigma_Id_score_matrix(sample, state_from_par(par), control = control)
  }
  sample_fn <- function(n_aux, par) {
    state <- state_from_par(par)
    rnormal_sigma_Id(n_aux, state$mu, state$sigma, control = control)
  }
  fisher_hat <- normal_sigma_Id_fisher_information(theta_hat, control)
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
  prepared$vhat_method <- "normal_sigma_Id_analytic_fisher"
  prepared$paper_Vhat <- -fisher_hat
  prepared$paper_Vhat_inverse <- solve(-fisher_hat)
  prepared$paper_Vhat_method <- "analytic_spherical_score_jacobian"
  prepared$paper_Vhat_diagnostics <- fast_multiplier_matrix_condition_diagnostics(-fisher_hat)
  prepared$score_parameterization <- "beta=(mu,sigma)"
  prepared$auxiliary_sampler <- "normal_sigma_Id"
  prepared
}

make_normal_sigma_Id_spec <- function() {
  # The general-Gaussian adapter already owns the Euclidean distance and the
  # noncentral-chi-square distance-profile evaluation.  Reusing it avoids any
  # second implementation of the same profile formula.
  spec <- make_mvnormal_spec(unknown_param = "both")
  base_profile_eval <- spec$profile_eval
  base_profile_matrix_eval <- spec$profile_matrix_eval
  base_sample_profile_matrix_eval <- spec$sample_profile_matrix_eval
  spec$name <- "normal_sigma_Id_euclidean"
  spec$fit_theta <- function(data, weights = NULL, null, control = list()) {
    fit_normal_sigma_Id_theta(data, weights, null, control)
  }
  spec$normalize_data <- normalize_normal_sigma_Id_data
  spec$profile_eval <- function(omega, t, theta, control = list()) {
    fitted <- normalize_normal_sigma_Id_theta(theta, control = control)
    base_profile_eval(omega, t, list(mu = fitted$mu, Sigma = fitted$Sigma), control)
  }
  spec$profile_matrix_eval <- function(omega_grid, t_grid, theta, control = list()) {
    fitted <- normalize_normal_sigma_Id_theta(theta, control = control)
    base_profile_matrix_eval(
      omega_grid, t_grid, list(mu = fitted$mu, Sigma = fitted$Sigma), control
    )
  }
  spec$sample_profile_matrix_eval <- function(data, distance_matrix, theta,
                                               control = list()) {
    fitted <- normalize_normal_sigma_Id_theta(theta, control = control)
    base_sample_profile_matrix_eval(
      data, distance_matrix, list(mu = fitted$mu, Sigma = fitted$Sigma), control
    )
  }
  spec$score_matrix <- normal_sigma_Id_score_matrix
  spec$fisher_information <- normal_sigma_Id_fisher_information
  spec$score_jacobian_V <- normal_sigma_Id_score_jacobian_V
  spec$weighted_mle <- TRUE
  spec$fast_multiplier_prepare <- function(data, theta_hat,
                                            ks_prep = NULL,
                                            cvm_prep = NULL,
                                            control = list()) {
    prepare_normal_sigma_Id_fast_multiplier(
      spec = make_normal_sigma_Id_spec(), data = data, theta_hat = theta_hat,
      ks_prep = ks_prep, cvm_prep = cvm_prep, control = control
    )
  }
  validate_model_spec(spec)
  spec
}
