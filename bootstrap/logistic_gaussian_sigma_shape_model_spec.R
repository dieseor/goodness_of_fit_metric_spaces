# Restricted Logistic-Gaussian model used only by the two Section 7 scenarios
# based on Maphosa et al. (2026):
#
#   ilr(X) ~ N_d(mu, sigma^2 A_d),
#   A_d = diag(1, ..., 1, 1/(d + 1)),
#   mu in R^d, sigma > 0.
#
# A_d is the covariance shape obtained when the logistic-normal model
#   Y ~ N_d(0, I_d),  X_k = exp(Y_k) / (1 + sum_j exp(Y_j))
# is expressed in the ilr coordinates used by this repository.
#
# This file is deliberately self-contained relative to the existing
# Logistic-Gaussian adapter. It does not modify the unrestricted adapter.

if (!exists("make_logistic_gaussian_spec", mode = "function")) {
  lg_sigma_shape_candidates <- c(
    file.path("bootstrap", "logistic_gaussian_model_spec.R"),
    "logistic_gaussian_model_spec.R",
    file.path("..", "bootstrap", "logistic_gaussian_model_spec.R"),
    file.path("..", "..", "bootstrap", "logistic_gaussian_model_spec.R")
  )
  lg_sigma_shape_path <- lg_sigma_shape_candidates[
    file.exists(lg_sigma_shape_candidates)
  ][1L]
  if (is.na(lg_sigma_shape_path)) {
    stop("Could not locate bootstrap/logistic_gaussian_model_spec.R.")
  }
  source(lg_sigma_shape_path)
}

logistic_gaussian_maphosa_shape <- function(ilr_dim) {
  d <- as.integer(ilr_dim)
  if (length(d) != 1L || !is.finite(d) || d < 1L) {
    stop("`ilr_dim` must be a strictly positive integer.")
  }
  diag(c(rep.int(1, max(0L, d - 1L)), 1 / (d + 1)), nrow = d, ncol = d)
}

logistic_gaussian_maphosa_shape_inverse <- function(ilr_dim) {
  d <- as.integer(ilr_dim)
  if (length(d) != 1L || !is.finite(d) || d < 1L) {
    stop("`ilr_dim` must be a strictly positive integer.")
  }
  diag(c(rep.int(1, max(0L, d - 1L)), d + 1), nrow = d, ncol = d)
}

normalize_logistic_gaussian_sigma_shape_theta <- function(theta,
                                                           ambient_dim = NULL,
                                                           control = list()) {
  if (!is.list(theta)) {
    stop("Logistic-Gaussian sigma-shape theta must be a list containing `mu_ilr` and `sigma`.")
  }

  mu_ilr <- as.numeric(theta$mu_ilr %||% theta$mu)
  sigma <- suppressWarnings(as.numeric(theta$sigma))
  if (!length(mu_ilr) || any(!is.finite(mu_ilr))) {
    stop("Logistic-Gaussian sigma-shape theta requires a finite, non-empty `mu_ilr` vector.")
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("Logistic-Gaussian sigma-shape `sigma` must be a strictly positive finite scalar.")
  }

  ambient_dim <- as.integer(ambient_dim %||% theta$ambient_dim %||% (length(mu_ilr) + 1L))
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 2L ||
      length(mu_ilr) != ambient_dim - 1L) {
    stop("Logistic-Gaussian sigma-shape theta has incompatible ambient dimension.")
  }

  d <- ambient_dim - 1L
  shape_ilr <- logistic_gaussian_maphosa_shape(d)
  Sigma_ilr <- sigma^2 * shape_ilr
  base <- normalize_logistic_gaussian_theta(
    list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr),
    ambient_dim = ambient_dim,
    control = control
  )

  c(base, list(
    sigma = sigma,
    sigma2 = sigma^2,
    shape_ilr = shape_ilr,
    shape_inverse_ilr = logistic_gaussian_maphosa_shape_inverse(d)
  ))
}

rlogistic_gaussian_sigma_shape <- function(n, mu_ilr, sigma = 1, control = list()) {
  theta <- normalize_logistic_gaussian_sigma_shape_theta(
    list(mu_ilr = mu_ilr, sigma = sigma),
    control = control
  )
  rlogistic_gaussian_simplex(
    n = n,
    mu_ilr = theta$mu_ilr,
    Sigma_ilr = theta$Sigma_ilr,
    control = control
  )
}

logistic_gaussian_sigma_shape_loglik <- function(data, theta, weights = NULL,
                                                  control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr
  theta <- normalize_logistic_gaussian_sigma_shape_theta(
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
  transformed <- centered %*% theta$shape_inverse_ilr
  quadratic <- rowSums(transformed * centered)
  d <- theta$ilr_dim
  logdet_shape <- as.numeric(determinant(theta$shape_ilr, logarithm = TRUE)$modulus)

  sum(probability_weights * (
    -0.5 * d * log(2 * pi) -
      d * log(theta$sigma) -
      0.5 * logdet_shape -
      quadratic / (2 * theta$sigma2)
  ))
}

logistic_gaussian_sigma_shape_score_matrix <- function(data, theta, control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr
  theta <- normalize_logistic_gaussian_sigma_shape_theta(
    theta,
    ambient_dim = normalized$ambient_dim,
    control = control
  )

  centered <- sweep(z, 2L, theta$mu_ilr, FUN = "-")
  shape_weighted <- centered %*% theta$shape_inverse_ilr
  quadratic <- rowSums(shape_weighted * centered)
  score_mu <- shape_weighted / theta$sigma2
  score_sigma <- -theta$ilr_dim / theta$sigma + quadratic / theta$sigma^3

  out <- cbind(score_mu, sigma = score_sigma)
  colnames(out)[seq_len(theta$ilr_dim)] <- paste0("mu", seq_len(theta$ilr_dim))
  out
}

logistic_gaussian_sigma_shape_fisher_information <- function(theta, control = list()) {
  theta <- normalize_logistic_gaussian_sigma_shape_theta(theta, control = control)
  d <- theta$ilr_dim
  information <- matrix(0, nrow = d + 1L, ncol = d + 1L)
  information[seq_len(d), seq_len(d)] <- theta$shape_inverse_ilr / theta$sigma2
  information[d + 1L, d + 1L] <- 2 * d / theta$sigma2
  information
}

# Paper convention: V = E[dot psi] = - Fisher information for the score.
logistic_gaussian_sigma_shape_score_jacobian_V <- function(theta, control = list()) {
  -logistic_gaussian_sigma_shape_fisher_information(theta, control = control)
}

fit_logistic_gaussian_sigma_shape_theta <- function(data, weights = NULL, null,
                                                      control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  z <- normalized$ilr

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_logistic_gaussian_sigma_shape_theta(
      null$theta,
      ambient_dim = normalized$ambient_dim,
      control = control
    ))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  probability_weights <- if (is.null(weights)) {
    rep.int(1 / nrow(z), nrow(z))
  } else {
    normalize_probability_weights(weights, nrow(z))
  }

  mu_hat <- colSums(z * probability_weights)
  centered <- sweep(z, 2L, mu_hat, FUN = "-")
  shape_inverse <- logistic_gaussian_maphosa_shape_inverse(normalized$ilr_dim)
  quadratic <- rowSums((centered %*% shape_inverse) * centered)
  sigma2_hat <- sum(probability_weights * quadratic) / normalized$ilr_dim

  if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
    stop("Logistic-Gaussian sigma-shape MLE stopped: the fitted scale is not strictly positive.")
  }

  fit <- normalize_logistic_gaussian_sigma_shape_theta(
    list(mu_ilr = mu_hat, sigma = sqrt(sigma2_hat)),
    ambient_dim = normalized$ambient_dim,
    control = control
  )
  fit$loglik <- logistic_gaussian_sigma_shape_loglik(
    normalized,
    fit,
    probability_weights,
    control
  )
  fit$fit_diagnostics <- list(
    weighted = !is.null(weights),
    optimizer = "closed_form_mle",
    covariance_shape = "maphosa_M2_ilr"
  )
  fit
}

prepare_logistic_gaussian_sigma_shape_fast_multiplier <- function(spec, data, theta_hat,
                                                                   ks_prep = NULL,
                                                                   cvm_prep = NULL,
                                                                   control = list()) {
  normalized <- normalize_logistic_gaussian_data(data, control)
  theta_hat <- normalize_logistic_gaussian_sigma_shape_theta(
    theta_hat,
    ambient_dim = normalized$ambient_dim,
    control = control
  )
  d <- theta_hat$ilr_dim
  par0 <- c(theta_hat$mu_ilr, theta_hat$sigma)

  state_from_par <- function(par) {
    par <- as.numeric(par)
    if (length(par) != d + 1L) {
      stop("Logistic-Gaussian sigma-shape fast-multiplier parameters have incompatible dimension.")
    }
    normalize_logistic_gaussian_sigma_shape_theta(
      list(mu_ilr = par[seq_len(d)], sigma = par[[d + 1L]]),
      ambient_dim = normalized$ambient_dim,
      control = control
    )
  }

  score_matrix_fn <- function(sample, par) {
    logistic_gaussian_sigma_shape_score_matrix(
      sample,
      state_from_par(par),
      control = control
    )
  }

  sample_fn <- function(n_aux, par) {
    state <- state_from_par(par)
    rlogistic_gaussian_sigma_shape(
      n = n_aux,
      mu_ilr = state$mu_ilr,
      sigma = state$sigma,
      control = control
    )
  }

  fisher_hat <- logistic_gaussian_sigma_shape_fisher_information(theta_hat, control)
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

  prepared$vhat_method <- "logistic_gaussian_sigma_shape_analytic_fisher"
  prepared$paper_Vhat <- -fisher_hat
  prepared$paper_Vhat_inverse <- solve(-fisher_hat)
  prepared$paper_Vhat_method <- "analytic_sigma_shape_score_jacobian"
  prepared$paper_Vhat_diagnostics <- fast_multiplier_matrix_condition_diagnostics(-fisher_hat)
  prepared$score_parameterization <- "beta=(mu_ilr,sigma)"
  prepared$auxiliary_sampler <- "logistic_gaussian_sigma_shape"
  prepared$covariance_shape <- theta_hat$shape_ilr
  prepared
}

make_logistic_gaussian_sigma_shape_spec <- function() {
  # Reuse all Aitchison-distance and Logistic-Gaussian profile machinery.
  # Only the fitted parameter space, score and Fisher information are replaced.
  spec <- make_logistic_gaussian_spec(unknown_param = "both")
  base_profile_eval <- spec$profile_eval
  base_profile_matrix_eval <- spec$profile_matrix_eval
  base_sample_profile_matrix_eval <- spec$sample_profile_matrix_eval

  spec$name <- "logistic_gaussian_sigma_shape_aitchison"
  spec$fit_theta <- function(data, weights = NULL, null, control = list()) {
    fit_logistic_gaussian_sigma_shape_theta(data, weights, null, control)
  }
  spec$normalize_data <- normalize_logistic_gaussian_data

  spec$profile_eval <- function(omega, t, theta, control = list()) {
    fitted <- normalize_logistic_gaussian_sigma_shape_theta(theta, control = control)
    base_profile_eval(
      omega,
      t,
      list(mu_ilr = fitted$mu_ilr, Sigma_ilr = fitted$Sigma_ilr),
      control
    )
  }

  if (is.function(base_profile_matrix_eval)) {
    spec$profile_matrix_eval <- function(omega_grid, t_grid, theta, control = list()) {
      fitted <- normalize_logistic_gaussian_sigma_shape_theta(theta, control = control)
      base_profile_matrix_eval(
        omega_grid,
        t_grid,
        list(mu_ilr = fitted$mu_ilr, Sigma_ilr = fitted$Sigma_ilr),
        control
      )
    }
  }

  if (is.function(base_sample_profile_matrix_eval)) {
    spec$sample_profile_matrix_eval <- function(data, distance_matrix, theta,
                                                 control = list()) {
      fitted <- normalize_logistic_gaussian_sigma_shape_theta(theta, control = control)
      base_sample_profile_matrix_eval(
        data,
        distance_matrix,
        list(mu_ilr = fitted$mu_ilr, Sigma_ilr = fitted$Sigma_ilr),
        control
      )
    }
  }

  spec$score_matrix <- logistic_gaussian_sigma_shape_score_matrix
  spec$fisher_information <- logistic_gaussian_sigma_shape_fisher_information
  spec$score_jacobian_V <- logistic_gaussian_sigma_shape_score_jacobian_V
  spec$weighted_mle <- TRUE
  spec$fast_multiplier_prepare <- function(data, theta_hat,
                                            ks_prep = NULL,
                                            cvm_prep = NULL,
                                            control = list()) {
    prepare_logistic_gaussian_sigma_shape_fast_multiplier(
      spec = make_logistic_gaussian_sigma_shape_spec(),
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
