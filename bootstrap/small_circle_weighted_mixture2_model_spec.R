# Weighted two-small-circles model adapter for multiplier bootstrap GOF tests

resolve_small_circle_weighted_mixture2_model_spec_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

model_specs_path_small_circle_weighted_mixture2 <- resolve_small_circle_weighted_mixture2_model_spec_path("bootstrap", "model_specs.R")
if (!exists("new_model_spec", mode = "function")) {
  source(model_specs_path_small_circle_weighted_mixture2)
}

normalize_small_circle_weighted_mixture2_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Weighted small-circle-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_small_circle_weighted_mixture2_theta <- function(theta, ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Weighted small-circle-mixture theta must be a list containing `mu`, `pi`, `kappa1`, `nu1`, `kappa2`, `nu2`.")
  }

  params <- small_circle_weighted_mixture2_normalize_theta(theta, ambient_dim = ambient_dim)
  list(
    mu = params$mu,
    pi = params$pi,
    kappa1 = params$kappa1,
    nu1 = params$nu1,
    kappa2 = params$kappa2,
    nu2 = params$nu2,
    ambient_dim = length(params$mu)
  )
}

fit_small_circle_weighted_mixture2_theta <- function(data,
                                                     weights = NULL,
                                                     null,
                                                     control = list()) {
  x <- normalize_small_circle_weighted_mixture2_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_small_circle_weighted_mixture2_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- small_circle_weighted_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  theta <- normalize_small_circle_weighted_mixture2_theta(fit, ambient_dim = ncol(x))
  c(theta, fit[setdiff(names(fit), names(theta))])
}

prepare_small_circle_weighted_mixture2_fast_multiplier <- function(spec,
                                                                   data,
                                                                   theta_hat,
                                                                   ks_prep = NULL,
                                                                   cvm_prep = NULL,
                                                                   control = list()) {
  x <- normalize_small_circle_weighted_mixture2_data(data, control)
  theta_hat <- normalize_small_circle_weighted_mixture2_theta(theta_hat, ambient_dim = ncol(x))
  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  nu_upper <- 1 - as.numeric(control$small_circle_weighted_mixture2_nu_eps %||% 1e-6)
  par0 <- c(
    0,
    0,
    stats::qlogis(theta_hat$pi),
    log(pmax(expm1(theta_hat$kappa1), .Machine$double.eps)),
    small_circle_inverse_logistic_bounded(theta_hat$nu1, upper = nu_upper),
    log(pmax(expm1(theta_hat$kappa2), .Machine$double.eps)),
    small_circle_inverse_logistic_bounded(theta_hat$nu2, upper = nu_upper)
  )
  kappa_min <- as.numeric(control$small_circle_weighted_mixture2_kappa_min %||% 1e-8)
  kappa_max <- as.numeric(control$small_circle_weighted_mixture2_kappa_max %||% 1e6)
  nu_eps <- as.numeric(control$small_circle_weighted_mixture2_nu_eps %||% 1e-6)
  weight_eps <- as.numeric(control$small_circle_weighted_mixture2_weight_eps %||% 1e-6)

  softplus <- function(x) log1p(exp(x))
  d_softplus <- function(x) stats::plogis(x)
  d_bounded_logistic <- function(raw) {
    prob <- stats::plogis(raw)
    (1 - nu_eps) * prob * (1 - prob)
  }

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_small_circle_weighted_mixture2_theta(
      list(
        mu = mapped$mu,
        pi = rotational_bounded_weight(par[[3L]], weight_eps = weight_eps),
        kappa1 = min(max(softplus(par[[4L]]), kappa_min), kappa_max),
        nu1 = small_circle_logistic_bounded(par[[5L]], upper = 1 - nu_eps),
        kappa2 = min(max(softplus(par[[6L]]), kappa_min), kappa_max),
        nu2 = small_circle_logistic_bounded(par[[7L]], upper = 1 - nu_eps)
      ),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_small_circle_weighted_mixture2_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    s1 <- z
    s2 <- -z
    log_f1 <- fast_multiplier_small_circle_component_log_density(
      s = s1,
      kappa = theta_state$kappa1,
      nu = theta_state$nu1
    )
    log_f2 <- fast_multiplier_small_circle_component_log_density(
      s = s2,
      kappa = theta_state$kappa2,
      nu = theta_state$nu2
    )
    log_mix <- rotational_logsumexp2(
      log(theta_state$pi) + log_f1,
      log1p(-theta_state$pi) + log_f2
    )
    r1 <- exp(log(theta_state$pi) + log_f1 - log_mix)
    r2 <- 1 - r1

    coeff_mu <- r1 * (-2 * theta_state$kappa1 * (s1 - theta_state$nu1)) +
      r2 * (2 * theta_state$kappa2 * (s2 - theta_state$nu2))
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))
    score_comp1 <- fast_multiplier_small_circle_component_scores_natural(
      s = s1,
      kappa = theta_state$kappa1,
      nu = theta_state$nu1
    )
    score_comp2 <- fast_multiplier_small_circle_component_scores_natural(
      s = s2,
      kappa = theta_state$kappa2,
      nu = theta_state$nu2
    )

    cbind(
      score_mu,
      r1 - theta_state$pi,
      r1 * score_comp1[, 1L] * d_softplus(par[[4L]]),
      r1 * score_comp1[, 2L] * d_bounded_logistic(par[[5L]]),
      r2 * score_comp2[, 1L] * d_softplus(par[[6L]]),
      r2 * score_comp2[, 2L] * d_bounded_logistic(par[[7L]])
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_small_circle_weighted_mixture2(
      n = n_aux,
      mu = theta_state$mu,
      pi = theta_state$pi,
      kappa1 = theta_state$kappa1,
      nu1 = theta_state$nu1,
      kappa2 = theta_state$kappa2,
      nu2 = theta_state$nu2
    )
  }

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = x,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = score_matrix_fn,
    sample_fn = sample_fn
  )
}

make_small_circle_weighted_mixture2_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("small_circle_weighted_mixture2_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_small_circle_weighted_mixture2_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_small_circle_weighted_mixture2_data(data, control)
      omega_matrix <- normalize_small_circle_weighted_mixture2_data(omega, control)

      if (ncol(x) != ncol(omega_matrix)) {
        stop("`data` and `omega` have incompatible ambient dimensions.")
      }

      dot_products <- pmin(pmax(x %*% t(omega_matrix), -1), 1)
      if (identical(distance_type, "chordal")) {
        sqrt(2 * (1 - dot_products))
      } else {
        acos(dot_products)
      }
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_small_circle_weighted_mixture2_theta(theta)
      distance_profile_small_circle_weighted_mixture2(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        pi = theta$pi,
        kappa1 = theta$kappa1,
        nu1 = theta$nu1,
        kappa2 = theta$kappa2,
        nu2 = theta$nu2,
        distance_type = distance_type,
        method = control$small_circle_weighted_mixture2_profile_method %||% "legendre",
        l_max = as.integer(control$small_circle_weighted_mixture2_L_max %||% 200L),
        quad_n = as.integer(control$small_circle_weighted_mixture2_quad_n %||% 400L),
        tol = as.numeric(control$small_circle_weighted_mixture2_tol %||% 1e-10),
        validate_against_integral = isTRUE(control$small_circle_weighted_mixture2_validate_against_integral),
        validation_tol = as.numeric(control$small_circle_weighted_mixture2_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_small_circle_weighted_mixture2_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_small_circle_weighted_mixture2_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_small_circle_weighted_mixture2_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_small_circle_weighted_mixture2_theta(theta)
        distance_profile_small_circle_weighted_mixture2_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          pi = theta$pi,
          kappa1 = theta$kappa1,
          nu1 = theta$nu1,
          kappa2 = theta$kappa2,
          nu2 = theta$nu2,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$small_circle_weighted_mixture2_profile_method %||% "legendre",
          l_max = as.integer(control$small_circle_weighted_mixture2_L_max %||% 200L),
          quad_n = as.integer(control$small_circle_weighted_mixture2_quad_n %||% 400L),
          tol = as.numeric(control$small_circle_weighted_mixture2_tol %||% 1e-10)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        prepare_small_circle_weighted_mixture2_fast_multiplier(
          spec = make_small_circle_weighted_mixture2_spec(distance_type = distance_type),
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
