# Two-component beta-mixture model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_beta_mixture2 <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_beta_mixture2 <- model_specs_candidates_beta_mixture2[file.exists(model_specs_candidates_beta_mixture2)][1L]
  if (is.na(model_specs_path_beta_mixture2)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_beta_mixture2)
}

prepare_beta_mixture2_fast_multiplier <- function(spec,
                                                  data,
                                                  theta_hat,
                                                  ks_prep = NULL,
                                                  cvm_prep = NULL,
                                                  control = list()) {
  x <- normalize_beta_mixture2_data(data, control)
  theta_hat <- normalize_beta_mixture2_theta(theta_hat, ambient_dim = ncol(x))
  fast_shape_regular_eps <- as.numeric(control$beta_mixture2_fast_shape_regular_eps %||% 0)
  shape_values <- c(theta_hat$alpha1, theta_hat$beta1, theta_hat$alpha2, theta_hat$beta2)
  if (any(!is.finite(shape_values))) {
    stop("The fitted beta_mixture2 shape parameters are not finite.")
  }
  if (min(shape_values) <= 1 + fast_shape_regular_eps) {
    return(list(
      fallback_to_reestimated = TRUE,
      fallback_reason = "beta_mixture2_shape_nonregular"
    ))
  }
  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(
    0,
    0,
    stats::qlogis(theta_hat$weight1),
    log(theta_hat$alpha1),
    log(theta_hat$beta1),
    log(theta_hat$alpha2),
    log(theta_hat$beta2)
  )
  weight_eps <- as.numeric(control$beta_mixture2_weight_eps %||% 0.01)
  shape_lower <- as.numeric(control$beta_mixture2_shape_lower %||% 0.05)
  shape_upper <- as.numeric(control$beta_mixture2_shape_upper %||% 1e3)
  y_eps <- as.numeric(control$beta_mixture2_eps %||% 1e-12)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_beta_mixture2_theta(
      list(
        mu = mapped$mu,
        weight1 = rotational_bounded_weight(par[[3L]], weight_eps = weight_eps),
        alpha1 = rotational_positive_parameter(par[[4L]], lower = shape_lower, upper = shape_upper),
        beta1 = rotational_positive_parameter(par[[5L]], lower = shape_lower, upper = shape_upper),
        alpha2 = rotational_positive_parameter(par[[6L]], lower = shape_lower, upper = shape_upper),
        beta2 = rotational_positive_parameter(par[[7L]], lower = shape_lower, upper = shape_upper)
      ),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_beta_mixture2_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    y <- rotational_clamp_unit_interval((z + 1) / 2, eps = y_eps)
    log_b1 <- stats::dbeta(y, theta_state$alpha1, theta_state$beta1, log = TRUE)
    log_b2 <- stats::dbeta(y, theta_state$alpha2, theta_state$beta2, log = TRUE)
    log_m <- rotational_logsumexp2(
      log(theta_state$weight1) + log_b1,
      log1p(-theta_state$weight1) + log_b2
    )
    r1 <- exp(log(theta_state$weight1) + log_b1 - log_m)
    r2 <- 1 - r1

    g1 <- (theta_state$alpha1 - 1) / y - (theta_state$beta1 - 1) / (1 - y)
    g2 <- (theta_state$alpha2 - 1) / y - (theta_state$beta2 - 1) / (1 - y)
    coeff_mu <- 0.5 * (r1 * g1 + r2 * g2)
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))

    cbind(
      score_mu,
      r1 - theta_state$weight1,
      theta_state$alpha1 * r1 * (
        log(y) - digamma(theta_state$alpha1) + digamma(theta_state$alpha1 + theta_state$beta1)
      ),
      theta_state$beta1 * r1 * (
        log1p(-y) - digamma(theta_state$beta1) + digamma(theta_state$alpha1 + theta_state$beta1)
      ),
      theta_state$alpha2 * r2 * (
        log(y) - digamma(theta_state$alpha2) + digamma(theta_state$alpha2 + theta_state$beta2)
      ),
      theta_state$beta2 * r2 * (
        log1p(-y) - digamma(theta_state$beta2) + digamma(theta_state$alpha2 + theta_state$beta2)
      )
    )
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_beta_mixture2(
      n = n_aux,
      mu = theta_state$mu,
      weight1 = theta_state$weight1,
      alpha1 = theta_state$alpha1,
      beta1 = theta_state$beta1,
      alpha2 = theta_state$alpha2,
      beta2 = theta_state$beta2
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
normalize_beta_mixture2_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Rotational beta-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_beta_mixture2_theta <- function(theta, ambient_dim = 3L) {
  beta_mixture2_normalize_theta(theta, ambient_dim = ambient_dim)
}

fit_beta_mixture2_theta <- function(data,
                                               weights = NULL,
                                               null,
                                               control = list()) {
  x <- normalize_beta_mixture2_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_beta_mixture2_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- beta_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  normalize_beta_mixture2_theta(fit, ambient_dim = ncol(x))
}

make_beta_mixture2_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("beta_mixture2_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_beta_mixture2_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_beta_mixture2_data(data, control)
      omega_matrix <- normalize_beta_mixture2_data(omega, control)

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
      theta <- normalize_beta_mixture2_theta(theta)
      distance_profile_beta_mixture2(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        weight1 = theta$weight1,
        alpha1 = theta$alpha1,
        beta1 = theta$beta1,
        alpha2 = theta$alpha2,
        beta2 = theta$beta2,
        distance_type = distance_type,
        method = control$beta_mixture2_profile_method %||% "legendre",
        l_max = as.integer(control$beta_mixture2_L_max %||% 150L),
        quad_n = as.integer(control$beta_mixture2_quad_n %||% 1000L),
        tol = as.numeric(control$beta_mixture2_tol %||% 1e-6),
        validate_against_integral = isTRUE(control$beta_mixture2_validate_against_integral),
        validation_tol = as.numeric(control$beta_mixture2_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_beta_mixture2_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_beta_mixture2_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_beta_mixture2_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_beta_mixture2_theta(theta)
        distance_profile_beta_mixture2_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$beta_mixture2_profile_method %||% "legendre",
          l_max = as.integer(control$beta_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$beta_mixture2_quad_n %||% 1000L),
          tol = as.numeric(control$beta_mixture2_tol %||% 1e-6)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_beta_mixture2_theta(theta)
        distance_profile_beta_mixture2_cvm_grid(
          X = data,
          mu = theta$mu,
          weight1 = theta$weight1,
          alpha1 = theta$alpha1,
          beta1 = theta$beta1,
          alpha2 = theta$alpha2,
          beta2 = theta$beta2,
          method = control$beta_mixture2_profile_method %||% "legendre",
          l_max = as.integer(control$beta_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$beta_mixture2_quad_n %||% 1000L),
          tol = as.numeric(control$beta_mixture2_tol %||% 1e-6)
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
        prepare_beta_mixture2_fast_multiplier(
          spec = make_beta_mixture2_spec(distance_type = distance_type),
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
