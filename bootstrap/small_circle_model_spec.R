# Small Circle (Bingham--Mardia) model adapter for multiplier bootstrap GOF tests

resolve_small_circle_model_spec_path <- function(...) {
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

model_specs_path_small_circle <- resolve_small_circle_model_spec_path("bootstrap", "model_specs.R")
if (!exists("new_model_spec", mode = "function")) {
  source(model_specs_path_small_circle)
}

normalize_small_circle_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Small Circle data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_small_circle_theta <- function(theta, ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Small Circle theta must be a list containing `mu`, `kappa` and `nu`.")
  }

  params <- small_circle_validate_parameters(
    mu = theta$mu,
    kappa = theta$kappa,
    nu = theta$nu,
    allow_negative_nu = TRUE
  )
  canonical <- small_circle_canonicalize_theta(params$mu, params$nu)
  params$mu <- canonical$mu
  params$nu <- canonical$nu

  if (length(params$mu) != ambient_dim) {
    stop("Small Circle theta has incompatible ambient dimension.")
  }

  list(
    mu = params$mu,
    kappa = params$kappa,
    nu = params$nu,
    ambient_dim = length(params$mu)
  )
}

fit_small_circle_theta <- function(data,
                                   weights = NULL,
                                   null,
                                   control = list()) {
  x <- normalize_small_circle_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_small_circle_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- small_circle_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  normalize_small_circle_theta(fit, ambient_dim = ncol(x))
}

prepare_small_circle_fast_multiplier <- function(spec,
                                                 data,
                                                 theta_hat,
                                                 ks_prep = NULL,
                                                 cvm_prep = NULL,
                                                 control = list()) {
  x <- normalize_small_circle_data(data, control)
  theta_hat <- normalize_small_circle_theta(theta_hat, ambient_dim = ncol(x))
  chart <- fast_multiplier_sphere_chart(theta_hat$mu)
  par0 <- c(0, 0, theta_hat$kappa, theta_hat$nu)

  state_from_par <- function(par) {
    mapped <- fast_multiplier_sphere_chart_map(chart, par[1:2])
    normalize_small_circle_theta(
      list(mu = mapped$mu, kappa = par[[3L]], nu = par[[4L]]),
      ambient_dim = ncol(x)
    )
  }

  score_matrix_fn <- function(sample, par) {
    theta_state <- state_from_par(par)
    sample <- normalize_small_circle_data(sample, control)
    jac_mu <- fast_multiplier_sphere_chart_jacobian(chart, par[1:2])
    z <- pmin(pmax(as.numeric(sample %*% theta_state$mu), -1), 1)
    coeff_mu <- -2 * theta_state$kappa * (z - theta_state$nu)
    score_mu <- t(vapply(seq_len(nrow(sample)), function(i) {
      drop(t(jac_mu) %*% (coeff_mu[[i]] * sample[i, ]))
    }, numeric(2L)))
    score_scalar <- fast_multiplier_small_circle_component_scores_natural(
      s = z,
      kappa = theta_state$kappa,
      nu = theta_state$nu
    )
    cbind(score_mu, score_scalar)
  }

  sample_fn <- function(n_aux, par) {
    theta_state <- state_from_par(par)
    r_sph_small_circle(
      n = n_aux,
      mu = theta_state$mu,
      kappa = theta_state$kappa,
      nu = theta_state$nu
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

make_small_circle_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("small_circle_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_small_circle_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_small_circle_data(data, control)
      omega_matrix <- normalize_small_circle_data(omega, control)

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
      theta <- normalize_small_circle_theta(theta)
      distance_profile_small_circle(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        kappa = theta$kappa,
        nu = theta$nu,
        distance_type = distance_type,
        method = control$small_circle_profile_method %||% "legendre",
        l_max = as.integer(control$small_circle_L_max %||% 200L),
        quad_n = as.integer(control$small_circle_quad_n %||% 400L),
        tol = as.numeric(control$small_circle_tol %||% 1e-10),
        validate_against_integral = isTRUE(control$small_circle_validate_against_integral),
        validation_tol = as.numeric(control$small_circle_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_small_circle_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_small_circle_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_small_circle_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_small_circle_theta(theta)
        distance_profile_small_circle_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          kappa = theta$kappa,
          nu = theta$nu,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$small_circle_profile_method %||% "legendre",
          l_max = as.integer(control$small_circle_L_max %||% 200L),
          quad_n = as.integer(control$small_circle_quad_n %||% 400L),
          tol = as.numeric(control$small_circle_tol %||% 1e-10)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_small_circle_theta(theta)
        distance_profile_small_circle_cvm_grid(
          X = data,
          mu = theta$mu,
          kappa = theta$kappa,
          nu = theta$nu,
          method = control$small_circle_profile_method %||% "legendre",
          l_max = as.integer(control$small_circle_L_max %||% 200L),
          quad_n = as.integer(control$small_circle_quad_n %||% 400L),
          tol = as.numeric(control$small_circle_tol %||% 1e-10)
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
        prepare_small_circle_fast_multiplier(
          spec = make_small_circle_spec(distance_type = distance_type),
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
