# Two-component logit-normal-mixture model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_logitnormal_mixture2 <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_logitnormal_mixture2 <- model_specs_candidates_logitnormal_mixture2[file.exists(model_specs_candidates_logitnormal_mixture2)][1L]
  if (is.na(model_specs_path_logitnormal_mixture2)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_logitnormal_mixture2)
}

normalize_logitnormal_mixture2_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Rotational logit-normal-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_logitnormal_mixture2_theta <- function(theta, ambient_dim = 3L) {
  logitnormal_mixture2_normalize_theta(theta, ambient_dim = ambient_dim)
}

fit_logitnormal_mixture2_theta <- function(data,
                                                      weights = NULL,
                                                      null,
                                                      control = list()) {
  x <- normalize_logitnormal_mixture2_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_logitnormal_mixture2_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- logitnormal_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )

  normalize_logitnormal_mixture2_theta(fit, ambient_dim = ncol(x))
}

make_logitnormal_mixture2_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("logitnormal_mixture2_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_logitnormal_mixture2_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_logitnormal_mixture2_data(data, control)
      omega_matrix <- normalize_logitnormal_mixture2_data(omega, control)

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
      theta <- normalize_logitnormal_mixture2_theta(theta)
      distance_profile_logitnormal_mixture2(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        weight1 = theta$weight1,
        mean1 = theta$mean1,
        sd1 = theta$sd1,
        mean2 = theta$mean2,
        sd2 = theta$sd2,
        distance_type = distance_type,
        method = control$logitnormal_mixture2_profile_method %||% "integral",
        l_max = as.integer(control$logitnormal_mixture2_L_max %||% 150L),
        quad_n = as.integer(control$logitnormal_mixture2_quad_n %||% 1000L),
        tol = as.numeric(control$logitnormal_mixture2_tol %||% 1e-6),
        eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12),
        validate_against_integral = isTRUE(control$logitnormal_mixture2_validate_against_integral),
        validation_tol = as.numeric(control$logitnormal_mixture2_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_logitnormal_mixture2_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_logitnormal_mixture2_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_logitnormal_mixture2_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_logitnormal_mixture2_theta(theta)
        distance_profile_logitnormal_mixture2_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$logitnormal_mixture2_profile_method %||% "integral",
          l_max = as.integer(control$logitnormal_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$logitnormal_mixture2_quad_n %||% 400L),
          tol = as.numeric(control$logitnormal_mixture2_tol %||% 1e-6),
          eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_logitnormal_mixture2_theta(theta)
        distance_profile_logitnormal_mixture2_cvm_grid(
          X = data,
          mu = theta$mu,
          weight1 = theta$weight1,
          mean1 = theta$mean1,
          sd1 = theta$sd1,
          mean2 = theta$mean2,
          sd2 = theta$sd2,
          method = control$logitnormal_mixture2_profile_method %||% "integral",
          l_max = as.integer(control$logitnormal_mixture2_L_max %||% 150L),
          quad_n = as.integer(control$logitnormal_mixture2_quad_n %||% 400L),
          tol = as.numeric(control$logitnormal_mixture2_tol %||% 1e-6),
          eps = as.numeric(control$logitnormal_mixture2_eps %||% 1e-12)
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE
    )
  )
}
