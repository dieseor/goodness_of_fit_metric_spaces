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
      weighted_mle = TRUE
    )
  )
}
