# Normal model adapter for multiplier bootstrap GOF tests

if (!exists("new_model_spec", mode = "function")) {
  model_specs_candidates_normal <- c(
    file.path("bootstrap", "model_specs.R"), "model_specs.R",
    file.path("..", "bootstrap", "model_specs.R"),
    file.path("..", "..", "bootstrap", "model_specs.R")
  )
  model_specs_path_normal <- model_specs_candidates_normal[file.exists(model_specs_candidates_normal)][1L]
  if (is.na(model_specs_path_normal)) stop("Could not locate bootstrap/model_specs.R.")
  source(model_specs_path_normal)
}

normalize_normal_data <- function(data, control = list()) {
  if (is.matrix(data) || is.data.frame(data)) {
    if (ncol(as.matrix(data)) != 1L) {
      stop("Normal data must be a numeric vector or a one-column matrix/data frame.")
    }
    data <- as.matrix(data)[, 1]
  }

  data <- as.numeric(data)
  if (length(data) == 0) {
    stop("`data` cannot be empty.")
  }
  if (any(!is.finite(data))) {
    stop("Normal data must be finite.")
  }

  data
}

normalize_normal_theta <- function(theta) {
  if (is.numeric(theta) && length(theta) == 2L && is.null(names(theta))) {
    theta <- list(mu = theta[[1]], sigma = theta[[2]])
  }

  if (!is.list(theta)) {
    stop("Normal theta must be a list with entries `mu` and `sigma`.")
  }

  mu <- as.numeric(theta$mu)
  sigma <- as.numeric(theta$sigma)

  if (length(mu) != 1L || !is.finite(mu)) {
    stop("Normal theta requires a finite scalar `mu`.")
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("Normal theta requires a strictly positive finite scalar `sigma`.")
  }

  list(mu = mu, sigma = sigma)
}

fit_normal_theta <- function(data,
                             weights = NULL,
                             null,
                             unknown_param = NULL,
                             control = list()) {
  x <- normalize_normal_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_normal_theta(null$theta))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (is.null(unknown_param)) {
    stop("Normal composite null requires `unknown_param` in the spec.")
  }

  fixed <- null$fixed %||% list()
  prob_weights <- if (is.null(weights)) {
    rep.int(1 / length(x), length(x))
  } else {
    normalize_probability_weights(weights, length(x))
  }

  weighted_mean <- function(values) {
    sum(prob_weights * values)
  }

  weighted_sigma <- function(values, center) {
    sqrt(sum(prob_weights * (values - center)^2))
  }

  if (identical(unknown_param, "mu")) {
    sigma_fixed <- as.numeric(fixed$sigma)
    if (length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0) {
      stop("Composite normal null with unknown `mu` requires `null$fixed$sigma`.")
    }
    return(list(mu = weighted_mean(x), sigma = sigma_fixed))
  }

  if (identical(unknown_param, "sigma")) {
    mu_fixed <- as.numeric(fixed$mu)
    if (length(mu_fixed) != 1L || !is.finite(mu_fixed)) {
      stop("Composite normal null with unknown `sigma` requires `null$fixed$mu`.")
    }
    return(list(mu = mu_fixed, sigma = weighted_sigma(x, mu_fixed)))
  }

  if (identical(unknown_param, "both")) {
    mu_hat <- weighted_mean(x)
    return(list(mu = mu_hat, sigma = weighted_sigma(x, mu_hat)))
  }

  stop("Unsupported `unknown_param` for the normal model.")
}

prepare_normal_fast_multiplier <- function(spec,
                                           data,
                                           theta_hat,
                                           ks_prep = NULL,
                                           cvm_prep = NULL,
                                           control = list(),
                                           unknown_param = "both") {
  x <- normalize_normal_data(data, control)
  theta_hat <- normalize_normal_theta(theta_hat)
  unknown_param <- tolower(as.character(unknown_param %||% "both"))

  if (!unknown_param %in% c("mu", "sigma", "both")) {
    stop("Unsupported `unknown_param` for the fast normal multiplier bootstrap.")
  }

  par0 <- switch(
    unknown_param,
    mu = c(theta_hat$mu),
    sigma = c(log(theta_hat$sigma)),
    both = c(theta_hat$mu, log(theta_hat$sigma))
  )

  state_from_par <- function(par) {
    if (identical(unknown_param, "mu")) {
      normalize_normal_theta(list(mu = par[[1L]], sigma = theta_hat$sigma))
    } else if (identical(unknown_param, "sigma")) {
      normalize_normal_theta(list(mu = theta_hat$mu, sigma = exp(par[[1L]])))
    } else {
      normalize_normal_theta(list(mu = par[[1L]], sigma = exp(par[[2L]])))
    }
  }

  score_matrix_fn <- function(sample, par) {
    theta <- state_from_par(par)
    sample <- normalize_normal_data(sample)
    centered <- sample - theta$mu
    sigma_sq <- theta$sigma^2

    if (identical(unknown_param, "mu")) {
      return(matrix(centered / sigma_sq, ncol = 1L))
    }
    if (identical(unknown_param, "sigma")) {
      return(matrix(-1 + (centered^2) / sigma_sq, ncol = 1L))
    }

    cbind(
      centered / sigma_sq,
      -1 + (centered^2) / sigma_sq
    )
  }

  sample_fn <- function(n_aux, par) {
    theta <- state_from_par(par)
    stats::rnorm(n_aux, mean = theta$mu, sd = theta$sigma)
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

make_normal_spec <- function(unknown_param = NULL) {
  new_model_spec(
    name = "normal",
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_normal_theta(
        data = data,
        weights = weights,
        null = null,
        unknown_param = unknown_param,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_normal_data(data, control)
      omega_vec <- as.numeric(omega)
      if (length(omega_vec) == 0) {
        stop("`omega` cannot be empty.")
      }
      abs(outer(x, omega_vec, FUN = "-"))
    },
    profile_eval = function(omega, t, theta, control = list()) {
      theta <- normalize_normal_theta(theta)
      theoretical_distance_profile_normal(
        omega = as.numeric(omega),
        mu = theta$mu,
        sigma = theta$sigma,
        t_values = as.numeric(t)
      )
    },
    normalize_data = normalize_normal_data,
    n_obs = function(data, control = list()) {
      length(normalize_normal_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_normal_data(data, control)[[idx]]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        if (!identical(distance_profile_backend_current(), "cpp")) return(NULL)
        theta <- normalize_normal_theta(theta)
        omega_vec <- as.numeric(omega_grid)
        t_matrix <- matrix(
          rep(as.numeric(t_grid), each = length(omega_vec)),
          nrow = length(omega_vec),
          ncol = length(t_grid)
        )
        distance_profile_cpp_call(
          "cpp_dp_normal_profile_matrix",
          omega_vec,
          t_matrix,
          theta$mu,
          theta$sigma
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_profile_backend_current(), "cpp")) return(NULL)
        theta <- normalize_normal_theta(theta)
        x <- normalize_normal_data(data, control)
        distance_profile_cpp_call(
          "cpp_dp_normal_profile_matrix",
          x,
          as.matrix(distance_matrix),
          theta$mu,
          theta$sigma
        )
      },
      unknown_param = unknown_param,
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        if (is.null(unknown_param)) {
          stop("The fast normal multiplier bootstrap requires `unknown_param` in the model spec.")
        }
        prepare_normal_fast_multiplier(
          spec = make_normal_spec(unknown_param = unknown_param),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control,
          unknown_param = unknown_param
        )
      }
    )
  )
}
