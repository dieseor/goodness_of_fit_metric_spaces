# Parsimonious temporal-model adapter for the joint sunspot GOF tests.

resolve_sunspots_joint_parsimonious_spec_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

if (!exists("make_sunspots_joint_time_space_spec", mode = "function")) {
  source(resolve_sunspots_joint_parsimonious_spec_path(
    "bootstrap", "sunspots_joint_time_space_model_spec.R"
  ))
}
if (!exists(
  "fit_sunspots_cycle23_joint_time_space_parsimonious",
  mode = "function"
)) {
  source(resolve_sunspots_joint_parsimonious_spec_path(
    "real_data",
    "sunspots",
    "sunspots_cycle23_joint_time_models_parsimonious.R"
  ))
}

require_sunspots_joint_parsimonious_spec_dependencies <- function() {
  required <- c(
    "normalize_sunspots_joint_spec_data",
    "sunspots_joint_spec_distance_matrix",
    "sunspots_joint_profile_block",
    "sunspots_joint_profile_block_sorted",
    "sunspots_joint_conditional_legendre_coefficients",
    "fit_sunspots_cycle23_joint_time_space_parsimonious",
    "sunspots_joint_parsimonious_time_quadrature",
    "sunspots_joint_parsimonious_pack_par",
    "sunspots_joint_parsimonious_score_matrix",
    "sample_sunspots_joint_time_space_parsimonious"
  )
  missing <- required[
    !vapply(required, exists, logical(1L), mode = "function")
  ]
  if (length(missing) > 0L) {
    stop(sprintf(
      "Missing parsimonious joint-model dependencies: %s.",
      paste(missing, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

fit_sunspots_joint_parsimonious_spec_theta <- function(
    data,
    weights = NULL,
    null,
    control = list()) {
  require_sunspots_joint_parsimonious_spec_dependencies()
  normalized <- normalize_sunspots_joint_spec_data(data, control)
  time_model <- normalize_sunspots_joint_parsimonious_time_model(
    control$time_model %||% "beta"
  )
  hemisphere_regression <-
    sunspots_time_varying_normalize_hemisphere_regression(
      control$hemisphere_regression %||% "shared"
    )

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least `type`.")
  }
  if (identical(null$type, "simple")) {
    theta_simple <- null$theta
    if (!is.list(theta_simple) ||
        is.null(theta_simple$eta_hat) ||
        is.null(theta_simple$theta_hat)) {
      stop(
        "For a simple null, `null$theta` must contain `eta_hat` and `theta_hat`."
      )
    }
    return(fit_sunspots_cycle23_joint_time_space_parsimonious(
      x = normalized$x,
      s = normalized$s,
      time_model = time_model,
      hemisphere_regression = hemisphere_regression,
      control = control,
      eta_hat = theta_simple$eta_hat,
      theta_hat = theta_simple$theta_hat
    ))
  }
  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }
  if (!is.null(weights)) {
    stop(
      paste(
        "The parsimonious joint sunspot specs only support",
        "`bootstrap_method = 'fast_multiplier'`;",
        "weighted refits are not implemented."
      )
    )
  }

  fit_sunspots_cycle23_joint_time_space_parsimonious(
    x = normalized$x,
    s = normalized$s,
    time_model = time_model,
    hemisphere_regression = hemisphere_regression,
    control = control
  )
}

sunspots_joint_parsimonious_spec_profile_eval <- function(
    omega,
    t,
    theta,
    control = list()) {
  require_sunspots_joint_parsimonious_spec_dependencies()
  omega <- as.numeric(omega)
  if (length(omega) != 4L || any(!is.finite(omega))) {
    stop("`omega` must be a finite vector (x1, x2, x3, s).")
  }
  center <- sunspots_joint_validate_data(
    x = matrix(omega[1:3], nrow = 1L),
    s = omega[[4L]]
  )
  fit <- theta
  if (!is.list(fit) ||
      is.null(fit$eta_hat) ||
      is.null(fit$theta_hat) ||
      is.null(fit$time_model)) {
    stop(
      "`theta` must be a fitted parsimonious joint object."
    )
  }

  quadrature <- sunspots_joint_parsimonious_time_quadrature(
    fit$eta_hat,
    n_nodes = as.integer(control$time_quad_n %||% 64L),
    time_model = fit$time_model,
    control = control
  )
  coefficients <- sunspots_joint_conditional_legendre_coefficients(
    fit$theta_hat,
    quadrature$nodes,
    l_max = as.integer(control$profile_l_max %||% 100L),
    quad_n = as.integer(control$profile_quad_n %||% 400L)
  )

  as.numeric(sunspots_joint_profile_block(
    radii = matrix(as.numeric(t), nrow = 1L),
    rho = center$x[1L, 3L],
    center_s = center$s,
    time_nodes = quadrature$nodes,
    time_weights = quadrature$weights,
    coefficients = coefficients,
    backend = control$distance_profile_backend %||% "auto"
  ))
}

sunspots_joint_parsimonious_spec_sorted_prepare <- function(
    data,
    sorted_distance_matrix,
    theta,
    control = list()) {
  require_sunspots_joint_parsimonious_spec_dependencies()
  normalized <- normalize_sunspots_joint_spec_data(data, control)
  fit <- theta
  if (!is.list(fit) ||
      is.null(fit$eta_hat) ||
      is.null(fit$theta_hat) ||
      is.null(fit$time_model)) {
    stop("`theta` must be a fitted parsimonious joint object.")
  }

  sorted_distance_matrix <- as.matrix(sorted_distance_matrix)
  n <- nrow(normalized$x)
  if (!identical(dim(sorted_distance_matrix), c(n, n))) {
    stop(
      "`sorted_distance_matrix` has incompatible dimensions."
    )
  }

  quadrature <- sunspots_joint_parsimonious_time_quadrature(
    fit$eta_hat,
    n_nodes = as.integer(control$time_quad_n %||% 64L),
    time_model = fit$time_model,
    control = control
  )
  coefficients <- sunspots_joint_conditional_legendre_coefficients(
    fit$theta_hat,
    quadrature$nodes,
    l_max = as.integer(control$profile_l_max %||% 100L),
    quad_n = as.integer(control$profile_quad_n %||% 400L)
  )
  backend <- sunspots_joint_effective_backend(
    control$distance_profile_backend %||% "r"
  )
  if (identical(backend, "cpp")) ensure_distance_profile_cpp_loaded()

  list(
    n = n,
    rho = as.numeric(normalized$x[, 3L]),
    center_s = as.numeric(normalized$s),
    time_nodes = quadrature$nodes,
    time_weights = quadrature$weights,
    coefficients = coefficients,
    backend = backend
  )
}

sunspots_joint_parsimonious_spec_sorted_block_eval <- function(
    data,
    sorted_distance_matrix,
    theta,
    row_indices,
    prepared = NULL,
    control = list()) {
  sorted_distance_matrix <- as.matrix(sorted_distance_matrix)
  n <- nrow(sorted_distance_matrix)
  if (ncol(sorted_distance_matrix) != n) {
    stop("`sorted_distance_matrix` must be square.")
  }
  row_indices <- as.integer(row_indices)
  if (length(row_indices) == 0L ||
      any(!is.finite(row_indices)) ||
      any(row_indices < 1L | row_indices > n)) {
    stop("`row_indices` contains invalid sample indices.")
  }
  if (is.null(prepared)) {
    prepared <- sunspots_joint_parsimonious_spec_sorted_prepare(
      data = data,
      sorted_distance_matrix = sorted_distance_matrix,
      theta = theta,
      control = control
    )
  }
  if (!is.list(prepared) ||
      !identical(as.integer(prepared$n), as.integer(n)) ||
      length(prepared$rho) != n ||
      length(prepared$center_s) != n) {
    stop("The prepared parsimonious profile context is incompatible.")
  }

  values <- sunspots_joint_profile_block_sorted(
    radii = sorted_distance_matrix[row_indices, , drop = FALSE],
    rho = prepared$rho[row_indices],
    center_s = prepared$center_s[row_indices],
    time_nodes = prepared$time_nodes,
    time_weights = prepared$time_weights,
    coefficients = prepared$coefficients,
    backend = prepared$backend
  )
  matrix(
    as.numeric(values),
    nrow = length(row_indices),
    ncol = n
  )
}

sunspots_joint_parsimonious_spec_profile_matrix_eval <- function(
    data,
    distance_matrix,
    theta,
    control = list()) {
  prepared <- sunspots_joint_parsimonious_spec_sorted_prepare(
    data = data,
    sorted_distance_matrix = distance_matrix,
    theta = theta,
    control = control
  )
  block_size <- as.integer(control$center_block_size %||% 8L)
  if (length(block_size) != 1L ||
      !is.finite(block_size) ||
      block_size <= 0L) {
    stop("`control$center_block_size` must be positive.")
  }

  n <- prepared$n
  out <- matrix(0, nrow = n, ncol = n)
  for (block_start in seq.int(1L, n, by = block_size)) {
    block_end <- min(block_start + block_size - 1L, n)
    row_indices <- block_start:block_end
    out[row_indices, ] <- sunspots_joint_profile_block(
      radii = distance_matrix[row_indices, , drop = FALSE],
      rho = prepared$rho[row_indices],
      center_s = prepared$center_s[row_indices],
      time_nodes = prepared$time_nodes,
      time_weights = prepared$time_weights,
      coefficients = prepared$coefficients,
      backend = prepared$backend
    )
  }
  out
}

prepare_sunspots_joint_parsimonious_fast_multiplier <- function(
    spec,
    data,
    theta_hat,
    ks_prep = NULL,
    cvm_prep = NULL,
    control = list()) {
  require_sunspots_joint_parsimonious_spec_dependencies()
  normalized <- normalize_sunspots_joint_spec_data(data, control)
  if (!is.list(theta_hat) ||
      is.null(theta_hat$eta_hat) ||
      is.null(theta_hat$theta_hat) ||
      is.null(theta_hat$time_model)) {
    stop("`theta_hat` must be a fitted parsimonious joint object.")
  }
  if (!sunspots_joint_parsimonious_time_fast_regular(
    theta_hat$eta_hat,
    time_model = theta_hat$time_model,
    control = control
  )) {
    stop(
      paste(
        "The parsimonious temporal MLE is on a boundary or is not",
        "locally identified; fast multiplier preparation is invalid."
      ),
      call. = FALSE
    )
  }

  par0 <- sunspots_joint_parsimonious_pack_par(
    theta_hat, control = control
  )
  prepare_fast_multiplier_score_model(
    spec = spec,
    data = normalized,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = function(sample, par) {
      sunspots_joint_parsimonious_score_matrix(
        sample, par, control = control
      )
    },
    sample_fn = function(n_aux, par) {
      sample_sunspots_joint_time_space_parsimonious(
        n_aux, par, control = control
      )
    }
  )
}

make_sunspots_joint_time_space_parsimonious_spec <- function(
    time_model = c("beta", "uniform_beta"),
    hemisphere_regression = c("shared", "asymmetric")) {
  time_model <- normalize_sunspots_joint_parsimonious_time_model(
    time_model
  )
  hemisphere_regression <-
    sunspots_time_varying_normalize_hemisphere_regression(
      match.arg(hemisphere_regression)
    )
  spec_name <- sprintf(
    "sunspots_joint_time_space_%s_%s",
    time_model,
    hemisphere_regression
  )

  add_control <- function(control) {
    utils::modifyList(control, list(
      time_model = time_model,
      hemisphere_regression = hemisphere_regression
    ))
  }

  new_model_spec(
    name = spec_name,
    fit_theta = function(
        data,
        weights = NULL,
        null,
        control = list()) {
      fit_sunspots_joint_parsimonious_spec_theta(
        data = data,
        weights = weights,
        null = null,
        control = add_control(control)
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      sunspots_joint_spec_distance_matrix(
        data, omega, control = add_control(control)
      )
    },
    profile_eval = function(omega, t, theta, control = list()) {
      sunspots_joint_parsimonious_spec_profile_eval(
        omega, t, theta, control = add_control(control)
      )
    },
    normalize_data = function(data, control = list()) {
      normalize_sunspots_joint_spec_data(
        data, control = add_control(control)
      )
    },
    n_obs = function(data, control = list()) {
      normalized <- normalize_sunspots_joint_spec_data(
        data, add_control(control)
      )
      nrow(normalized$x)
    },
    observation_at = function(data, idx, control = list()) {
      normalized <- normalize_sunspots_joint_spec_data(
        data, add_control(control)
      )
      sunspots_joint_spec_obs_row(normalized, idx)
    },
    extras = list(
      sample_profile_matrix_eval = function(
          data,
          distance_matrix,
          theta,
          control = list()) {
        sunspots_joint_parsimonious_spec_profile_matrix_eval(
          data,
          distance_matrix,
          theta,
          control = add_control(control)
        )
      },
      sample_profile_sorted_prepare = function(
          data,
          sorted_distance_matrix,
          theta,
          control = list()) {
        sunspots_joint_parsimonious_spec_sorted_prepare(
          data,
          sorted_distance_matrix,
          theta,
          control = add_control(control)
        )
      },
      sample_profile_sorted_block_eval = function(
          data,
          sorted_distance_matrix,
          theta,
          row_indices,
          prepared = NULL,
          control = list()) {
        sunspots_joint_parsimonious_spec_sorted_block_eval(
          data,
          sorted_distance_matrix,
          theta,
          row_indices,
          prepared = prepared,
          control = add_control(control)
        )
      },
      fast_multiplier_prepare = function(
          data,
          theta_hat,
          ks_prep = NULL,
          cvm_prep = NULL,
          control = list()) {
        control <- add_control(control)
        prepare_sunspots_joint_parsimonious_fast_multiplier(
          spec = make_sunspots_joint_time_space_parsimonious_spec(
            time_model = time_model,
            hemisphere_regression = hemisphere_regression
          ),
          data = data,
          theta_hat = theta_hat,
          ks_prep = ks_prep,
          cvm_prep = cvm_prep,
          control = control
        )
      },
      weighted_mle = FALSE
    )
  )
}
