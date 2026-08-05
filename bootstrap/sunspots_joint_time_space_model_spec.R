# Sunspots joint time-space model adapter for multiplier bootstrap GOF tests

resolve_sunspots_joint_spec_path <- function(...) {
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

if (!exists("new_model_spec", mode = "function")) {
  source(resolve_sunspots_joint_spec_path("bootstrap", "model_specs.R"))
}

require_sunspots_joint_spec_dependencies <- function() {
  required <- c(
    "sunspots_joint_validate_data",
    "fit_sunspots_cycle23_joint_time_space",
    "sunspots_joint_distance",
    "sunspots_joint_time_quadrature",
    "sunspots_joint_conditional_legendre_coefficients",
    "sunspots_joint_profile_block",
    "sunspots_joint_pack_par",
    "sunspots_joint_score_matrix",
    "sample_sunspots_joint_time_space"
  )
  missing <- required[!vapply(required, exists, logical(1L), mode = "function")]
  if (length(missing) > 0L) {
    stop(
      sprintf(
        paste(
          "Sunspots joint model dependencies are missing:",
          "%s.",
          "Source real_data/sunspots/sunspots_cycle23_joint_time_space.R before constructing this spec."
        ),
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

normalize_sunspots_joint_spec_data <- function(data, control = list()) {
  require_sunspots_joint_spec_dependencies()
  if (is.list(data) && !is.null(data$x) && !is.null(data$s)) {
    return(sunspots_joint_validate_data(data$x, data$s))
  }

  mat <- as.matrix(data)
  if (!is.matrix(mat) || ncol(mat) != 4L || nrow(mat) == 0L || any(!is.finite(mat))) {
    stop("Joint sunspots data must be a non-empty finite matrix with four columns (x1, x2, x3, s), or a list with fields `x` and `s`.")
  }

  sunspots_joint_validate_data(
    x = mat[, 1:3, drop = FALSE],
    s = mat[, 4L]
  )
}

sunspots_joint_spec_obs_row <- function(data, idx) {
  c(data$x[idx, , drop = TRUE], data$s[idx])
}

fit_sunspots_joint_spec_theta <- function(data,
                                          weights = NULL,
                                          null,
                                          control = list()) {
  require_sunspots_joint_spec_dependencies()
  normalized <- normalize_sunspots_joint_spec_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least `type`.")
  }

  if (identical(null$type, "simple")) {
    theta_simple <- null$theta
    if (!is.list(theta_simple) || is.null(theta_simple$eta_hat) || is.null(theta_simple$theta_hat)) {
      stop("For `null$type = 'simple'`, `null$theta` must contain `eta_hat` and `theta_hat`.")
    }
    return(fit_sunspots_cycle23_joint_time_space(
      x = normalized$x,
      s = normalized$s,
      hemisphere_regression = control$hemisphere_regression %||% "asymmetric",
      control = control,
      eta_hat = theta_simple$eta_hat,
      theta_hat = theta_simple$theta_hat
    ))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  if (!is.null(weights)) {
    stop("The sunspots joint model spec only supports `bootstrap_method = 'fast_multiplier'`; weighted composite refits are not implemented.")
  }

  fit_sunspots_cycle23_joint_time_space(
    x = normalized$x,
    s = normalized$s,
    hemisphere_regression = control$hemisphere_regression %||% "asymmetric",
    control = control
  )
}

sunspots_joint_spec_distance_matrix <- function(data, omega, control = list()) {
  require_sunspots_joint_spec_dependencies()
  left <- normalize_sunspots_joint_spec_data(data, control)
  right <- normalize_sunspots_joint_spec_data(omega, control)

  out <- matrix(0, nrow = nrow(left$x), ncol = nrow(right$x))
  for (j in seq_len(nrow(right$x))) {
    out[, j] <- sunspots_joint_distance(
      x = left$x,
      s = left$s,
      omega = right$x[j, , drop = TRUE],
      center_s = right$s[j]
    )
  }
  out
}

sunspots_joint_spec_profile_eval <- function(omega, t, theta, control = list()) {
  require_sunspots_joint_spec_dependencies()
  omega <- as.numeric(omega)
  if (length(omega) != 4L || any(!is.finite(omega))) {
    stop("`omega` must be a finite vector with four entries (x1, x2, x3, s).")
  }
  center <- sunspots_joint_validate_data(
    x = matrix(omega[1:3], nrow = 1L),
    s = omega[[4L]]
  )
  fit <- theta
  if (!is.list(fit) || is.null(fit$eta_hat) || is.null(fit$theta_hat)) {
    stop("`theta` must be the fitted joint object containing `eta_hat` and `theta_hat`.")
  }

  t <- as.numeric(t)
  quadrature <- sunspots_joint_time_quadrature(
    fit$eta_hat,
    n_nodes = as.integer(control$time_quad_n %||% 64L),
    control = control
  )
  coefficients <- sunspots_joint_conditional_legendre_coefficients(
    fit$theta_hat,
    quadrature$nodes,
    l_max = as.integer(control$profile_l_max %||% 100L),
    quad_n = as.integer(control$profile_quad_n %||% 400L)
  )

  as.numeric(sunspots_joint_profile_block(
    radii = matrix(t, nrow = 1L),
    rho = center$x[1L, 3L],
    center_s = center$s,
    time_nodes = quadrature$nodes,
    time_weights = quadrature$weights,
    coefficients = coefficients,
    backend = control$distance_profile_backend %||% "auto"
  ))
}

sunspots_joint_spec_sample_profile_matrix_eval <- function(data,
                                                           distance_matrix,
                                                           theta,
                                                           control = list()) {
  require_sunspots_joint_spec_dependencies()
  normalized <- normalize_sunspots_joint_spec_data(data, control)
  fit <- theta
  if (!is.list(fit) || is.null(fit$eta_hat) || is.null(fit$theta_hat)) {
    stop("`theta` must be the fitted joint object containing `eta_hat` and `theta_hat`.")
  }

  distance_matrix <- as.matrix(distance_matrix)
  if (!identical(dim(distance_matrix), c(nrow(normalized$x), nrow(normalized$x)))) {
    stop("`distance_matrix` has incompatible dimensions for sample-profile evaluation.")
  }

  quadrature <- sunspots_joint_time_quadrature(
    fit$eta_hat,
    n_nodes = as.integer(control$time_quad_n %||% 64L),
    control = control
  )
  coefficients <- sunspots_joint_conditional_legendre_coefficients(
    fit$theta_hat,
    quadrature$nodes,
    l_max = as.integer(control$profile_l_max %||% 100L),
    quad_n = as.integer(control$profile_quad_n %||% 400L)
  )

  block_size <- as.integer(control$center_block_size %||% 8L)
  if (!is.finite(block_size) || block_size <= 0L) {
    stop("`control$center_block_size` must be a strictly positive integer.")
  }

  n <- nrow(normalized$x)
  out <- matrix(0, nrow = n, ncol = n)
  for (block_start in seq.int(1L, n, by = block_size)) {
    block_end <- min(block_start + block_size - 1L, n)
    idx <- block_start:block_end
    out[idx, ] <- sunspots_joint_profile_block(
      radii = distance_matrix[idx, , drop = FALSE],
      rho = normalized$x[idx, 3L],
      center_s = normalized$s[idx],
      time_nodes = quadrature$nodes,
      time_weights = quadrature$weights,
      coefficients = coefficients,
      backend = control$distance_profile_backend %||% "auto"
    )
  }

  out
}

prepare_sunspots_joint_fast_multiplier <- function(spec,
                                                   data,
                                                   theta_hat,
                                                   ks_prep = NULL,
                                                   cvm_prep = NULL,
                                                   control = list()) {
  require_sunspots_joint_spec_dependencies()
  normalized <- normalize_sunspots_joint_spec_data(data, control)

  if (!is.list(theta_hat) || is.null(theta_hat$eta_hat) || is.null(theta_hat$theta_hat)) {
    stop("`theta_hat` must be a fitted joint object containing `eta_hat` and `theta_hat`.")
  }

  boundary_flags <- theta_hat$eta_hat$boundary_flags %||% list()
  if (isTRUE(boundary_flags$weight) ||
      isTRUE(boundary_flags$shape_lower) ||
      isTRUE(boundary_flags$shape_upper)) {
    stop(
      paste(
        "Temporal beta-mixture MLE is on an admissible boundary;",
        "fast multiplier preparation is invalid for the joint model."
      ),
      call. = FALSE
    )
  }

  par0 <- sunspots_joint_pack_par(theta_hat, control = control)

  prepare_fast_multiplier_score_model(
    spec = spec,
    data = normalized,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = cvm_prep,
    control = control,
    par0 = par0,
    score_matrix_fn = function(data, par) {
      sunspots_joint_score_matrix(data, par, control = control)
    },
    sample_fn = function(n_aux, par) {
      sample_sunspots_joint_time_space(n_aux, par, control = control)
    }
  )
}

make_sunspots_joint_time_space_spec <- function(
    hemisphere_regression = c("asymmetric", "shared")) {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(
    match.arg(hemisphere_regression)
  )

  new_model_spec(
    name = sprintf("sunspots_joint_time_space_%s", hemisphere_regression),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      control <- utils::modifyList(control, list(
        hemisphere_regression = hemisphere_regression
      ))
      fit_sunspots_joint_spec_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      control <- utils::modifyList(control, list(
        hemisphere_regression = hemisphere_regression
      ))
      sunspots_joint_spec_distance_matrix(data, omega, control = control)
    },
    profile_eval = function(omega, t, theta, control = list()) {
      control <- utils::modifyList(control, list(
        hemisphere_regression = hemisphere_regression
      ))
      sunspots_joint_spec_profile_eval(omega, t, theta, control = control)
    },
    normalize_data = function(data, control = list()) {
      control <- utils::modifyList(control, list(
        hemisphere_regression = hemisphere_regression
      ))
      normalize_sunspots_joint_spec_data(data, control = control)
    },
    n_obs = function(data, control = list()) {
      normalized <- normalize_sunspots_joint_spec_data(data, control)
      nrow(normalized$x)
    },
    observation_at = function(data, idx, control = list()) {
      normalized <- normalize_sunspots_joint_spec_data(data, control)
      sunspots_joint_spec_obs_row(normalized, idx)
    },
    extras = list(
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        control <- utils::modifyList(control, list(
          hemisphere_regression = hemisphere_regression
        ))
        sunspots_joint_spec_sample_profile_matrix_eval(
          data = data,
          distance_matrix = distance_matrix,
          theta = theta,
          control = control
        )
      },
      fast_multiplier_prepare = function(data,
                                         theta_hat,
                                         ks_prep = NULL,
                                         cvm_prep = NULL,
                                         control = list()) {
        control <- utils::modifyList(control, list(
          hemisphere_regression = hemisphere_regression
        ))
        prepare_sunspots_joint_fast_multiplier(
          spec = make_sunspots_joint_time_space_spec(
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
