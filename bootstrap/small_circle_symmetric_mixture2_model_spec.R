# Symmetric two-small-circles model adapter for multiplier bootstrap GOF tests

resolve_small_circle_symmetric_mixture2_model_spec_path <- function(...) {
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

model_specs_path_small_circle_symmetric_mixture2 <- resolve_small_circle_symmetric_mixture2_model_spec_path("bootstrap", "model_specs.R")
if (!exists("new_model_spec", mode = "function")) {
  source(model_specs_path_small_circle_symmetric_mixture2)
}

normalize_small_circle_symmetric_mixture2_data <- function(data, control = list()) {
  x <- jp_normalize_unit_matrix(data, arg_name = "`data`", min_ncol = 3L)
  if (ncol(x) != 3L) {
    stop("Symmetric small-circle-mixture data must be a non-empty S^2 sample with three columns.")
  }
  x
}

normalize_small_circle_symmetric_mixture2_theta <- function(theta, ambient_dim = 3L) {
  if (!is.list(theta)) {
    stop("Symmetric small-circle-mixture theta must be a list containing `mu`, `kappa` and `nu`.")
  }

  params <- small_circle_symmetric_mixture2_normalize_theta(theta, ambient_dim = ambient_dim)
  list(
    mu = params$mu,
    kappa = params$kappa,
    nu = params$nu,
    ambient_dim = length(params$mu)
  )
}

small_circle_symmetric_mixture2_range_string <- function(x) {
  x <- as.numeric(x)
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0L) {
    return("no finite values")
  }
  sprintf("[%.16g, %.16g]", min(finite_x), max(finite_x))
}

small_circle_symmetric_mixture2_has_nonfinite_theta <- function(theta) {
  if (is.null(theta) || !is.list(theta)) {
    return(TRUE)
  }
  components <- c(theta$mu, theta$kappa, theta$nu)
  any(!is.finite(as.numeric(components)))
}

emit_small_circle_symmetric_mixture2_bootstrap_diag <- function(context,
                                                                replicate_index = NA_integer_,
                                                                block_label = NULL,
                                                                theta_star = NULL,
                                                                mle_loglik = NA_real_,
                                                                mle_convergence = NA_integer_,
                                                                x = NULL,
                                                                dot_matrix = NULL,
                                                                dot_block = NULL,
                                                                empirical_block = NULL,
                                                                theoretical_block = NULL,
                                                                empirical_obs_block = NULL,
                                                                theoretical_obs_block = NULL,
                                                                process_block = NULL,
                                                                warnings = character(),
                                                                extra = NULL) {
  header <- sprintf(
    "[small_circle_symmetric_mixture2][bootstrap diag] context=%s | replicate=%s%s",
    context,
    if (is.finite(replicate_index)) as.integer(replicate_index) else "NA",
    if (!is.null(block_label)) paste0(" | block=", block_label) else ""
  )
  message(header)
  if (!is.null(theta_star)) {
    message(sprintf(
      "  theta_star: mu=(%.16g, %.16g, %.16g), kappa=%.16g, nu=%.16g",
      as.numeric(theta_star$mu[[1L]]),
      as.numeric(theta_star$mu[[2L]]),
      as.numeric(theta_star$mu[[3L]]),
      as.numeric(theta_star$kappa),
      as.numeric(theta_star$nu)
    ))
  } else {
    message("  theta_star: NULL")
  }
  message(sprintf(
    "  mle: loglik=%s, convergence=%s",
    format(mle_loglik, digits = 16),
    if (is.na(mle_convergence)) "NA" else as.character(as.integer(mle_convergence))
  ))
  if (!is.null(x)) {
    message(sprintf("  X_star range: not available in multiplier bootstrap; original x range=%s", small_circle_symmetric_mixture2_range_string(x)))
  }
  if (!is.null(dot_matrix)) {
    dot_matrix_text <- if (is.character(dot_matrix) && length(dot_matrix) == 1L) {
      dot_matrix
    } else {
      small_circle_symmetric_mixture2_range_string(dot_matrix)
    }
    message(sprintf("  dot_matrix range: %s", dot_matrix_text))
  }
  if (!is.null(dot_block)) {
    message(sprintf("  dot_block range: %s", small_circle_symmetric_mixture2_range_string(dot_block)))
  }
  if (!is.null(empirical_block)) {
    message(sprintf("  profile empirical range: %s", small_circle_symmetric_mixture2_range_string(empirical_block)))
  }
  if (!is.null(theoretical_block)) {
    message(sprintf("  profile theoretical range: %s", small_circle_symmetric_mixture2_range_string(theoretical_block)))
  }
  if (!is.null(empirical_obs_block)) {
    message(sprintf("  profile empirical_obs range: %s", small_circle_symmetric_mixture2_range_string(empirical_obs_block)))
  }
  if (!is.null(theoretical_obs_block)) {
    message(sprintf("  profile theoretical_obs range: %s", small_circle_symmetric_mixture2_range_string(theoretical_obs_block)))
  }
  if (!is.null(process_block)) {
    message(sprintf("  process block range: %s", small_circle_symmetric_mixture2_range_string(process_block)))
  }
  if (length(warnings) > 0L) {
    for (warning_text in warnings) {
      message(sprintf("  warning: %s", warning_text))
    }
  }
  if (!is.null(extra) && length(extra) > 0L) {
    for (line in as.character(extra)) {
      message(sprintf("  %s", line))
    }
  }
}

small_circle_symmetric_mixture2_full_dot_range <- function(x) {
  dot_matrix <- pmin(pmax(x %*% t(x), -1), 1)
  small_circle_symmetric_mixture2_range_string(dot_matrix)
}

small_circle_symmetric_mixture2_empirical_profile_block <- function(dot_block,
                                                                    weights = NULL,
                                                                    return_order = FALSE) {
  dot_block <- as.matrix(dot_block)
  block_rows <- nrow(dot_block)
  n <- ncol(dot_block)

  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    weights <- as.numeric(weights)
  }

  if (length(weights) != n) {
    stop("`weights` must have length equal to ncol(`dot_block`).")
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and nonnegative.")
  }

  total_weight <- sum(weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("`weights` must have strictly positive finite sum.")
  }

  order_matrix <- t(vapply(seq_len(block_rows), function(i) {
    as.integer(order(dot_block[i, ], decreasing = TRUE))
  }, integer(n)))
  rank_matrix <- t(vapply(seq_len(block_rows), function(i) {
    as.integer(rank(-dot_block[i, ], ties.method = "max"))
  }, integer(n)))

  ordered_weights_matrix <- matrix(weights[order_matrix], nrow = block_rows, ncol = n)
  cumulative_weights_matrix <- ordered_weights_matrix
  if (n >= 2L) {
    for (j in 2:n) {
      cumulative_weights_matrix[, j] <- cumulative_weights_matrix[, j] +
        cumulative_weights_matrix[, j - 1L]
    }
  }

  row_index_matrix <- matrix(rep.int(seq_len(block_rows), n), nrow = block_rows, ncol = n)
  local_linear_index <- row_index_matrix + (rank_matrix - 1L) * block_rows
  profile_block <- matrix(cumulative_weights_matrix[local_linear_index] / total_weight, nrow = block_rows, ncol = n)

  if (any(!is.finite(profile_block)) || any(profile_block < -1e-12) || any(profile_block > 1 + 1e-12)) {
    stop("Empirical profile block produced values outside [0, 1].")
  }

  if (isTRUE(return_order)) {
    list(profile = profile_block, order_matrix = order_matrix, rank_matrix = rank_matrix)
  } else {
    profile_block
  }
}

small_circle_symmetric_mixture2_cvm_prepare <- function(data,
                                                        theta_hat,
                                                        distance_type = "geodesic",
                                                        control = list()) {
  x <- normalize_small_circle_symmetric_mixture2_data(data, control)
  theta_hat <- normalize_small_circle_symmetric_mixture2_theta(theta_hat, ambient_dim = ncol(x))
  n <- nrow(x)
  block_size <- as.integer(control$small_circle_symmetric_mixture2_cvm_block_size %||% 128L)
  block_starts <- seq.int(1L, n, by = max(block_size, 1L))
  cvm_sum <- 0
  for (block_start in block_starts) {
    block_end <- min(block_start + block_size - 1L, n)
    block_idx <- block_start:block_end
    x_block <- x[block_idx, , drop = FALSE]
    dot_block <- pmin(pmax(x_block %*% t(x), -1), 1)
    empirical_block <- small_circle_symmetric_mixture2_empirical_profile_block(dot_block)
    theoretical_block <- small_circle_symmetric_mixture2_cvm_profile_block(
      X_block = x_block,
      dot_threshold_block = dot_block,
      mu = theta_hat$mu,
      kappa = theta_hat$kappa,
      nu = theta_hat$nu,
      distance_type = distance_type,
      method = control$small_circle_symmetric_mixture2_profile_method %||% "legendre",
      l_max = as.integer(control$small_circle_symmetric_mixture2_L_max %||% 200L),
      quad_n = as.integer(control$small_circle_symmetric_mixture2_quad_n %||% 400L),
      tol = as.numeric(control$small_circle_symmetric_mixture2_tol %||% 1e-10),
      control = control
    )
    process_block <- sqrt(n) * (empirical_block - theoretical_block)
    cvm_sum <- cvm_sum + sum(process_block^2)
  }

  list(
    x = x,
    block_size = block_size,
    theta_hat = theta_hat,
    statistic = cvm_sum / (n * n),
    observed_process_available = FALSE
  )
}

small_circle_symmetric_mixture2_cvm_bootstrap_stat <- function(data,
                                                               normalized_weights,
                                                               theta_star,
                                                               cvm_prep,
                                                               null,
                                                               distance_type = "geodesic",
                                                               control = list(),
                                                               scale_factor = 1) {
  x <- normalize_small_circle_symmetric_mixture2_data(data, control)
  replicate_index <- as.integer(control$small_circle_symmetric_mixture2_bootstrap_replicate_index %||% NA_integer_)
  bootstrap_warnings <- as.character(control$small_circle_symmetric_mixture2_bootstrap_warnings %||% character())
  mle_loglik <- as.numeric(control$small_circle_symmetric_mixture2_bootstrap_loglik %||% NA_real_)
  mle_convergence <- as.integer(control$small_circle_symmetric_mixture2_bootstrap_convergence %||% NA_integer_)
  theta_star <- if (!is.null(theta_star)) {
    normalize_small_circle_symmetric_mixture2_theta(theta_star, ambient_dim = ncol(x))
  } else {
    NULL
  }
  if (!identical(null$type, "simple") && small_circle_symmetric_mixture2_has_nonfinite_theta(theta_star)) {
    emit_small_circle_symmetric_mixture2_bootstrap_diag(
      context = "theta_star_nonfinite",
      replicate_index = replicate_index,
      theta_star = theta_star,
      mle_loglik = mle_loglik,
      mle_convergence = mle_convergence,
      x = x,
      dot_matrix = small_circle_symmetric_mixture2_full_dot_range(x),
      warnings = bootstrap_warnings
    )
    return(NA_real_)
  }
  if (any(!is.finite(as.numeric(x)))) {
    emit_small_circle_symmetric_mixture2_bootstrap_diag(
      context = "x_nonfinite",
      replicate_index = replicate_index,
      theta_star = theta_star,
      mle_loglik = mle_loglik,
      mle_convergence = mle_convergence,
      x = x,
      dot_matrix = small_circle_symmetric_mixture2_full_dot_range(x),
      warnings = bootstrap_warnings
    )
    return(NA_real_)
  }
  theta_hat <- normalize_small_circle_symmetric_mixture2_theta(cvm_prep$theta_hat, ambient_dim = ncol(x))
  n <- nrow(x)
  block_size <- as.integer(cvm_prep$block_size)
  block_starts <- seq.int(1L, n, by = max(block_size, 1L))
  cvm_sum <- 0

  for (block_start in block_starts) {
    block_end <- min(block_start + block_size - 1L, n)
    block_idx <- block_start:block_end
    block_label <- sprintf("%d:%d", block_start, block_end)
    x_block <- x[block_idx, , drop = FALSE]
    dot_block <- pmin(pmax(x_block %*% t(x), -1), 1)
    empirical_star_block <- small_circle_symmetric_mixture2_empirical_profile_block(dot_block, weights = normalized_weights)
    empirical_obs_block <- small_circle_symmetric_mixture2_empirical_profile_block(dot_block)

    if (any(!is.finite(empirical_star_block)) ||
        any(empirical_star_block < -1e-12) ||
        any(empirical_star_block > 1 + 1e-12)) {
      emit_small_circle_symmetric_mixture2_bootstrap_diag(
        context = "empirical_star_block_invalid",
        replicate_index = replicate_index,
        block_label = block_label,
        theta_star = theta_star,
        mle_loglik = mle_loglik,
        mle_convergence = mle_convergence,
        x = x,
        dot_matrix = small_circle_symmetric_mixture2_full_dot_range(x),
        dot_block = dot_block,
        empirical_block = empirical_star_block,
        empirical_obs_block = empirical_obs_block,
        warnings = bootstrap_warnings
      )
      return(NA_real_)
    }

    if (identical(null$type, "simple")) {
      process_block <- scale_factor * sqrt(n) * (empirical_star_block - empirical_obs_block)
      theoretical_star_block <- NULL
      theoretical_obs_block <- NULL
    } else {
      theoretical_star_block <- small_circle_symmetric_mixture2_cvm_profile_block(
        X_block = x[block_idx, , drop = FALSE],
        dot_threshold_block = dot_block,
        mu = theta_star$mu,
        kappa = theta_star$kappa,
        nu = theta_star$nu,
        distance_type = distance_type,
        method = control$small_circle_symmetric_mixture2_profile_method %||% "legendre",
        l_max = as.integer(control$small_circle_symmetric_mixture2_L_max %||% 200L),
        quad_n = as.integer(control$small_circle_symmetric_mixture2_quad_n %||% 400L),
        tol = as.numeric(control$small_circle_symmetric_mixture2_tol %||% 1e-10),
        control = control
      )
      theoretical_obs_block <- small_circle_symmetric_mixture2_cvm_profile_block(
        X_block = x[block_idx, , drop = FALSE],
        dot_threshold_block = dot_block,
        mu = theta_hat$mu,
        kappa = theta_hat$kappa,
        nu = theta_hat$nu,
        distance_type = distance_type,
        method = control$small_circle_symmetric_mixture2_profile_method %||% "legendre",
        l_max = as.integer(control$small_circle_symmetric_mixture2_L_max %||% 200L),
        quad_n = as.integer(control$small_circle_symmetric_mixture2_quad_n %||% 400L),
        tol = as.numeric(control$small_circle_symmetric_mixture2_tol %||% 1e-10),
        control = control
      )
      process_block <- scale_factor * sqrt(n) * (
        (empirical_star_block - theoretical_star_block) -
          (empirical_obs_block - theoretical_obs_block)
      )
    }

    if ((!is.null(theoretical_star_block) &&
         (any(!is.finite(theoretical_star_block)) ||
            any(theoretical_star_block < -1e-12) ||
            any(theoretical_star_block > 1 + 1e-12))) ||
        any(!is.finite(empirical_obs_block)) ||
        any(empirical_obs_block < -1e-12) ||
        any(empirical_obs_block > 1 + 1e-12) ||
        (!is.null(theoretical_obs_block) &&
         (any(!is.finite(theoretical_obs_block)) ||
            any(theoretical_obs_block < -1e-12) ||
            any(theoretical_obs_block > 1 + 1e-12))) ||
        any(!is.finite(process_block))) {
      emit_small_circle_symmetric_mixture2_bootstrap_diag(
        context = "block_invalid",
        replicate_index = replicate_index,
        block_label = block_label,
        theta_star = theta_star,
        mle_loglik = mle_loglik,
        mle_convergence = mle_convergence,
        x = x,
        dot_matrix = small_circle_symmetric_mixture2_full_dot_range(x),
        dot_block = dot_block,
        empirical_block = empirical_star_block,
        theoretical_block = theoretical_star_block,
        empirical_obs_block = empirical_obs_block,
        theoretical_obs_block = theoretical_obs_block,
        process_block = process_block,
        warnings = bootstrap_warnings
      )
      return(NA_real_)
    }

    cvm_sum <- cvm_sum + sum(process_block^2)
  }
  bootstrap_stat <- cvm_sum / (n * n)
  if (!is.finite(bootstrap_stat)) {
    emit_small_circle_symmetric_mixture2_bootstrap_diag(
      context = "bootstrap_stat_nonfinite",
      replicate_index = replicate_index,
      theta_star = theta_star,
      mle_loglik = mle_loglik,
      mle_convergence = mle_convergence,
      x = x,
      dot_matrix = small_circle_symmetric_mixture2_full_dot_range(x),
      warnings = bootstrap_warnings,
      extra = sprintf("bootstrap_stat=%s", format(bootstrap_stat, digits = 16))
    )
    return(NA_real_)
  }

  bootstrap_stat
}

fit_small_circle_symmetric_mixture2_theta <- function(data,
                                                      weights = NULL,
                                                      null,
                                                      control = list()) {
  x <- normalize_small_circle_symmetric_mixture2_data(data, control)

  if (!is.list(null) || is.null(null$type)) {
    stop("`null` must be a list containing at least the field `type`.")
  }

  if (identical(null$type, "simple")) {
    return(normalize_small_circle_symmetric_mixture2_theta(null$theta, ambient_dim = ncol(x)))
  }

  if (!identical(null$type, "composite")) {
    stop("`null$type` must be either `simple` or `composite`.")
  }

  fit <- small_circle_symmetric_mixture2_mle_s2_weighted(
    x = x,
    weights = weights,
    control = control
  )
  theta <- normalize_small_circle_symmetric_mixture2_theta(fit, ambient_dim = ncol(x))
  c(theta, fit[setdiff(names(fit), names(theta))])
}

make_small_circle_symmetric_mixture2_spec <- function(distance_type = c("chordal", "geodesic")) {
  distance_type <- match.arg(distance_type)

  new_model_spec(
    name = sprintf("small_circle_symmetric_mixture2_%s", distance_type),
    fit_theta = function(data, weights = NULL, null, control = list()) {
      fit_small_circle_symmetric_mixture2_theta(
        data = data,
        weights = weights,
        null = null,
        control = control
      )
    },
    distance_matrix = function(data, omega, control = list()) {
      x <- normalize_small_circle_symmetric_mixture2_data(data, control)
      omega_matrix <- normalize_small_circle_symmetric_mixture2_data(omega, control)

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
      theta <- normalize_small_circle_symmetric_mixture2_theta(theta)
      distance_profile_small_circle_symmetric_mixture2(
        omega = omega,
        t_values = as.numeric(t),
        mu = theta$mu,
        kappa = theta$kappa,
        nu = theta$nu,
        distance_type = distance_type,
        method = control$small_circle_symmetric_mixture2_profile_method %||% "legendre",
        l_max = as.integer(control$small_circle_symmetric_mixture2_L_max %||% 200L),
        quad_n = as.integer(control$small_circle_symmetric_mixture2_quad_n %||% 400L),
        tol = as.numeric(control$small_circle_symmetric_mixture2_tol %||% 1e-10),
        validate_against_integral = isTRUE(control$small_circle_symmetric_mixture2_validate_against_integral),
        validation_tol = as.numeric(control$small_circle_symmetric_mixture2_validation_tol %||% 5e-6)
      )
    },
    normalize_data = normalize_small_circle_symmetric_mixture2_data,
    n_obs = function(data, control = list()) {
      nrow(normalize_small_circle_symmetric_mixture2_data(data, control))
    },
    observation_at = function(data, idx, control = list()) {
      normalize_small_circle_symmetric_mixture2_data(data, control)[idx, , drop = TRUE]
    },
    extras = list(
      profile_matrix_eval = function(omega_grid, t_grid, theta, control = list()) {
        theta <- normalize_small_circle_symmetric_mixture2_theta(theta)
        distance_profile_small_circle_symmetric_mixture2_grid(
          omega_grid = omega_grid,
          mu = theta$mu,
          kappa = theta$kappa,
          nu = theta$nu,
          t_grid = t_grid,
          distance_type = distance_type,
          method = control$small_circle_symmetric_mixture2_profile_method %||% "legendre",
          l_max = as.integer(control$small_circle_symmetric_mixture2_L_max %||% 200L),
          quad_n = as.integer(control$small_circle_symmetric_mixture2_quad_n %||% 400L),
          tol = as.numeric(control$small_circle_symmetric_mixture2_tol %||% 1e-10)
        )
      },
      sample_profile_matrix_eval = function(data, distance_matrix, theta, control = list()) {
        if (!identical(distance_type, "geodesic")) {
          return(NULL)
        }

        theta <- normalize_small_circle_symmetric_mixture2_theta(theta)
        distance_profile_small_circle_symmetric_mixture2_cvm_grid(
          X = data,
          mu = theta$mu,
          kappa = theta$kappa,
          nu = theta$nu,
          distance_matrix = distance_matrix,
          distance_type = distance_type,
          control = control,
          method = control$small_circle_symmetric_mixture2_profile_method %||% "legendre",
          l_max = as.integer(control$small_circle_symmetric_mixture2_L_max %||% 200L),
          quad_n = as.integer(control$small_circle_symmetric_mixture2_quad_n %||% 400L),
          tol = as.numeric(control$small_circle_symmetric_mixture2_tol %||% 1e-10)
        )
      },
      cvm_prepare = function(data, theta_hat, control = list()) {
        small_circle_symmetric_mixture2_cvm_prepare(
          data = data,
          theta_hat = theta_hat,
          distance_type = distance_type,
          control = control
        )
      },
      cvm_bootstrap_stat = function(data, normalized_weights, theta_star, cvm_prep, null, control = list(), scale_factor = 1) {
        small_circle_symmetric_mixture2_cvm_bootstrap_stat(
          data = data,
          normalized_weights = normalized_weights,
          theta_star = theta_star,
          cvm_prep = cvm_prep,
          null = null,
          distance_type = distance_type,
          control = control,
          scale_factor = scale_factor
        )
      },
      distance_type = distance_type,
      unknown_param = "theta",
      weighted_mle = TRUE
    )
  )
}
