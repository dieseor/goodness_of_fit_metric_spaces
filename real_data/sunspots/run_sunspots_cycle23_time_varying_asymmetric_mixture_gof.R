#!/usr/bin/env Rscript

# Fast global KS test for the time-varying cycle-23 small-circle model.
# The temporal ranks U are treated as a fixed design. The null profile is the
# global law induced by averaging P_theta{d(X, omega) <= t | U_i} over all
# retained design points, deliberately avoiding a time-local/windowed test.

resolve_sunspots_time_gof_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

time_model_path_sunspots_gof <- resolve_sunspots_time_gof_path(
  "real_data", "sunspots", "run_sunspots_cycle23_time_varying_asymmetric_mixture.R"
)
model_specs_path_sunspots_gof <- resolve_sunspots_time_gof_path("bootstrap", "model_specs.R")
bootstrap_path_sunspots_gof <- resolve_sunspots_time_gof_path("bootstrap", "multiplier_bootstrap.R")
source(time_model_path_sunspots_gof)
source(model_specs_path_sunspots_gof)
source(bootstrap_path_sunspots_gof)

normalize_sunspots_time_gof_data <- function(x, u) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  if (ncol(x) != 3L) stop("The conditional sunspots GOF supports only S^2 data.")
  u <- as.numeric(u)
  if (length(u) != nrow(x) || any(!is.finite(u)) || any(u <= 0) || any(u >= 1)) {
    stop("`u` must have one finite value in (0, 1) for each row of `x`.")
  }
  list(x = x, u = u)
}

sunspots_time_gof_path_state <- function(start_raw, end_raw, u, nu_eps) {
  prob_eps <- sqrt(.Machine$double.eps)
  p_start <- pmin(pmax(stats::plogis(start_raw), prob_eps), 1 - prob_eps)
  p_end <- pmin(pmax(stats::plogis(end_raw), prob_eps), 1 - prob_eps)
  range_nu <- 1 - 2 * nu_eps
  start <- nu_eps + range_nu * p_start
  end <- nu_eps + (start - nu_eps) * p_end
  d_start <- range_nu * p_start * (1 - p_start)
  d_end_start <- p_end * d_start
  d_end_end <- (start - nu_eps) * p_end * (1 - p_end)

  list(
    nu = (1 - u) * start + u * end,
    d_start_raw = (1 - u) * d_start + u * d_end_start,
    d_end_raw = u * d_end_end
  )
}

sunspots_time_gof_state_from_par <- function(par, control = list()) {
  nu_eps <- as.numeric(control$nu_eps %||% 1e-6)
  c_min <- as.numeric(control$c_min %||% 1e-8)
  c_max <- as.numeric(control$c_max %||% 1e6)
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(
    control$hemisphere_regression %||% "asymmetric"
  )
  theta <- sunspots_time_varying_unpack_par(
    par, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
  list(
    theta = theta,
    nu_eps = nu_eps,
    c_min = c_min,
    c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
}

sunspots_time_gof_score_matrix <- function(data, par, control = list()) {
  data <- normalize_sunspots_time_gof_data(data$x, data$u)
  state <- sunspots_time_gof_state_from_par(par, control)
  theta <- state$theta
  north_path <- sunspots_time_gof_path_state(par[[1L]], par[[2L]], data$u, state$nu_eps)
  south_path <- if (identical(state$hemisphere_regression, "shared")) {
    north_path
  } else {
    sunspots_time_gof_path_state(par[[3L]], par[[4L]], data$u, state$nu_eps)
  }
  z <- pmin(pmax(data$x[, 3L], -1), 1)

  log_north <- sunspots_time_varying_component_log_axis_density(z, theta$c, north_path$nu)
  log_south <- sunspots_time_varying_component_log_axis_density(-z, theta$c, south_path$nu)
  log_mix <- rotational_logsumexp2(log(0.5) + log_north, log(0.5) + log_south)
  responsibility_north <- exp(log(0.5) + log_north - log_mix)
  responsibility_south <- 1 - responsibility_north

  score_north <- fast_multiplier_small_circle_component_scores_natural(
    s = z, kappa = theta$c, nu = north_path$nu
  )
  score_south <- fast_multiplier_small_circle_component_scores_natural(
    s = -z, kappa = theta$c, nu = south_path$nu
  )
  c_derivative <- stats::plogis(par[[if (identical(state$hemisphere_regression, "shared")) 3L else 5L]])

  score_c <- (responsibility_north * score_north[, 1L] +
    responsibility_south * score_south[, 1L]) * c_derivative
  if (identical(state$hemisphere_regression, "shared")) {
    return(cbind(
      (responsibility_north * score_north[, 2L] + responsibility_south * score_south[, 2L]) *
        north_path$d_start_raw,
      (responsibility_north * score_north[, 2L] + responsibility_south * score_south[, 2L]) *
        north_path$d_end_raw,
      score_c
    ))
  }
  cbind(
    responsibility_north * score_north[, 2L] * north_path$d_start_raw,
    responsibility_north * score_north[, 2L] * north_path$d_end_raw,
    responsibility_south * score_south[, 2L] * south_path$d_start_raw,
    responsibility_south * score_south[, 2L] * south_path$d_end_raw,
    score_c
  )
}

sample_sunspots_time_gof_conditional <- function(n, u_reference, par, control = list()) {
  n <- as.integer(n)
  if (!is.finite(n) || n <= 0L) stop("`n` must be a strictly positive integer.")
  u_reference <- as.numeric(u_reference)
  if (length(u_reference) == 0L) stop("`u_reference` cannot be empty.")
  state <- sunspots_time_gof_state_from_par(par, control)
  theta <- state$theta
  u <- sample(u_reference, size = n, replace = TRUE)
  nu <- sunspots_time_varying_nu(u, theta)
  north <- stats::runif(n) <= 0.5
  x <- matrix(0, nrow = n, ncol = 3L)

  for (i in seq_len(n)) {
    x[i, ] <- r_sph_small_circle(
      n = 1L,
      mu = if (north[[i]]) c(0, 0, 1) else c(0, 0, -1),
      kappa = theta$c,
      nu = if (north[[i]]) nu$north[[i]] else nu$south[[i]],
      check = FALSE
    )
  }
  list(x = x, u = u)
}

sunspots_time_gof_average_legendre_coefficients <- function(theta,
                                                            u_reference,
                                                            l_max = 100L,
                                                            quad_n = 400L) {
  theta <- sunspots_time_varying_validate_theta(theta)
  u_reference <- as.numeric(u_reference)
  if (length(u_reference) == 0L || any(!is.finite(u_reference))) {
    stop("`u_reference` must contain finite values.")
  }

  quadrature <- small_circle_gauss_legendre(as.integer(quad_n))
  z_nodes <- quadrature$nodes
  nu <- sunspots_time_varying_nu(u_reference, theta)
  log_norm_north <- small_circle_log_norm_constant(theta$c, nu$north)
  log_norm_south <- small_circle_log_norm_constant(theta$c, nu$south)
  north_log_density <- sweep(
    -theta$c * outer(z_nodes, nu$north, FUN = "-")^2,
    2L,
    log(2) + log_norm_north,
    FUN = "-"
  )
  south_log_density <- sweep(
    -theta$c * outer(-z_nodes, nu$south, FUN = "-")^2,
    2L,
    log(2) + log_norm_south,
    FUN = "-"
  )
  average_density <- 0.5 * (rowMeans(exp(north_log_density)) + rowMeans(exp(south_log_density)))
  # The Legendre helper uses the rotational density with integral two, not the
  # axial probability density (which integrates to one).
  rotational_density <- 2 * average_density
  legendre <- small_circle_legendre_matrix(z_nodes, l_max = as.integer(l_max))
  ell <- 0:as.integer(l_max)
  coefficients <- ((2 * ell + 1) / 2) * as.numeric(crossprod(
    legendre,
    quadrature$weights * rotational_density
  ))
  coefficients[[1L]] <- 1
  coefficients
}

sunspots_time_gof_profile_row <- function(omega,
                                          sorted_distances,
                                          coefficients) {
  omega <- jp_normalize_unit_vector(omega, arg_name = "`omega`", min_length = 3L)
  sorted_distances <- as.numeric(sorted_distances)
  thresholds <- sphere_distance_to_dot_threshold(sorted_distances, distance_type = "geodesic")
  profile <- 1 - small_circle_projection_cdf_legendre(
    x = thresholds,
    r = omega[[3L]],
    coefficients = coefficients
  )
  small_circle_monotone_clip(sorted_distances, profile, upper_bound = pi)
}

sunspots_time_gof_select_centers <- function(n, n_sample_centers = 100L, seed = 20260711L) {
  if (is.infinite(n_sample_centers) || n_sample_centers >= n) return(seq_len(n))
  n_sample_centers <- as.integer(n_sample_centers)
  if (!is.finite(n_sample_centers) || n_sample_centers <= 0L) {
    stop("`n_sample_centers` must be a positive integer or Inf.")
  }
  set.seed(as.integer(seed))
  sort(sample.int(n, size = n_sample_centers, replace = FALSE))
}

sunspots_time_gof_prepare_centers <- function(data,
                                              theta,
                                              center_indices,
                                              l_max = 100L,
                                              quad_n = 400L) {
  coefficients <- sunspots_time_gof_average_legendre_coefficients(
    theta = theta, u_reference = data$u, l_max = l_max, quad_n = quad_n
  )
  n <- nrow(data$x)
  centers <- lapply(center_indices, function(center_index) {
    omega <- data$x[center_index, ]
    distances <- acos(pmin(pmax(as.numeric(data$x %*% omega), -1), 1))
    order_index <- order(distances)
    sorted_distances <- distances[order_index]
    tie_end_positions <- sorted_tie_end_positions(sorted_distances)
    empirical <- tie_end_positions / n
    theoretical <- sunspots_time_gof_profile_row(omega, sorted_distances, coefficients)
    observed_process <- sqrt(n) * (empirical - theoretical)
    list(
      center_index = center_index,
      omega = omega,
      order_index = order_index,
      sorted_distances = sorted_distances,
      tie_end_positions = tie_end_positions,
      ks_observed_statistic = max(abs(observed_process)),
      cvm_observed_sum = sum(observed_process^2)
    )
  })
  n_centers <- length(centers)
  list(
    coefficients = coefficients,
    centers = centers,
    ks_statistic = max(vapply(centers, `[[`, numeric(1L), "ks_observed_statistic")),
    cvm_statistic = sum(vapply(centers, `[[`, numeric(1L), "cvm_observed_sum")) / (n_centers * n),
    # Retained for callers that used the original KS-only helper.
    statistic = max(vapply(centers, `[[`, numeric(1L), "ks_observed_statistic"))
  )
}

sunspots_time_gof_prepare_fast_corrections <- function(data,
                                                       theta_hat,
                                                       centers,
                                                       derivative_mc_size,
                                                       seed,
                                                       control = list()) {
  par_hat <- sunspots_time_varying_pack_theta(
    theta_hat,
    nu_eps = as.numeric(control$nu_eps %||% 1e-6),
    c_min = as.numeric(control$c_min %||% 1e-8),
    c_max = as.numeric(control$c_max %||% 1e6),
    hemisphere_regression = control$hemisphere_regression %||% "asymmetric"
  )
  set.seed(as.integer(seed))
  auxiliary <- sample_sunspots_time_gof_conditional(
    n = derivative_mc_size, u_reference = data$u, par = par_hat, control = control
  )
  score_observed <- sunspots_time_gof_score_matrix(data, par_hat, control)
  score_auxiliary <- sunspots_time_gof_score_matrix(auxiliary, par_hat, control)
  vhat <- crossprod(score_auxiliary) / nrow(score_auxiliary)
  diagnostics <- fast_multiplier_vhat_diagnostics(
    S_obs = score_observed, Psi_aux = score_auxiliary, Vhat = vhat, par0 = par_hat
  )
  fast_multiplier_validate_vhat(
    Vhat = vhat,
    diagnostics = diagnostics,
    label = "conditional sunspots fast multiplier preparation",
    rcond_tol = as.numeric(control$fast_multiplier_vhat_rcond_tol %||% 1e-12)
  )
  auxiliary_scores_solved <- score_auxiliary %*% t(fast_multiplier_solve_vhat(
    vhat, diag(ncol(vhat)), label = "the conditional sunspots score correction"
  ))

  for (i in seq_along(centers)) {
    auxiliary_distances <- acos(pmin(pmax(as.numeric(auxiliary$x %*% centers[[i]]$omega), -1), 1))
    auxiliary_order <- order(auxiliary_distances)
    auxiliary_sorted <- auxiliary_distances[auxiliary_order]
    cumulative_scores <- col_cumsums_base(auxiliary_scores_solved[auxiliary_order, , drop = FALSE]) /
      nrow(auxiliary_scores_solved)
    score_basis <- rbind(rep(0, ncol(cumulative_scores)), cumulative_scores)
    selected_counts <- findInterval(centers[[i]]$sorted_distances, auxiliary_sorted)
    centers[[i]]$correction <- score_basis[selected_counts + 1L, , drop = FALSE]
  }

  list(
    centers = centers,
    score_observed = score_observed,
    vhat = vhat,
    diagnostics = diagnostics,
    derivative_mc_size = nrow(auxiliary$x)
  )
}

sunspots_time_gof_fast_statistics <- function(score_observed,
                                              centers,
                                              statistics = c("ks", "cvm"),
                                              B,
                                              seed,
                                              n_cores = 1L,
                                              bootstrap_block_size = 25L) {
  statistics <- unique(tolower(as.character(statistics)))
  if (length(statistics) == 0L || any(!statistics %in% c("ks", "cvm"))) {
    stop("`statistics` must contain one or both of 'ks' and 'cvm'.")
  }
  want_ks <- "ks" %in% statistics
  want_cvm <- "cvm" %in% statistics
  n <- nrow(score_observed)
  n_centers <- length(centers)
  if (n_centers == 0L) stop("`centers` must not be empty.")
  multiplier_spec <- resolve_multiplier_spec(NULL)
  raw_weights <- generate_multiplier_matrix(B = as.integer(B), n = n, multiplier_spec = multiplier_spec,
                                             seed = as.integer(seed))
  centered_weights <- raw_weights / rowMeans(raw_weights) - 1
  scale_factor <- multiplier_spec$mean / multiplier_spec$sd
  score_sums <- centered_weights %*% score_observed
  replicate_blocks <- split(seq_len(nrow(centered_weights)),
                            ceiling(seq_len(nrow(centered_weights)) / as.integer(bootstrap_block_size)))

  one_block <- function(indices) {
    maxima <- if (want_ks) numeric(length(indices)) else NULL
    cvm_sums <- if (want_cvm) numeric(length(indices)) else NULL
    for (center in centers) {
      ordered_weights <- centered_weights[indices, center$order_index, drop = FALSE]
      empirical <- row_cumsums_base(ordered_weights)
      empirical <- empirical[, center$tie_end_positions, drop = FALSE]
      correction <- score_sums[indices, , drop = FALSE] %*% t(center$correction)
      process <- scale_factor * (empirical - correction) / sqrt(n)
      if (want_ks) maxima <- pmax(maxima, apply(abs(process), 1L, max))
      if (want_cvm) cvm_sums <- cvm_sums + rowSums(process^2)
    }
    list(indices = indices, ks = maxima, cvm_sum = cvm_sums)
  }

  n_cores <- min(max(1L, as.integer(n_cores)), length(replicate_blocks))
  block_results <- if (n_cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(replicate_blocks, one_block, mc.cores = n_cores, mc.preschedule = TRUE)
  } else {
    lapply(replicate_blocks, one_block)
  }
  output <- list()
  if (want_ks) output$ks <- numeric(B)
  if (want_cvm) output$cvm <- numeric(B)
  for (result in block_results) {
    if (want_ks) output$ks[result$indices] <- result$ks
    if (want_cvm) output$cvm[result$indices] <- result$cvm_sum / (n_centers * n)
  }
  output
}

sunspots_time_gof_fast_ks_statistics <- function(score_observed,
                                                 centers,
                                                 B,
                                                 seed,
                                                 n_cores = 1L,
                                                 bootstrap_block_size = 25L) {
  sunspots_time_gof_fast_statistics(
    score_observed = score_observed,
    centers = centers,
    statistics = "ks",
    B = B,
    seed = seed,
    n_cores = n_cores,
    bootstrap_block_size = bootstrap_block_size
  )$ks
}

run_sunspots_cycle23_time_varying_asymmetric_mixture_gof <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle23_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle23_time_varying_asymmetric_mixture_gof_sampleks_fast"),
    start_date = "1997-06-01",
    end_date = "2006-01-01",
    statistics = "ks",
    hemisphere_regression = "asymmetric",
    B = 1000L,
    n_cores = 6L,
    seed = 20260711L,
    n_sample_centers = 100L,
    derivative_mc_size = 2000L,
    bootstrap_block_size = 25L,
    control = list()) {
  statistics <- unique(tolower(as.character(statistics)))
  if (length(statistics) == 0L || any(!statistics %in% c("ks", "cvm"))) {
    stop("`statistics` must contain one or both of 'ks' and 'cvm'.")
  }
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  control <- utils::modifyList(control, list(hemisphere_regression = hemisphere_regression))
  total_start <- proc.time()[["elapsed"]]
  retained <- prepare_sunspots_cycle23_time_varying_data(input_csv, start_date, end_date)
  data <- normalize_sunspots_time_gof_data(
    x = as.matrix(retained[, c("x1", "x2", "x3")]), u = retained$u
  )
  mle_start <- proc.time()[["elapsed"]]
  theta_hat <- fit_sunspots_time_varying_asymmetric_mixture(
    data$x, data$u, hemisphere_regression = hemisphere_regression, control = control
  )
  mle_seconds <- proc.time()[["elapsed"]] - mle_start
  center_indices <- sunspots_time_gof_select_centers(
    nrow(data$x), n_sample_centers = n_sample_centers, seed = seed
  )
  profile_l_max <- as.integer(control$profile_l_max %||% 100L)
  profile_quad_n <- as.integer(control$profile_quad_n %||% 400L)
  observed_start <- proc.time()[["elapsed"]]
  observed <- sunspots_time_gof_prepare_centers(
    data = data, theta = theta_hat, center_indices = center_indices,
    l_max = profile_l_max, quad_n = profile_quad_n
  )
  observed_seconds <- proc.time()[["elapsed"]] - observed_start
  fast_prep_start <- proc.time()[["elapsed"]]
  fast_prep <- sunspots_time_gof_prepare_fast_corrections(
    data = data, theta_hat = theta_hat, centers = observed$centers,
    derivative_mc_size = as.integer(derivative_mc_size), seed = as.integer(seed) + 1L,
    control = control
  )
  fast_prep_seconds <- proc.time()[["elapsed"]] - fast_prep_start
  bootstrap_start <- proc.time()[["elapsed"]]
  bootstrap_statistics <- sunspots_time_gof_fast_statistics(
    score_observed = fast_prep$score_observed,
    centers = fast_prep$centers,
    statistics = statistics,
    B = as.integer(B), seed = as.integer(seed) + 2L, n_cores = n_cores,
    bootstrap_block_size = bootstrap_block_size
  )
  bootstrap_seconds <- proc.time()[["elapsed"]] - bootstrap_start
  total_seconds <- proc.time()[["elapsed"]] - total_start
  observed_statistics <- c(ks = observed$ks_statistic, cvm = observed$cvm_statistic)[statistics]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- do.call(rbind, lapply(statistics, function(statistic_name) {
    statistic_values <- bootstrap_statistics[[statistic_name]]
    data.frame(
      model = sprintf("time_varying_%s_small_circle_mixture", hemisphere_regression),
      hemisphere_regression = hemisphere_regression,
      n_parameters = theta_hat$n_parameters,
      statistic_type = if (identical(statistic_name, "ks")) {
        "ks_sample_centers_unique_distances_global"
      } else {
        "cvm_sample_centers_unique_distances_global"
      },
      bootstrap_method = "fast_multiplier",
      n = nrow(data$x), n_sample_centers = length(center_indices), B = as.integer(B),
      derivative_mc_size = fast_prep$derivative_mc_size, n_cores = as.integer(n_cores),
      start_date = start_date, end_date_exclusive = end_date,
      profile_scope = "global_average_over_fixed_time_design",
      statistic = observed_statistics[[statistic_name]],
      critical_value_0.95 = as.numeric(stats::quantile(statistic_values, probs = 0.95, names = FALSE, type = 8)),
      p_value = (1 + sum(statistic_values >= observed_statistics[[statistic_name]])) / (length(statistic_values) + 1),
      a_N = theta_hat$a_N, b_N = theta_hat$b_N, a_S = theta_hat$a_S, b_S = theta_hat$b_S,
      c = theta_hat$c, loglik = theta_hat$loglik, convergence = theta_hat$opt$convergence,
      vhat_rcond = fast_prep$diagnostics$Vhat_rcond,
      mle_seconds = mle_seconds, observed_profile_seconds = observed_seconds,
      fast_preparation_seconds = fast_prep_seconds, bootstrap_seconds = bootstrap_seconds,
      total_seconds = total_seconds,
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(summary, file.path(output_dir, "cycle23_time_varying_asymmetric_samplegof_fast_results.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    center_rank = seq_along(center_indices), center_index = center_indices,
    date = retained$date[center_indices], u = retained$u[center_indices],
    x1 = data$x[center_indices, 1L], x2 = data$x[center_indices, 2L], x3 = data$x[center_indices, 3L]
  ), file.path(output_dir, "cycle23_time_varying_asymmetric_samplegof_centers.csv"), row.names = FALSE)
  for (statistic_name in statistics) {
    utils::write.csv(
      data.frame(replicate = seq_along(bootstrap_statistics[[statistic_name]]), statistic = bootstrap_statistics[[statistic_name]]),
      file.path(output_dir, sprintf("cycle23_time_varying_asymmetric_samplegof_fast_bootstrap_%s.csv", statistic_name)),
      row.names = FALSE
    )
  }
  utils::write.csv(data.frame(
    mle_seconds = mle_seconds,
    observed_profile_seconds = observed_seconds,
    fast_preparation_seconds = fast_prep_seconds,
    bootstrap_seconds = bootstrap_seconds,
    total_seconds = total_seconds,
    stringsAsFactors = FALSE
  ), file.path(output_dir, "cycle23_time_varying_asymmetric_samplegof_fast_timing.csv"), row.names = FALSE)
  if (identical(statistics, "ks")) {
    # Preserve the filenames emitted by the original KS-only runner.
    utils::write.csv(summary, file.path(output_dir, "cycle23_time_varying_asymmetric_sampleks_fast_results.csv"), row.names = FALSE)
    utils::write.csv(data.frame(
      center_rank = seq_along(center_indices), center_index = center_indices,
      date = retained$date[center_indices], u = retained$u[center_indices],
      x1 = data$x[center_indices, 1L], x2 = data$x[center_indices, 2L], x3 = data$x[center_indices, 3L]
    ), file.path(output_dir, "cycle23_time_varying_asymmetric_sampleks_centers.csv"), row.names = FALSE)
    utils::write.csv(data.frame(replicate = seq_along(bootstrap_statistics$ks), statistic = bootstrap_statistics$ks),
                     file.path(output_dir, "cycle23_time_varying_asymmetric_sampleks_fast_bootstrap.csv"), row.names = FALSE)
  }
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  invisible(list(summary = summary, theta_hat = theta_hat, observed_statistics = observed_statistics,
                 bootstrap_statistics = bootstrap_statistics, centers = center_indices,
                 timing = list(mle_seconds = mle_seconds, observed_profile_seconds = observed_seconds,
                               fast_preparation_seconds = fast_prep_seconds,
                               bootstrap_seconds = bootstrap_seconds, total_seconds = total_seconds),
                 output_dir = output_dir))
}

parse_sunspots_time_gof_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  out <- list()
  integer_keys <- c("B", "n_cores", "seed", "derivative_mc_size", "bootstrap_block_size")
  character_keys <- c("input_csv", "output_dir", "start_date", "end_date", "hemisphere_regression")
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat("Options: --statistics=ks,cvm --hemisphere_regression=asymmetric|shared --B=INTEGER --n_cores=INTEGER --n_sample_centers=INTEGER|all --derivative_mc_size=INTEGER --bootstrap_block_size=INTEGER --start_date=YYYY-MM-DD --end_date=YYYY-MM-DD --output_dir=PATH\n")
      quit(save = "no", status = 0L)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop(sprintf("Invalid option: %s", arg))
    if (parts[[1L]] %in% integer_keys) out[[parts[[1L]]]] <- as.integer(parts[[2L]])
    if (parts[[1L]] == "n_sample_centers") {
      out[[parts[[1L]]]] <- if (tolower(parts[[2L]]) %in% c("all", "inf", "infinity")) Inf else as.integer(parts[[2L]])
    }
    if (parts[[1L]] == "statistics") {
      out[[parts[[1L]]]] <- strsplit(tolower(parts[[2L]]), ",", fixed = TRUE)[[1L]]
    }
    if (parts[[1L]] %in% character_keys) out[[parts[[1L]]]] <- parts[[2L]]
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(run_sunspots_cycle23_time_varying_asymmetric_mixture_gof, parse_sunspots_time_gof_args())
}
