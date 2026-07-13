#!/usr/bin/env Rscript

resolve_cycle24_weighted_worst_ks_path <- function(...) {
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

utils_path_cycle24_weighted_worst_ks <- resolve_cycle24_weighted_worst_ks_path("utils.R")
bootstrap_path_cycle24_weighted_worst_ks <- resolve_cycle24_weighted_worst_ks_path("bootstrap", "multiplier_bootstrap.R")
spec_path_cycle24_weighted_worst_ks <- resolve_cycle24_weighted_worst_ks_path("bootstrap", "small_circle_weighted_mixture2_model_spec.R")

source(utils_path_cycle24_weighted_worst_ks)
source(spec_path_cycle24_weighted_worst_ks)
source(bootstrap_path_cycle24_weighted_worst_ks)

build_component_projected_cdf_cycle24_weighted <- function(a_values,
                                                           omega,
                                                           mu,
                                                           kappa,
                                                           nu,
                                                           l_max = 200L,
                                                           quad_n = 400L,
                                                           tol = 1e-10) {
  a_values <- pmin(pmax(as.numeric(a_values), -1), 1)
  t_values <- acos(a_values)
  1 - distance_profile_small_circle(
    omega = omega,
    t_values = t_values,
    mu = mu,
    kappa = kappa,
    nu = nu,
    distance_type = "geodesic",
    method = "legendre",
    l_max = as.integer(l_max),
    quad_n = as.integer(quad_n),
    tol = as.numeric(tol)
  )
}

build_weighted_projected_cdf_cycle24 <- function(a_values,
                                                 omega,
                                                 theta,
                                                 l_max = 200L,
                                                 quad_n = 400L,
                                                 tol = 1e-10) {
  theta$pi * build_component_projected_cdf_cycle24_weighted(
    a_values = a_values,
    omega = omega,
    mu = theta$mu,
    kappa = theta$kappa1,
    nu = theta$nu1,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  ) + (1 - theta$pi) * build_component_projected_cdf_cycle24_weighted(
    a_values = a_values,
    omega = omega,
    mu = -theta$mu,
    kappa = theta$kappa2,
    nu = theta$nu2,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )
}

validate_projected_cdf_cycle24_weighted <- function(name, grid_values, fitted_cdf, eps = 1e-10) {
  warning_log <- character()
  if (any(!is.finite(fitted_cdf))) {
    stop(sprintf("%s projected CDF contains non-finite values.", name))
  }

  cdf_checked <- as.numeric(fitted_cdf)
  if (any(cdf_checked < -eps | cdf_checked > 1 + eps)) {
    warning_log <- c(
      warning_log,
      sprintf(
        "%s projected CDF left [0,1]: min=%.16f max=%.16f",
        name,
        min(cdf_checked),
        max(cdf_checked)
      )
    )
  }
  if (any(cdf_checked < 0 | cdf_checked > 1)) {
    warning_log <- c(warning_log, sprintf("%s projected CDF clipped into [0,1].", name))
    cdf_checked <- pmin(pmax(cdf_checked, 0), 1)
  }
  if (any(diff(cdf_checked) < -eps)) {
    warning_log <- c(
      warning_log,
      sprintf(
        "%s projected CDF is not monotone before isotonic repair: minimum increment %.16e",
        name,
        min(diff(cdf_checked))
      )
    )
    cdf_checked <- stats::isoreg(grid_values, cdf_checked)$yf
    cdf_checked <- pmin(pmax(cdf_checked, 0), 1)
  }

  list(cdf = cdf_checked, warnings = warning_log)
}

plot_cycle24_weighted_mixture_worst_ks_projection <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle24_weighted_mixture_worst_ks_projection"),
    include_simulated_ecdf = TRUE,
    simulation_seed = 20260604L,
    simulation_size = NULL,
    M_value = 60L,
    ks_t_points = 200L,
    grid_size = 1001L,
    l_max = 200L,
    quad_n = 400L,
    tol = 1e-10) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  theta_hat <- small_circle_weighted_mixture2_normalize_theta(list(
    pi = 0.529731,
    mu = c(-0.001131, -0.004108, 0.999991),
    kappa1 = 26.806931,
    nu1 = 0.237492,
    kappa2 = 24.109693,
    nu2 = 0.268405
  ))

  sunspots_df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  x <- as.matrix(sunspots_df[, c("x1", "x2", "x3")])
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  n <- nrow(x)
  if (is.null(simulation_size)) {
    simulation_size <- n
  }

  spec <- make_small_circle_weighted_mixture2_spec(distance_type = "geodesic")
  ks_grid <- make_sample_unique_distance_ks_grid()
  ks_prep <- prepare_ks_observed_data(
    data = x,
    spec = spec,
    theta_hat = theta_hat,
    ks_grid = ks_grid,
    control = list(
      small_circle_weighted_mixture2_profile_method = "legendre",
      small_circle_weighted_mixture2_L_max = as.integer(l_max),
      small_circle_weighted_mixture2_quad_n = as.integer(quad_n),
      small_circle_weighted_mixture2_tol = as.numeric(tol)
    )
  )

  process_abs <- abs(ks_prep$process_matrix)
  max_index <- which.max(process_abs)
  dims <- dim(process_abs)
  omega_index <- ((max_index - 1L) %% dims[[1L]]) + 1L
  t_index <- ((max_index - 1L) %/% dims[[1L]]) + 1L

  omega_star <- as.numeric(grid_point_at(ks_prep$omega_grid, omega_index))
  t_star <- as.numeric(ks_prep$t_grid[[t_index]])
  dot_threshold <- cos(t_star)
  observed_profile_at_star <- as.numeric(ks_prep$empirical_profile[omega_index, t_index])
  theoretical_profile_at_star <- as.numeric(ks_prep$theoretical_profile[omega_index, t_index])
  observed_ecdf_at_threshold <- 1 - observed_profile_at_star
  fitted_cdf_at_threshold <- 1 - theoretical_profile_at_star
  diff_observed_minus_fitted <- observed_ecdf_at_threshold - fitted_cdf_at_threshold
  abs_diff <- abs(diff_observed_minus_fitted)
  sqrt_n_abs_diff <- sqrt(n) * abs_diff

  y_observed <- pmin(pmax(as.numeric(x %*% omega_star), -1), 1)
  ecdf_y <- stats::ecdf(y_observed)
  y_grid <- sort(unique(c(seq(-1, 1, length.out = as.integer(grid_size)), y_observed, dot_threshold)))
  fitted_cdf <- build_weighted_projected_cdf_cycle24(
    a_values = y_grid,
    omega = omega_star,
    theta = theta_hat,
    l_max = l_max,
    quad_n = quad_n,
    tol = tol
  )
  validation <- validate_projected_cdf_cycle24_weighted(
    name = "weighted mixture worst-KS projection",
    grid_values = y_grid,
    fitted_cdf = fitted_cdf
  )
  fitted_cdf <- validation$cdf
  observed_ecdf <- as.numeric(ecdf_y(y_grid))

  cdf_grid_df <- data.frame(
    y = y_grid,
    fitted_cdf = fitted_cdf,
    observed_ecdf = observed_ecdf,
    diff_observed_minus_fitted = observed_ecdf - fitted_cdf,
    abs_diff = abs(observed_ecdf - fitted_cdf),
    stringsAsFactors = FALSE
  )

  max_projected_idx <- which.max(cdf_grid_df$abs_diff)
  projected_data_df <- data.frame(
    y_observed = sort(y_observed),
    ecdf_observed = seq_len(n) / n,
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_simulated_ecdf)) {
    set.seed(as.integer(simulation_seed))
    x_simulated <- r_sph_small_circle_weighted_mixture2(
      n = as.integer(simulation_size),
      mu = theta_hat$mu,
      pi = theta_hat$pi,
      kappa1 = theta_hat$kappa1,
      nu1 = theta_hat$nu1,
      kappa2 = theta_hat$kappa2,
      nu2 = theta_hat$nu2
    )
    y_simulated <- pmin(pmax(as.numeric(x_simulated %*% omega_star), -1), 1)
    projected_data_df$y_simulated <- sort(y_simulated)
    projected_data_df$ecdf_simulated <- seq_len(simulation_size) / simulation_size
  }

  summary_df <- data.frame(
    n = n,
    omega_1 = omega_star[[1L]],
    omega_2 = omega_star[[2L]],
    omega_3 = omega_star[[3L]],
    t_star = t_star,
    dot_threshold = dot_threshold,
    observed_ecdf_at_threshold = observed_ecdf_at_threshold,
    fitted_cdf_at_threshold = fitted_cdf_at_threshold,
    diff_observed_minus_fitted = diff_observed_minus_fitted,
    abs_diff = abs_diff,
    sqrt_n_abs_diff = sqrt_n_abs_diff,
    max_abs_diff_projected = cdf_grid_df$abs_diff[[max_projected_idx]],
    where_max_diff_y = cdf_grid_df$y[[max_projected_idx]],
    stringsAsFactors = FALSE
  )

  png_path <- file.path(output_dir, "cycle24_weighted_mixture_worst_ks_projected_cdf.png")
  grDevices::png(png_path, width = 1400, height = 1000, res = 160)
  plot(
    NA,
    NA,
    xlim = c(-1, 1),
    ylim = c(0, 1),
    xlab = "y = omega*' x",
    ylab = "CDF",
    main = "Cycle 24 sunspots: projected CDF on the worst KS direction"
  )
  y_emp <- sort(unique(c(-1, y_observed, 1)))
  lines(y_emp, stats::ecdf(y_observed)(y_emp), type = "s", lwd = 2.2, col = "black")
  lines(cdf_grid_df$y, cdf_grid_df$fitted_cdf, lwd = 2.4, col = "#2166ac")
  abline(v = dot_threshold, lty = 2, lwd = 1.4, col = "#b2182b")
  if (isTRUE(include_simulated_ecdf)) {
    y_sim_emp <- sort(unique(c(-1, projected_data_df$y_simulated, 1)))
    lines(y_sim_emp, stats::ecdf(projected_data_df$y_simulated)(y_sim_emp), type = "s", lwd = 1.4, lty = 3, col = "#1b9e77")
  }
  grid(col = "#d9d9d9")
  legend(
    "topleft",
    legend = c(
      "Observed ECDF",
      "Fitted weighted mixture",
      "Worst-cap threshold",
      if (isTRUE(include_simulated_ecdf)) "Simulated ECDF" else NULL
    ),
    col = c("black", "#2166ac", "#b2182b", if (isTRUE(include_simulated_ecdf)) "#1b9e77" else NULL),
    lwd = c(2.2, 2.4, 1.4, if (isTRUE(include_simulated_ecdf)) 1.4 else NULL),
    lty = c(1, 1, 2, if (isTRUE(include_simulated_ecdf)) 3 else NULL),
    bty = "n"
  )
  grDevices::dev.off()

  utils::write.csv(
    cdf_grid_df,
    file = file.path(output_dir, "sunspots_cycle24_weighted_mixture_worst_ks_projected_cdf_grid.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    projected_data_df,
    file = file.path(output_dir, "sunspots_cycle24_weighted_mixture_worst_ks_projected_data.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary_df,
    file = file.path(output_dir, "sunspots_cycle24_weighted_mixture_worst_ks_summary.csv"),
    row.names = FALSE
  )
  warning_path <- file.path(output_dir, "sunspots_cycle24_weighted_mixture_worst_ks_warnings.txt")
  if (length(validation$warnings) == 0L) {
    writeLines("No warnings.", con = warning_path)
  } else {
    writeLines(validation$warnings, con = warning_path)
  }

  invisible(list(
    summary = summary_df,
    cdf_grid = cdf_grid_df,
    projected_data = projected_data_df,
    warnings = validation$warnings,
    png_path = png_path
  ))
}

if (sys.nframe() == 0L) {
  plot_cycle24_weighted_mixture_worst_ks_projection()
}
