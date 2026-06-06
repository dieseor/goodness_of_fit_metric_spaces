#!/usr/bin/env Rscript

resolve_cycle24_weighted_projection_path <- function(...) {
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

utils_path_cycle24_weighted_projection <- resolve_cycle24_weighted_projection_path("utils.R")
source(utils_path_cycle24_weighted_projection)

sunspots_projected_gof_cycle24_weighted <- function(z, cdf_fun, grid_size = 1001L) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  if (length(z) == 0L) {
    stop("`z` must contain at least one finite projected value.")
  }

  ecdf_z <- stats::ecdf(z)
  z_grid <- sort(unique(c(seq(-1, 1, length.out = as.integer(grid_size)), z)))
  fitted_cdf <- as.numeric(cdf_fun(z_grid))
  empirical_cdf <- as.numeric(ecdf_z(z_grid))

  list(
    z_grid = z_grid,
    fitted_cdf = fitted_cdf,
    empirical_cdf = empirical_cdf,
    ks = max(abs(empirical_cdf - fitted_cdf)),
    cvm = mean((empirical_cdf - fitted_cdf)^2)
  )
}

build_projected_component_cdf <- function(z_grid,
                                          omega,
                                          mu,
                                          kappa,
                                          nu,
                                          method = "legendre",
                                          l_max = 200L,
                                          quad_n = 400L,
                                          tol = 1e-10) {
  z_grid <- pmin(pmax(as.numeric(z_grid), -1), 1)
  t_grid <- acos(z_grid)
  1 - distance_profile_small_circle(
    omega = omega,
    t_values = t_grid,
    mu = mu,
    kappa = kappa,
    nu = nu,
    distance_type = "geodesic",
    method = method,
    l_max = as.integer(l_max),
    quad_n = as.integer(quad_n),
    tol = as.numeric(tol)
  )
}

validate_projected_cdf <- function(name, z_grid, fitted_cdf, clip_eps = 1e-10) {
  warnings <- character()
  if (length(fitted_cdf) != length(z_grid)) {
    stop(sprintf("Length mismatch in %s projected CDF.", name))
  }
  if (any(!is.finite(fitted_cdf))) {
    stop(sprintf("%s projected CDF contains non-finite values.", name))
  }

  cdf_checked <- fitted_cdf
  if (any(cdf_checked < -clip_eps | cdf_checked > 1 + clip_eps)) {
    warnings <- c(
      warnings,
      sprintf(
        "%s projected CDF left [0,1]: min=%.16f, max=%.16f",
        name,
        min(cdf_checked),
        max(cdf_checked)
      )
    )
  }
  if (any(cdf_checked < 0 | cdf_checked > 1)) {
    warnings <- c(
      warnings,
      sprintf("%s projected CDF clipped into [0,1].", name)
    )
    cdf_checked <- pmin(pmax(cdf_checked, 0), 1)
  }

  diffs <- diff(cdf_checked)
  if (any(diffs < -clip_eps)) {
    warnings <- c(
      warnings,
      sprintf(
        "%s projected CDF is not monotone before isotonic repair: minimum increment %.16e",
        name,
        min(diffs)
      )
    )
    cdf_checked <- stats::isoreg(z_grid, cdf_checked)$yf
    cdf_checked <- pmin(pmax(cdf_checked, 0), 1)
  }

  if (cdf_checked[[1L]] > 1e-4) {
    warnings <- c(
      warnings,
      sprintf("%s projected CDF at z=-1 is %.16f, not close to 0.", name, cdf_checked[[1L]])
    )
  }
  if (cdf_checked[[length(cdf_checked)]] < 1 - 1e-4) {
    warnings <- c(
      warnings,
      sprintf("%s projected CDF at z=1 is %.16f, not close to 1.", name, cdf_checked[[length(cdf_checked)]])
    )
  }

  list(cdf = cdf_checked, warnings = warnings)
}

build_projected_data_frame <- function(z_values) {
  z_sorted <- sort(as.numeric(z_values))
  n <- length(z_sorted)
  data.frame(
    z_observed = z_sorted,
    ecdf_observed = seq_len(n) / n,
    stringsAsFactors = FALSE
  )
}

add_simulated_projection_columns <- function(data_df, z_simulated) {
  z_sim_sorted <- sort(as.numeric(z_simulated))
  n_sim <- length(z_sim_sorted)
  data_df$z_simulated <- NA_real_
  data_df$ecdf_simulated <- NA_real_
  data_df$z_simulated[seq_len(n_sim)] <- z_sim_sorted
  data_df$ecdf_simulated[seq_len(n_sim)] <- seq_len(n_sim) / n_sim
  data_df
}

plot_cycle24_weighted_mixture_projected_cdf <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle24_weighted_mixture_projected_cdf"),
    include_simulated_ecdf = TRUE,
    simulation_seed = 20260604L,
    simulation_size = NULL,
    grid_size = 1001L,
    l_max = 200L,
    quad_n = 400L,
    tol = 1e-10) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  weighted_theta <- small_circle_weighted_mixture2_normalize_theta(list(
    pi = 0.529731,
    mu = c(-0.001131, -0.004108, 0.999991),
    kappa1 = 26.806931,
    nu1 = 0.237492,
    kappa2 = 24.109693,
    nu2 = 0.268405
  ))
  symmetric_theta <- small_circle_symmetric_mixture2_normalize_theta(list(
    mu = c(0.002746, 0.002837, -0.999992),
    kappa = 25.045825,
    nu = 0.251825
  ))

  sunspots_df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  x <- as.matrix(sunspots_df[, c("x1", "x2", "x3")])
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  north_pole <- c(0, 0, 1)
  z_observed <- pmin(pmax(as.numeric(x %*% north_pole), -1), 1)
  n <- nrow(x)
  if (is.null(simulation_size)) {
    simulation_size <- n
  }

  warning_log <- character()

  weighted_cdf_fun <- function(z_grid) {
    weighted_theta$pi * build_projected_component_cdf(
      z_grid = z_grid,
      omega = north_pole,
      mu = weighted_theta$mu,
      kappa = weighted_theta$kappa1,
      nu = weighted_theta$nu1,
      method = "legendre",
      l_max = l_max,
      quad_n = quad_n,
      tol = tol
    ) + (1 - weighted_theta$pi) * build_projected_component_cdf(
      z_grid = z_grid,
      omega = north_pole,
      mu = -weighted_theta$mu,
      kappa = weighted_theta$kappa2,
      nu = weighted_theta$nu2,
      method = "legendre",
      l_max = l_max,
      quad_n = quad_n,
      tol = tol
    )
  }

  symmetric_cdf_fun <- function(z_grid) {
    1 - distance_profile_small_circle_symmetric_mixture2(
      omega = north_pole,
      t_values = acos(pmin(pmax(as.numeric(z_grid), -1), 1)),
      mu = symmetric_theta$mu,
      kappa = symmetric_theta$kappa,
      nu = symmetric_theta$nu,
      distance_type = "geodesic",
      method = "legendre",
      l_max = as.integer(l_max),
      quad_n = as.integer(quad_n),
      tol = as.numeric(tol)
    )
  }

  projected_fit <- sunspots_projected_gof_cycle24_weighted(
    z = z_observed,
    cdf_fun = weighted_cdf_fun,
    grid_size = grid_size
  )
  projected_fit_symmetric <- sunspots_projected_gof_cycle24_weighted(
    z = z_observed,
    cdf_fun = symmetric_cdf_fun,
    grid_size = grid_size
  )

  weighted_validation <- validate_projected_cdf(
    name = "weighted mixture",
    z_grid = projected_fit$z_grid,
    fitted_cdf = projected_fit$fitted_cdf
  )
  symmetric_validation <- validate_projected_cdf(
    name = "symmetric mixture",
    z_grid = projected_fit_symmetric$z_grid,
    fitted_cdf = projected_fit_symmetric$fitted_cdf
  )
  projected_fit$fitted_cdf <- weighted_validation$cdf
  projected_fit_symmetric$fitted_cdf <- symmetric_validation$cdf
  warning_log <- c(warning_log, weighted_validation$warnings, symmetric_validation$warnings)

  cdf_grid_df <- data.frame(
    z = projected_fit$z_grid,
    fitted_cdf = projected_fit$fitted_cdf,
    observed_ecdf = projected_fit$empirical_cdf,
    diff_observed_minus_fitted = projected_fit$empirical_cdf - projected_fit$fitted_cdf,
    abs_diff = abs(projected_fit$empirical_cdf - projected_fit$fitted_cdf),
    symmetric_fitted_cdf = projected_fit_symmetric$fitted_cdf,
    diff_observed_minus_symmetric = projected_fit$empirical_cdf - projected_fit_symmetric$fitted_cdf,
    abs_diff_symmetric = abs(projected_fit$empirical_cdf - projected_fit_symmetric$fitted_cdf),
    stringsAsFactors = FALSE
  )

  max_idx <- which.max(cdf_grid_df$abs_diff)
  z_projected_df <- build_projected_data_frame(z_observed)

  z_simulated <- NULL
  if (isTRUE(include_simulated_ecdf)) {
    set.seed(as.integer(simulation_seed))
    x_simulated <- r_sph_small_circle_weighted_mixture2(
      n = as.integer(simulation_size),
      mu = weighted_theta$mu,
      pi = weighted_theta$pi,
      kappa1 = weighted_theta$kappa1,
      nu1 = weighted_theta$nu1,
      kappa2 = weighted_theta$kappa2,
      nu2 = weighted_theta$nu2
    )
    z_simulated <- pmin(pmax(as.numeric(x_simulated %*% north_pole), -1), 1)
    z_projected_df <- add_simulated_projection_columns(z_projected_df, z_simulated)

    simulated_fit <- sunspots_projected_gof_cycle24_weighted(
      z = z_simulated,
      cdf_fun = weighted_cdf_fun,
      grid_size = grid_size
    )
    sim_validation <- validate_projected_cdf(
      name = "weighted mixture (simulation check)",
      z_grid = simulated_fit$z_grid,
      fitted_cdf = simulated_fit$fitted_cdf
    )
    warning_log <- c(warning_log, sim_validation$warnings)
  }

  equator_bands <- c(0.05, 0.075, 0.10, 0.15, 0.20, 0.25)
  observed_equator_mass <- vapply(equator_bands, function(eps) {
    mean(abs(z_observed) <= eps)
  }, numeric(1L))
  fitted_equator_mass <- vapply(equator_bands, function(eps) {
    weighted_cdf_fun(eps) - weighted_cdf_fun(-eps)
  }, numeric(1L))

  summary_df <- data.frame(
    n = n,
    pi_hat = weighted_theta$pi,
    mu_1 = weighted_theta$mu[[1L]],
    mu_2 = weighted_theta$mu[[2L]],
    mu_3 = weighted_theta$mu[[3L]],
    kappa1_hat = weighted_theta$kappa1,
    nu1_hat = weighted_theta$nu1,
    kappa2_hat = weighted_theta$kappa2,
    nu2_hat = weighted_theta$nu2,
    max_abs_diff_projected = cdf_grid_df$abs_diff[[max_idx]],
    where_max_diff_z = cdf_grid_df$z[[max_idx]],
    observed_P_z_positive = mean(z_observed > 0),
    fitted_P_z_positive = 1 - weighted_cdf_fun(0),
    symmetric_P_z_positive = 1 - symmetric_cdf_fun(0),
    weighted_minus_symmetric_max_abs_diff = max(abs(cdf_grid_df$fitted_cdf - cdf_grid_df$symmetric_fitted_cdf)),
    stringsAsFactors = FALSE
  )

  for (j in seq_along(equator_bands)) {
    suffix <- gsub("\\.", "p", formatC(equator_bands[[j]], format = "f", digits = 3))
    summary_df[[paste0("observed_mass_near_equator_le_", suffix)]] <- observed_equator_mass[[j]]
    summary_df[[paste0("fitted_mass_near_equator_le_", suffix)]] <- fitted_equator_mass[[j]]
    summary_df[[paste0("observed_minus_fitted_mass_near_equator_le_", suffix)]] <-
      observed_equator_mass[[j]] - fitted_equator_mass[[j]]
  }

  png_path <- file.path(output_dir, "cycle24_weighted_mixture_north_pole_projected_cdf.png")
  grDevices::png(png_path, width = 1400, height = 1000, res = 160)
  plot(
    NA,
    NA,
    xlim = c(-1, 1),
    ylim = c(0, 1),
    xlab = "z = e3' x = x3",
    ylab = "CDF",
    main = "Cycle 24 sunspots: north-pole projected CDF"
  )
  z_emp <- sort(unique(c(-1, z_observed, 1)))
  lines(z_emp, stats::ecdf(z_observed)(z_emp), type = "s", lwd = 2.2, col = "black")
  lines(cdf_grid_df$z, cdf_grid_df$fitted_cdf, lwd = 2.4, col = "#2166ac")
  lines(cdf_grid_df$z, cdf_grid_df$symmetric_fitted_cdf, lwd = 2.2, lty = 2, col = "#b2182b")
  if (!is.null(z_simulated)) {
    z_sim_emp <- sort(unique(c(-1, z_simulated, 1)))
    lines(z_sim_emp, stats::ecdf(z_simulated)(z_sim_emp), type = "s", lwd = 1.4, lty = 3, col = "#1b9e77")
  }
  grid(col = "#d9d9d9")
  legend(
    "topleft",
    legend = c(
      "Observed ECDF",
      "Fitted weighted mixture",
      "Fitted symmetric mixture",
      if (!is.null(z_simulated)) "Simulated ECDF" else NULL
    ),
    col = c("black", "#2166ac", "#b2182b", if (!is.null(z_simulated)) "#1b9e77" else NULL),
    lwd = c(2.2, 2.4, 2.2, if (!is.null(z_simulated)) 1.4 else NULL),
    lty = c(1, 1, 2, if (!is.null(z_simulated)) 3 else NULL),
    bty = "n"
  )
  grDevices::dev.off()

  utils::write.csv(
    cdf_grid_df[, c("z", "fitted_cdf", "observed_ecdf", "diff_observed_minus_fitted", "abs_diff")],
    file = file.path(output_dir, "sunspots_cycle24_weighted_mixture_north_pole_projected_cdf_grid.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    z_projected_df,
    file = file.path(output_dir, "sunspots_cycle24_weighted_mixture_north_pole_projected_data.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary_df,
    file = file.path(output_dir, "sunspots_cycle24_weighted_mixture_projected_summary.csv"),
    row.names = FALSE
  )

  warning_log <- unique(warning_log[nzchar(warning_log)])
  warning_path <- file.path(output_dir, "sunspots_cycle24_weighted_mixture_projected_warnings.txt")
  if (length(warning_log) == 0L) {
    writeLines("No warnings.", con = warning_path)
  } else {
    writeLines(warning_log, con = warning_path)
  }

  invisible(list(
    cdf_grid = cdf_grid_df,
    projected_data = z_projected_df,
    summary = summary_df,
    warnings = warning_log,
    png_path = png_path
  ))
}

if (sys.nframe() == 0L) {
  plot_cycle24_weighted_mixture_projected_cdf()
}
