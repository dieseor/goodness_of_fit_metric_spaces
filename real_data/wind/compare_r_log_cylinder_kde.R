#!/usr/bin/env Rscript

# Compare the reflected-r KDE with the Gaussian log(r) KDE against the same
# fitted HvMF density.  All discrepancies are integrated with respect to
# dtheta dz after transforming the reflected-r estimator by q_z(theta,z) =
# q_r(theta, exp(z)) exp(z).  Hence the numerical comparisons use one common
# dominating measure and are not artifacts of the coordinate change.

script_path <- function() {
  file_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_args) == 1L) {
    return(normalizePath(sub("^--file=", "", file_args), mustWork = TRUE))
  }
  normalizePath(file.path("real_data", "wind", "compare_r_log_cylinder_kde.R"),
                mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(script_path()), "..", ".."), mustWork = TRUE)
source(file.path(repo_root, "real_data", "wind", "plot_risoe_cylinder_density_contours.R"))

old_summary_path <- file.path(
  repo_root, "real_data", "wind", "cylinder_density_contours",
  "risoe_125m_cylinder_density_contours_summary.csv"
)
log_summary_path <- file.path(
  repo_root, "real_data", "wind", "log_cylinder_density_contours",
  "risoe_125m_log_cylinder_density_contours_summary.csv"
)
output_path <- file.path(
  repo_root, "real_data", "wind", "log_cylinder_density_contours",
  "risoe_125m_r_vs_log_kde_fit_comparison.csv"
)

log_reflected_gaussian_kernel <- function(r, center, bandwidth) {
  a <- stats::dnorm((r - center) / bandwidth, log = TRUE)
  b <- stats::dnorm((r + center) / bandwidth, log = TRUE)
  pmax(a, b) + log1p(exp(-abs(a - b))) - log(bandwidth)
}

old_reflected_r_kde_on_log_grid <- function(z_grid, theta_grid, case, h_r, kappa_theta) {
  r_grid <- exp(z_grid)
  density_r <- numeric(length(z_grid))
  for (j in seq_len(nrow(case))) {
    density_r <- density_r + exp(
      log_von_mises_kernel(theta_grid - case$theta[[j]], kappa_theta) +
        log_reflected_gaussian_kernel(r_grid, case$r[[j]], h_r)
    )
  }
  density_r / nrow(case) * r_grid
}

fit_from_summary <- function(row) {
  chi <- row$mu_r_hat_log[[1L]]
  theta <- row$mu_theta_deg_hat_log[[1L]] * pi / 180
  list(
    mu = c(cosh(chi), sinh(chi) * cos(theta), sinh(chi) * sin(theta)),
    kappa = row$kappa_hat_log[[1L]],
    sinh_chi = sinh(chi)
  )
}

integrated_discrepancies <- function(parametric, estimate, theta_values, z_values) {
  n_theta <- length(theta_values)
  n_z <- length(z_values)
  to_matrix <- function(x) t(matrix(x, nrow = n_theta, ncol = n_z))
  p <- to_matrix(parametric)
  q <- to_matrix(estimate)
  integrate <- function(x) trapezoid_2d_mass(x, theta_values, z_values)

  mass_p <- integrate(p)
  mass_q <- integrate(q)
  p_normalized <- p / mass_p
  q_normalized <- q / mass_q

  c(
    mass_parametric_on_grid = mass_p,
    mass_kde_on_grid = mass_q,
    l1 = integrate(abs(p_normalized - q_normalized)),
    total_variation = 0.5 * integrate(abs(p_normalized - q_normalized)),
    hellinger_squared = 0.5 * integrate((sqrt(p_normalized) - sqrt(q_normalized))^2),
    integrated_squared_error = integrate((p_normalized - q_normalized)^2)
  )
}

compare_estimators <- function(n_theta = 241L, n_z = 401L) {
  if (!file.exists(old_summary_path) || !file.exists(log_summary_path)) {
    stop("Both old-r and log-r summary CSV files are required.", call. = FALSE)
  }
  old_summary <- utils::read.csv(old_summary_path, stringsAsFactors = FALSE)
  log_summary <- utils::read.csv(log_summary_path, stringsAsFactors = FALSE)
  key <- c("window", "pattern")
  summaries <- merge(old_summary, log_summary, by = key, suffixes = c("_old", "_log"))
  if (nrow(summaries) != 12L) {
    stop("Expected 12 matched window-pattern cases in the summary CSV files.", call. = FALSE)
  }

  all_data <- load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"
  )
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")
  theta_values <- seq(0, 2 * pi, length.out = n_theta)
  z_values <- seq(-10, log(10), length.out = n_z)
  grid <- expand.grid(theta = theta_values, z = z_values)

  rows <- lapply(seq_len(nrow(summaries)), function(i) {
    summary_row <- summaries[i, , drop = FALSE]
    case <- build_case(selected, summary_row$window, summary_row$pattern)
    fit <- fit_from_summary(summary_row)

    if (abs(fit$kappa - summary_row$kappa_hat_log) > 1e-10 ||
        abs(summary_row$kappa_hat_old - summary_row$kappa_hat_log) > 1e-10) {
      stop("Parametric fits differ between the two summary CSV files.", call. = FALSE)
    }

    parametric <- hvmf_log_cylinder_density(grid$z, grid$theta, fit)
    old_kde <- old_reflected_r_kde_on_log_grid(
      grid$z, grid$theta, case,
      h_r = summary_row$kde_h_r,
      kappa_theta = summary_row$kde_kappa_theta_old
    )
    log_kde <- kde_log_cylinder_density(
      grid$z, grid$theta, case$z, case$theta,
      bandwidths = c(
        g_z = summary_row$kde_g_z,
        kappa_theta = summary_row$kde_kappa_theta_log
      )
    )

    old_metrics <- integrated_discrepancies(parametric, old_kde, theta_values, z_values)
    log_metrics <- integrated_discrepancies(parametric, log_kde, theta_values, z_values)
    data.frame(
      window = summary_row$window,
      pattern = summary_row$pattern,
      n = nrow(case),
      old_l1 = old_metrics[["l1"]],
      log_l1 = log_metrics[["l1"]],
      log_minus_old_l1 = log_metrics[["l1"]] - old_metrics[["l1"]],
      old_total_variation = old_metrics[["total_variation"]],
      log_total_variation = log_metrics[["total_variation"]],
      old_hellinger_squared = old_metrics[["hellinger_squared"]],
      log_hellinger_squared = log_metrics[["hellinger_squared"]],
      old_ise = old_metrics[["integrated_squared_error"]],
      log_ise = log_metrics[["integrated_squared_error"]],
      old_mass_on_comparison_grid = old_metrics[["mass_kde_on_grid"]],
      log_mass_on_comparison_grid = log_metrics[["mass_kde_on_grid"]],
      parametric_mass_on_comparison_grid = log_metrics[["mass_parametric_on_grid"]]
    )
  })
  result <- do.call(rbind, rows)
  result <- result[order(result$pattern, match(result$window, names(month_windows))), ]
  rownames(result) <- NULL
  utils::write.csv(result, output_path, row.names = FALSE)
  result
}

if (sys.nframe() == 0L) {
  result <- compare_estimators()
  print(result, row.names = FALSE, digits = 6)
  cat("\nOutput CSV:\n", output_path, "\n", sep = "")
}
