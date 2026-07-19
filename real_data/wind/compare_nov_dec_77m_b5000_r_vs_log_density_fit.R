#!/usr/bin/env Rscript

# Coordinate-invariant comparison of fitted HvMF and nonparametric KDE for
# the strict-step-four Risoe 77 m November-December B=5000 samples.

source(file.path(
  "real_data", "wind", "plot_risoe_77m_nov_dec_b5000_r_cylinder_density_contours.R"
))

r_summary_path <- file.path(
  results_dir, "r_cylinder_density_contours",
  "risoe_77m_r_cylinder_density_contours_summary.csv"
)
log_summary_path <- file.path(
  results_dir, "log_cylinder_density_contours",
  "risoe_77m_log_cylinder_density_contours_summary.csv"
)
output_path <- file.path(results_dir, "nov_dec_77m_b5000_r_vs_log_density_fit_comparison.csv")

integrated_metrics <- function(parametric, nonparametric, theta_values, linear_values) {
  n_theta <- length(theta_values)
  n_linear <- length(linear_values)
  to_matrix <- function(x) t(matrix(x, nrow = n_theta, ncol = n_linear))
  p <- to_matrix(parametric)
  q <- to_matrix(nonparametric)
  integrate <- function(x) trapezoid_2d_mass(x, theta_values, linear_values)
  mass_p <- integrate(p)
  mass_q <- integrate(q)
  p <- p / mass_p
  q <- q / mass_q
  l1 <- integrate(abs(p - q))
  h2 <- 0.5 * integrate((sqrt(p) - sqrt(q))^2)
  c(
    parametric_mass_on_grid = mass_p,
    nonparametric_mass_on_grid = mass_q,
    l1 = l1,
    total_variation = l1 / 2,
    overlap = 1 - l1 / 2,
    hellinger_squared = h2,
    hellinger_affinity = 1 - h2
  )
}

compare_nov_dec_77m_r_vs_log <- function(n_theta = 241L, n_linear = 401L) {
  r_summary <- utils::read.csv(r_summary_path, stringsAsFactors = FALSE)
  log_summary <- utils::read.csv(log_summary_path, stringsAsFactors = FALSE)
  summaries <- merge(
    r_summary, log_summary, by = c("window", "pattern"), suffixes = c("_r", "_log")
  )
  if (nrow(summaries) != 4L) {
    stop("Expected four matched starts in the 77 m r/log summaries.", call. = FALSE)
  }

  all_data <- load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"
  )
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")
  theta_values <- seq(0, 2 * pi, length.out = n_theta)
  r_values <- seq(0, 10, length.out = n_linear)
  z_values <- seq(-10, log(10), length.out = n_linear)
  r_grid <- expand.grid(theta = theta_values, r = r_values)
  z_grid <- expand.grid(theta = theta_values, z = z_values)

  rows <- lapply(seq_len(nrow(summaries)), function(i) {
    row <- summaries[i, , drop = FALSE]
    case <- build_case(selected, "nov_dec", row$pattern)
    fit <- hvmf_mle_h2(as.matrix(case[, c("x0", "x1", "x2")]))
    if (nrow(case) != row$n_r || nrow(case) != row$n_log) {
      stop("Sample-size mismatch in the 77 m density comparison.", call. = FALSE)
    }
    p_r <- hvmf_r_cylinder_density(r_grid$r, r_grid$theta, fit)
    q_r <- kde_r_cylinder_density(
      r_grid$r, r_grid$theta, case$r, case$theta,
      c(h_r = row$kde_h_r, kappa_theta = row$kde_kappa_theta_r)
    )
    p_log <- hvmf_log_cylinder_density(z_grid$z, z_grid$theta, fit)
    q_log <- kde_log_cylinder_density(
      z_grid$z, z_grid$theta, case$z, case$theta,
      c(g_z = row$kde_g_z, kappa_theta = row$kde_kappa_theta_log)
    )
    metrics_r <- integrated_metrics(p_r, q_r, theta_values, r_values)
    metrics_log <- integrated_metrics(p_log, q_log, theta_values, z_values)
    data.frame(
      pattern = row$pattern, n = nrow(case),
      r_l1 = metrics_r[["l1"]], log_l1 = metrics_log[["l1"]],
      r_total_variation = metrics_r[["total_variation"]],
      log_total_variation = metrics_log[["total_variation"]],
      r_overlap = metrics_r[["overlap"]], log_overlap = metrics_log[["overlap"]],
      r_hellinger_squared = metrics_r[["hellinger_squared"]],
      log_hellinger_squared = metrics_log[["hellinger_squared"]],
      r_hellinger_affinity = metrics_r[["hellinger_affinity"]],
      log_hellinger_affinity = metrics_log[["hellinger_affinity"]],
      r_parametric_mass_on_grid = metrics_r[["parametric_mass_on_grid"]],
      r_nonparametric_mass_on_grid = metrics_r[["nonparametric_mass_on_grid"]],
      log_parametric_mass_on_grid = metrics_log[["parametric_mass_on_grid"]],
      log_nonparametric_mass_on_grid = metrics_log[["nonparametric_mass_on_grid"]]
    )
  })
  result <- do.call(rbind, rows)
  result <- result[order(match(result$pattern, paste0("start", 1:4))), ]
  rownames(result) <- NULL
  utils::write.csv(result, output_path, row.names = FALSE)
  result
}

if (sys.nframe() == 0L) {
  result <- compare_nov_dec_77m_r_vs_log()
  print(result, row.names = FALSE, digits = 6)
  cat("\nOutput CSV:\n", output_path, "\n", sep = "")
}
