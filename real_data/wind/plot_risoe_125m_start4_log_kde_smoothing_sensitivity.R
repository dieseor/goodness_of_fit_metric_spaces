#!/usr/bin/env Rscript

# Sensitivity plot for Risoe 125 m, November-December, start4.
# Relative to the automatically LCV-selected product KDE, both effective
# bandwidths are multiplied by 1.25: g_z <- 1.25 g_z and
# h_theta <- 1.25 h_theta, equivalently kappa_theta <- kappa_theta / 1.25^2.

source(file.path(
  "real_data", "wind", "plot_risoe_125m_nov_dec_b5000_log_cylinder_density_contours.R"
))

smoothing_factor <- 1.25
automatic_bandwidth_selector <- select_cylinder_bandwidths
select_cylinder_bandwidths <- function(z, theta) {
  bandwidths <- automatic_bandwidth_selector(z, theta)
  bandwidths[["g_z"]] <- smoothing_factor * bandwidths[["g_z"]]
  bandwidths[["kappa_theta"]] <-
    bandwidths[["kappa_theta"]] / smoothing_factor^2
  bandwidths
}
plot_bandwidth_note <- sprintf(
  "; sensitivity KDE bandwidths × %.2f", smoothing_factor
)

run_start4_smoothing_sensitivity <- function(
    output_dir = file.path(
      results_dir, "log_cylinder_density_contours_sensitivity_bw1p25"
    )) {
  validate_plot_samples()
  run_cylinder_plots(output_dir = output_dir, patterns = "start4")
}

if (sys.nframe() == 0L) {
  result <- run_start4_smoothing_sensitivity()
  print(result, row.names = FALSE, digits = 6)
}
