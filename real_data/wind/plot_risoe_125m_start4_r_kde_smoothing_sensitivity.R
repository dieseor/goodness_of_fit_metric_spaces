#!/usr/bin/env Rscript

# Sensitivity plot in the original r coordinate for Risoe 125 m,
# November-December, start4. Relative to the automatically LCV-selected
# reflected product KDE, both effective bandwidths are multiplied by 1.25:
# h_r <- 1.25 h_r and h_theta <- 1.25 h_theta, equivalently
# kappa_theta <- kappa_theta / 1.25^2.

source(file.path(
  "real_data", "wind", "plot_risoe_125m_nov_dec_b5000_r_cylinder_density_contours.R"
))

smoothing_factor <- 1.25
automatic_r_bandwidth_selector <- select_r_bandwidths
select_r_bandwidths <- function(r, theta) {
  bandwidths <- automatic_r_bandwidth_selector(r, theta)
  bandwidths[["h_r"]] <- smoothing_factor * bandwidths[["h_r"]]
  bandwidths[["kappa_theta"]] <-
    bandwidths[["kappa_theta"]] / smoothing_factor^2
  bandwidths
}
plot_bandwidth_note <- sprintf(
  "; KDE bw × %.2f", smoothing_factor
)

run_start4_r_smoothing_sensitivity <- function(
    output_dir = file.path(
      results_dir, "r_cylinder_density_contours_sensitivity_bw1p25"
    )) {
  validate_plot_samples()
  all_data <- load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"), fixed_tz = "UTC"
  )
  selected <- select_noon_all_months(all_data, fixed_tz = "UTC")
  case <- build_case(selected, "nov_dec", "start4")
  summary <- plot_r_pattern(case, "start4", output_dir)
  utils::write.csv(
    summary,
    file.path(output_dir, "risoe_125m_r_cylinder_density_contours_summary.csv"),
    row.names = FALSE
  )
  invisible(summary)
}

if (sys.nframe() == 0L) {
  result <- run_start4_r_smoothing_sensitivity()
  print(result, row.names = FALSE, digits = 6)
}
