#!/usr/bin/env Rscript

# Original-r cylinder density plots for the strict-step-four Risoe 77 m
# November-December B=5000 analysis. Both estimates are densities with
# respect to dtheta dr; the fitted HvMF therefore includes sinh(r).

source(file.path(
  "real_data", "wind", "plot_risoe_77m_nov_dec_b5000_log_cylinder_density_contours.R"
))
source(file.path(
  "real_data", "wind", "plot_risoe_125m_nov_dec_b5000_r_cylinder_density_contours.R"
))

run_nov_dec_77m_b5000_r_density_plots <- function(
    patterns = names(day_patterns),
    output_dir = file.path(results_dir, "r_cylinder_density_contours")) {
  run_nov_dec_b5000_r_density_plots(patterns = patterns, output_dir = output_dir)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  requested_patterns <- if (length(args) == 0L) names(day_patterns) else args
  result <- run_nov_dec_77m_b5000_r_density_plots(patterns = requested_patterns)
  print(result, row.names = FALSE, digits = 5)
}
