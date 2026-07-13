#!/usr/bin/env Rscript

# Final power experiment for the selected nearby HvMF alternative.  It does
# not alter manuscript sources.

source("scripts/run_power_mixtures_paper.R")

hvmf_point_at_angular_offset <- function(offset_deg) {
  chi <- acosh(sqrt(2))
  theta <- pi / 4 + offset_deg * pi / 180
  hyperboloid_point_h2(chi, theta)
}

screening_scenarios <- list(
  hvmf_k200_d12 = list(kappa = 200, offset_deg = 12)
)

mu0_screening <- c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2))

paper_scenario_catalog <- function() {
  lapply(names(screening_scenarios), function(name) {
    setup <- screening_scenarios[[name]]
    list(
      model = "hvmf",
      alternative = sprintf("nearby_mixture_kappa_%d_offset_%d", setup$kappa, setup$offset_deg),
      mu0 = mu0_screening,
      mu1 = hvmf_point_at_angular_offset(setup$offset_deg),
      kappa = setup$kappa
    )
  }) |> stats::setNames(names(screening_scenarios))
}

paper_block_scenarios <- function(block) {
  if (!identical(block, "hvmf")) {
    stop("This script only defines the selected HvMF power scenario.")
  }
  names(screening_scenarios)
}

output_dir <- file.path("simulation_results", "power_hvmf_kappa200_offset12_B5000")

result <- run_power_mixtures_paper_block(
  block = "hvmf",
  output_dir = output_dir,
  M = 1000L,
  B = 5000L,
  n_values = c(50L, 100L, 200L),
  beta_values = c(0.25, 0.5, 1),
  statistics = c("ks", "cvm"),
  bootstrap_method = "fast_multiplier",
  bootstrap_n_cores = 1L,
  derivative_mc_size = 1000L,
  fast_multiplier_cvm_block_size = 50L,
  n_cores_outer = 4L,
  show_progress = TRUE
)

summary_df <- result$summary
summary_df$offset_deg <- vapply(summary_df$scenario, function(name) {
  screening_scenarios[[name]]$offset_deg
}, numeric(1))
summary_df$kappa <- vapply(summary_df$scenario, function(name) {
  screening_scenarios[[name]]$kappa
}, numeric(1))
summary_df$geodesic_separation <- acosh(2 - cos(summary_df$offset_deg * pi / 180))
summary_df <- summary_df[order(summary_df$kappa, summary_df$offset_deg, summary_df$n, summary_df$beta), ]

utils::write.csv(
  summary_df,
  file.path(output_dir, "screening_summary_with_geometry.csv"),
  row.names = FALSE
)

print(summary_df[, c("scenario", "kappa", "offset_deg", "geodesic_separation", "n", "beta", "power_ks_005", "power_cvm_005")])
