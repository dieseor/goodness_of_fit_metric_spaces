#!/usr/bin/env Rscript

# Regenerate the vMF convergence figures currently included in the manuscript.
# Both the empirical and limiting processes are evaluated on the same fixed
# 10-by-10 grid; this is intentionally not the sample-based KS supremum.

script_path <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
if (length(script_path) != 1L || !nzchar(script_path)) {
  stop("Run this script with Rscript.")
}

script_dir <- dirname(normalizePath(script_path))
project_root <- dirname(script_dir)
setwd(project_root)

suppressPackageStartupMessages({
  library(ggplot2)
  library(rotasym)
})
invisible(capture.output(source("utils.R")))
invisible(capture.output(source("convergence_empirical_process/gaussian_process_vmf.R")))

set.seed(42)

M <- 10000L
N_CORES <- 10L
OMEGA_POINTS <- 10L
T_POINTS <- 10L
DISTANCE_TYPE <- "geodesic"
MU <- c(1, 0, 0)
KAPPAS <- c(0.5, 1, 5)
OUTPUT_DIR <- "output/convergence/gaussian_process/vmf_paper_grid10x10"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# These grids are shared by every n, and by the empirical and limiting laws.
omega_grid <- generate_canonical_lattice(OMEGA_POINTS, dim = length(MU))
t_grid <- seq(1e-8, pi - 1e-8, length.out = T_POINTS)
mu_label <- "mu1p0_0p0_0p0"

save_figure_pair <- function(result, prefix) {
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0(prefix, ".png")),
    plot = result$plot,
    width = 12,
    height = 8,
    dpi = 300
  )
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0("qq_", prefix, ".png")),
    plot = result$qq_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

run_setting <- function(kappa, h0) {
  label <- sprintf("%s, kappa = %g", h0, kappa)
  message("Generating ", label, "...")

  result <- NULL
  invisible(capture.output({
    result <- visualize_convergence_to_limit_vmf(
      n_values = c(50, 100, 500),
      mu = MU,
      kappa = kappa,
      distance_type = DISTANCE_TYPE,
      omega_grid = omega_grid,
      t_grid = t_grid,
      M = M,
      n_cores = N_CORES,
      seed = 42,
      n50_adjust_multiplier = if (h0 == "simple") 6 else 4.5,
      h0 = h0,
      unknown_param = if (h0 == "composite") "xi" else NULL,
      empirical_ks_mode = "grid",
      cov_method = "integral_s2_simple",
      qqplot = TRUE
    )
  }))

  filename_prefix <- if (h0 == "simple") {
    sprintf(
      "simple_vmf_kappa%g_%s_%s_M%d_grid%dx%d",
      kappa, mu_label, DISTANCE_TYPE, M, OMEGA_POINTS, T_POINTS
    )
  } else {
    sprintf(
      "comp_vmf_xi_kappa%g_%s_%s_M%d_grid%dx%d",
      kappa, mu_label, DISTANCE_TYPE, M, OMEGA_POINTS, T_POINTS
    )
  }
  save_figure_pair(result, filename_prefix)
}

for (kappa in KAPPAS) run_setting(kappa, "simple")
for (kappa in KAPPAS) run_setting(kappa, "composite")

message("Created 12 PNG files in: ", normalizePath(OUTPUT_DIR))
