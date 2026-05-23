#!/usr/bin/env Rscript
# Example runner for QQ-plot helpers in tests/test_utils.R
source(file.path("tests", "test_utils.R"))
source(file.path("convergence_empirical_process", "gaussian_process_normal.R"))
source(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))
library(ggplot2)

cat("Running QQ plot examples (small M for speed)...\n")

# Normal example
omega_grid_norm <- seq(-1, 1, length.out = 5)
t_grid_norm <- seq(0, 1, length.out = 5)
mu_norm <- 0
sigma_norm <- 1
n_values <- c(10, 50, 100)
M <- 200

plot_norm <- qqplot_empirical_vs_limit_normal(omega_grid_norm, t_grid_norm, mu_norm, sigma_norm, n_values = n_values, M = M, n_cores = 1, seed = 123)
ggplot2::ggsave(file.path("output", "qqplot_normal_example.png"), plot = plot_norm, width = 6, height = 5)
cat("Saved QQ plot for Normal to output/qqplot_normal_example.png\n")

# vMF example
omega_grid_vmf <- generate_canonical_lattice(5, dim = 3)
t_grid_vmf <- seq(0 + 1e-8, 2 - 1e-8, length.out = 5)
mu_vmf <- c(0, 0, 1)
kappa_vmf <- 2
n_values <- c(10, 50, 100)
M <- 200

plot_vmf <- qqplot_empirical_vs_limit_vmf(omega_grid_vmf, t_grid_vmf, mu_vmf, kappa_vmf, n_values = n_values, M = M, n_mc_samples = 500, n_cores = 1, seed = 123)
ggplot2::ggsave(file.path("output", "qqplot_vmf_example.png"), plot = plot_vmf, width = 6, height = 5)
cat("Saved QQ plot for vMF to output/qqplot_vmf_example.png\n")

cat("Finished QQ plot examples\n")
