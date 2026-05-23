# Regenerate a selected subset of normal-distribution plots with stronger
# smoothing for n=50 only.

library(ggplot2)

project_root <- normalizePath(file.path(getwd()))
source(file.path(project_root, "utils.R"))
source(file.path(project_root, "convergence_empirical_process", "gaussian_process_normal.R"))

output_dir <- file.path(project_root, "output", "gaussian_process_normal")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

M <- 10000
omega_points <- 10
t_points <- 10
n_cores <- 10
n50_adjust_multiplier <- 4

# Keep the same global grid used in variants so these files are directly comparable
# with prior grid10x10 outputs.
mu_values_global <- c(-3, 0, 3)
sigma_values_global <- c(1, 2, 5)
max_sigma <- max(sigma_values_global)
min_mu <- min(mu_values_global)
max_mu <- max(mu_values_global)
t_max_global <- max_sigma * qnorm(0.995)
omega_min_global <- min_mu - t_max_global
omega_max_global <- max_mu + t_max_global
omega_grid_global <- seq(omega_min_global, omega_max_global, length.out = omega_points)
t_grid_global <- seq(0, t_max_global, length.out = t_points)

run_and_save <- function(mu, sigma, h0, unknown_param = NULL, filename) {
  cat("\nRunning:", filename, "\n")

  res <- visualize_convergence_to_limit_normal(
    n_values = c(50, 100, 500),
    mu = mu,
    sigma = sigma,
    omega_grid = omega_grid_global,
    t_grid = t_grid_global,
    M = M,
    n_cores = n_cores,
    h0 = h0,
    unknown_param = unknown_param,
    n50_adjust_multiplier = n50_adjust_multiplier,
    qqplot = FALSE
  )

  out_path <- file.path(output_dir, filename)
  ggsave(out_path, res$plot, width = 12, height = 8, dpi = 300)
  cat("Saved:", out_path, "\n")
}

# 1) comp_unk_sigma_mu_3_sigma_2_M10000_grid10x10
run_and_save(
  mu = 3,
  sigma = 2,
  h0 = "composite",
  unknown_param = "sigma",
  filename = "comp_unk_sigma_mu_3_sigma_2_M10000_grid10x10.png"
)

# 2) comp_unk_sigma_mu_-3_sigma_2_M10000_grid10x10
run_and_save(
  mu = -3,
  sigma = 2,
  h0 = "composite",
  unknown_param = "sigma",
  filename = "comp_unk_sigma_mu_-3_sigma_2_M10000_grid10x10.png"
)

# 3) simple_mu3_sigma5_M10000_grid10x10
run_and_save(
  mu = 3,
  sigma = 5,
  h0 = "simple",
  filename = "simple_mu3_sigma5_M10000_grid10x10.png"
)

# 4) simple_mu-3_sigma5_M10000_grid10x10
run_and_save(
  mu = -3,
  sigma = 5,
  h0 = "simple",
  filename = "simple_mu-3_sigma5_M10000_grid10x10.png"
)

cat("\nDone: regenerated selected 4 plots with n50_adjust_multiplier=", n50_adjust_multiplier, "\n", sep = "")
