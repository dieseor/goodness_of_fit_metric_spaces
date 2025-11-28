# ============================================================================
# GAUSSIAN PROCESS LIMIT THEORY - NORMAL DISTRIBUTION EXECUTION
# ============================================================================
# This script runs the comparison between empirical and limit Gaussian
# processes for goodness-of-fit testing with Normal distribution.
#
# Main function: visualize_convergence_to_limit_normal()
# Output: Histogram comparison plots and KS test results
# ============================================================================

# Load necessary libraries
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(parallel)

# Load generic utilities and Normal-specific functions
source(file.path("utils.R"))
source(file.path("convergence_empirical_process", "gaussian_process_normal.R"))

# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================

# Reproducibility
set.seed(42)  # Fixed seed for reproducible results across runs

# Grid parameters
OMEGA_MIN <- -2.5      # Minimum omega value
OMEGA_MAX <- 2.5       # Maximum omega value
OMEGA_POINTS <- 5       # Number of omega points in grid

T_MIN <- 0             # Minimum t value
T_MAX <- 5           # Maximum t value
T_POINTS <- 5         # Number of t points in grid

# Simulation parameters
M_SIMULATIONS <- 10000 # Number of Monte Carlo simulations

# Distribution parameters
MU_TRUE <- 0           # True mean of X ~ N(mu, sigma^2)
SIGMA_TRUE <- 1        # True standard deviation of X

# Parallelization
N_CORES <- 10          # Number of cores for parallel computation

# H0 settings: 'simple' or 'composite' with unknown_param in c(NULL, 'mu', 'sigma', 'both')
H0_TYPE <- 'simple'
UNKNOWN_PARAM <- NULL

# Density smoothing adjustment for density curves (similar to vMF)
DENSITY_ADJUST <- 1.5

# ============================================================================
# MAIN EXECUTION
# ============================================================================

run_simple_normal <- function(scenarios = NULL, output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_cores = N_CORES, density_adjust = DENSITY_ADJUST, omega_grid = NULL, t_grid = NULL) {
  if (is.null(output_dir)) {
    output_dir <- "output/gaussian_process_normal"
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(scenarios)) {
    scenarios <- list(
      list(mu = 0, sigma = 1, label = "N(0,1)"),
      list(mu = 3, sigma = 1, label = "N(3,1)"),
      list(mu = 0, sigma = 5, label = "N(0,25)")
    )
  }
  all_results <- list()
  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    cat("\n")
    cat("====================================================================\n")
    cat("  SCENARIO", i, ":", scenario$label, "\n")
    cat("====================================================================\n")
    cat("\n")
    convergence_result <- visualize_convergence_to_limit_normal(
      n_values = c(50, 100, 500),
      mu = scenario$mu,
      sigma = scenario$sigma,
      omega_points = omega_points,
      t_points = t_points,
      omega_grid = omega_grid,
      t_grid = t_grid,
      M = M,
      n_cores = n_cores,
      h0 = H0_TYPE, density_adjust = density_adjust
    )
    convergence_filename <- sprintf("simple_mu%g_sigma%g_M%d_grid%dx%d.png",
                    scenario$mu,
                    scenario$sigma,
                    as.integer(M),
                    as.integer(omega_points),
                    as.integer(t_points))
    convergence_path <- file.path(output_dir, convergence_filename)
    ggsave(convergence_path, convergence_result$plot, width = 12, height = 8, dpi = 300)
    cat("\nSaved plot for", scenario$label, "to:", convergence_path, "\n")
    all_results[[scenario$label]] <- convergence_result
  }
  return(all_results)
}

run_composite_normal <- function(scenarios = NULL, output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_cores = N_CORES, density_adjust = DENSITY_ADJUST, unknown_param = c('mu','sigma','both'), omega_grid = NULL, t_grid = NULL) {
  if (is.null(output_dir)) {
    output_dir <- "output/gaussian_process_normal"
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(scenarios)) scenarios <- list(
    list(mu = 0, sigma = 1, label = "N(0,1)"),
    list(mu = 3, sigma = 1, label = "N(3,1)"),
    list(mu = 0, sigma = 5, label = "N(0,25)")
  )
  # Use all provided scenarios; we want to run every scenario with each unknown_param case
  composite_scenarios <- scenarios
  all_results <- list()
  for (i in seq_along(composite_scenarios)) {
    scenario <- composite_scenarios[[i]]
    cat("\n====================================================================\n")
    cat("  COMPOSITE H0:", scenario$label, "\n")
    cat("====================================================================\n")
    cat("\n")
    for (unk in unknown_param) {
      cat("\n  COMPOSITE H0: scenario =", scenario$label, "unknown_param =", unk, "\n")
      convergence_result_comp <- visualize_convergence_to_limit_normal(
        n_values = c(50, 100, 500),
      mu = scenario$mu,
      sigma = scenario$sigma,
      omega_points = omega_points,
      t_points = t_points,
      omega_grid = omega_grid,
      t_grid = t_grid,
      M = M,
      n_cores = n_cores,
      h0 = 'composite',
      unknown_param = unk,
      density_adjust = density_adjust
    )
      convergence_filename_comp <- sprintf("comp_unk_%s_mu_%g_sigma_%g_M%d_grid%dx%d.png",
                        unk,
                                          scenario$mu,
                                          scenario$sigma,
                                        M,
                                        omega_points,
                                        t_points)
    convergence_path_comp <- file.path(output_dir, convergence_filename_comp)
    ggsave(convergence_path_comp, convergence_result_comp$plot, width = 12, height = 8, dpi = 300)
    cat("\nSaved composite-null plot for", scenario$label, "to:", convergence_path_comp, "\n")
    cat("  limit mean:", round(mean(as.numeric(convergence_result_comp$limit_values), na.rm=TRUE), 4), "sd:", round(sd(as.numeric(convergence_result_comp$limit_values), na.rm=TRUE), 4), " (unknown_param=", unk, ")\n")
    all_results[[paste0(scenario$label, "_comp_unk_", unk)]] <- convergence_result_comp
    }
  }
  return(all_results)
}


# Mimic run_composite_variants for the simple null case
run_simple_variants <- function(mu_values = c(0, -1, 1), sigma_values = c(1, 2, 5), output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_cores = N_CORES, density_adjust = DENSITY_ADJUST) {
  scenarios <- list()
  for (mu_val in mu_values) {
    for (sigma_val in sigma_values) {
      scenarios[[length(scenarios) + 1]] <- list(mu = mu_val, sigma = sigma_val, label = sprintf("N(mu=%g,sigma=%g)", mu_val, sigma_val))
    }
  }
  results <- list()
  summary_rows <- list()
  max_sigma <- max(sigma_values)
  min_mu <- min(mu_values)
  max_mu <- max(mu_values)
  t_max_global <- max_sigma * qnorm(0.995)
  omega_min_global <- min_mu - t_max_global
  omega_max_global <- max_mu + t_max_global
  omega_grid_global <- seq(omega_min_global, omega_max_global, length.out = omega_points)
  t_grid_global <- seq(0, t_max_global, length.out = t_points)

  # For simple case, use run_simple_normal to process all scenarios in one batch
  res <- run_simple_normal(
    scenarios = scenarios,
    output_dir = output_dir,
    M = M,
    omega_points = omega_points,
    t_points = t_points,
    n_cores = n_cores,
    density_adjust = density_adjust,
    omega_grid = omega_grid_global,
    t_grid = t_grid_global
  )
  results <- res
  # Summarize KS test p-values to CSV rows
  for (scn_name in names(res)) {
    res_scn <- res[[scn_name]]
    mu_val <- scenarios[[which(sapply(scenarios, function(s) s$label == scn_name))]]$mu
    sigma_val <- scenarios[[which(sapply(scenarios, function(s) s$label == scn_name))]]$sigma
    ks <- res_scn$ks_results
    for (nchar in names(ks)) {
      ksr <- ks[[nchar]]
      pval <- if (!is.null(ksr$p.value)) ksr$p.value else NA_real_
      stat <- if (!is.null(ksr$statistic)) as.numeric(ksr$statistic) else NA_real_
      summary_rows[[length(summary_rows) + 1]] <- list(
        scenario = scn_name,
        mu = mu_val,
        sigma = sigma_val,
        n = as.numeric(nchar),
        p_value = pval,
        D_stat = stat,
        M = M
      )
    }
  }
  # Write summary CSV if output_dir provided
  if (is.null(output_dir)) output_dir <- "output/gaussian_process_normal"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, lapply(summary_rows, as.data.frame, stringsAsFactors = FALSE))
    csv_name <- file.path(output_dir, sprintf("simple_variants_summary_M%d.csv", M))
    write.csv(summary_df, csv_name, row.names = FALSE)
    cat("Saved simple variants summary to", csv_name, "\n")
  }
  return(results)
}

run_composite_variants <- function(mu_values = c(0, -1, 1), sigma_values = c(1, 2, 5), unknown_params = c('mu', 'sigma', 'both'), nominal_mu = 0, output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_cores = N_CORES, density_adjust = DENSITY_ADJUST) {
  scenarios <- list()
  for (mu_val in mu_values) {
    for (sigma_val in sigma_values) {
      scenarios[[length(scenarios) + 1]] <- list(mu = mu_val, sigma = sigma_val, label = sprintf("N(mu=%g,sigma=%g)", mu_val, sigma_val))
    }
  }
  # Run composite for all scenarios (do not filter to mu==0 so we get the full grid)
  results <- list()
  summary_rows <- list()
  # Create a fixed global grid spanning all scenarios so results are comparable across sigma
  max_sigma <- max(sigma_values)
  min_mu <- min(mu_values)
  max_mu <- max(mu_values)
  t_max_global <- max_sigma * qnorm(0.995)
  omega_min_global <- min_mu - t_max_global
  omega_max_global <- max_mu + t_max_global
  omega_grid_global <- seq(omega_min_global, omega_max_global, length.out = omega_points)
  t_grid_global <- seq(0, t_max_global, length.out = t_points)

  for (unk in unknown_params) {
    cat('\nRunning composite variants with unknown_param = ', unk, '\n')
    results[[unk]] <- run_composite_normal(scenarios = scenarios, output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_cores = n_cores, density_adjust = density_adjust, unknown_param = unk, omega_grid = omega_grid_global, t_grid = t_grid_global)
    # Summarize KS test p-values to CSV rows
    for (scn_name in names(results[[unk]])) {
      res_scn <- results[[unk]][[scn_name]]
      # scn_name is like 'N(mu=0,sigma=1)_comp_unk_mu' ; strip the suffix to find original scenario label
      orig_label <- sub("_comp_unk_.*$", "", scn_name)
      scn_info <- Filter(function(s) s$label == orig_label, scenarios)[[1]]
      mu_val <- scn_info$mu
      sigma_val <- scn_info$sigma
      ks <- res_scn$ks_results
      for (nchar in names(ks)) {
        ksr <- ks[[nchar]]
        pval <- if (!is.null(ksr$p.value)) ksr$p.value else NA_real_
        stat <- if (!is.null(ksr$statistic)) as.numeric(ksr$statistic) else NA_real_
        summary_rows[[length(summary_rows) + 1]] <- list(
          unknown_param = unk,
          scenario = scn_name,
          mu = mu_val,
          sigma = sigma_val,
          n = as.numeric(nchar),
          p_value = pval,
          D_stat = stat,
          M = M
        )
      }
    }
  }
  # Write summary CSV if output_dir provided
  if (is.null(output_dir)) output_dir <- "output/gaussian_process_normal"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, lapply(summary_rows, as.data.frame, stringsAsFactors = FALSE))
  csv_name <- file.path(output_dir, sprintf("comp_variants_summary_M%d_%s.csv", M, paste(unknown_params, collapse = ",")))
    write.csv(summary_df, csv_name, row.names = FALSE)
    cat("Saved composite variants summary to", csv_name, "\n")
  }
  return(results)
}

run_all_normal <- function(output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_cores = N_CORES, density_adjust = DENSITY_ADJUST) {
  simple_res <- run_simple_normal(output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_cores = n_cores, density_adjust = density_adjust)
  composite_res <- run_composite_normal(output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_cores = n_cores, density_adjust = density_adjust)
  simple_variants_res <- run_simple_variants(mu_values = c(0, -3, 3), sigma_values = c(1, 2, 5), output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_cores = n_cores, density_adjust = density_adjust)
  composite_variants_res <- run_composite_variants(mu_values = c(0, -3, 3), sigma_values = c(1, 2, 5), unknown_params = c("mu", "sigma", "both"), output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_cores = n_cores, density_adjust = density_adjust)
  return(list(
    simple = simple_res,
    composite = composite_res,
    simple_variants = simple_variants_res,
    composite_variants = composite_variants_res
  ))
}


# Only run if executed directly (not when sourced by tests)
if (sys.nframe() == 0) {
  # Allow first command-line arg to choose 'simple', 'composite', 'simple_variants', 'composite_variants', or 'all'
  args <- commandArgs(trailingOnly = TRUE)
  action <- if (length(args) > 0) args[1] else 'all'
  # Allow second arg to override M
  if (length(args) > 1 && !is.na(as.numeric(args[2]))) M_SIMULATIONS <- as.numeric(args[2])
  # Map action
  if (action == 'simple') {
    run_simple_normal(M = M_SIMULATIONS)
  } else if (action == 'composite') {
    run_composite_normal(M = M_SIMULATIONS)
  } else if (action == 'simple_variants') {
    # Run a grid of mu/sigma combos for simple testing
    run_simple_variants(mu_values = c(0, -3, 3), sigma_values = c(1, 2, 5), M = M_SIMULATIONS)
  } else if (action == 'composite_variants') {
    # Run a grid of mu/sigma combos for composite testing
    run_composite_variants(mu_values = c(0, -3, 3), sigma_values = c(1, 2, 5), unknown_params = c("mu", "sigma", "both"),  M = M_SIMULATIONS)
  } else {
    run_all_normal(M = M_SIMULATIONS)
  }
  cat("\nAction executed:", action, "\n")
}

