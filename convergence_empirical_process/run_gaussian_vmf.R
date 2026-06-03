# ============================================================================
# GAUSSIAN PROCESS LIMIT THEORY - VON MISES-FISHER DISTRIBUTION EXECUTION
# ============================================================================
# This script runs the comparison between empirical and limit Gaussian
# processes for goodness-of-fit testing with von Mises-Fisher distribution
# on the unit sphere.
#
# Main function: visualize_convergence_to_limit_vmf()
# Output: Histogram comparison plots and KS test results
# ============================================================================

# Load necessary libraries
library(rotasym) 
library(ggplot2)
library(parallel)

# Get script directory for proper sourcing
script_path <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
if (length(script_path) == 0 || script_path == "") {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
project_root <- dirname(script_dir)

# Load generic utilities and vMF-specific functions
source(file.path(project_root, "utils.R"))
source(file.path(script_dir, "gaussian_process_vmf.R"))

# Get current working directory (equivalent to os.getcwd())

# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================

# Reproducibility
set.seed(42)  # Fixed seed for reproducible results across runs

# Grid parameters - 10x10 GRID
OMEGA_POINTS <- 10     # Number of omega points (canonical lattice on sphere)
T_POINTS <- 10         # Number of t points in grid

# Simulation parameters - FULL SIMULATIONS
M_SIMULATIONS <- 10000     # Number of Monte Carlo simulations for supremum
N_MC_COVARIANCE <- 100000   # Number of MC samples for covariance (50K+ for large grids!)

# Density smoothing adjustment: increase this >1.0 to smooth density curves for large M
DENSITY_ADJUST <- 1.5

# Distribution parameters (vMF on unit sphere S^2 in R^3)
MU_VMF <- c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3))   # Mean direction (north pole)
KAPPA_VMF <- 0.5       # Concentration parameter

# Distance type default - will be possibly overwritten by CLI args at bottom
DISTANCE_TYPE <- "geodesic"  # Either "chordal" or "geodesic"

# Parallelization
N_CORES <- 10          # Number of cores for parallel computation

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Only run if executed directly (not when sourced by tests)
## ---------------------------------------------------------------------------
## Execution dispatch - mirror normal runner structure
## ---------------------------------------------------------------------------
run_all_vmf <- function(output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_cores = N_CORES, density_adjust = DENSITY_ADJUST, distance_type = 'geodesic') {
  simple_res <- run_simple_vmf(output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_mc_samples = N_MC_COVARIANCE, n_cores = n_cores, distance_type = distance_type, density_adjust = density_adjust)
  comp_res <- run_composite_vmf(output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_mc_samples = N_MC_COVARIANCE, n_cores = n_cores, distance_type = distance_type, density_adjust = density_adjust)
  # Also run variants for the simple null (mirror composite variants)
  simple_variants_res <- run_simple_variants_vmf(mu_values = list(c(1,0,0), 
  c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), c(-1/sqrt(2), 0, -1/sqrt(2))), kappa_values = c(0.5, 
  1, 5), output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_mc_samples = N_MC_COVARIANCE, n_cores = n_cores, distance_type = distance_type, density_adjust = density_adjust)
  # Also run composite variants as part of 'all'
  comp_variants_res <- run_composite_variants_vmf(mu_values = list(c(1,0,0), 
  c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), c(-1/sqrt(2), 0, -1/sqrt(2))), kappa_values = c(0.5, 
  1, 5), unknown_params = c('xi'), output_dir = output_dir, M = M, omega_points = omega_points, t_points = t_points, n_mc_samples = N_MC_COVARIANCE, n_cores = n_cores, distance_type = distance_type, density_adjust = density_adjust)
  return(c(simple_res, comp_res, simple_variants_res))
}


# Helper: composite variants for vMF
run_composite_variants_vmf <- function(mu_values = list(c(1,0,0), 
c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), c(-1/sqrt(2), 0, -1/sqrt(2))), kappa_values = c(0.5, 
1, 5), unknown_params = c('xi'), output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_mc_samples = N_MC_COVARIANCE, n_cores = N_CORES, distance_type = 'geodesic', density_adjust = DENSITY_ADJUST, n50_multiplier_overrides = NULL) {
  if (is.null(output_dir)) output_dir <- 'output/convergence/gaussian_process/vmf'
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  scenarios <- list()
  for (mu_val in mu_values) {
    # Build a compact mu label safe for file names like mu_1p000_0p000_0p000
    mu_label <- paste0('mu', paste0(sapply(mu_val, function(x) gsub("\\.", 'p', sprintf('%.1f', x))), collapse = '_'))
    for (k in kappa_values) {
      scenarios[[length(scenarios) + 1]] <- list(mu = mu_val, kappa = k, label = sprintf('vMF_kappa%g_%s', k, mu_label), mu_label = mu_label)
    }
  }
  results <- list()
  # Compute a shared absolute grid spanning possible mu/kappa ranges (rough heuristic)
  # For the sphere, we may reuse canonical lattice and same omega_grid across scenarios
  omega_grid_global <- generate_canonical_lattice(omega_points, dim = length(MU_VMF))
  t_grid_global <- seq(0 + 1e-8, ifelse(distance_type == 'chordal', 2 - 1e-8, pi - 1e-8), length.out = t_points)
  for (unk in unknown_params) {
    cat('Running composite variants for vMF unknown param =', unk, '\n')
    results[[unk]] <- list()
    for (scn in scenarios) {
      # Check if there's a custom multiplier for this mu-kappa combination
      n50_mult <- 4.5  # Default multiplier for n=50 adjustment.
      if (!is.null(n50_multiplier_overrides)) {
        # Check for a match: list(mu_idx = 2, kappa_idx = 1) -> multiplier = 5.5
        for (override in n50_multiplier_overrides) {
          # Find mu index
          mu_idx <- which(sapply(mu_values, function(m) all(m == scn$mu)))[1]
          kappa_idx <- which(kappa_values == scn$kappa)[1]
          if (!is.na(mu_idx) && !is.na(kappa_idx) && override$mu_idx == mu_idx && override$kappa_idx == kappa_idx) {
            n50_mult <- override$multiplier
            cat(sprintf("  Using custom n50_multiplier = %.1f for mu_idx=%d, kappa_idx=%d\n", n50_mult, mu_idx, kappa_idx))
            break
          }
        }
      }
      res <- visualize_convergence_to_limit_vmf(n_values = c(50, 100, 500), mu = scn$mu, kappa = scn$kappa, distance_type = distance_type, omega_grid = omega_grid_global, t_grid = t_grid_global, M = M, n_mc_samples = n_mc_samples, n_cores = n_cores, seed = 42, n50_adjust_multiplier = n50_mult, h0 = 'composite', unknown_param = unk, qqplot = TRUE)
      fname <- file.path(output_dir, sprintf('comp_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', unk, scn$kappa, scn$mu_label, distance_type, M, omega_points, t_points))
      if (!is.null(res$qq_plot)) {
        qq_fname <- file.path(output_dir, sprintf('qq_comp_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', unk, scn$kappa, scn$mu_label, distance_type, M, omega_points, t_points))
        ggsave(qq_fname, res$qq_plot, width = 8, height = 6, dpi = 300)
      }
      ggsave(fname, res$plot, width = 12, height = 8, dpi = 300)
      results[[unk]][[scn$label]] <- res
    }
  }
  return(results)
}


# Helper: simple variants for vMF (mirror composite variants behavior)
run_simple_variants_vmf <- function(mu_values = list(c(1,0,0), 
c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), c(-1/sqrt(2), 0, -1/sqrt(2))), kappa_values = c(0.5, 
1, 5), output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_mc_samples = N_MC_COVARIANCE, n_cores = N_CORES, distance_type = 'geodesic', density_adjust = DENSITY_ADJUST, cov_method = c("mc", "exact_s1_simple", "integral_s2_simple")) {
  cov_method <- match.arg(cov_method)
  if (is.null(output_dir)) output_dir <- 'output/convergence/gaussian_process/vmf'
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  scenarios <- list()
  for (mu_val in mu_values) {
    # Build a compact mu label safe for file names like mu_1p000_0p000_0p000
    mu_label <- paste0('mu', paste0(sapply(mu_val, function(x) gsub("\\.", 'p', sprintf('%.1f', x))), collapse = '_'))
    for (k in kappa_values) {
      scenarios[[length(scenarios) + 1]] <- list(mu = mu_val, kappa = k, label = sprintf('vMF_kappa%g_%s', k, mu_label), mu_label = mu_label)
    }
  }
  results <- list()
  omega_grid_global <- generate_canonical_lattice(omega_points, dim = length(MU_VMF))
  t_grid_global <- seq(0 + 1e-8, ifelse(distance_type == 'chordal', 2 - 1e-8, pi - 1e-8), length.out = t_points)
  cat('Running simple variants for vMF (h0 = simple)\n')
  for (scn in scenarios) {
    res <- visualize_convergence_to_limit_vmf(n_values = c(50, 100, 500), mu = scn$mu, kappa = scn$kappa, distance_type = distance_type, omega_grid = omega_grid_global, t_grid = t_grid_global, M = M, n_mc_samples = n_mc_samples, n_cores = n_cores, seed = 42, n50_adjust_multiplier = 6, h0 = 'simple', cov_method = cov_method, qqplot = TRUE)
    fname <- file.path(output_dir, sprintf('simple_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', cov_method, scn$kappa, scn$mu_label, distance_type, M, omega_points, t_points))
    if (!is.null(res$qq_plot)) {
      qq_fname <- file.path(output_dir, sprintf('qq_simple_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', cov_method, scn$kappa, scn$mu_label, distance_type, M, omega_points, t_points))
      ggsave(qq_fname, res$qq_plot, width = 8, height = 6, dpi = 300)
    }
    ggsave(fname, res$plot, width = 12, height = 8, dpi = 300)
    results[[scn$label]] <- res
  }
  return(results)
}


# ---------------------------------------------------------------------------
# New: simple/composite run helpers for vMF (mirror Normal structure)
# ---------------------------------------------------------------------------
run_simple_vmf <- function(scenarios = NULL, output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_mc_samples = N_MC_COVARIANCE, n_cores = N_CORES, density_adjust = DENSITY_ADJUST, distance_type = 'geodesic', cov_method = c("mc", "exact_s1_simple", "integral_s2_simple")) {
  cov_method <- match.arg(cov_method)
  if (is.null(output_dir)) output_dir <- 'output/convergence/gaussian_process/vmf'
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(scenarios)) {
    mu_label <- paste0('mu', paste0(sapply(MU_VMF, function(x) gsub("\\.", 'p', sprintf('%.1f', x))), collapse = '_'))
    scenarios <- list(list(mu = MU_VMF, kappa = KAPPA_VMF, label = sprintf('vMF_kappa%g_%s', KAPPA_VMF, mu_label), mu_label = mu_label))
  }
  results <- list()
  for (s in scenarios) {
    res <- visualize_convergence_to_limit_vmf(n_values = c(50, 100, 500), mu = s$mu, kappa = s$kappa, distance_type = distance_type, omega_points = omega_points, t_points = t_points, M = M, n_mc_samples = n_mc_samples, n_cores = n_cores, seed = 42, n50_adjust_multiplier = 6, h0 = 'simple', cov_method = cov_method, qqplot = TRUE)
    fname <- file.path(output_dir, sprintf('simple_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', cov_method, s$kappa, s$mu_label, distance_type, M, omega_points, t_points))
    if (!is.null(res$qq_plot)) {
      qq_fname <- file.path(output_dir, sprintf('qq_simple_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', cov_method, s$kappa, s$mu_label, distance_type, M, omega_points, t_points))
      ggsave(qq_fname, res$qq_plot, width = 8, height = 6, dpi = 300)
    }
    ggsave(fname, res$plot, width = 12, height = 8, dpi = 300)
    results[[s$label]] <- res
  }
  return(results)
}

run_simple_integral_vmf <- function(scenarios = NULL,
                                    output_dir = NULL,
                                    M = M_SIMULATIONS,
                                    omega_points = OMEGA_POINTS,
                                    t_points = T_POINTS,
                                    n_mc_samples = N_MC_COVARIANCE,
                                    n_cores = N_CORES,
                                    density_adjust = DENSITY_ADJUST,
                                    distance_type = 'geodesic') {
  run_simple_vmf(
    scenarios = scenarios,
    output_dir = output_dir,
    M = M,
    omega_points = omega_points,
    t_points = t_points,
    n_mc_samples = n_mc_samples,
    n_cores = n_cores,
    density_adjust = density_adjust,
    distance_type = distance_type,
    cov_method = "integral_s2_simple"
  )
}

run_simple_variants_integral_vmf <- function(mu_values = list(c(1,0,0),
                                                              c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)),
                                                              c(-1/sqrt(2), 0, -1/sqrt(2))),
                                             kappa_values = c(0.5, 1, 5),
                                             output_dir = NULL,
                                             M = M_SIMULATIONS,
                                             omega_points = OMEGA_POINTS,
                                             t_points = T_POINTS,
                                             n_mc_samples = N_MC_COVARIANCE,
                                             n_cores = N_CORES,
                                             distance_type = 'geodesic',
                                             density_adjust = DENSITY_ADJUST) {
  run_simple_variants_vmf(
    mu_values = mu_values,
    kappa_values = kappa_values,
    output_dir = output_dir,
    M = M,
    omega_points = omega_points,
    t_points = t_points,
    n_mc_samples = n_mc_samples,
    n_cores = n_cores,
    distance_type = distance_type,
    density_adjust = density_adjust,
    cov_method = "integral_s2_simple"
  )
}

run_composite_vmf <- function(scenarios = NULL, output_dir = NULL, M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, n_mc_samples = N_MC_COVARIANCE, n_cores = N_CORES, distance_type = 'geodesic', unknown_param = c('xi')) {
  if (is.null(output_dir)) output_dir <- 'output/convergence/gaussian_process/vmf'
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(scenarios)) {
    mu_label <- paste0('mu', paste0(sapply(MU_VMF, function(x) gsub("\\.", 'p', sprintf('%.2f', x))), collapse = '_'))
    scenarios <- list(list(mu = MU_VMF, kappa = KAPPA_VMF, label = sprintf('vMF_kappa%g_%s', KAPPA_VMF, mu_label), mu_label = mu_label))
  }
  results <- list()
  for (s in scenarios) {
    for (unk in unknown_param) {
      res <- visualize_convergence_to_limit_vmf(n_values = c(50, 100, 500), mu = s$mu, kappa = s$kappa, distance_type = distance_type, omega_points = omega_points, t_points = t_points, M = M, n_mc_samples = n_mc_samples, n_cores = n_cores, seed = 42, n50_adjust_multiplier = 4.5, h0 = 'composite', unknown_param = unk, qqplot = TRUE)
      fname <- file.path(output_dir, sprintf('comp_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', unk, s$kappa, s$mu_label, distance_type, M, omega_points, t_points))
      if (!is.null(res$qq_plot)) {
        qq_fname <- file.path(output_dir, sprintf('qq_comp_vmf_%s_kappa%g_%s_%s_M%d_grid%dx%d.png', unk, s$kappa, s$mu_label, distance_type, M, omega_points, t_points))
        ggsave(qq_fname, res$qq_plot, width = 8, height = 6, dpi = 300)
      }
      ggsave(fname, res$plot, width = 12, height = 8, dpi = 300)
      results[[paste0(s$label, '_comp_', unk)]] <- res
    }
  }
  return(results)
}


## ---------------------------------------------------------------------------
## CLI dispatch: only run this after all functions are defined
## ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  cat("  GAUSSIAN PROCESS LIMIT - vMF\n")
  args <- commandArgs(trailingOnly = TRUE)
  action <- if (length(args) > 0) args[1] else 'all'
  # Allow overriding M on command line as second arg
  if (length(args) > 1 && !is.na(as.numeric(args[2]))) M_SIMULATIONS <- as.numeric(args[2])
  # Optionally allow overriding the distance type via the third arg
  if (length(args) > 2 && args[3] %in% c('geodesic', 'chordal')) DISTANCE_TYPE <- args[3]

  if (action == 'simple') {
    run_simple_vmf(M = M_SIMULATIONS, distance_type = DISTANCE_TYPE)
  } else if (action == 'simple_variants') {
    run_simple_variants_vmf(mu_values = list(c(1,0,0), 
    c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), c(-1/sqrt(2), 0, -1/sqrt(2))), kappa_values = c(0.5, 
    1, 5), M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, distance_type = DISTANCE_TYPE)
  } else if (action == 'composite') {
    run_composite_vmf(M = M_SIMULATIONS, distance_type = DISTANCE_TYPE)
  } else if (action == 'composite_variants') {
    run_composite_variants_vmf(mu_values = list(c(1,0,0), 
    c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), c(-1/sqrt(2), 0, -1/sqrt(2))), kappa_values = c(0.5, 
    1, 5), unknown_params = c('xi'), M = M_SIMULATIONS, omega_points = OMEGA_POINTS, t_points = T_POINTS, distance_type = DISTANCE_TYPE)
  } else if (action == 'all') {
    run_all_vmf(M = M_SIMULATIONS, distance_type = DISTANCE_TYPE)
  } else {
    run_simple_vmf(M = M_SIMULATIONS, distance_type = DISTANCE_TYPE)
  }
  cat('\nAction executed:', action, '\n')
}
