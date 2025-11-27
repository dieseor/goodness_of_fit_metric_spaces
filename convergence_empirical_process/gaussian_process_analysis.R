# POSSIBLE DUPLICATE FILE?

# ============================================================================
# GAUSSIAN PROCESS LIMIT THEORY - MAIN ANALYSIS
# ============================================================================
# This script implements the comparison between empirical and limit Gaussian
# processes for goodness-of-fit testing in metric spaces.
#
# Main function: compare_processes()
# Output: Histogram comparison plot and KS test results
# ============================================================================

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(gridExtra)

# ============================================================================
# CONFIGURATION PARAMETERS - EDIT HERE TO CHANGE GRID SIZES AND SIMULATIONS
# ============================================================================

# Reproducibility
set.seed(42)  # Fixed seed for reproducible results across runs

# Grid parameters
OMEGA_MIN <- -2.5      # Minimum omega value
OMEGA_MAX <- 2.5       # Maximum omega value
OMEGA_POINTS <- 50    # Number of omega points in grid

T_MIN <- 0           # Minimum t value
T_MAX <- 2.5           # Maximum t value
T_POINTS <- 50        # Number of t points in grid

# Simulation parameters
N_SAMPLE <- 1000       # Sample size for empirical process
M_SIMULATIONS <- 10000  # Number of Monte Carlo simulations

# Distribution parameters
MU_TRUE <- 0         # True mean of X ~ N(mu, sigma^2)
SIGMA_TRUE <- 1      # True standard deviation of X

# Parallelization
N_CORES <- 10       # Number of cores for parallel computation

# ============================================================================
# AUXILIARY FUNCTIONS (Internal)
# ============================================================================

#' Compute theoretical distance profile for univariate normal distribution
#' @param omega Location parameter(s) (can be vector)
#' @param mu Mean of X ~ N(mu, sigma^2)
#' @param sigma Standard deviation of X
#' @param t_values Distance threshold(s) (can be vector)
#' @return Probability that distance from X to omega is <= t
compute_distance_profile <- function(omega, mu, sigma, t_values) {
  # Vectorize to handle different length inputs
  n <- max(length(omega), length(t_values))
  omega <- rep(omega, length.out = n)
  t_values <- rep(t_values, length.out = n)
  
  result <- numeric(n)
  valid_t <- t_values > 0
  
  if (any(valid_t)) {
    upper_bound <- (omega[valid_t] + t_values[valid_t] - mu) / sigma
    lower_bound <- (omega[valid_t] - t_values[valid_t] - mu) / sigma
    result[valid_t] <- pnorm(upper_bound) - pnorm(lower_bound)
  }
  
  return(result)
}

#' Compute joint probability for covariance calculation
#' @param omega1 First location parameter
#' @param t1 First distance threshold
#' @param omega2 Second location parameter
#' @param t2 Second distance threshold
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @return Joint probability P(|X - omega1| <= t1, |X - omega2| <= t2)
compute_joint_probability <- function(omega1, t1, omega2, t2, mu, sigma) {
  if (t1 <= 0 || t2 <= 0) {
    return(0)
  }
  
  # Intervals
  lower1 <- omega1 - t1
  upper1 <- omega1 + t1
  lower2 <- omega2 - t2
  upper2 <- omega2 + t2
  
  # Intersection
  intersection_lower <- max(lower1, lower2)
  intersection_upper <- min(upper1, upper2)
  
  if (intersection_lower >= intersection_upper) {
    return(0)
  }
  
  # P(intersection_lower < X < intersection_upper)
  return(pnorm((intersection_upper - mu) / sigma) - 
         pnorm((intersection_lower - mu) / sigma))
}

#' Compute one row of the covariance matrix (vectorized)
#' @param omega1 Location parameter for row
#' @param t1 Distance threshold for row
#' @param omega2_vec Vector of location parameters for columns
#' @param t2_vec Vector of distance thresholds for columns
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @return Vector of covariances for the entire row
compute_covariance_row <- function(omega1, t1, omega2_vec, t2_vec, mu, sigma) {
  n <- length(omega2_vec)
  
  # Vectorized computation of interval intersections
  lower1 <- omega1 - t1
  upper1 <- omega1 + t1
  lower2_vec <- omega2_vec - t2_vec
  upper2_vec <- omega2_vec + t2_vec
  
  intersection_lower_vec <- pmax(lower1, lower2_vec)
  intersection_upper_vec <- pmin(upper1, upper2_vec)
  
  valid_intersections <- intersection_lower_vec < intersection_upper_vec
  
  # Joint probabilities (vectorized)
  joint_prob_vec <- numeric(n)
  if (any(valid_intersections)) {
    valid_lower <- intersection_lower_vec[valid_intersections]
    valid_upper <- intersection_upper_vec[valid_intersections]
    joint_prob_vec[valid_intersections] <- 
      pnorm((valid_upper - mu) / sigma) - pnorm((valid_lower - mu) / sigma)
  }
  
  # Marginal probabilities
  f_omega1_t1 <- compute_distance_profile(omega1, mu, sigma, t1)
  f_omega2_t2_vec <- compute_distance_profile(omega2_vec, mu, sigma, t2_vec)
  
  # Covariance = Joint - Product of marginals
  cov_vec <- joint_prob_vec - f_omega1_t1 * f_omega2_t2_vec
  
  return(cov_vec)
}

# ============================================================================
# MAIN FUNCTIONS (Public)
# ============================================================================

#' Create covariance matrix for Gaussian process (vectorized + parallelized)
#' @param omega_grid Vector of omega values
#' @param t_grid Vector of t values
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @param n_cores Number of cores for parallel computation (default: 10)
#' @return Covariance matrix of size (n_omega * n_t) x (n_omega * n_t)
create_covariance_matrix <- function(omega_grid, t_grid, mu, sigma, n_cores = 10) {
  n_omega <- length(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  
  # Create all combinations of (omega, t) pairs
  grid_combinations <- expand.grid(omega = omega_grid, t = t_grid)
  
  cat("Creating covariance matrix of size", n_total, "x", n_total, "\n")
  cat("Using", n_cores, "cores\n")

  # Extract vectors for vectorized operations
  omega_vec <- grid_combinations$omega
  t_vec <- grid_combinations$t
  
  # Setup parallel cluster
  library(parallel)
  cl <- makeCluster(n_cores)
  
  # Set seed for reproducibility in parallel workers
  clusterSetRNGStream(cl, iseed = 42)
  
  # Export helper function to workers
  clusterEvalQ(cl, {
    compute_distance_profile <- function(omega, mu, sigma, t_values) {
      n <- max(length(omega), length(t_values))
      omega <- rep(omega, length.out = n)
      t_values <- rep(t_values, length.out = n)
      result <- numeric(n)
      valid_t <- t_values > 0
      if (any(valid_t)) {
        upper_bound <- (omega[valid_t] + t_values[valid_t] - mu) / sigma
        lower_bound <- (omega[valid_t] - t_values[valid_t] - mu) / sigma
        result[valid_t] <- pnorm(upper_bound) - pnorm(lower_bound)
      }
      return(result)
    }
    
    compute_covariance_row <- function(omega1, t1, omega2_vec, t2_vec, mu, sigma) {
      n <- length(omega2_vec)
      lower1 <- omega1 - t1
      upper1 <- omega1 + t1
      lower2_vec <- omega2_vec - t2_vec
      upper2_vec <- omega2_vec + t2_vec
      intersection_lower_vec <- pmax(lower1, lower2_vec)
      intersection_upper_vec <- pmin(upper1, upper2_vec)
      valid_intersections <- intersection_lower_vec < intersection_upper_vec
      joint_prob_vec <- numeric(n)
      if (any(valid_intersections)) {
        valid_lower <- intersection_lower_vec[valid_intersections]
        valid_upper <- intersection_upper_vec[valid_intersections]
        joint_prob_vec[valid_intersections] <- 
          pnorm((valid_upper - mu) / sigma) - pnorm((valid_lower - mu) / sigma)
      }
      f_omega1_t1 <- compute_distance_profile(omega1, mu, sigma, t1)
      f_omega2_t2_vec <- compute_distance_profile(omega2_vec, mu, sigma, t2_vec)
      cov_vec <- joint_prob_vec - f_omega1_t1 * f_omega2_t2_vec
      return(cov_vec)
    }
  })
  
  clusterExport(cl, c("omega_vec", "t_vec", "mu", "sigma"), envir = environment())
  
  start_time <- Sys.time()
  
  # Parallel computation of covariance matrix
  cov_matrix <- tryCatch({
    # Distribute rows across cores (round-robin for load balancing)
    row_chunks <- vector("list", n_cores)
    for (i in 1:n_cores) {
      row_chunks[[i]] <- seq(from = i, to = n_total, by = n_cores)
    }
    
    cat("Computing", sum(sapply(row_chunks, length)), "rows in parallel\n")
    
    # Compute rows in parallel
    results <- parLapply(cl, row_chunks, function(row_indices) {
      local_rows <- list()
      for (i in row_indices) {
        row_vec <- compute_covariance_row(
          omega_vec[i], t_vec[i], omega_vec, t_vec, mu, sigma
        )
        local_rows[[length(local_rows) + 1]] <- list(i = i, row = row_vec)
      }
      return(local_rows)
    })
    
    # Assemble matrix
    cat("Assembling matrix from parallel results...\n")
    full_matrix <- matrix(0, n_total, n_total)
    for (worker_results in results) {
      for (result in worker_results) {
        full_matrix[result$i, ] <- result$row
      }
    }
    
    full_matrix
  }, finally = {
    stopCluster(cl)
  })
  
  end_time <- Sys.time()
  time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (time_elapsed < 1) {
    cat("Vectorized + parallel covariance matrix created in", round(time_elapsed, 1), "seconds!\n")
  } else {
    cat("Vectorized + parallel covariance matrix created in", round(time_elapsed, 1), "seconds!\n")
  }
  
  return(cov_matrix)
}

#' Simulate the limiting Gaussian process
#' @param omega_grid Vector of omega values
#' @param t_grid Vector of t values
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @param M Number of Monte Carlo simulations
#' @param n_cores Number of cores for covariance matrix computation
#' @return Vector of M supremum values from the Gaussian process
simulate_limit_gaussian <- function(omega_grid, t_grid, mu, sigma, M = 10000, n_cores = 10) {
  
  cat("=== Simulating Gaussian Process Limit (VECTORIZED) ===\n")
  cat("Grid size: ", length(omega_grid), "omegas x", length(t_grid), "t values\n")
  cat("Monte Carlo simulations:", M, "\n\n")
  
  # Create covariance matrix
  cov_matrix <- create_covariance_matrix(omega_grid, t_grid, mu, sigma, n_cores)
  n_total <- nrow(cov_matrix)
  
  # Memory estimate
  memory_gb <- (n_total^2 * 8) / (1024^3)
  cat("Generating", M, "multivariate normal samples...\n")
  
  # Generate M samples from multivariate normal
  gaussian_samples <- rmvnorm(M, mean = rep(0, n_total), sigma = cov_matrix)
  
  # Compute supremum for each sample
  supremum_values <- apply(gaussian_samples, 1, function(row) max(abs(row)))
  
  cat("Supremum statistics, limiting Gaussian process:\n")
  cat("  Mean:", round(mean(supremum_values), 4), "\n")
  cat("  Median:", round(median(supremum_values), 4), "\n")
  cat("  Max:", round(max(supremum_values), 4), "\n")
  cat("  Min:", round(min(supremum_values), 4), "\n\n")
  
  return(supremum_values)
}

#' Simulate the empirical process (VECTORIZED + PARALLELIZED)
#' @param omega_grid Vector of omega values
#' @param t_grid Vector of t values
#' @param n Sample size
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @param M Number of Monte Carlo simulations
#' @return Vector of M supremum values from the empirical process
simulate_empirical_process <- function(omega_grid, t_grid, n, mu, sigma, M = 10000, n_cores = 10) {
  
  cat("=== Simulating Empirical Process ===\n")
  cat("Grid size:", length(omega_grid), "omegas x", length(t_grid), "t values\n")
  cat("Sample size n =", n, "\n")
  cat("Monte Carlo simulations:", M, "\n")
  cat("Using", n_cores, "cores\n\n")
  
  n_omega <- length(omega_grid)
  n_t <- length(t_grid)
  
  # Pre-compute theoretical values for all (omega, t) pairs (vectorized)
  grid_df <- expand.grid(omega = omega_grid, t = t_grid)
  F_theoretical <- compute_distance_profile(grid_df$omega, mu, sigma, grid_df$t)
  F_theoretical_matrix <- matrix(F_theoretical, nrow = n_omega, ncol = n_t)
  # Dimensión: n_omega × n_t
  
  # Setup parallel cluster
  library(parallel)
  cl <- makeCluster(n_cores)
  
  # Set seed for reproducibility in parallel workers
  clusterSetRNGStream(cl, iseed = 123)
  
  # Export helper function to workers
  clusterEvalQ(cl, {
    compute_distance_profile <- function(omega, mu, sigma, t_values) {
      n <- max(length(omega), length(t_values))
      omega <- rep(omega, length.out = n)
      t_values <- rep(t_values, length.out = n)
      result <- numeric(n)
      valid_t <- t_values > 0
      if (any(valid_t)) {
        upper_bound <- (omega[valid_t] + t_values[valid_t] - mu) / sigma
        lower_bound <- (omega[valid_t] - t_values[valid_t] - mu) / sigma
        result[valid_t] <- pnorm(upper_bound) - pnorm(lower_bound)
      }
      return(result)
    }
  })
  
  # Export necessary variables
  clusterExport(cl, c("omega_grid", "t_grid", "n", "mu", "sigma", 
                      "F_theoretical_matrix", "n_omega", "n_t"), 
                envir = environment())
  
  start_time <- Sys.time()
  
  # Parallel computation of simulations
  supremum_values <- tryCatch({
    # Distribute simulations across cores (round-robin)
    sim_chunks <- vector("list", n_cores)
    for (i in 1:n_cores) {
      sim_chunks[[i]] <- seq(from = i, to = M, by = n_cores)
    }
    
    cat("Computing", M, "simulations in parallel across", n_cores, "cores...\n")
    
    # Run simulations in parallel
    results <- parLapply(cl, sim_chunks, function(sim_indices) {
      local_supremums <- numeric(length(sim_indices))
      
      for (idx in seq_along(sim_indices)) {
        # Generate sample
        sample_x <- rnorm(n, mean = mu, sd = sigma)
        
        # VECTORIZED: Compute distance matrix (n × n_omega)
        distance_matrix <- outer(sample_x, omega_grid, function(x, w) abs(x - w))
        
        # VECTORIZED: Compute F_hat for ALL (omega, t) pairs
        F_hat_matrix <- sapply(t_grid, function(t_val) {
          colMeans(distance_matrix <= t_val)
        })
        
        # Compute scaled differences
        scaled_diff_matrix <- sqrt(n) * abs(F_hat_matrix - F_theoretical_matrix)
        
        # Take supremum
        local_supremums[idx] <- max(scaled_diff_matrix)
      }
      
      return(local_supremums)
    })
    
    # Combine results from all workers
    cat("Assembling results from parallel workers...\n")
    all_supremums <- numeric(M)
    for (worker_idx in 1:n_cores) {
      indices <- sim_chunks[[worker_idx]]
      all_supremums[indices] <- results[[worker_idx]]
    }
    
    all_supremums
  }, finally = {
    stopCluster(cl)
  })
  
  end_time <- Sys.time()
  time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("Empirical process simulated in", round(time_elapsed, 1), "seconds\n\n")
  
  cat("Supremum statistics, empirical process:\n")
  cat("  Mean:", round(mean(supremum_values), 4), "\n")
  cat("  Median:", round(median(supremum_values), 4), "\n")
  cat("  Max:", round(max(supremum_values), 4), "\n")
  cat("  Min:", round(min(supremum_values), 4), "\n\n")
  
  return(supremum_values)
}

#' Compare empirical and limit processes side by side
#' @return List with limit_values, empirical_values, and ggplot object
compare_processes <- function() {
  cat("=== Comparing Empirical vs Limit Processes ===\n")
  
  # Use global configuration parameters
  omega_grid <- seq(OMEGA_MIN, OMEGA_MAX, length.out = OMEGA_POINTS)
  t_grid <- seq(T_MIN, T_MAX, length.out = T_POINTS)
  
  cat("Omega grid:", OMEGA_POINTS, "points from", OMEGA_MIN, "to", OMEGA_MAX, "\n")
  cat("t grid:", T_POINTS, "points from", T_MIN, "to", T_MAX, "\n")
  cat("Total grid points:", OMEGA_POINTS * T_POINTS, "\n")
  cat("Monte Carlo simulations:", M_SIMULATIONS, "each process\n")
  cat("Using", N_CORES, "cores for parallelization\n\n")
  
  # Memory warning
  n_grid <- OMEGA_POINTS * T_POINTS
  memory_gb <- (n_grid^2 * 8) / (1024^3)
  cat("WARNING: This will create a", n_grid, "x", n_grid, "covariance matrix\n")
  cat("Estimated memory usage: ~", round(memory_gb, 2), "GB\n\n")
  
  # Simulate limiting Gaussian process
  limit_values <- simulate_limit_gaussian(
    omega_grid, t_grid, MU_TRUE, SIGMA_TRUE, M_SIMULATIONS, N_CORES
  )
  
  # Simulate empirical process
  empirical_values <- simulate_empirical_process(
    omega_grid, t_grid, N_SAMPLE, MU_TRUE, SIGMA_TRUE, M_SIMULATIONS, N_CORES
  )
  
  # Statistical comparison
  cat("=== Statistical Comparison ===\n")
  ks_result <- safe_ks_test(limit_values, empirical_values)
  pval <- if (!is.null(ks_result$p.value)) ks_result$p.value else NA_real_
  if (is.na(pval)) {
    cat("Kolmogorov-Smirnov test p-value: NA (test not run / invalid input)\n\n")
  } else {
    cat("Kolmogorov-Smirnov test p-value:", round(pval, 4), "\n")
    if (pval > 0.05) {
      cat("-> Distributions are NOT significantly different :) \n\n")
    } else {
      cat("-> Distributions ARE significantly different :( \n\n")
    }
  }
  
  # Create comparison plot
  comparison_data <- data.frame(
    values = c(limit_values, empirical_values),
    process = rep(c("Limit Gaussian", paste0("Empirical (n=", N_SAMPLE, ")")), each = M_SIMULATIONS)
  )
  
  # Colors for histograms
  limit_color <- "#0066CC"  # Dark blue for limit
  empirical_color <- "#FF3333"  # Red for empirical
  
  p_comparison <- ggplot(comparison_data, aes(x = values, fill = process)) +
    geom_histogram(
      aes(y = after_stat(density)),
      alpha = 0.3, 
      position = "identity", 
      bins = 50
    ) +
    geom_density(
      aes(color = process, linetype = process),
      linewidth = 1.0,
      adjust = 1.0,
      key_glyph = draw_key_path
    ) +
    scale_fill_manual(
      values = setNames(
        c(limit_color, empirical_color),
        c("Limit Gaussian", paste0("Empirical (n=", N_SAMPLE, ")"))
      ),
      name = "Histogram"
    ) +
    scale_color_manual(
      values = setNames(
        c(limit_color, empirical_color),
        c("Limit Gaussian", paste0("Empirical (n=", N_SAMPLE, ")"))
      ),
      name = "Histogram"
    ) +
    scale_linetype_manual(
      values = c("dashed", "solid"),
      name = "Histogram"
    ) +
    labs(
      x = "Supremum of the process",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19),
      axis.title.y = element_text(size = 19)
    )
  
  return(list(
    limit_values = limit_values,
    empirical_values = empirical_values,
    plot = p_comparison,
    ks_test = ks_result
  ))
}

#' Visualize convergence of empirical process to limit process as n grows
#' Shows how density of empirical process approaches the Gaussian limit
#' @param n_values Vector of sample sizes to compare (default: c(10, 50, 100, 500))
#' @param mu Mean of the normal distribution (default: uses MU_TRUE)
#' @param sigma Standard deviation of the normal distribution (default: uses SIGMA_TRUE)
#' @return List with all simulated values and convergence plot
visualize_convergence_to_limit <- function(n_values = c(10, 50, 100, 500), 
                                          mu = NULL, 
                                          sigma = NULL) {
  # Use provided mu/sigma or fall back to global constants
  if (is.null(mu)) mu <- MU_TRUE
  if (is.null(sigma)) sigma <- SIGMA_TRUE
  
  cat("=== Visualizing Convergence to Gaussian Limit ===\n")
  cat("Sample sizes to compare:", paste(n_values, collapse = ", "), "\n")
  cat("Distribution: N(μ =", mu, ", σ² =", sigma^2, ")\n")
  
  # Automatic grid computation based on mu and sigma
  # For X ~ N(mu, sigma^2), we want P(|X - mu| < T_MAX) >= 0.99
  # This gives T_MAX = sigma * qnorm(0.995) ≈ 2.576 * sigma
  t_max_auto <- sigma * qnorm(0.995)
  omega_max_auto <- t_max_auto
  omega_min_auto <- -t_max_auto
  t_min_auto <- 0
  
  omega_grid <- seq(omega_min_auto, omega_max_auto, length.out = OMEGA_POINTS)
  t_grid <- seq(t_min_auto, t_max_auto, length.out = T_POINTS)
  
  cat("  Omega grid:", OMEGA_POINTS, "points from", round(omega_min_auto, 2), 
      "to", round(omega_max_auto, 2), "\n")
  cat("  t grid:", T_POINTS, "points from", round(t_min_auto, 2), 
      "to", round(t_max_auto, 2), "\n")
  cat("Monte Carlo simulations:", M_SIMULATIONS, "each process\n")
  cat("Using", N_CORES, "cores for parallelization\n\n")
  
  # Simulate limiting Gaussian process (once)
  cat("Simulating limit Gaussian process...\n")
  limit_values <- simulate_limit_gaussian(
    omega_grid, t_grid, mu, sigma, M_SIMULATIONS, N_CORES
  )
  cat("  Mean:", round(mean(limit_values), 4), "\n")
  cat("  SD:", round(sd(limit_values), 4), "\n\n")
  
  # Simulate empirical process for each n
  empirical_data <- list()
  for (n in n_values) {
    cat("Simulating empirical process with n =", n, "...\n")
    empirical_values <- simulate_empirical_process(
      omega_grid, t_grid, n, mu, sigma, M_SIMULATIONS, N_CORES
    )
    empirical_data[[as.character(n)]] <- empirical_values
  }
  
  # Prepare data for plotting
  # Histograms: Limit process and largest n empirical process
  n_max <- max(n_values)
  histogram_data <- data.frame(
    values = c(limit_values, empirical_data[[as.character(n_max)]]),
    process = rep(c("Limit Gaussian", paste0("Empirical (n=", n_max, ")")), 
                  each = M_SIMULATIONS)
  )
  
  # Density curves: All n values
  density_data <- data.frame()
  for (n in n_values) {
    n_char <- as.character(n)
    density_data <- rbind(
      density_data,
      data.frame(
        values = empirical_data[[n_char]],
        n = n,
        label = paste0("n=", n)
      )
    )
  }
  
  # Add limit process to density data
  density_data_full <- rbind(
    density_data,
    data.frame(
      values = limit_values,
      n = Inf,
      label = "Limit Gaussian"
    )
  )
  
  # Create color palette using rainbow
  # Limit process: dark blue (#0066CC)
  # Empirical processes: rainbow colors (reversed so red is at largest n)
  limit_color <- "#0066CC"  # Dark blue for limit
  n_colors <- length(n_values)
  
  # Use rainbow colors for empirical processes (reversed)
  empirical_colors <- rev(rainbow(n_colors))
  
  # Combine colors: empirical + limit
  all_colors <- c(empirical_colors, limit_color)
  names(all_colors) <- c(paste0("n=", n_values), "Limit Gaussian")
  
  # Create plot WITHOUT automatic legend
  p_convergence <- ggplot() +
    geom_histogram(
      data = histogram_data,
      aes(x = values, y = after_stat(density), fill = process),
      alpha = 0.3,
      position = "identity",
      bins = 50,
      show.legend = FALSE
    ) +
    geom_density(
      data = density_data_full,
      aes(x = values, color = label, linetype = label),
      linewidth = 1.0,
      adjust = 1.0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = all_colors) +
    scale_linetype_manual(values = rep("solid", n_colors + 1)) +
    scale_fill_manual(
      values = setNames(
        c(limit_color, "#FF3333"),
        c("Limit Gaussian", paste0("Empirical (n=", n_max, ")"))
      )
    ) +
    labs(
      x = "Supremum of the process",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19),
      axis.title.y = element_text(size = 19)
    ) +
    coord_cartesian(xlim = c(NA, 2))
  
  # Manual legend positioning: compute positions relative to x-axis limits or data
  if (!is.null(xlim) && length(xlim) == 2 && !any(is.na(xlim))) {
    x_upper <- xlim[2]
  } else {
    x_candidates <- histogram_data$values
    if (exists('limit_values') && !is.null(limit_values)) x_candidates <- c(x_candidates, limit_values)
    x_data_max <- if (length(x_candidates) > 0) max(x_candidates, na.rm = TRUE) else 2.0
    x_upper <- max(2.0, x_data_max * 1.05)
  }
  legend_x <- x_upper * 0.76
  legend_y_start <- 1.55
  legend_spacing <- 0.06
  legend_line_length <- x_upper * 0.03
  legend_rect_height <- 0.023
  legend_rect_half_width <- max(0.2, x_upper * 0.12)
  legend_line_x_start <- legend_x - legend_rect_half_width * 0.35
  legend_text_x <- legend_x - legend_rect_half_width * 0.2
  
  # Add white background box for legend
  p_convergence <- p_convergence +
    annotate("rect",
             xmin = legend_x - legend_rect_half_width, xmax = legend_x + legend_rect_half_width,
             ymin = legend_y_start - (n_colors + 2) * legend_spacing - 0.03,
             ymax = legend_y_start + 0.02,
             fill = "transparent")
  
  # Show all empirical density curves
  for (i in 1:n_colors) {
    y_pos <- legend_y_start - (i-1) * legend_spacing
    
    # Draw line (all solid)
    p_convergence <- p_convergence +
      annotate("segment",
           x = legend_line_x_start, xend = legend_line_x_start + legend_line_length,
               y = y_pos, yend = y_pos,
               color = all_colors[i],
               linetype = "solid",
               linewidth = 1.0) +
      annotate("text",
           x = legend_text_x,
               y = y_pos,
               label = paste0("italic(n)==", n_values[i]),
               parse = TRUE,
               hjust = 0,
               size = 5.5,
               color = "black")
  }
  
  # Add Limit Gaussian density line entry (solid, not dashed)
  y_pos_limit <- legend_y_start - n_colors * legend_spacing
  p_convergence <- p_convergence +
    annotate("segment",
         x = legend_line_x_start, xend = legend_line_x_start + legend_line_length,
             y = y_pos_limit, yend = y_pos_limit,
             color = limit_color,
             linetype = "solid",
             linewidth = 1.0) +
    annotate("text",
         x = legend_text_x,
             y = y_pos_limit,
             label = "'𝔾'[mu[theta[0]]]",
             parse = TRUE,
             hjust = 0,
             size = 5.5,
             color = "black")
  
  # Histogram boxes at the end (continuing with same spacing)
  # Limit Gaussian histogram - position after all density lines
  y_pos_hist1 <- legend_y_start - (n_colors + 1) * legend_spacing
  p_convergence <- p_convergence +
    annotate("rect",
         xmin = legend_line_x_start, xmax = legend_line_x_start + legend_line_length,
             ymin = y_pos_hist1 - legend_rect_height, ymax = y_pos_hist1 + legend_rect_height,
             fill = limit_color, alpha = 0.3, color = NA) +
    annotate("text",
         x = legend_text_x,
             y = y_pos_hist1,
             label = "'𝔾'[mu[theta[0]]]",
             parse = TRUE,
             hjust = 0,
             size = 5.5,
             color = "black")
  
  # Largest n histogram - next position
  y_pos_hist2 <- legend_y_start - (n_colors + 2) * legend_spacing
  p_convergence <- p_convergence +
    annotate("rect",
         xmin = legend_line_x_start, xmax = legend_line_x_start + legend_line_length,
             ymin = y_pos_hist2 - legend_rect_height, ymax = y_pos_hist2 + legend_rect_height,
             fill = "#FF3333", alpha = 0.3, color = NA) +
    annotate("text",
         x = legend_text_x,
             y = y_pos_hist2,
             label = paste0("italic(n)==", n_max),
             parse = TRUE,
             hjust = 0,
             size = 5.5,
             color = "black")
  
  # Statistical comparison for each n
  cat("=== Kolmogorov-Smirnov Tests vs Limit ===\n")
  ks_results <- list()
  for (n in n_values) {
    n_char <- as.character(n)
    ks_test <- safe_ks_test(limit_values, empirical_data[[n_char]])
    ks_results[[n_char]] <- ks_test
    
    cat(sprintf("n=%5d: p-value = %.4f", n, ks_test$p.value))
    if (ks_test$p.value > 0.05) {
      cat(" -> Not significantly different ✓\n")
    } else {
      cat(" -> Significantly different\n")
    }
  }
  cat("\n")
  
  return(list(
    limit_values = limit_values,
    empirical_data = empirical_data,
    n_values = n_values,
    plot = p_convergence,
    ks_results = ks_results
  ))
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Only run if executed directly (not when sourced by tests)
if (sys.nframe() == 0) {
  cat("\n")
  cat("====================================================================\n")
  cat("  GAUSSIAN PROCESS LIMIT ANALYSIS - GOODNESS-OF-FIT TESTING\n")
  cat("====================================================================\n")
  cat("\n")
  
  # Create output directory
  output_dir <- "output/gaussian_process"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define the three scenarios to test
  scenarios <- list(
    list(mu = 0, sigma = 1, label = "N(0,1)"),
    list(mu = 3, sigma = 1, label = "N(3,1)"),
    list(mu = 0, sigma = 5, label = "N(0,25)")
  )
  
  # Run convergence analysis for each scenario
  all_results <- list()
  
  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    
    cat("\n")
    cat("====================================================================\n")
    cat("  SCENARIO", i, ":", scenario$label, "\n")
    cat("====================================================================\n")
    cat("\n")
    
    # Run convergence visualization with specific mu and sigma
    convergence_result <- visualize_convergence_to_limit(
      n_values = c(10, 50, 100, 500),
      mu = scenario$mu,
      sigma = scenario$sigma
    )
    
    # Save convergence plot
    convergence_filename <- sprintf("convergence_mu%g_sigma%g_M%d_grid%dx%d.png", 
                                    scenario$mu,
                                    scenario$sigma,
                                    M_SIMULATIONS, 
                                    OMEGA_POINTS, 
                                    T_POINTS)
    convergence_path <- file.path(output_dir, convergence_filename)
    
    ggsave(convergence_path, convergence_result$plot, 
           width = 12, height = 8, dpi = 300)
    cat("\nSaved plot for", scenario$label, "to:", convergence_path, "\n")
    
    # Store results
    all_results[[scenario$label]] <- convergence_result
  }
  
  cat("\n")
  cat("====================================================================\n")
  cat("  ALL SCENARIOS COMPLETED\n")
  cat("====================================================================\n")
  cat("  Sample sizes tested:", paste(c(10, 50, 100, 500), collapse = ", "), "\n")
  cat("  Simulations per n: M =", M_SIMULATIONS, "\n")
  cat("  Grid points: ω =", OMEGA_POINTS, ", t =", T_POINTS, "(auto-computed)\n")
  cat("  Cores:", N_CORES, "\n")
  cat("  Scenarios tested:\n")
  for (scenario in scenarios) {
    cat("    -", scenario$label, "\n")
  }
  cat("  All plots saved to", output_dir, "directory\n")
  cat("\n")
}