# ============================================================================
# GAUSSIAN PROCESS LIMIT THEORY - PERFORMANCE TESTS
# ============================================================================
# This script tests the performance improvement from vectorization
# by comparing the original nested-loop implementation with the
# vectorized version.
#
# Main function: test_vectorization_performance()
# ============================================================================

library(mvtnorm)

# Source the main analysis file for vectorized functions
source("gaussian_process_analysis.R")

# ============================================================================
# BASELINE IMPLEMENTATIONS (for comparison only)
# ============================================================================

#' Simulate empirical process WITHOUT vectorization (ORIGINAL - SLOW)
#' @param omega_grid Vector of omega values
#' @param t_grid Vector of t values
#' @param n Sample size
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @param M Number of Monte Carlo simulations
#' @return Vector of M supremum values
simulate_empirical_process_loops <- function(omega_grid, t_grid, n, mu, sigma, M = 1000) {
  supremum_values <- numeric(M)
  
  for (i in 1:M) {
    # Generate sample
    sample_x <- rnorm(n, mean = mu, sd = sigma)
    
    # Nested loops over grid (SLOW)
    max_difference <- 0
    for (omega in omega_grid) {
      for (t in t_grid) {
        # Empirical distance profile
        empirical_distances <- abs(sample_x - omega)
        F_hat_omega_t <- mean(empirical_distances <= t)
        
        # Theoretical distance profile
        F_omega_t <- compute_distance_profile(omega, mu, sigma, t)
        
        # Scaled difference
        scaled_difference <- sqrt(n) * abs(F_hat_omega_t - F_omega_t)
        
        # Update maximum
        max_difference <- max(max_difference, scaled_difference)
      }
    }
    
    supremum_values[i] <- max_difference
  }
  
  return(supremum_values)
}

#' Compute covariance between two grid points (SCALAR VERSION - SLOW)
#' @param omega1 First location parameter
#' @param t1 First distance threshold
#' @param omega2 Second location parameter
#' @param t2 Second distance threshold
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @return Covariance value
covariance_gaussian_process <- function(omega1, t1, omega2, t2, mu, sigma) {
  if (t1 <= 0 || t2 <= 0) {
    return(0)
  }
  
  # Joint probability
  joint_prob <- compute_joint_probability(omega1, t1, omega2, t2, mu, sigma)
  
  # Marginal probabilities
  f_omega1_t1 <- compute_distance_profile(omega1, mu, sigma, t1)
  f_omega2_t2 <- compute_distance_profile(omega2, mu, sigma, t2)
  
  # Covariance
  return(joint_prob - f_omega1_t1 * f_omega2_t2)
}

#' Create covariance matrix using nested loops (ORIGINAL VERSION - SLOW)
#' @param omega_grid Vector of omega values
#' @param t_grid Vector of t values
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @return Covariance matrix
create_covariance_matrix_loops <- function(omega_grid, t_grid, mu, sigma) {
  n_omega <- length(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  
  grid_combinations <- expand.grid(omega = omega_grid, t = t_grid)
  cov_matrix <- matrix(0, n_total, n_total)
  
  total_iterations <- n_total * n_total
  progress_interval <- max(1, floor(total_iterations / 5))
  iteration_count <- 0
  
  for (i in 1:n_total) {
    for (j in 1:n_total) {
      iteration_count <- iteration_count + 1
      
      if (iteration_count %% progress_interval == 0) {
        progress_pct <- round(100 * iteration_count / total_iterations, 1)
        cat("  Progress:", progress_pct, "%\n")
      }
      
      cov_matrix[i, j] <- covariance_gaussian_process(
        grid_combinations$omega[i],
        grid_combinations$t[i],
        grid_combinations$omega[j],
        grid_combinations$t[j],
        mu, sigma
      )
    }
  }
  
  cat("Covariance matrix created successfully\n")
  return(cov_matrix)
}

# ============================================================================
# PERFORMANCE TESTING
# ============================================================================

#' Test covariance matrix vectorization performance across different grid sizes
#' @return Data frame with performance comparison results
test_covariance_vectorization <- function() {
  cat("\n")
  cat("====================================================================\n")
  cat("  COVARIANCE MATRIX VECTORIZATION TEST\n")
  cat("====================================================================\n")
  cat("\n")
  cat("This test compares the original nested-loop implementation\n")
  cat("with the vectorized version (using single core to isolate\n")
  cat("the impact of vectorization without parallelization overhead).\n")
  cat("\n")
  
  mu <- 0
  sigma <- 1
  
  # Test different grid sizes
  test_sizes <- list(
    small = list(omega = seq(-1, 1, length.out = 5), t = seq(0, 1, length.out = 5)),
    medium = list(omega = seq(-2, 2, length.out = 10), t = seq(0, 5, length.out = 10)),
    large = list(omega = seq(-5, 5, length.out = 20), t = seq(0, 10, length.out = 20)),
    xlarge = list(omega = seq(-5, 5, length.out = 50), t = seq(0, 10, length.out = 50))
  )
  
  results <- data.frame(
    size = character(),
    method = character(),
    grid_points = numeric(),
    time_seconds = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size_name in names(test_sizes)) {
    omega_grid <- test_sizes[[size_name]]$omega
    t_grid <- test_sizes[[size_name]]$t
    grid_points <- length(omega_grid) * length(t_grid)
    
    cat("Testing", size_name, "grid:", length(omega_grid), "x", length(t_grid), "=", grid_points, "points\n")
    
    # Test original method (nested loops - element by element)
    cat("  Running original method (nested loops)...\n")
    start_time <- Sys.time()
    cov_matrix_original <- create_covariance_matrix_loops(omega_grid, t_grid, mu, sigma)
    end_time <- Sys.time()
    time_original <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results <- rbind(results, data.frame(
      size = size_name,
      method = "Original (Nested Loops)",
      grid_points = grid_points,
      time_seconds = time_original
    ))
    
    cat("    Time:", round(time_original, 3), "seconds\n")
    
    # Test vectorized method WITHOUT parallelization (n_cores = 1)
    cat("  Running vectorized method (single core, no parallelization)...\n")
    start_time <- Sys.time()
    cov_matrix_vectorized <- create_covariance_matrix(omega_grid, t_grid, mu, sigma, n_cores = 1)
    end_time <- Sys.time()
    time_vectorized <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results <- rbind(results, data.frame(
      size = size_name,
      method = "Vectorized (Single Core)",
      grid_points = grid_points,
      time_seconds = time_vectorized
    ))
    
    cat("    Time:", round(time_vectorized, 3), "seconds\n")
    
    # Compare accuracy
    max_diff <- max(abs(cov_matrix_original - cov_matrix_vectorized))
    cat("    Maximum difference between methods:", formatC(max_diff, format = "e", digits = 0), "\n")
    
    if (max_diff < 1e-10) {
      cat("    ✅ Results are numerically identical!\n")
    } else {
      cat("    ⚠️  Results differ by more than 1e-10\n")
    }
    
    # Calculate speedup (covariance matrix)
    speedup <- time_original / time_vectorized
    time_reduction <- (1 - time_vectorized/time_original) * 100
    
    cat("    🚀 Vectorization speedup:", round(speedup, 2), "x faster\n")
    if (time_reduction > 0) {
      cat("    💡 Vectorization reduced time by", round(time_reduction, 1), "%\n")
    } else {
      cat("    ⚠️  Vectorization increased time by", round(abs(time_reduction), 1), "% (overhead dominates)\n")
    }
    
    cat("\n")
  }
  
  # Summary table
  cat("=== PERFORMANCE SUMMARY ===\n")
  print(results)
  cat("\n")
  
  # Calculate average speedup
  original_times <- results[results$method == "Original (Nested Loops)", "time_seconds"]
  vectorized_times <- results[results$method == "Vectorized (Single Core)", "time_seconds"]
  
  speedups <- original_times / vectorized_times
  avg_speedup <- mean(speedups)
  avg_reduction <- mean((1 - vectorized_times/original_times) * 100)
  
  cat("📊 VECTORIZATION IMPACT:\n")
  cat("Average speedup from vectorization:", round(avg_speedup, 1), "x\n")
  cat("Average time reduction:", round(avg_reduction, 1), "%\n")
  cat("\n")
  
  return(results)
}

#' Test empirical process vectorization
#' @return Data frame with comparison results
test_empirical_vectorization <- function() {
  cat("\n")
  cat("====================================================================\n")
  cat("  EMPIRICAL PROCESS VECTORIZATION TEST\n")
  cat("====================================================================\n")
  cat("\n")
  cat("This test compares the original nested-loop empirical process\n")
  cat("with the vectorized version (using 1 core to isolate vectorization).\n")
  cat("\n")
  
  mu <- 0
  sigma <- 1
  n <- 1000
  
  # Test different grid sizes
  test_sizes <- list(
    small = list(omega = seq(-1, 1, length.out = 5), t = seq(0, 1, length.out = 5), M = 100),
    medium = list(omega = seq(-2, 2, length.out = 10), t = seq(0, 5, length.out = 10), M = 100),
    large = list(omega = seq(-5, 5, length.out = 20), t = seq(0, 10, length.out = 20), M = 100),
    xlarge = list(omega = seq(-5, 5, length.out = 50), t = seq(0, 10, length.out = 50), M = 100)
  )
  
  results <- data.frame(
    size = character(),
    method = character(),
    grid_points = numeric(),
    simulations = numeric(),
    time_seconds = numeric(),
    mean_supremum = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size_name in names(test_sizes)) {
    omega_grid <- test_sizes[[size_name]]$omega
    t_grid <- test_sizes[[size_name]]$t
    M <- test_sizes[[size_name]]$M
    grid_points <- length(omega_grid) * length(t_grid)
    
    cat("Testing", size_name, "grid:", length(omega_grid), "x", length(t_grid), "=", grid_points, "points,", M, "simulations\n")
    
    # Test original method (nested loops)
    cat("  Running original method (nested loops)...\n")
    start_time <- Sys.time()
    supremum_original <- simulate_empirical_process_loops(omega_grid, t_grid, n, mu, sigma, M)
    end_time <- Sys.time()
    time_original <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results <- rbind(results, data.frame(
      size = size_name,
      method = "Original (Nested Loops)",
      grid_points = grid_points,
      simulations = M,
      time_seconds = time_original,
      mean_supremum = mean(supremum_original)
    ))
    
    cat("    Time:", round(time_original, 3), "seconds\n")
    cat("    Mean supremum:", round(mean(supremum_original), 4), "\n")
    
    # Test vectorized method (1 core = no parallelization)
    cat("  Running vectorized method (1 core, no parallelization)...\n")
    start_time <- Sys.time()
    supremum_vectorized <- simulate_empirical_process(omega_grid, t_grid, n, mu, sigma, M, n_cores = 1)
    end_time <- Sys.time()
    time_vectorized <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results <- rbind(results, data.frame(
      size = size_name,
      method = "Vectorized (Single Core)",
      grid_points = grid_points,
      simulations = M,
      time_seconds = time_vectorized,
      mean_supremum = mean(supremum_vectorized)
    ))
    
    cat("    Time:", round(time_vectorized, 3), "seconds\n")
    cat("    Mean supremum:", round(mean(supremum_vectorized), 4), "\n")
    
    # Compare distributions (not exact values due to different RNG in parallel)
    cat("    Comparing distributions (means should be similar):\n")
    cat("      Original mean:", round(mean(supremum_original), 4), "\n")
    cat("      Vectorized mean:", round(mean(supremum_vectorized), 4), "\n")
    cat("      Difference:", round(abs(mean(supremum_original) - mean(supremum_vectorized)), 4), "\n")
    
    if (abs(mean(supremum_original) - mean(supremum_vectorized)) < 0.05) {
      cat("    ✅ Distributions have similar means (difference < 0.05)\n")
    } else {
      cat("    ⚠️  Distributions differ in mean\n")
    }
    
    # Calculate speedup (empirical process)
    speedup <- time_original / time_vectorized
    time_reduction <- (1 - time_vectorized/time_original) * 100
    
    cat("    🚀 Vectorization speedup:", round(speedup, 2), "x faster\n")
    if (time_reduction > 0) {
      cat("    💡 Vectorization reduced time by", round(time_reduction, 1), "%\n")
    } else {
      cat("    ⚠️  Vectorization increased time by", round(abs(time_reduction), 1), "% (overhead dominates)\n")
    }
    
    cat("\n")
  }
  
  # Summary table
  cat("=== EMPIRICAL PROCESS PERFORMANCE SUMMARY ===\n")
  print(results)
  cat("\n")
  
  # Calculate average speedup
  original_times <- results[results$method == "Original (Nested Loops)", "time_seconds"]
  vectorized_times <- results[results$method == "Vectorized (Single Core)", "time_seconds"]
  
  speedups <- original_times / vectorized_times
  avg_speedup <- mean(speedups)
  avg_reduction <- mean((1 - vectorized_times/original_times) * 100)
  
  cat("📊 VECTORIZATION IMPACT ON EMPIRICAL PROCESS:\n")
  cat("Average speedup from vectorization:", round(avg_speedup, 1), "x\n")
  cat("Average time reduction:", round(avg_reduction, 1), "%\n")
  cat("\n")
  
  return(results)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("\n")
cat("====================================================================\n")
cat("  RUNNING VECTORIZATION PERFORMANCE TESTS\n")
cat("====================================================================\n")
cat("\n")

# Run covariance matrix tests
cat("### PART 1: COVARIANCE MATRIX TESTS ###\n\n")
covariance_results <- test_covariance_vectorization()

cat("\n### PART 2: EMPIRICAL PROCESS TESTS ###\n\n")
empirical_results <- test_empirical_vectorization()

cat("=== ALL VECTORIZATION TESTING COMPLETED ===\n")
cat("Results show the impact of vectorization on computation time.\n")
cat("For small grids, parallelization overhead may dominate.\n")
cat("For large grids, vectorization provides exponential speedup.\n")
