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
source(file.path("tests", "gaussian_process_analysis.R"))
# Source canonical normal functions for testing
source(file.path("convergence_empirical_process", "gaussian_process_normal.R"))
# Source test helpers
source(file.path("tests", "test_utils.R"))

# ============================================================================
# BASELINE IMPLEMENTATIONS (for comparison only)
# ============================================================================

## NOTE: The baseline looped implementation moved into tests/test_utils.R

## NOTE: The scalar covariance helper moved into tests/test_utils.R

## NOTE: The nested loop covariance matrix helper moved into tests/test_utils.R

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
    cov_matrix_vectorized <- cov_normal(omega_grid, t_grid, mu, sigma, n_cores = 1)
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
    supremum_vectorized <- simulate_empirical_process_normal(omega_grid, t_grid, n, mu, sigma, M, n_cores = 1)
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
