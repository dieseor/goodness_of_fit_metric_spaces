# Test to measure the impact of triangular inequality optimization
# on covariance matrix computation time

library(rotasym)
library(parallel)

source(file.path("utils.R"))
source(file.path("tests", "test_utils.R"))
source(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))

cat("\n")
cat("=========================================================\n")
cat("  TESTING TRIANGULAR INEQUALITY OPTIMIZATION IMPACT\n")
cat("=========================================================\n")
cat("\n")

# Test parameters
OMEGA_POINTS <- 30
T_POINTS <- 30
N_MC_SAMPLES <- 10000
MU_VMF <- c(0, 0, 1)
KAPPA_VMF <- 2.0
N_CORES <- 10

# Test both distance types
for (DISTANCE_TYPE in c("chordal", "geodesic")) {
  
  cat("\n")
  cat("=========================================================\n")
  cat("  Testing distance type:", DISTANCE_TYPE, "\n")
  cat("=========================================================\n")
  cat("\n")
  
  # Generate grid
  omega_grid <- generate_canonical_lattice(OMEGA_POINTS)
  if (DISTANCE_TYPE == "chordal") {
    t_grid <- seq(0, 2, length.out = T_POINTS)
  } else {
    t_grid <- seq(0, pi, length.out = T_POINTS)
  }
  
  n_omega <- nrow(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  q <- 3
  
  cat("Grid size:", n_total, "(", n_omega, "omegas x", n_t, "t values)\n")
  cat("MC samples:", N_MC_SAMPLES, "\n")
  cat("Cores:", N_CORES, "\n\n")
  
  # Pre-compute constants
  A_q_kappa <- besselI(KAPPA_VMF, nu = (q + 1) / 2, expon.scaled = TRUE) / 
               besselI(KAPPA_VMF, nu = (q - 1) / 2, expon.scaled = TRUE)
  var_X <- diag(q) - A_q_kappa^2 * diag(q) - (A_q_kappa / KAPPA_VMF) * outer(MU_VMF, MU_VMF)
  
  # Generate MC samples
  cat("Generating", N_MC_SAMPLES, "MC samples...\n")
  mc_samples <- rotasym::r_vMF(n = N_MC_SAMPLES, mu = MU_VMF, kappa = KAPPA_VMF)
  cat("Done.\n\n")
  
  # Test WITH optimization (current version)
  cat("--- Testing WITH triangular inequality optimization ---\n")
  
  start_time_with <- Sys.time()
  
  # Just compute first 100 entries as a sample
  n_test <- 100
  count_skipped <- 0
  
  for (i in 1:n_test) {
    omega1_idx <- ((i - 1) %% n_omega) + 1
    t1_idx <- ((i - 1) %/% n_omega) + 1
    omega1 <- omega_grid[omega1_idx, ]
    t1 <- t_grid[t1_idx]
    
    for (j in 1:n_test) {
      omega2_idx <- ((j - 1) %% n_omega) + 1
      t2_idx <- ((j - 1) %/% n_omega) + 1
      omega2 <- omega_grid[omega2_idx, ]
      t2 <- t_grid[t2_idx]
      
      # Check if optimization would skip
      d_omega <- sphere_distance(omega1, omega2, DISTANCE_TYPE)
      if (d_omega > t1 + t2) {
        count_skipped <- count_skipped + 1
      }
      
      # Compute (with optimization inside the function)
      P_joint <- compute_joint_probability_vmf(omega1, t1, omega2, t2, MU_VMF, KAPPA_VMF,
                                              DISTANCE_TYPE, mc_samples)
    }
  }
  
  end_time_with <- Sys.time()
  time_with <- as.numeric(difftime(end_time_with, start_time_with, units = "secs"))
  
  cat("Time WITH optimization:", round(time_with, 2), "seconds\n")
  cat("Entries skipped:", count_skipped, "out of", n_test * n_test, 
      sprintf("(%.1f%%)\n", 100 * count_skipped / (n_test * n_test)))
  
  # Estimate total time for full matrix
  total_entries <- n_total * n_total
  estimated_full_time <- time_with * total_entries / (n_test * n_test)
  cat("Estimated time for full", n_total, "x", n_total, "matrix:", 
      round(estimated_full_time / 60, 1), "minutes\n\n")
  
  # Test WITHOUT optimization
  cat("--- Testing WITHOUT triangular inequality optimization ---\n")
  
  
  start_time_without <- Sys.time()
  
  for (i in 1:n_test) {
    omega1_idx <- ((i - 1) %% n_omega) + 1
    t1_idx <- ((i - 1) %/% n_omega) + 1
    omega1 <- omega_grid[omega1_idx, ]
    t1 <- t_grid[t1_idx]
    
    for (j in 1:n_test) {
      omega2_idx <- ((j - 1) %% n_omega) + 1
      t2_idx <- ((j - 1) %/% n_omega) + 1
      omega2 <- omega_grid[omega2_idx, ]
      t2 <- t_grid[t2_idx]
      
      # Compute WITHOUT optimization
      P_joint <- compute_joint_probability_vmf_no_opt(omega1, t1, omega2, t2, MU_VMF, KAPPA_VMF,
                                                      DISTANCE_TYPE, mc_samples)
    }
  }
  
  end_time_without <- Sys.time()
  time_without <- as.numeric(difftime(end_time_without, start_time_without, units = "secs"))
  
  cat("Time WITHOUT optimization:", round(time_without, 2), "seconds\n")
  
  # Estimate total time for full matrix
  estimated_full_time_no_opt <- time_without * total_entries / (n_test * n_test)
  cat("Estimated time for full", n_total, "x", n_total, "matrix:", 
      round(estimated_full_time_no_opt / 60, 1), "minutes\n\n")
  
  # Summary
  speedup <- time_without / time_with
  time_saved <- estimated_full_time_no_opt - estimated_full_time
  
  cat("=== SUMMARY ===\n")
  cat("Speedup factor:", round(speedup, 2), "x\n")
  cat("Time saved (estimated for full matrix):", round(time_saved / 60, 1), "minutes\n")
  cat("Percentage improvement:", round(100 * (1 - 1/speedup), 1), "%\n")
  cat("\n")
}

cat("\n")
cat("=========================================================\n")
cat("  TEST COMPLETED\n")
cat("=========================================================\n")
