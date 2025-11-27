# Test vectorization and edge case fixes for vMF Gaussian process

library(rotasym)
source("R/gaussian_process_vmf.R")
source("R/utils.R")

cat("=== Testing vMF Edge Cases and Vectorization ===\n\n")

# Setup
mu <- c(0, 0, 1)
kappa <- 2.0

# ============================================================================
# TEST 1: Edge cases in distance profile
# ============================================================================
cat("TEST 1: Edge case handling for distance profiles\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

omega <- c(0, 0, 1)  # Same as mu

# Test chordal distance
cat("\nChordal distance:\n")
t_values_chordal <- c(1e-5, 0.01, 0.5, 1.0, 1.5, 1.99, 2.0 - 1e-6)
F_chordal <- theoretical_distance_profile_vmf(omega, mu, kappa, t_values_chordal, "chordal")
for (i in seq_along(t_values_chordal)) {
  cat(sprintf("  t = %.6f: F = %.6f\n", t_values_chordal[i], F_chordal[i]))
}

# Test geodesic distance  
cat("\nGeodesic distance:\n")
t_values_geodesic <- c(1e-5, 0.01, 0.5, 1.0, 2.0, 3.0, pi - 1e-6)
F_geodesic <- theoretical_distance_profile_vmf(omega, mu, kappa, t_values_geodesic, "geodesic")
for (i in seq_along(t_values_geodesic)) {
  cat(sprintf("  t = %.6f: F = %.6f\n", t_values_geodesic[i], F_geodesic[i]))
}

# Check edge cases are correct
cat("\nEdge case checks:\n")
if (F_chordal[1] < 0.01) {
  cat("  ✓ Chordal: F(t≈0) ≈ 0\n")
} else {
  cat("  ✗ FAIL: Chordal F(t≈0) = ", F_chordal[1], " should be ≈ 0\n")
}
if (F_chordal[length(F_chordal)] > 0.99) {
  cat("  ✓ Chordal: F(t≈2) ≈ 1\n")
} else {
  cat("  ✗ FAIL: Chordal F(t≈2) = ", F_chordal[length(F_chordal)], " should be ≈ 1\n")
}
if (F_geodesic[1] < 0.01) {
  cat("  ✓ Geodesic: F(t≈0) ≈ 0\n")
} else {
  cat("  ✗ FAIL: Geodesic F(t≈0) = ", F_geodesic[1], " should be ≈ 0\n")
}
if (F_geodesic[length(F_geodesic)] > 0.99) {
  cat("  ✓ Geodesic: F(t≈π) ≈ 1\n")
} else {
  cat("  ✗ FAIL: Geodesic F(t≈π) = ", F_geodesic[length(F_geodesic)], " should be ≈ 1\n")
}

# ============================================================================
# TEST 2: Vectorization correctness (small grid)
# ============================================================================
cat("\n\nTEST 2: Vectorization correctness\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Small grid for testing
omega_grid <- generate_canonical_lattice(5)  # 5 omega points
t_grid <- seq(0.1, 1.0, length.out = 5)      # 5 t points
n_mc_samples <- 5000

cat("Grid size: 5 omegas × 5 t values = 25 entries\n")
cat("MC samples:", n_mc_samples, "\n\n")

# Generate MC samples
set.seed(42)
mc_samples <- rotasym::r_vMF(n = n_mc_samples, mu = mu, kappa = kappa)

# Pre-compute constants
q_sphere <- length(mu) - 1
A_q_kappa <- besselI(kappa, nu = (q_sphere + 1) / 2, expon.scaled = TRUE) / 
             besselI(kappa, nu = (q_sphere - 1) / 2, expon.scaled = TRUE)
q_ambient <- length(mu)
scalar_coef <- 1 - A_q_kappa^2 - ((q_ambient) * A_q_kappa / kappa)
var_X <- (A_q_kappa / kappa) * diag(q_ambient) + scalar_coef * outer(mu, mu)

# Test ONE row with vectorized function
cat("Computing row 1 with vectorized function...\n")
start_time <- Sys.time()
row_vectorized <- row_cov_vmf(
  h0 = 'simple', unknown_param = NULL, idx = 1, 
  omega_grid = omega_grid, 
  t_grid = t_grid, 
  mu = mu, 
  kappa = kappa,
  distance_type = "chordal", 
  mc_samples = mc_samples,
  A_q_kappa = A_q_kappa, 
  var_X = var_X
)
time_vectorized <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("  Time:", round(time_vectorized, 3), "seconds\n")
cat("  Row length:", length(row_vectorized), "\n")
cat("  Row range: [", round(min(row_vectorized), 6), ",", round(max(row_vectorized), 6), "]\n")
cat("  Row mean:", round(mean(row_vectorized), 6), "\n\n")

# Verify it's the correct length
if (length(row_vectorized) == 25) {
  cat("  ✓ Row has correct length (25)\n")
} else {
  cat("  ✗ FAIL: Row length is", length(row_vectorized), "instead of 25\n")
}

# Check for NaN/Inf
if (all(is.finite(row_vectorized))) {
  cat("  ✓ All values are finite\n")
} else {
  cat("  ✗ FAIL: Some values are NaN or Inf\n")
}

# ============================================================================
# TEST 3: Spot check specific covariance values
# ============================================================================
cat("\n\nTEST 3: Spot check covariance values\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Check diagonal entry (should be positive)
diag_entry <- row_vectorized[1]  # Covariance of (omega1, t1) with itself
cat("Diagonal entry (var): ", round(diag_entry, 6), "\n")
if (diag_entry > 0) {
  cat("  ✓ Diagonal is positive\n")
} else {
  cat("  ✗ FAIL: Diagonal should be positive\n")
}

# Check symmetry property (compute entry (2,1) and compare with (1,2))
cat("\nComputing row 2 to check symmetry...\n")
row2_vectorized <- row_cov_vmf(
  h0 = 'simple', unknown_param = NULL, idx = 2, 
  omega_grid = omega_grid, 
  t_grid = t_grid, 
  mu = mu, 
  kappa = kappa,
  distance_type = "chordal", 
  mc_samples = mc_samples,
  A_q_kappa = A_q_kappa, 
  var_X = var_X
)

entry_1_2 <- row_vectorized[2]  # (1,2)
entry_2_1 <- row2_vectorized[1]  # (2,1)
diff <- abs(entry_1_2 - entry_2_1)
cat("  Entry (1,2):", round(entry_1_2, 8), "\n")
cat("  Entry (2,1):", round(entry_2_1, 8), "\n")
cat("  Difference:", round(diff, 10), "\n")

if (diff < 1e-6) {
  cat("  ✓ Matrix is symmetric\n")
} else {
  cat("  ⚠ Warning: Difference is", diff, "(Monte Carlo variation)\n")
}

cat("\n=== All tests completed ===\n")
