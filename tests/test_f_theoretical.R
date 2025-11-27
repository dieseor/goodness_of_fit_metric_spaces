# Sanity check: verify F_theoretical values make sense

library(rotasym)

source(file.path("R", "utils.R"))
source(file.path("R", "gaussian_process_vmf.R"))

omega_grid <- generate_canonical_lattice(10)
t_grid <- seq(0 + 1e-8, 2 - 1e-8, length.out = 10)
mu <- c(0, 0, 1)
kappa <- 2.0

# Compute F_theoretical
F_theoretical_matrix <- matrix(0, nrow = 10, ncol = 10)
for (i in 1:10) {
  F_theoretical_matrix[i, ] <- compute_distance_profile_vmf(
    omega_grid[i, ], mu, kappa, t_grid, "chordal"
  )
}

cat("F_theoretical matrix:\n")
cat("First row (omega at north pole):\n")
print(round(F_theoretical_matrix[1, ], 4))
cat("\nLast row (omega somewhere else):\n")
print(round(F_theoretical_matrix[10, ], 4))

cat("\nExpected behavior:\n")
cat("- For t close to 0: F should be ~0\n")
cat("- For t close to 2: F should be ~1\n")
cat("- F should be monotone increasing in t\n")

cat("\nChecking monotonicity for first row:\n")
diffs <- diff(F_theoretical_matrix[1, ])
cat("All differences positive?", all(diffs > -1e-10), "\n")
cat("Min diff:", min(diffs), "Max diff:", max(diffs), "\n")
