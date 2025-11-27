# Compare theoretical Var(X) formula with empirical sample variance
# for vMF distribution

library(rotasym)

# Parameters
MU <- c(-1/sqrt(2), 1/sqrt(2), 0)
KAPPA <- 2.0
N_SAMPLES <- 1000000
q <- 2  # Sphere dimension (S²)
q_ambient <- q + 1  # Ambient dimension (ℝ³)

samples <- rotasym::r_vMF(n = N_SAMPLES, mu = MU, kappa = KAPPA)

# Compute A_q(κ)
A_q <- besselI(KAPPA, nu = (q + 1) / 2, expon.scaled = TRUE) / 
       besselI(KAPPA, nu = (q - 1) / 2, expon.scaled = TRUE)

# ============================================================================
# EMPIRICAL VARIANCE
# ============================================================================

# Sample mean
sample_mean <- colMeans(samples)
cat("Sample mean:\n")
print(round(sample_mean, 6))

cat("Theoretical mean = A_q(κ) * μ:\n")
print(round(A_q * MU, 6))

# Sample variance-covariance matrix
sample_var <- cov(samples)
cat("Sample Var(X) = cov(samples):\n")
print(round(sample_var, 6))

# ============================================================================
# THEORETICAL FORMULA (CURRENT IMPLEMENTATION)
# ============================================================================


scalar_coef <- 1 - A_q^2 - ((q + 1) * A_q / KAPPA)
var_theory_current <- (A_q / KAPPA) * diag(q_ambient) + scalar_coef * outer(MU, MU)

cat("Theoretical Var(X):\n")
print(round(var_theory_current, 6))

# ============================================================================
# COMPARISON
# ============================================================================

cat("Difference (Theoretical - Sample):\n")
diff_matrix <- var_theory_current - sample_var
print(round(diff_matrix, 6))

cat("Maximum absolute difference:", round(max(abs(diff_matrix)), 6), "\n")
cat("Relative error (max):", round(max(abs(diff_matrix)) / max(abs(sample_var)), 4), "\n\n")