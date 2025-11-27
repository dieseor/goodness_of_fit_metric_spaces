# Generate HvMF test data for visualization
# Load the HvMF functions
source("hvmf_distance_profile_analysis.R")

# Set parameters
set.seed(123)
n_samples <- 1000
kappa <- 10  # Use moderate concentration for better visualization
q <- 2  # H^2 hyperboloid

# Define mu (mean direction on hyperboloid)
mu_spatial <- c(1/sqrt(2), 1/sqrt(2))
mu <- c(cosh(0.5), sinh(0.5) * mu_spatial[1], sinh(0.5) * mu_spatial[2])

cat("Generating HvMF samples...\n")
cat("Parameters:\n")
cat("  n =", n_samples, "\n")
cat("  κ =", kappa, "\n")
cat("  μ =", paste(round(mu, 4), collapse = ", "), "\n")

# Verify mu is on hyperboloid
constraint <- -mu[1]^2 + mu[2]^2 + mu[3]^2
cat("  μ constraint (-x1² + x2² + x3²) =", round(constraint, 8), "\n")

# Generate samples
samples <- generate_hvmf_samples(n_samples, mu, kappa)

cat("Generated", nrow(samples), "samples\n")

# Check constraints for generated samples
constraints <- apply(samples, 1, function(x) -x[1]^2 + x[2]^2 + x[3]^2)
cat("Sample constraint check:\n")
cat("  Mean constraint:", round(mean(constraints), 8), "\n")
cat("  Max deviation from -1:", round(max(abs(constraints + 1)), 8), "\n")

# Create data frame for export
hvmf_data <- data.frame(
  V1 = samples[, 1],
  V2 = samples[, 2], 
  V3 = samples[, 3]
)

# Add mu as the first row for reference
mu_data <- data.frame(V1 = mu[1], V2 = mu[2], V3 = mu[3])
hvmf_data_with_mu <- rbind(mu_data, hvmf_data)
hvmf_data_with_mu$type <- c("mu", rep("sample", n_samples))

# Save to CSV
output_file <- "hvmf_test_data.csv"
write.csv(hvmf_data_with_mu, output_file, row.names = FALSE)

cat("Data saved to:", output_file, "\n")
cat("First few rows:\n")
print(head(hvmf_data_with_mu))