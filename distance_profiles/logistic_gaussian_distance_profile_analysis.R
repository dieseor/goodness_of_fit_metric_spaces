# Logistic Gaussian Distance Profile Analysis on the Simplex
# Implementation for validating distance profiles on the simplex using Aitchison distance

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(gridExtra)
  library(Matrix)
  library(MASS)        # For mvrnorm
  library(pracma)      # For special functions
})

# Load the main framework
source(file.path("distance_profiles", "distance_profile_analysis.R"))

# Logistic Gaussian and Aitchison Distance Functions

#' Softmax transformation
#' @param y Vector in R^q
#' @return Vector on simplex Delta^{q-1}
softmax <- function(y) {
  exp_y <- exp(y)  # Exact formula as in the theory
  return(exp_y / sum(exp_y))
}

#' Centered log-ratio transform
#' @param x Vector on simplex Delta^{q-1}
#' @return Vector in subspace S = {z : 1^T z = 0}
clr <- function(x) {
  log_x <- log(x)
  mean_log <- mean(log_x)
  return(log_x - mean_log)
}

#' Aitchison distance between two simplex vectors
#' @param x First simplex vector
#' @param omega Second simplex vector (reference)
#' @return Aitchison distance
aitchison_distance <- function(x, omega) {
  clr_x <- clr(x)
  clr_omega <- clr(omega)
  return(sqrt(sum((clr_x - clr_omega)^2)))
}

#' Generate samples from Logistic Gaussian distribution
#' @param n Number of samples
#' @param mu Mean vector in R^q
#' @param Sigma Covariance matrix (q x q SPD)
#' @return Matrix with samples on simplex (n x q)
generate_logistic_gaussian_samples <- function(n, mu, Sigma) {
  # Generate multivariate normal samples
  Y <- mvrnorm(n, mu, Sigma)
  
  # Apply softmax transformation row-wise
  X <- t(apply(Y, 1, softmax))
  
  return(X)
}

#' Compute Aitchison distances from omega to all samples
#' @param omega Reference simplex vector
#' @param data Matrix with simplex samples (n x q)
#' @return Vector of Aitchison distances
compute_aitchison_distances <- function(omega, data) {
  n <- nrow(data)
  distances <- numeric(n)
  
  for (i in 1:n) {
    distances[i] <- aitchison_distance(data[i, ], omega)
  }
  
  return(distances)
}

#' Theoretical distance profile for Logistic Gaussian distribution
#' @param omega Reference simplex vector
#' @param mu Mean vector of Y ~ N(mu, Sigma)
#' @param Sigma Covariance matrix of Y ~ N(mu, Sigma)
#' @param t_values Vector of t values to evaluate
#' @return Vector of theoretical probabilities
theoretical_distance_profile_logistic_gaussian <- function(omega, mu, Sigma, t_values) {
  q <- length(mu)
  
  # Centering matrix H = I - (1/q) * 1 * 1^T
  H <- diag(q) - (1/q) * matrix(1, q, q)
  
  # Mean and covariance of clr(X)
  clr_omega <- clr(omega)
  m <- H %*% mu - clr_omega  # Mean of Z = H*Y - clr(omega)
  C <- H %*% Sigma %*% H  # Covariance of Z
  
  # Eigendecomposition of C
  eigen_C <- eigen(C, symmetric = TRUE)
  lambda <- eigen_C$values
  U <- eigen_C$vectors
  
  # Keep only positive eigenvalues (numerical tolerance)
  tol <- 1e-10
  positive_idx <- lambda > tol
  lambda_pos <- lambda[positive_idx]
  U_pos <- U[, positive_idx, drop = FALSE]
  r <- sum(positive_idx)  # Effective rank
  
  if (r == 0) {
    stop("Degenerate covariance: all eigenvalues are (near) zero.")
  }
  
  # Transformed mean: mu_prime = U^T * m
  nu <- as.vector(t(U_pos) %*% m)
  
  # Noncentrality parameters: delta_i = (mu_prime_i)^2 / lambda_i
  delta <- nu^2 / lambda_pos
  
  # Compute distance profile using weighted chi-squared distribution
  sapply(t_values, function(t) {
    # We need P(d_A <= t) = P(d_A^2 <= t^2)
    threshold <- t^2
    
    # Use sphunif for weighted sum of noncentral chi-squared variables
    # sum(lambda_pos[j] * chi^2_1(delta[j]))
    sphunif::p_wschisq(threshold, weights = lambda_pos, df = rep(1, r), ncp = delta)
  })
}

#' Special case for q=2: analytical distance profile
#' @param omega Reference simplex vector (length 2)
#' @param mu Mean vector of Y ~ N(mu, Sigma) (length 2)
#' @param Sigma Covariance matrix of Y ~ N(mu, Sigma) (2x2)
#' @param t_values Vector of t values to evaluate
#' @return Vector of theoretical probabilities

#' Run comprehensive Logistic Gaussian distance profile analysis
#' @param output_suffix Suffix for output files
run_logistic_gaussian_distance_profile_analysis <- function(output_suffix = "") {
  set.seed(42)  # For reproducibility
  
  cat("=== Logistic Gaussian Distance Profile Analysis ===\n\n")
  
  # Set parameters for Logistic Gaussian distribution
  q <- 3  # Dimension (3-simplex)
  
  # Mean vector mu
  mu <- c(1, 0.5, -0.5)
  
  # Covariance matrix Sigma (positive definite)
  Sigma <- matrix(c(1.0, 0.3, 0.1,
                    0.3, 0.8, 0.2,
                    0.1, 0.2, 0.6), 3, 3)
  
  cat("Distribution parameters:\n")
  cat("Simplex dimension q =", q, "\n")
  cat("Mean vector μ =", paste(round(mu, 3), collapse = ", "), "\n")
  cat("Covariance matrix Σ:\n")
  print(round(Sigma, 3))
  cat("Distance type: Aitchison\n\n")
  
  # Define three different omega values on the simplex
  omega_values <- list()
  
  # Omega 1: Uniform on simplex
  omega_values[[1]] <- c(1/3, 1/3, 1/3)
  
  # Omega 2: Skewed towards first component
  omega_values[[2]] <- c(0.6, 0.3, 0.1)
  
  # Omega 3: Skewed towards second component
  omega_values[[3]] <- c(0.2, 0.7, 0.1)
  
  cat("Omega values (simplex vectors):\n")
  for (i in seq_along(omega_values)) {
    cat("ω", i, "=", paste(round(omega_values[[i]], 3), collapse = ", "), "\n")
    # Verify simplex constraint
    omega <- omega_values[[i]]
    sum_check <- sum(omega)
    min_check <- min(omega)
    cat("Sum:", round(sum_check, 6), "(should be 1)")
    cat(", Min:", round(min_check, 6), "(should be >= 0)\n")
  }
  cat("\n")
  
  # Create data generator function
  data_generator <- function(N) {
    generate_logistic_gaussian_samples(N, mu, Sigma)
  }
  
  # Set reasonable t_max for Aitchison distances
  t_max <- 3  # Aitchison distances are typically moderate
  
  # Create theoretical profile function
  theoretical_profile <- function(omega, t_values) {
    theoretical_distance_profile_logistic_gaussian(omega, mu, Sigma, t_values)
  }
  
  output_dir <- file.path("output", paste0("logistic_gaussian", output_suffix))
  file_prefix <- "lg_aitchison_dp"
  
  # Set legend positions
  legend_positions <- c("bottom_right", "bottom_right", "top_left")
  
  plots <- create_distance_profile_analysis(
    omega_values = omega_values,
    data_generator = data_generator,
    theoretical_profile = theoretical_profile,
    distance = aitchison_distance,
    sample_sizes = c(50, 200),
    n_simulations = 10,
    t_max = t_max,
    save_plots = TRUE,
    output_dir = output_dir,
    file_prefix = file_prefix,
    legend_positions = legend_positions
  )
    
  cat("\nLogistic Gaussian Aitchison distance analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(plots)
}

cat("Logistic Gaussian analysis functions loaded successfully!\n")
cat("Main function: run_logistic_gaussian_distance_profile_analysis()\n")
cat("Test function: test_logistic_gaussian_sampling()\n")

#' Test function to verify Logistic Gaussian sampling and distance computation
test_logistic_gaussian_sampling <- function(N = 10) {
  cat("Testing Logistic Gaussian sampling with n =", N, "...\n")
  
  # Set up test parameters
  q <- 3
  mu <- c(1, 0, -0.5)
  Sigma <- diag(q)
  
  # Generate samples
  tryCatch({
    samples <- generate_logistic_gaussian_samples(N, mu, Sigma)
    cat("✓ Successfully generated", N, "Logistic Gaussian samples\n")
    cat("Sample matrix dimensions:", dim(samples), "\n")
    
    # Check simplex constraints
    sums <- rowSums(samples)
    mins <- apply(samples, 1, min)
    cat("Row sums (should be 1):", round(range(sums), 6), "\n")
    cat("Minimum values (should be >= 0):", round(range(mins), 6), "\n")
    
    # Test distance computation
    cat("\nTesting distance computations...\n")
    omega <- c(1/3, 1/3, 1/3)
    
    distances <- compute_aitchison_distances(omega, samples)
    cat("✓ Distance computation successful\n")
    cat("Sample distances:", round(head(distances, 5), 3), "...\n")
    
    # Test theoretical computation
    cat("\nTesting theoretical distance profile...\n")
    t_vals <- c(0.5, 1.0, 1.5)
    theory <- theoretical_distance_profile_logistic_gaussian(omega, mu, Sigma, t_vals)
    cat("Theoretical values at t =", t_vals, ":", round(theory, 3), "\n")
    
    return(list(samples = samples, distances = distances, theory = theory))
    
  }, error = function(e) {
    cat("✗ Error in Logistic Gaussian analysis:", e$message, "\n")
    return(NULL)
  })
}

# Run the analysis
run_logistic_gaussian_distance_profile_analysis()

# Test basic functionality  
# test_logistic_gaussian_sampling(10)
# test_logistic_gaussian_q2()
