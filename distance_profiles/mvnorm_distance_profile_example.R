# Multivariate Normal Distance Profile Analysis Example
# Complete implementation for validating distance profiles

# Load the analysis framework
source(file.path("distance_profiles", "distance_profile_analysis.R"))

# Set random seed for reproducibility
set.seed(12345)


#' Euclidean distance between two points
#' @param x1 First point (vector)
#' @param x2 Second point (vector)
#' @return Euclidean distance
euclidean_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}

#' Compute distances from omega to all points in a dataset
#' @param omega Reference point (vector)
#' @param data Matrix with observations in rows
#' @return Vector of distances
compute_euclidean_distances <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  apply(data, 1, function(x) euclidean_distance(omega, x))
}



#' Theoretical distance profile for multivariate normal distribution
#' @param omega Reference point (vector)
#' @param mu Mean vector of the distribution
#' @param Sigma Covariance matrix of the distribution
#' @param t_values Vector of t values to evaluate
#' @return Vector of theoretical probabilities
theoretical_distance_profile_mvnorm <- function(omega, mu, Sigma, t_values) {
  # Compute Z = X - omega, so Z ~ N(mu - omega, Sigma)
  z_mean <- mu - omega
  
  # Eigendecomposition of Sigma
  eigen_decomp <- eigen(Sigma)
  U <- eigen_decomp$vectors
  lambda <- eigen_decomp$values
  
  # Compute nu = U^T * (mu - omega)
  nu <- as.vector(t(U) %*% z_mean)
  
  # For each t, compute P(||Z||^2 <= t^2)
  sapply(t_values, function(t) {
    if (t <= 0) return(0)
    
    # We need P(sum(T_i^2) <= t^2) where T_i ~ N(nu_i, lambda_i)
    # This is equivalent to a weighted sum of non-central chi-squared variables
    
    # Use sphunif::p_wschisq for weighted sum of chi-squared
    # Convert to sum of chi-squared form: sum(lambda_i * X_i) where X_i ~ chi^2_1(delta_i)
    # with non-centrality parameter delta_i = nu_i^2 / lambda_i
    
    weights <- lambda
    ncp <- nu^2 / lambda  # non-centrality parameters
    
    # Use weighted chi-squared distribution
      sphunif::p_wschisq(t^2, weights = weights, df = rep(1, length(weights)), ncp = ncp)
    })
}

#' Generate samples from multivariate normal distribution
#' @param n Sample size
#' @param mu Mean vector
#' @param Sigma Covariance matrix
#' @return Matrix with samples in rows
generate_mvnorm_samples <- function(n, mu, Sigma) {
  mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
}


#' Example implementation for multivariate normal distribution
run_mvnorm_distance_profile_analysis <- function() {
  
  cat("=== Multivariate Normal Distance Profile Analysis ===\n\n")
  
  # Define parameters for the multivariate normal distribution
  q <- 3  # dimension
  mu <- c(3, -2, 1)  # mean vector

  Sigma <- matrix(c(1, 0.3, 0.3, 1, -2, 1.4, -3.1, 2.1, -2), nrow = 3)  # covariance matrix
  Sigma <- Sigma %*% t(Sigma)

  cat("Distribution parameters:\n")
  cat("Dimension q =", q, "\n")
  cat("Mean μ =", paste(mu, collapse = ", "), "\n")
  cat("Covariance matrix Σ:\n")
  print(Sigma)
  cat("\n")
  
  # Define three different omega values (reduced from 4 for better visualization)
  omega_values <- list(
    c(3, -2, 1),      # At the mean
    c(1, 1, 1),      # Displaced from mean
    c(-1, 0.5, 0.5)    # Another displacement
  )

  
  cat("Omega values:\n")
  for (i in 1:length(omega_values)) {
    cat("ω", i, "=", paste(omega_values[[i]], collapse = ", "), "\n")
  }
  cat("\n")
  
  # Create data generator function
  data_generator <- function(n) {
    generate_mvnorm_samples(n, mu, Sigma)
  }
  
  # Create theoretical profile function
  theoretical_profile <- function(omega, t_values) {
    theoretical_distance_profile_mvnorm(omega, mu, Sigma, t_values)
  }  

  # Run the complete analysis
  cat("Starting analysis...\n")
  cat("This will generate 6 plots (3 omegas × 2 sample sizes)\n\n")
  
  plots <- create_distance_profile_analysis(
    omega_values = omega_values,
    data_generator = data_generator,
    theoretical_profile = theoretical_profile,
    distance = euclidean_distance,
    sample_sizes = c(50, 200),
    n_simulations = 10,
    t_max = 18,  # Reasonable range for this example
    save_plots = TRUE,
    output_dir = "output/mvnorm",
    file_prefix = "mvnorm_dp"
  )
  
  return(plots)
}  

# If running this script directly, execute the analysis
cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")

# Run full analysis
plots <- run_mvnorm_distance_profile_analysis()
