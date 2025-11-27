# Hyperbolic von Mises-Fisher Distance Profile Analysis on the Hyperboloid
# Implementation for validating distance profiles on the unit hyperboloid

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sphunif)
  library(rotasym)
  library(viridis)
  library(gridExtra)
  library(pracma)
  library(gbutils)
  library(parallel)
})

# Load the main framework
source(file.path("R", "distance_profile_analysis.R"))

# HvMF Sampling Functions (from your implementation)

# Density function 
f_u <- function(u, kappa, chi) {
  exp(log(kappa) + kappa - log(2) + u + log1p(-exp(-2*u)) - 
        kappa/2 * (exp(chi-u) + exp(-chi + u)) + 
        log(besselI(kappa * sinh(chi) * sinh(u), nu = 0, expon.scaled = TRUE)) 
  )
}

# Define the CDF function
cdf <- function(ub_ls, kappa, chi) {
  results <- numeric(length(ub_ls))  # Initialize a vector to store results
  # Loop through all upper bounds
  for (i in seq_along(ub_ls)) {
    ub <- ub_ls[i]
    results[i] <- integrate(f_u, lower = 0, upper = ub, kappa = kappa, 
                            chi = chi)$value
  }
  return(results)
}

# create dataframe with the quantile values to interpolate from
calculate_qvals <- function(kappa, chi, upper) {
  base = seq(0, 1-1e-3, by = 0.01)
  quantile_values <- numeric(length(base))
  for (i in seq_along(base)) {
    quantile_values[i] <- gbutils::cdf2quantile(p = base[i], lower = 0, 
                                                upper = upper, cdf = cdf, 
                                                kappa = kappa, chi = chi
    ) 
  }
  df = data.frame(base, quantile_values)
  return(df)
}

# create function that estimates the quantle by interpolating 
# from the values of df
quantile_interpolation <- function(p, df) {
  return(approx(x = df$base, y = df$quantile_values, xout = p)$y)
}

r_u = function(n, quantile_values) {
  aux = runif(n, min=0, max = 0.99)
  return(quantile_interpolation(p = aux, df = quantile_values))
}


# GEN SAMPLES FUNCTION
r_hvmf = function(i, kappa, chi, theta, upper, quantile_values) {
  # quantile_values <- calculate_qvals(kappa = kappa, chi = chi[i], upper = upper)
  u = r_u(1, quantile_values)
  e_w = r_vMF(1, mu = c(cos(theta[i]), sin(theta[i])), kappa = 
                kappa*sinh(chi[i])*sinh(u)
  )
  return(c(cosh(u), sinh(u) * e_w))
} 


#' Generate samples from hyperbolic von Mises-Fisher distribution on hyperboloid
#' @param n Sample size
#' @param mu Mean direction (hyperboloid vector)
#' @param kappa Concentration parameter
#' @return Matrix with samples in rows (each row is a hyperboloid vector)
generate_hvmf_samples <- function(N, mu, kappa) {
  # Parameters for sampling
  upper = 3  # Upper bound for sampling
  
  # Extract spatial direction from the input mu parameter
  # mu should be on the hyperboloid: -mu[1]^2 + sum(mu[2:end]^2) = -1
  q <- length(mu) - 1  # Dimension of hyperboloid
  
  # Convert to hyperboloid coordinates using the input mu's spatial direction
  
  # Compute chi and theta based on mu's spatial direction
  chi <- acosh(mu[1])
  
  quantile_values <- calculate_qvals(kappa = 50, chi = chi, upper = upper)

  chi = rep(chi, N)  # Repeat chi for all samples

  theta <- acos(mu[2] / sinh(chi))
  
  # ---- PARALLELIZE THE SAMPLING BELOW ----
  num_cores <- parallel::detectCores() - 2  # Leave 2 cores free
  cl <- parallel::makeCluster(num_cores)
  
  # Load required libraries and export variables/functions
  parallel::clusterEvalQ(cl, {
    library(rotasym)
    library(pracma)
    library(gbutils)
  })
  
  parallel::clusterExport(cl, varlist = c("f_u", "cdf", "r_u", "quantile_interpolation",
                                          "calculate_qvals", "mu", "r_vMF", "r_hvmf", 
                                          "upper", "N", "kappa", "chi"), envir = environment())
  
  # Parallel loop: each core computes one r_hvmf sample for index i
  hyp_points <- parallel::parLapply(cl, seq_along(chi), function(i) {
    r_hvmf(i, kappa, chi, theta, upper, quantile_values)
  })
  
  parallel::stopCluster(cl)
  
  # Combine and save
  hyp_points <- do.call(rbind, hyp_points)

  return(hyp_points)
  
}

#' Compute Minkowski pseudo-inner product
#' @param x First hyperboloid vector
#' @param y Second hyperboloid vector  
#' @return Minkowski inner product -x1*y1 + sum(xi*yi for i=2:q+1)
minkowski_inner_product <- function(x, y) {
  -x[1] * y[1] + sum(x[2:length(x)] * y[2:length(y)])
}

#' Compute geodesic distance on hyperboloid
#' @param x1 First point (hyperboloid vector)
#' @param x2 Second point (hyperboloid vector)
#' @return Geodesic distance on hyperboloid
hyperbolic_geodesic_distance <- function(x1, x2) {
  inner_prod <- minkowski_inner_product(x1, x2)
  acosh(-inner_prod)
}

#' Project vector to unit hyperboloid
#' @param x Vector in R^{q+1}
#' @return Vector projected to hyperboloid with x[1] > 0
project_to_hyperboloid <- function(x) {
  # For unit hyperboloid: -x1^2 + sum(xi^2 for i=2:q+1) = -1
  # This gives: x1^2 = 1 + sum(xi^2 for i=2:q+1)
  
  spatial_part <- x[2:length(x)]
  spatial_norm_sq <- sum(spatial_part^2)
  
  # Compute required first coordinate
  x1 <- sqrt(1 + spatial_norm_sq)
  
  return(c(x1, spatial_part))
}

#' Compute distances from omega to all points in hyperboloid data
#' @param omega Reference point (hyperboloid vector)
#' @param data Matrix with observations in rows (hyperboloid vectors)
#' @return Vector of geodesic distances on hyperboloid
compute_hyperboloid_distances <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  apply(data, 1, function(x) hyperbolic_geodesic_distance(omega, x))
}

#' Theoretical distance profile for hyperbolic von Mises-Fisher using density of R
#' @param omega Reference point (hyperboloid vector)
#' @param mu Mean direction of HvMF
#' @param kappa Concentration parameter
#' @param t_values Vector of t values to evaluate
#' @return Vector of theoretical probabilities
theoretical_distance_profile_hvmf <- function(omega, mu, kappa, t_values) {
  # Compute mu_omega = (mu, omega) using Minkowski inner product
  mu_omega <- minkowski_inner_product(mu, omega)  # This should be <= -1
  q <- length(mu) - 1  # Dimension of hyperboloid
  
    sapply(t_values, function(t) {
    # For hyperboloid: F(t) = P(R <= t) = ∫[0 to t] f_R(r) dr
    # where f_R(r) is the density of the geodesic distance R
    
      # Density of R for r >= 0 (equation from your derivation):
      # f_R(r) = [c_q^HvMF(κ) * (sinh(r))^{q-1} * exp(κ(μ,ω)cosh(r))] / 
      #          [c_q^vMF(κ*sinh(r)*√((μ,ω)^2-1))]
      
      density_R <- function(r) {
        
        # Compute everything in log scale for numerical stability
        nu <- (q - 1) / 2  # Order of modified Bessel function
        
        # Use expon.scaled = TRUE to avoid underflow
        # besselK(κ, nu, expon.scaled = TRUE) returns e^κ * K_nu(κ)
        # So log(K_nu(κ)) = log(besselK(κ, nu, expon.scaled = TRUE)) - κ
        log_bessel_k <- log(besselK(kappa, nu = nu, expon.scaled = TRUE)) - kappa
        
        # Log numerator: normalization + power term + exponential term
        log_sinh_power_term <- (q - 1) * log(sinh(r))
        log_exp_term <- kappa * mu_omega * cosh(r)  # Note: mu_omega <= -1, so this is negative
        log_numerator <- nu * log(kappa) - nu * log(2 * pi) - log(2) - log_bessel_k + 
                        log_sinh_power_term + log_exp_term
        
        # Log denominator: log(c_q^vMF(κ*sinh(r)*√((μ,ω)^2-1)))
        # Handle numerical precision for sqrt((μ,ω)^2-1)
        mu_omega_sq_minus_1 <- mu_omega^2 - 1
        if (mu_omega_sq_minus_1 < 0) {
          if (abs(mu_omega_sq_minus_1) < 1e-12) {
            # Small numerical error - set to 0
            mu_omega_sq_term <- 0
          } else {
            # Large negative value indicates a real mathematical error
            stop(paste("Mathematical error: mu_omega^2 - 1 =", mu_omega_sq_minus_1, 
                      "unexpected negative value inside square root."))
          }
        } else {
          mu_omega_sq_term <- mu_omega_sq_minus_1
        }
        
        kappa_vMF <- kappa * sinh(r) * sqrt(mu_omega_sq_term)
        # Be careful: do nnot confuse extrinsic/intrinsic dimensions
        log_denominator <- rotasym::c_vMF(p = q, kappa = kappa_vMF, log = TRUE)
        
        # Return exp(log_numerator - log_denominator)
        result <- exp(log_numerator - log_denominator)
        
        return(result)
      }

      # Integrate from 0 to t: F(t) = ∫[0 to t] f_R(r) dr
      cdf_result <- integrate(density_R, 
                             lower = 0, 
                             upper = t, 
                             rel.tol = 1e-6, 
                             abs.tol = 1e-8,
                             stop.on.error = FALSE)
      
      # Return F(t) = P(R <= t)
      return(cdf_result$value)
      
  })
}


#' Run comprehensive hyperbolic von Mises-Fisher analysis
#' @param output_suffix Suffix for output files
run_hvmf_distance_profile_analysis <- function(output_suffix = "") {
    set.seed(1)
  cat("=== Hyperbolic von Mises-Fisher Distance Profile Analysis ===\n\n")
  
  # Set parameters for HvMF distribution
  q <- 2  # dimension of hyperboloid H^2 in R^3
  # Mean direction on hyperboloid - using a pattern consistent with your sampling
  # For H^2, we work with 3D vectors (x1, x2, x3) where -x1^2 + x2^2 + x3^2 = -1
  mu_spatial <- c(1/sqrt(2), 1/sqrt(2))  # 2D direction as in your example
  mu <- c(sqrt(2), mu_spatial[1], mu_spatial[2])  # Point on hyperboloid
  kappa <- 50      # concentration parameter (using values from your example)
  
  cat("Distribution parameters:\n")
  cat("Intrinsic dimension q =", q, "\n")
  cat("Mean direction μ =", paste(round(mu, 3), collapse = ", "), "\n")
  cat("Concentration κ =", kappa, "\n")
  cat("Distance type: geodesic (hyperboloid)\n\n")
  
  # Define three different omega values on the hyperboloid
  # Using different chi and theta values to create diverse reference points
  omega_values <- list()

  omega_values[[1]] <- mu
  # Omega 2: Different chi, same theta
  omega_values[[2]] <- c(sqrt(2), -1/sqrt(2), 1/sqrt(2))
  
  # Omega 3: Same chi, different theta
  omega_values[[3]] <- c(sqrt(2), -1/sqrt(2), -1/sqrt(2))
  
  cat("Omega values (hyperboloid vectors):\n")
  for (i in seq_along(omega_values)) {
    cat("ω", i, "=", paste(round(omega_values[[i]], 3), collapse = ", "), "\n")
    # Verify hyperboloid constraint
    omega <- omega_values[[i]]
    constraint <- -omega[1]^2 + sum(omega[2:length(omega)]^2)
    cat("Minkowski norm:", round(constraint, 6), "\n")
  }
  cat("\n")
  
  # Create data generator function
  data_generator <- function(N) {
    generate_hvmf_samples(N, mu, kappa)
  }
  
  # Set reasonable t_max for hyperboloid
  t_max <- pi - 1e-3  # Geodesic distances on hyperboloid can be large
  
  # Create theoretical profile function
  theoretical_profile <- function(omega, t_values) {
    theoretical_distance_profile_hvmf(omega, mu, kappa, t_values)
  }
  
  output_dir <- file.path("output", paste0("hvmf_geodesic", output_suffix))
  file_prefix <- "hvmf_geodesic_dp"
  
  # Set legend positions (can be customized based on empirical observation)
  # First and Second omega: bottom_right, Third omega: top_left 
  legend_positions <- c("bottom_right", "bottom_right", "top_left")
  
  plots <- create_distance_profile_analysis(
    omega_values = omega_values,
    data_generator = data_generator,
    theoretical_profile = theoretical_profile,
    distance = hyperbolic_geodesic_distance,
    sample_sizes = c(50, 200),
    n_simulations = 10,
    t_max = t_max,
    save_plots = TRUE,
    output_dir = output_dir,
    file_prefix = file_prefix,
    legend_positions = legend_positions
  )
  
  
  cat("\nHvMF geodesic distance analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(plots)
}

cat("Hyperbolic von Mises-Fisher analysis functions loaded successfully!\n")
cat("Main function: run_hvmf_distance_profile_analysis()\n")
cat("Test function: test_hvmf_sampling()\n")

#' Test function to verify HvMF sampling works
test_hvmf_sampling <- function(N = 10) {
  cat("Testing HvMF sampling with n =", N, "...\n")
  
  # Set up test parameters
  mu_spatial <- c(1/sqrt(2), 1/sqrt(2))
  mu <- c(sqrt(2), mu_spatial[1], mu_spatial[2])
  kappa <- 10  # Smaller kappa for testing
  
  # Generate samples
  tryCatch({
    samples <- generate_hvmf_samples(N, mu, kappa)
    cat("✓ Successfully generated", N, "HvMF samples\n")
    
    # Check hyperboloid constraints
    constraints <- apply(samples, 1, function(x) -x[1]^2 + x[2]^2 + x[3]^2)
    cat("Constraint values (should be ≈ -1):\n")
    print(round(constraints, 4))
    
    # Check distances
    cat("\nTesting distance computations...\n")
    omega <- c(cosh(1.0), sinh(1.0) * 0.8, sinh(1.0) * 0.6)
    omega <- project_to_hyperboloid(omega)
    
    distances <- compute_hyperboloid_distances(omega, samples)
    cat("Sample distances:", round(head(distances, 5), 3), "...\n")
    
    return(list(samples = samples, distances = distances))
    
  }, error = function(e) {
    cat("✗ Error in HvMF sampling:", e$message, "\n")
    return(NULL)
  })
}

# Run the analysis (commented out by default)
run_hvmf_distance_profile_analysis()

#test_hvmf_sampling(10)


