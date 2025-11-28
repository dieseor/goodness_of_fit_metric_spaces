# von Mises-Fisher Distance Profile Analysis on the Sphere
# Implementation for validating distance profiles on the unit sphere

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sphunif)
  library(rotasym)
  library(viridis)
  library(gridExtra)
})

# Load the main framework
source(file.path("distance_profiles", "distance_profile_analysis.R"))
# Load utilities (centralized implementations)
source(file.path("utils.R"))

#' Generate samples from von Mises-Fisher distribution on sphere
#' @param n Sample size
#' @param mu Mean direction (unit vector)
#' @param kappa Concentration parameter
#' @return Matrix with samples in rows (each row is a unit vector)
generate_vmf_samples <- function(n, mu, kappa) {
  rotasym::r_vMF(n, mu = mu, kappa = kappa)
}

#' Compute chordal distance on sphere
#' @param x1 First point (unit vector)
#' @param x2 Second point (unit vector)
#' @return Chordal distance
chordal_distance <- function(x1, x2) {
  sqrt(2 * (1 - sum(x1 * x2)))
}

#' Compute geodesic distance on sphere
#' @param x1 First point (unit vector)
#' @param x2 Second point (unit vector)
#' @return Geodesic distance
geodesic_distance <- function(x1, x2) {
  acos(sum(x1 * x2))
}

#' Compute distances from omega to all points in spherical data
#' @param omega result[valid_idx] <- numerator / denominator
#' @param data Data matrix (rows are the observations)
#' @param distance_type Either "chordal" or "geodesic"
#' @return Vector of distances
compute_chordal_distances <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
    apply(data, 1, function(x) chordal_distance(omega, x))
  }

compute_geodesic_distances <- function(omega, data) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  apply(data, 1, function(x) geodesic_distance(omega, x))
}

## theoretical_distance_profile_vmf moved to R/utils.R (centralized). Use that implementation.

#' Run comprehensive von Mises-Fisher analysis for both distance types
#' @param distance_type Either "chordal" or "geodesic"
#' @param output_suffix Suffix for output files
run_vmf_distance_profile_analysis <- function(distance_type = "chordal", output_suffix = "") {
  
  cat("=== von Mises-Fisher Distance Profile Analysis ===\n")
  cat("Distance type:", distance_type, "\n\n")
  
  # Set parameters for vMF distribution
  q <- 3  # dimension (unit sphere S^2 in R^3)
  mu <- c(0, 0, 1)  # mean direction (north pole)
  kappa <- 2.0      # concentration parameter
  
  cat("Distribution parameters:\n")
  cat("Dimension q =", q, "\n")
  cat("Mean direction μ =", paste(mu, collapse = ", "), "\n")
  cat("Concentration κ =", kappa, "\n")
  cat("Distance type:", distance_type, "\n\n")
  
  # Define three different omega values on the sphere
  omega_values <- list(
    c(0, 0, 1),      # North pole (same as mean)
    c(1, 0, 0),      # On equator
    c(0, 0, -1)      # South pole (opposite to mean)
  )
  
  # Normalize to ensure unit vectors
  omega_values <- lapply(omega_values, function(w) w / sqrt(sum(w^2)))
  
  cat("Omega values (unit vectors):\n")
  for (i in 1:length(omega_values)) {
    cat("ω", i, "=", paste(round(omega_values[[i]], 3), collapse = ", "), "\n")
  }
  cat("\n")
  
  # Create data generator function
  data_generator <- function(n) {
    generate_vmf_samples(n, mu, kappa)
  }

    t_max <- ifelse(distance_type == "chordal", 2 - 1e-8, pi - 1e-8)

    # Create theoretical profile function
    theoretical_profile <- function(omega, t_values) {
      theoretical_distance_profile_vmf(omega, mu, kappa, t_values, distance_type)
    }
  

  output_dir <- file.path("output", paste0("vmf_", distance_type, output_suffix))
  file_prefix <- paste0("vmf_", distance_type, "_dp")
  
  # Set legend positions based on omega values and distance type
  legend_positions <- c("bottom_right", "bottom_right", "bottom_right")  # Default
  
  if (distance_type == "chordal") {
    # For chordal: omega = (1, 0, 0) -> top_left, omega = (0, 0, -1) -> top_left
    legend_positions <- c("bottom_right", "top_left", "top_left")
  } else if (distance_type == "geodesic") {
    # For geodesic: omega = (1, 0, 0) -> top_left (note: you mentioned (-1,0,0) but omega_2 is (1,0,0))
    legend_positions <- c("bottom_right", "top_left", "bottom_right")
  }
  


if (distance_type == "chordal") {
  distance = compute_chordal_distances
} else if (distance_type == "geodesic") {
  distance = compute_geodesic_distances
} else {
  stop("distance_type must be either 'chordal' or 'geodesic'")
}

  plots <- create_distance_profile_analysis(
    omega_values = omega_values,
    data_generator = data_generator,
    theoretical_profile = theoretical_profile,
    distance = distance,
    sample_sizes = c(50, 200),
    n_simulations = 10,
    t_max = t_max,
    save_plots = TRUE,
    output_dir = output_dir,
    file_prefix = file_prefix,
    legend_positions = legend_positions
  )


  cat("\nvMF", distance_type, "distance analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(plots)
}

#' Run complete comparison of multivariate normal, von Mises-Fisher, and hyperbolic von Mises-Fisher
run_complete_distance_profile_comparison <- function() {
  
  cat("========================================\n")
  cat("COMPLETE DISTANCE PROFILE ANALYSIS\n")
  cat("========================================\n\n")
  
  # 1. Multivariate Normal Analysis
  cat("1. MULTIVARIATE NORMAL ANALYSIS\n")
  cat("--------------------------------\n")
  source("mvnorm_distance_profile_example.R")
  mvn_plots <- run_mvnorm_distance_profile_analysis()
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")
  
  # 2. von Mises-Fisher Chordal Distance Analysis
  cat("2. VON MISES-FISHER (CHORDAL DISTANCE)\n")
  cat("-------------------------------------\n")
  vmf_chordal_plots <- run_vmf_distance_profile_analysis("chordal")
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")
  
  # 3. von Mises-Fisher Geodesic Distance Analysis
  cat("3. VON MISES-FISHER (GEODESIC DISTANCE)\n")
  cat("--------------------------------------\n")
  vmf_geodesic_plots <- run_vmf_distance_profile_analysis("geodesic")
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")
  
  # 4. Hyperbolic von Mises-Fisher Analysis
  cat("4. HYPERBOLIC VON MISES-FISHER (GEODESIC DISTANCE)\n")
  cat("-------------------------------------------------\n")
  source("hvmf_distance_profile_analysis.R")
  hvmf_plots <- run_hvmf_distance_profile_analysis()
  
  cat("\n========================================\n")
  cat("ALL ANALYSES COMPLETE!\n")
  cat("========================================\n")
  cat("Check the following directories for results:\n")
  cat("- output/mvnorm/ (Multivariate Normal)\n")
  cat("- output/vmf_chordal/ (vMF Chordal Distance)\n")  
  cat("- output/vmf_geodesic/ (vMF Geodesic Distance)\n")
  cat("- output/hvmf_geodesic/ (HvMF Geodesic Distance)\n")
  
  return(list(
    mvn = mvn_plots,
    vmf_chordal = vmf_chordal_plots,
    vmf_geodesic = vmf_geodesic_plots,
    hvmf_geodesic = hvmf_plots
  ))
}

# For chordal distance analysis:
run_vmf_distance_profile_analysis("chordal")

# For geodesic distance analysis:  
run_vmf_distance_profile_analysis("geodesic")

# Or run both analyses:
# all_plots <- run_complete_distance_profile_comparison()
