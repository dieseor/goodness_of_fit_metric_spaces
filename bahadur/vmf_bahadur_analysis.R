# Bahadur Representation Analysis for vMF Distribution
# Validation of the MLE estimator asymptotic behavior
# Based on theoretical derivation for vMF(mu, kappa) with xi = mu * kappa

# Load required libraries
suppressPackageStartupMessages({
  library(movMF)       # For vMF MLE estimation
  library(sphunif)     # For vMF sampling and utilities
  library(ggplot2)     # For plotting
  library(dplyr)       # For data manipulation
  library(viridis)     # For colors
  library(pracma)      # For special functions (bessel)
  library(rotasym)
  library(parallel)    # For parallel processing
})

# Load utils and Bahadur analysis functions
source(file.path("utils.R"))
source("bahadur/bahadur_functions.R")

#' Run complete Bahadur representation analysis
#' @param mu_true True mean direction (unit vector on S^q)
#' @param kappa_true True concentration parameter (default: 1)
#' @param n_trajectories Number of trajectories to simulate (default: 20)
#' @param n_cores Number of CPU cores to use for parallel processing (default: 10)
#' @param min_n Minimum sample size (default: 50) - currently fixed
#' @param max_n Maximum sample size (default: 100000) - currently fixed
#' @param plot_suffix Suffix for plot filename (default: "" for default name)
#' @param seed Random seed for reproducibility (default: 123)
run_bahadur_analysis <- function(mu_true, kappa_true = 1, n_trajectories = 20, 
                                n_cores = 10, min_n = 50, max_n = 100000, plot_suffix = "", seed = 123) {
  
  cat("=== Bahadur Representation Analysis for vMF Distribution ===\n\n")
  
  q <- length(mu_true) - 1
  xi_true <- kappa_true * mu_true
  
  cat("True parameters:\n")
  cat("Sphere dimension: S^", q, "\n", sep = "")
  cat("True mu:", paste(round(mu_true, 4), collapse = ", "), "\n")
  cat("True kappa:", kappa_true, "\n")
  cat("True xi:", paste(round(xi_true, 4), collapse = ", "), "\n\n")
  
  # Create sample sizes with variable step sizes
  sample_sizes <- c(
    seq(50, 1000, by = 100),       # Step 100 until n=1000 (includes 1000)
    seq(1250, 5000, by = 250),     # Step 250 until n=5000 
    seq(5500, 10000, by = 500),    # Step 500 until n=10000  
    seq(11000, 50000, by = 1000),  # Step 1000 until n=50000
    seq(55000, 100000, by = 5000)  # Step 5000 until n=100000
  )
  
  cat("Total sample points:", length(sample_sizes), "\n")
  cat("Number of trajectories:", n_trajectories, "\n")
  cat("Number of cores:", n_cores, "\n")  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Run analysis for all trajectories in parallel
  cat("Setting up parallel cluster with", n_cores, "cores...\n")
  cl <- makeCluster(n_cores)
  
  # Use tryCatch for robust cluster management
  trajectory_results <- tryCatch({
    # Load required libraries and source script on each worker
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(movMF)
        library(sphunif)
        library(pracma)
        library(rotasym)
      })
      # Source utils and functions file so compute_mle_xi is available on workers
      source(file.path("R", "utils.R"))
      source("bahadur_functions.R")
    })
    
    cat("Running", n_trajectories, "trajectories in parallel...\n")
    
    # Run trajectories in parallel with direct parameter passing and deterministic seeds
    results <- parLapply(cl, 1:n_trajectories, function(traj, mu, kappa, sizes, base_seed) {
      # Set deterministic seed for this trajectory
      set.seed(base_seed + traj)
      analyze_single_trajectory(mu, kappa, sizes, traj)
    }, mu_true, kappa_true, sample_sizes, seed)
        
    results
  }, finally = {
    # Always stop cluster, even if error occurs
    if (exists("cl") && inherits(cl, "cluster")) {
      stopCluster(cl)
      cat("Cluster stopped successfully.\n")
    }
  })
  
  # Combine results (trajectory_results is a list of dataframes)
  all_results <- do.call(rbind, trajectory_results)
  
  # Check for problematic values - STOP if any issues found
  cat("\n=== Data Quality Check ===\n")
  cat("Total observations:", nrow(all_results), "\n")
  
  # Check each column
  na_diff <- sum(is.na(all_results$difference_norm))
  inf_diff <- sum(is.infinite(all_results$difference_norm))
  
  cat("difference_norm - NA count:", na_diff, "\n")
  cat("difference_norm - Inf count:", inf_diff, "\n")
  
  if (na_diff > 0 || inf_diff > 0) {
    cat("FOUND PROBLEMATIC VALUES!\n")
    problematic_rows <- which(is.na(all_results$difference_norm) | is.infinite(all_results$difference_norm))
    cat("Problematic rows:", head(problematic_rows, 10), "\n")
    print(all_results[head(problematic_rows, 5), ])
    stop("Data contains NA or Inf values")
  }
  
  cat("\nCreating plots...\n")
  
  # Calculate pointwise means
  pointwise_means_single <- all_results %>%
    group_by(n) %>%
    summarise(mean_difference_norm = mean(difference_norm), .groups = "drop")
  
  # Main plot: trajectories of difference norms
  # Create theoretical curves data for legend
  x_range <- range(all_results$n)
  theoretical_data <- data.frame(
    n = seq(x_range[1], x_range[2], length.out = 100)
  )
  theoretical_data$inv_sqrt_n <- 1/sqrt(theoretical_data$n)
  
  p1 <- ggplot(all_results, aes(x = n, y = difference_norm)) +
    # Individual trajectories
    geom_line(aes(group = trajectory, color = "Individual trajectories"), alpha = 0.6) +
    # Pointwise means
    geom_line(data = pointwise_means_single, aes(x = n, y = mean_difference_norm, color = "Mean"), 
              linetype = "dashed", linewidth = 0.5) +
    # Theoretical curves
    geom_line(data = theoretical_data, aes(x = n, y = inv_sqrt_n, color = "theoretical"), 
              linetype = "solid", linewidth = 0.5) +
    # Custom colors and legend
    scale_color_manual(
      values = c(
        "Individual trajectories" = "steelblue",
        "Mean" = "red",
        "theoretical" = "black"
      ),
      labels = c(
        "Individual trajectories" = expression("||LHS - RHS||"[2]),
        "Mean" = "Mean",
        "theoretical" = expression(frac(1, sqrt(n)))
      ),
      name = NULL  # Remove legend title
    ) +
    labs(
      x = expression("Sample size (" * italic(n) * ")"),
      y = NULL  # Remove y-axis label
    ) +
    # Set y-axis limit for better visualization
    ylim(0, 0.25) +
    theme_minimal() +
    theme(legend.position = c(0.85, 0.85),  # Move legend more inside
          legend.justification = c("right", "top"),
          legend.text = element_text(size = 19),  # Larger legend text (16+3)
          legend.key.size = unit(1.8, "lines"),  # Larger legend key size (1.5+0.3)
          axis.text.x = element_text(size = 19),  # Larger x-axis tick labels (14+5)
          axis.text.y = element_text(size = 19),  # Larger y-axis tick labels (14+5)
          axis.title.x = element_text(size = 19))  # Larger x-axis title (16+3)
  
  
  # Save plots
  plot_filename <- if (plot_suffix == "") {
    "output/bahadur_vmf/vmf_bahadur_trajectories.png"
  } else {
    paste0("output/bahadur_vmf/vmf_bahadur_mu", plot_suffix, ".png")
  }
  ggsave(plot_filename, p1, width = 12, height = 8, dpi = 600)
  
  cat("Plot saved to", plot_filename, "\n")
  
  # Summary statistics
  final_results <- all_results[all_results$n == max(sample_sizes), ]
  
  cat("\n=== Summary Statistics (n =", max(sample_sizes), ") ===\n")
  cat("Mean difference norm:", round(mean(final_results$difference_norm), 6), "\n")
  cat("Median difference norm:", round(median(final_results$difference_norm), 6), "\n")
  cat("SD difference norm:", round(sd(final_results$difference_norm), 6), "\n")
  cat("Max difference norm:", round(max(final_results$difference_norm), 6), "\n")
  
  return(all_results)
}

#' Run Bahadur analysis for multiple mu values and create comparison plot
#' @param mu_list List of mu vectors to compare
#' @param kappa_true True concentration parameter (default: 1)
#' @param n_trajectories Number of trajectories per mu (default: 20)
#' @param n_cores Number of CPU cores to use (default: 10)
run_bahadur_comparison <- function(mu_list = list(
                                    c(1, 0, 0),
                                    c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)),
                                    c(-1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
                                   ),
                                   kappa_true = 1, n_trajectories = 20, n_cores = 10) {
  
  cat("=== Bahadur Analysis Comparison for Multiple mu Values ===\n\n")
  
  # Define proper mathematical labels for each mu (hardcoded, no approximations)
  mu_labels <- c(
    expression("μ" * " = (1, 0, 0)"),
    expression("μ" * " = " * bgroup("(", list(frac(1, sqrt(3)), frac(1, sqrt(3)), frac(1, sqrt(3))), ")")),
    expression("μ" * " = " * bgroup("(", list(-frac(1, sqrt(3)), frac(1, sqrt(3)), frac(1, sqrt(3))), ")"))
  )
  
  # Run analysis for each mu
  all_comparison_results <- list()
  
  for (i in 1:length(mu_list)) {
    mu_current <- mu_list[[i]]
    cat("Running analysis for mu", i, ":", paste(round(mu_current, 4), collapse = ", "), "\n")
    
    # Use different seed for each mu to ensure different trajectories but reproducible results
    mu_seed <- 123 + i * 1000  # Seeds: 1123, 2123, 3123, etc.
    results <- run_bahadur_analysis(mu_current, kappa_true, n_trajectories, n_cores, 
                                   plot_suffix = i, seed = mu_seed)
    results$mu_index <- i
    results$mu_label <- mu_labels[i]  # Use proper mathematical expression
    
    all_comparison_results[[i]] <- results
    
    cat("Individual plot saved for mu", i, "\n")
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_comparison_results)
  
  # Since we already checked for NA/Inf values and would have stopped, no need to filter
  # Calculate pointwise means for comparison plot
  pointwise_means_comparison <- combined_results %>%
    group_by(n, mu_index) %>%
    summarise(mean_difference_norm = mean(difference_norm), .groups = "drop")
  
  # Create theoretical data for comparison plot
  x_range <- range(combined_results$n)
  theoretical_data <- data.frame(
    n = seq(x_range[1], x_range[2], length.out = 100)
  )
  theoretical_data$inv_sqrt_n <- 1/sqrt(theoretical_data$n)
  
  # Create comparison plot with facets
  # Add factor levels for proper ordering and mathematical expressions
  combined_results$mu_label_factor <- factor(
    combined_results$mu_index,
    levels = 1:3,
    labels = mu_labels
  )
  
  pointwise_means_comparison$mu_label_factor <- factor(
    pointwise_means_comparison$mu_index,
    levels = 1:3,
    labels = mu_labels
  )
  
  p_comparison <- ggplot(combined_results, aes(x = n, y = difference_norm)) +
    # Individual trajectories
    geom_line(aes(group = trajectory, color = "Individual trajectories"), alpha = 0.6) +
    # Pointwise means
    geom_line(data = pointwise_means_comparison, aes(x = n, y = mean_difference_norm, color = "Mean"), 
              linetype = "dashed", linewidth = 0.5) +
    # Theoretical curves
    geom_line(data = theoretical_data, aes(x = n, y = inv_sqrt_n, color = "theoretical"), 
              linetype = "solid", linewidth = 0.5) +
    # Custom colors and legend
    scale_color_manual(
      values = c(
        "Individual trajectories" = "steelblue",
        "Mean" = "red",
        "theoretical" = "black"
      ),
      labels = c(
        "Individual trajectories" = expression("||LHS - RHS||"[2]),
        "Mean" = "Mean",
        "theoretical" = expression(frac(1, sqrt(n)))
      ),
      name = NULL  # Remove legend title
    ) +
    labs(
      x = expression("Sample size (" * italic(n) * ")"),
      y = NULL  # Remove y-axis label
    ) +
    # Set y-axis limit for better visualization
    ylim(0, 0.25) +
    # Create facets for each mu with mathematical expressions
    facet_wrap(~ mu_label_factor, nrow = 1, labeller = label_parsed) +
    theme_minimal() +
    theme(
      # Position legend in the first panel (upper right)
      legend.position = c(0.16, 0.85),  # Positioned in first panel
      legend.justification = c("right", "top"),
      legend.text = element_text(size = 19),  # Larger legend text (16+3)
      legend.key.size = unit(1.8, "lines"),  # Larger legend key size (1.5+0.3)
      strip.text = element_text(size = 17),  # Larger facet labels (14+3)
      axis.text.x = element_text(size = 19),  # Larger x-axis tick labels (14+5)
      axis.text.y = element_text(size = 19),  # Larger y-axis tick labels (14+5)
      axis.title.x = element_text(size = 19)  # Larger x-axis title (16+3)
    )
  
  # Save comparison plot
  ggsave("output/bahadur_vmf/vmf_bahadur_comparison.png", p_comparison, width = 18, height = 6, dpi = 600)
  
  cat("\nComparison plot saved to output/bahadur_vmf/vmf_bahadur_comparison.png\n")
  cat("Individual plots saved to output/bahadur_vmf/vmf_bahadur_mu1.png, vmf_bahadur_mu2.png, vmf_bahadur_mu3.png\n")
  
  return(combined_results)
}

# Example usage with default parameters for S^2 (3D sphere)
# mu_example <- c(1, 0, 0)  # Point on S^2

# Run example analysis (commented out by default)
# results <- run_bahadur_analysis(mu_example, kappa_true = 1)

# Run comparison analysis with default 20 trajectories
comparison_results <- run_bahadur_comparison()
