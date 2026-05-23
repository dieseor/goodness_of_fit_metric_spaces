# Distance Profile Plotting and Analysis Functions
# Main functions for generating comparison plots and analysis

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sphunif)
  library(MASS)
  library(mvtnorm)
  library(viridis)
  library(gridExtra)
})

#' Compute distances from omega to all points in a dataset
#' @param omega Reference point (vector)
#' @param data Matrix with observations in rows
#' @return Vector of distances
compute_distances <- function(omega, data, distance) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  apply(data, 1, function(x) distance(omega, x))
}

#' Compute empirical distance profile
#' @param distances Vector of distances
#' @param t_values Vector of t values to evaluate the profile at
#' @return Vector of empirical probabilities F_hat(t)
empirical_distance_profile <- function(distances, t_values) {
  sapply(t_values, function(t) mean(distances <= t))
}

#' Generate confidence bands for empirical distance profile
#' @param n Sample size
#' @param p_values Theoretical probabilities at t_values
#' @param alpha Significance level (default 0.05 for 95% confidence)
#' @return List with lower and upper bounds
confidence_bands <- function(n, p_values, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  se <- sqrt(p_values * (1 - p_values) / n)
  
  lower <- pmax(0, p_values - z_alpha * se)
  upper <- pmin(1, p_values + z_alpha * se)
  
  list(lower = lower, upper = upper)
}

#' Analyze distance profiles for a single omega and sample size
#' @param omega Reference point
#' @param n_samples Sample size
#' @param n_simulations Number of simulation runs (default 10)
#' @param data_generator Function to generate samples
#' @param theoretical_profile Function to compute theoretical profile
#' @param t_max Maximum t value for evaluation
#' @param n_points Number of points to evaluate profile at
#' @return List with empirical profiles, theoretical profile, and confidence bands
analyze_distance_profile <- function(omega, n_samples, n_simulations = 10, distance,
                                    data_generator, theoretical_profile, 
                                    t_max = 5, n_points = 100) {
  
  # Create t values for evaluation
  t_values <- seq(0, t_max, length.out = n_points)
  
  # Generate theoretical profile
  theoretical <- theoretical_profile(omega, t_values)
  
  # Generate empirical profiles from multiple simulations
  empirical_profiles <- matrix(0, nrow = n_simulations, ncol = n_points)
  
  for (i in 1:n_simulations) {
    # Generate sample
    sample_data <- data_generator(n_samples)
    
    # Compute distances
    distances <- compute_distances(omega, sample_data, distance)
    
    # Compute empirical profile
    empirical_profiles[i, ] <- empirical_distance_profile(distances, t_values)
  }
  
  # Compute confidence bands based on theoretical profile
  conf_bands <- confidence_bands(n_samples, theoretical)
  
  list(
    t_values = t_values,
    theoretical = theoretical,
    empirical_profiles = empirical_profiles,
    conf_bands = conf_bands,
    omega = omega,
    n_samples = n_samples
  )
}

#' Create a single distance profile plot
#' @param analysis_result Result from analyze_distance_profile()
#' @param title Plot title (unused; validation plots are generated without titles)
#' @param legend_position Position for legend: "bottom_right" or "top_left"
#' @return ggplot object
plot_distance_profile <- function(analysis_result, title = NULL, legend_position = "bottom_right") {
  
  t_values <- analysis_result$t_values
  theoretical <- analysis_result$theoretical
  empirical_profiles <- analysis_result$empirical_profiles
  conf_bands <- analysis_result$conf_bands
  omega <- analysis_result$omega
  n_samples <- analysis_result$n_samples
  
  # Create data frame for empirical profiles
  empirical_df <- data.frame()
  for (i in seq_len(nrow(empirical_profiles))) {
    temp_df <- data.frame(
      t = t_values,
      F_hat = empirical_profiles[i, ],
      simulation = factor(i)
    )
    empirical_df <- rbind(empirical_df, temp_df)
  }
  
  # Create data frame for theoretical profile and confidence bands
  theoretical_df <- data.frame(
    t = t_values,
    F_theoretical = theoretical,
    lower = conf_bands$lower,
    upper = conf_bands$upper
  )
  
  # Create y-axis label with mathematical notation
  y_label <- bquote(F[omega]^mu * (t))
  
  # Calculate legend positioning
  if (legend_position == "top_left") {
    # Top left positioning
    legend_x_base <- max(t_values) * 0.04
    legend_y_top <- 0.97
    legend_y_mid <- 0.92
    legend_y_bot <- 0.87
    hjust_val <- 0
  } else {
    # Bottom right positioning (default)
    legend_x_base <- max(t_values) * 0.50
    legend_y_top <- 0.17
    legend_y_mid <- 0.12
    legend_y_bot <- 0.07
    hjust_val <- 0
  }

  # Create the plot without title and with proper legend
  p <- ggplot() +
    # Confidence bands
    geom_ribbon(data = theoretical_df, 
                aes(x = t, ymin = lower, ymax = upper),
                fill = "#6b6b9e", alpha = 0.3) +
    
    # Empirical profiles (using step function for ECDF)
    geom_step(data = empirical_df, 
              aes(x = t, y = F_hat, group = simulation),
              color = "#252323", alpha = 0.6, linewidth = 0.5, direction = "hv") +
    
    # Theoretical profile
    geom_line(data = theoretical_df, 
              aes(x = t, y = F_theoretical),
              color = "red", linewidth = 1) +
    
    
    # Create manual legend using annotate (no gaps)
    labs(
      x = "t",
      y = y_label
    ) +
    
    # Add legend annotations
    annotate("text", x = legend_x_base + max(t_values) * 0.04, y = legend_y_top, 
             label = "Distance profile", hjust = hjust_val, size = 5.0, color = "black") +
    annotate("segment", x = legend_x_base, xend = legend_x_base + max(t_values) * 0.03, 
             y = legend_y_top, yend = legend_y_top, color = "red", linewidth = 0.8) +
    
    annotate("text", x = legend_x_base + max(t_values) * 0.04, y = legend_y_mid, 
             label = "Empirical distance profile", hjust = hjust_val, size = 5.0, color = "black") +
    annotate("segment", x = legend_x_base, xend = legend_x_base + max(t_values) * 0.03, 
             y = legend_y_mid, yend = legend_y_mid, color = "#252323", linewidth = 0.8) +
    
    annotate("text", x = legend_x_base + max(t_values) * 0.04, y = legend_y_bot, 
             label = "Confidence band", hjust = hjust_val, size = 5.0, color = "black") +
    annotate("rect", xmin = legend_x_base, xmax = legend_x_base + max(t_values) * 0.03, 
             ymin = legend_y_bot - 0.005, ymax = legend_y_bot + 0.005, fill = "#6b6b9e", alpha = 0.3) +
    
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14, hjust = 0.5),
      axis.text = element_text(size = 13),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    
    coord_cartesian(xlim = c(0, max(t_values)), ylim = c(0, 1))
  
  return(p)
}

#' Main function to create comprehensive distance profile analysis
#' @param omega_values List of omega vectors (3 values)
#' @param data_generator Function to generate samples: function(n) -> matrix
#' @param theoretical_profile Function: function(omega, t_values) -> probabilities
#' @param sample_sizes Vector of sample sizes (default c(50, 200))
#' @param n_simulations Number of simulations per case (default 10)
#' @param t_max Maximum t value for plots (default 5)
#' @param save_plots Whether to save plots to files (default TRUE)
#' @param output_dir Output directory for plots (default "output")
#' @param file_prefix Prefix for saved files (default "dp")
#' @param legend_positions Vector of legend positions for each omega (default all "bottom_right")
#' @return List of all plots
create_distance_profile_analysis <- function(omega_values, 
                                           data_generator,
                                           theoretical_profile,
                                           distance,
                                           sample_sizes = c(50, 200),
                                           n_simulations = 10,
                                           t_max = 5,
                                           save_plots = TRUE,
                                           output_dir = "output",
                                           file_prefix = "dp",
                                           legend_positions = NULL) {
  
  # Allow flexible number of omega values (not just 3)
  n_omegas <- length(omega_values)
  
  # Set default legend positions if not provided
  if (is.null(legend_positions)) {
    legend_positions <- rep("bottom_right", n_omegas)
  }
  
  # Ensure output directory exists
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots <- list()
  plot_counter <- 1
  
  cat("Creating distance profile analysis...\n")
  cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
  cat("Number of simulations per case:", n_simulations, "\n")
  cat("Number of omega values:", length(omega_values), "\n\n")
  
  for (i in 1:length(omega_values)) {
    omega <- omega_values[[i]]
    omega_name <- paste0("omega_", i)
    
    cat("Processing omega", i, ":", paste(round(omega, 3), collapse = ", "), "\n")
    
    for (j in 1:length(sample_sizes)) {
      n <- sample_sizes[j]
      
      cat("  Sample size:", n, "... ")
      
      # Analyze distance profile
      analysis <- analyze_distance_profile(
        omega = omega,
        n_samples = n,
        n_simulations = n_simulations,
        data_generator = data_generator,
        theoretical_profile = theoretical_profile,
        distance = distance,
        t_max = t_max
      )
      
      # Create plot
      legend_pos <- if (i <= length(legend_positions)) legend_positions[i] else "bottom_right"
      plot <- plot_distance_profile(analysis, legend_position = legend_pos)
      
      # Store plot
      plot_name <- paste0(omega_name, "_n", n)
      plots[[plot_name]] <- plot
      
      # Save individual plot
      if (save_plots) {
        filename <- file.path(output_dir, paste0(file_prefix, "_", plot_name, ".png"))
        ggsave(filename, plot, width = 8, height = 6, dpi = 300)
      }
      
      cat("Done\n")
      plot_counter <- plot_counter + 1
    }
  }
  
  # Create combined plots with flexible layout
  if (save_plots) {
    cat("\nCreating combined plots...\n")
    
    # Determine optimal grid layout
    n_cols <- min(5, n_omegas)  # Max 5 columns
    n_rows <- ceiling(n_omegas / n_cols)
    
    # Combined plot for n=50
    plots_n50 <- plots[grepl("_n50", names(plots))]
    if (length(plots_n50) > 0) {
      combined_n50 <- grid.arrange(grobs = plots_n50, ncol = n_cols, nrow = n_rows)
      ggsave(file.path(output_dir, paste0(file_prefix, "_combined_n50.png")), 
             combined_n50, width = n_cols * 6, height = n_rows * 4, dpi = 300)
    }
    
    # Combined plot for n=200
    plots_n200 <- plots[grepl("_n200", names(plots))]
    if (length(plots_n200) > 0) {
      combined_n200 <- grid.arrange(grobs = plots_n200, ncol = n_cols, nrow = n_rows)
      ggsave(file.path(output_dir, paste0(file_prefix, "_combined_n200.png")), 
             combined_n200, width = n_cols * 6, height = n_rows * 4, dpi = 300)
    }
    
    cat("All plots saved to:", output_dir, "\n")
  }
  
  cat("\nAnalysis complete!\n")
  cat("Generated", length(plots), "individual plots\n")
  if (save_plots) cat("Saved", length(plots) + 2, "files total\n")
  
  return(plots)
}

cat("Distance profile analysis functions loaded successfully!\n")
cat("Main function: create_distance_profile_analysis()\n")
