# ============================================================================
# GAUSSIAN PROCESS FUNCTIONS - NORMAL DISTRIBUTION
# ============================================================================
# Functions specific to the univariate Normal distribution for Gaussian
# process limit theory in goodness-of-fit testing.
# ============================================================================

# Load generic utilities
utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else if (file.exists(file.path("..", "utils.R"))) {
  file.path("..", "utils.R")
} else {
  stop("Could not find utils.R in current directory or parent directory.")
}
source(utils_path)

# ============================================================================
# COVARIANCE MATRIX COMPUTATION FOR NORMAL MODEL
# Supports:
#   - simple null
#   - composite null with unknown_param in {"mu", "sigma", "both"}


# ============================================================================
# Calculate covariance row
# ============================================================================

row_cov_normal <- function(
    h0,
    unknown_param = NULL,  # "mu", "sigma", "both"
    omega1, t1,
    omega2_vec, t2_vec,
    mu, sigma
) {
  
  n <- length(omega2_vec)
  
  # ------------------------------
  # Base C matrix: C_{12}
  # ------------------------------
  lower1 <- omega1 - t1
  upper1 <- omega1 + t1
  
  lower2_vec <- omega2_vec - t2_vec
  upper2_vec <- omega2_vec + t2_vec
  
  intersect_low <- pmax(lower1, lower2_vec)
  intersect_up  <- pmin(upper1, upper2_vec)
  
  valid <- (intersect_low < intersect_up)
  joint <- numeric(n)
  
  if (any(valid)) {
    joint[valid] <- pnorm((intersect_up[valid] - mu)/sigma) -
                    pnorm((intersect_low[valid] - mu)/sigma)
  }
  
  f1 <- theoretical_distance_profile_normal(omega1, mu, sigma, t1)
  f2_vec <- theoretical_distance_profile_normal(omega2_vec, mu, sigma, t2_vec)
  
  C_base <- joint - f1 * f2_vec
  
  # If simple null, you're done
  if (h0 == "simple") {
    return(C_base)
  }
  
  # =====================================================================
  # Composite null
  # =====================================================================
  
  # Initialize with base covariance
  cov_vec <- C_base
  
  # ---------------------------------------------------------------------
  # Case 1: mu unknown
  # ---------------------------------------------------------------------
if (unknown_param == "mu") {
  dF1_mu <- dotF_mu_normal(omega1, mu, sigma, t1)
  dF2_mu <- dotF_mu_normal(omega2_vec, mu, sigma, t2_vec)
  cov_vec <- cov_vec -
    (sigma^2) * (dF1_mu * dF2_mu)
}

if (unknown_param == "sigma") {
  dF1_sigma <- dotF_sigma_normal(omega1, mu, sigma, t1)
  dF2_sigma <- dotF_sigma_normal(omega2_vec, mu, sigma, t2_vec)
  cov_vec <- cov_vec -
    (sigma^2 / 2) * (dF1_sigma * dF2_sigma)
}

if (unknown_param == "both") {
  dF1_mu <- dotF_mu_normal(omega1, mu, sigma, t1)
  dF2_mu <- dotF_mu_normal(omega2_vec, mu, sigma, t2_vec)
  dF1_sigma <- dotF_sigma_normal(omega1, mu, sigma, t1)
  dF2_sigma <- dotF_sigma_normal(omega2_vec, mu, sigma, t2_vec)
  
  cov_vec <- cov_vec -
    (sigma^2) * (dF1_mu * dF2_mu) -
    (sigma^2 / 2) * (dF1_sigma * dF2_sigma)
}
  return(cov_vec)
}

#' Create covariance matrix for Gaussian process (vectorized + parallelized)
#' @param omega_grid Vector of omega values
#' @param t_grid Vector of t values
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @param n_cores Number of cores for parallel computation (default: 7)
#' @return Covariance matrix of size (n_omega * n_t) x (n_omega * n_t)
cov_normal <- function(omega_grid, t_grid, mu, sigma, n_cores = 7, h0 = c("simple","composite"), unknown_param = NULL) {
  h0 <- match.arg(h0)
  n_omega <- length(omega_grid)
  n_t <- length(t_grid)
  n_total <- n_omega * n_t
  
  # Create all combinations of (omega, t) pairs
  grid_combinations <- expand.grid(omega = omega_grid, t = t_grid)
  
  cat("Creating covariance matrix of size", n_total, "x", n_total, "\n")
  cat("Using", n_cores, "cores\n")

  # Extract vectors for vectorized operations
  omega_vec <- grid_combinations$omega
  t_vec <- grid_combinations$t
  
  # Setup parallel cluster
  library(parallel)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)  # Ensure cluster is stopped even if interrupted
  
  # Set seed for reproducibility in parallel workers
  clusterSetRNGStream(cl, iseed = 42)
  
  # Export helper function to workers
  clusterEvalQ(cl, {
    # Load base namespace if needed; the workers will get functions via clusterExport
    NULL
  })
  
  clusterExport(cl, c("omega_vec", "t_vec", "mu", "sigma", "row_cov_normal", "theoretical_distance_profile_normal", "dotF_mu_normal", "dotF_sigma_normal", "h0", "unknown_param"), envir = environment())
  
  start_time <- Sys.time()
  
  # Parallel computation of covariance matrix
  # Distribute rows across cores (round-robin for load balancing)
  row_chunks <- vector("list", n_cores)
  for (i in 1:n_cores) {
    row_chunks[[i]] <- seq(from = i, to = n_total, by = n_cores)
  }
  
  cat("Computing", sum(sapply(row_chunks, length)), "rows in parallel\n")
  
  # Compute rows in parallel
  # pass h0 and unknown_param into row function so we can compute composite null
  results <- parLapply(cl, row_chunks, function(row_indices) {
    local_rows <- list()
    for (i in row_indices) {
      row_vec <- row_cov_normal(
        h0,
        unknown_param,
        omega_vec[i], t_vec[i], omega_vec, t_vec, mu, sigma
      )
      local_rows[[length(local_rows) + 1]] <- list(i = i, row = row_vec)
    }
    return(local_rows)
  })
  # Assemble matrix
  cat("Assembling matrix from parallel results...\n")
  cov_matrix <- matrix(0, n_total, n_total)
  for (worker_results in results) {
    for (result in worker_results) {
      cov_matrix[result$i, ] <- result$row
    }
  }
  
  end_time <- Sys.time()
  time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("Covariance matrix created in", round(time_elapsed, 1), "seconds!\n")
  
  return(cov_matrix)
}

# ============================================================================
# PROCESS SIMULATION
# ============================================================================

#' Simulate the limiting Gaussian process for Normal
#' @param omega_grid Vector of omega values
#' @param t_grid Vector of t values
#' @param mu Mean of X
#' @param sigma Standard deviation of X
#' @param M Number of Monte Carlo simulations
#' @param n_cores Number of cores for covariance matrix computation
#' @return Vector of M supremum values from the Gaussian process
simulate_limit_gaussian_normal <- function(omega_grid, t_grid, mu, sigma, M = 10000, n_cores = 7, h0 = c("simple","composite"), unknown_param = NULL) {
  
  cat("=== Simulating Gaussian Process Limit for Normal Distribution ===\n")
  
  # Create covariance matrix
  cov_matrix <- cov_normal(omega_grid, t_grid, mu, sigma, n_cores, h0 = h0, unknown_param = unknown_param)
  
  # Use generic function to simulate
  supremum_values <- simulate_limit_gaussian(cov_matrix, M)
  
  return(supremum_values)
}

compute_theoretical_sample_profile_normal <- function(center, radii, mu, sigma) {
  theoretical_distance_profile_normal(center, mu, sigma, radii)
}

compute_sample_ks_sup_normal <- function(sample_x,
                                         mu,
                                         sigma,
                                         h0 = c("simple", "composite"),
                                         unknown_param = c(NULL, "mu", "sigma", "both")) {
  h0 <- match.arg(h0)
  if (h0 == "composite") {
    unknown_param <- match.arg(unknown_param, choices = c("mu", "sigma", "both"))
    if (unknown_param == "mu") {
      mu_hat <- mean(sample_x)
      sigma_hat <- sigma
    } else if (unknown_param == "sigma") {
      mu_hat <- mu
      sigma_hat <- sqrt(mean((sample_x - mu_hat)^2))
    } else {
      mu_hat <- mean(sample_x)
      sigma_hat <- sqrt(mean((sample_x - mu_hat)^2))
    }
  } else {
    mu_hat <- mu
    sigma_hat <- sigma
  }

  n <- length(sample_x)
  max_diff <- 0
  for (i in seq_len(n)) {
    distances_i <- abs(sample_x - sample_x[[i]])
    order_i <- order(distances_i)
    sorted_distances_i <- distances_i[order_i]
    empirical_i <- seq_len(n) / n
    theoretical_i <- compute_theoretical_sample_profile_normal(
      center = sample_x[[i]],
      radii = sorted_distances_i,
      mu = mu_hat,
      sigma = sigma_hat
    )
    max_diff <- max(max_diff, max(abs(empirical_i - theoretical_i)))
  }

  sqrt(n) * max_diff
}

# ============================================================================
# SIMULATE EMPIRICAL PROCESS (VECTORIZED + PARALLELIZED)
# Returns a numeric vector of length M with suprema sqrt(n) * sup_{(omega,t)} |F_hat - F_theta|
# - Supports simple null (fixed mu,sigma) and composite null (unknown_param in "mu","sigma","both")
# - Requires theoretical_distance_profile_normal(omega_vector, mu, sigma, t_vector) to be defined
# ============================================================================

simulate_empirical_process_normal <- function(omega_grid,
                                              t_grid,
                                              n,
                                              mu,
                                              sigma,
                                              M = 10000,
                                              n_cores = 7,
                                              h0 = c("simple", "composite"),
                                              unknown_param = c(NULL, "mu", "sigma", "both"),
                                              empirical_ks_mode = c("sample", "grid")
) {
  # --- args and sanity checks
  h0 <- match.arg(h0)
  empirical_ks_mode <- match.arg(empirical_ks_mode)
  if (h0 == "composite") {
    unknown_param <- match.arg(unknown_param, choices = c("mu","sigma","both"))
  } else {
    unknown_param <- NULL
  }

  cat("=== Simulating Empirical Process for Normal Distribution ===\n")
  cat("Empirical KS mode:", empirical_ks_mode, "\n")

  n_omega <- length(omega_grid)
  n_t <- length(t_grid)
  
  # grid for later convenience (row-major: omega x t)
  grid_df <- expand.grid(omega = omega_grid, t = t_grid)
  
  # Precompute theoretical F matrix only for simple null
  if (h0 == "simple") {
    F_theoretical <- theoretical_distance_profile_normal(grid_df$omega, mu, sigma, grid_df$t)
    # arrange into n_omega x n_t matrix (rows = omega, cols = t)
    F_theoretical_matrix_global <- matrix(F_theoretical, nrow = n_omega, ncol = n_t)
  } else {
    # placeholder (not used on master) - workers will compute per-sim MLE
    F_theoretical_matrix_global <- NULL
  }
  
  # ---- setup parallel cluster
  library(parallel)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  
  # reproducible RNG streams
  clusterSetRNGStream(cl, iseed = 12345)
  
  # export required objects and functions to workers
  clusterExport(cl, c("omega_grid", "t_grid", "n", "mu", "sigma",
                      "grid_df", "n_omega", "n_t",
                      "F_theoretical_matrix_global",
                      "theoretical_distance_profile_normal",
                      "compute_theoretical_sample_profile_normal",
                      "compute_sample_ks_sup_normal",
                      "h0", "unknown_param", "empirical_ks_mode"),
                envir = environment())
  
  # worker-side load of packages if needed (none required here, but keep safe)
  clusterEvalQ(cl, {
    NULL
  })
  
  start_time <- Sys.time()
  cat("Computing", M, "simulations in parallel across", n_cores, "cores...\n")
  
  # split simulation indices across cores
  sim_chunks <- vector("list", n_cores)
  for (i in seq_len(n_cores)) sim_chunks[[i]] <- seq(from = i, to = M, by = n_cores)
  
  # parallel simulation: each worker gets a vector 'sim_indices' and returns a list with indices and values
  results_list <- parLapply(cl, sim_chunks, function(sim_indices) {
    # local copies of exported objects are available
    local_supremums <- numeric(length(sim_indices))
    
    # iterate simulations assigned to this worker
    for (j in seq_along(sim_indices)) {
      # generate sample
      sample_x <- rnorm(n, mean = mu, sd = sigma)

      if (identical(empirical_ks_mode, "sample")) {
        local_supremums[j] <- compute_sample_ks_sup_normal(
          sample_x = sample_x,
          mu = mu,
          sigma = sigma,
          h0 = h0,
          unknown_param = unknown_param
        )
      } else {
        # compute distance matrix: n x n_omega (rows: sample, cols: omega)
        distance_matrix <- outer(sample_x, omega_grid, FUN = function(x,w) abs(x - w))

        # compute F_hat for all (omega,t) pairs: result n_omega x n_t
        # For each t_val, proportion of sample points with distance <= t_val (col-wise)
        F_hat_matrix <- sapply(t_grid, function(t_val) {
          colMeans(distance_matrix <= t_val)
        })
        # ensure dimensions: n_omega x n_t
        if (!is.matrix(F_hat_matrix)) F_hat_matrix <- matrix(F_hat_matrix, nrow = n_omega, ncol = n_t)

        # compute theoretical F matrix for this simulation
        if (h0 == "simple") {
          F_theoretical_matrix <- F_theoretical_matrix_global
        } else {
          # estimate MLE from this sample according to unknown_param
          if (unknown_param == "mu") {
            mu_hat <- mean(sample_x)
            sigma_hat <- sigma
          } else if (unknown_param == "sigma") {
            mu_hat <- mu
            sigma_hat <- sqrt(mean((sample_x - mu_hat)^2))  # population-style estimator (consistent)
          } else if (unknown_param == "both") {
            mu_hat <- mean(sample_x)
            sigma_hat <- sqrt(mean((sample_x - mu_hat)^2))
          } else {
            stop("unknown_param must be one of 'mu','sigma','both' when h0 == 'composite'")
          }
          # compute theoretical distance profile at estimated params
          F_theoretical_vec <- theoretical_distance_profile_normal(grid_df$omega, mu_hat, sigma_hat, grid_df$t)
          F_theoretical_matrix <- matrix(F_theoretical_vec, nrow = n_omega, ncol = n_t)
        }

        # scaled absolute differences
        scaled_diff_matrix <- sqrt(n) * abs(F_hat_matrix - F_theoretical_matrix)

        # supremum over grid
        local_supremums[j] <- max(scaled_diff_matrix)
      }
    }
    
    # return indices + values so master can place them in correct order
    list(indices = sim_indices, values = local_supremums)
  })
  
  # combine results into a vector of length M in the correct order
  sup_vec <- numeric(M)
  for (res in results_list) {
    sup_vec[res$indices] <- res$values
  }
  
  end_time <- Sys.time()
  dur <- difftime(end_time, start_time, units = "secs")
  cat("Done. Elapsed (secs):", round(as.numeric(dur),2), "\n")
  
  return(sup_vec)
}

# ============================================================================
# VISUALIZATION AND ANALYSIS
# ============================================================================

#' Visualize convergence of empirical process to limit process as n grows
#' @param n_values Vector of sample sizes to compare
#' @param mu Mean of the normal distribution
#' @param sigma Standard deviation of the normal distribution
#' @param omega_points Number of omega grid points
#' @param t_points Number of t grid points
#' @param M Number of Monte Carlo simulations
#' @param n_cores Number of cores for parallelization
#' @return List with all simulated values and convergence plot
visualize_convergence_to_limit_normal <- function(n_values = c(50, 100, 500), 
                                                  mu = 0, 
                                                  sigma = 1,
                                                  omega_points = 50,
                                                  t_points = 50,
                                                  omega_grid = NULL,
                                                  t_grid = NULL,
                                                  M = 10000,
                                                  n_cores = 7,
                                                  h0 = c("simple","composite"),
                                                  unknown_param = NULL,
                                                  empirical_ks_mode = c("sample", "grid"),
                                                  n50_adjust_multiplier = 2,
                                                  xlim = NULL,
                                                  qqplot = FALSE,
                                                  qqplot_save = NULL) {
  h0 <- match.arg(h0)
  empirical_ks_mode <- match.arg(empirical_ks_mode)
  
  cat("=== Visualizing Convergence to Gaussian Limit (Normal Distribution) ===\n")
  cat("Sample sizes to compare:", paste(n_values, collapse = ", "), "\n")
  cat("Distribution: N(μ =", mu, ", σ² =", sigma^2, ")\n")
  cat("Empirical KS mode:", empirical_ks_mode, "\n")
  
  # Use provided grid if present, otherwise automatic grid computation based on mu and sigma
  if (is.null(omega_grid) || is.null(t_grid)) {
    t_max_auto <- sigma * qnorm(0.995)
    omega_max_auto <- mu + t_max_auto
    omega_min_auto <- mu - t_max_auto
    t_min_auto <- 0
    omega_grid <- seq(omega_min_auto, omega_max_auto, length.out = omega_points)
    t_grid <- seq(t_min_auto, t_max_auto, length.out = t_points)
  } else {
    # If provided, override points length
    omega_points <- length(omega_grid)
    t_points <- length(t_grid)
  }
  
    cat("Omega grid:", omega_points, "points from", round(min(omega_grid), 2), 
      "to", round(max(omega_grid), 2), "\n")
    cat("t grid:", t_points, "points from", round(min(t_grid), 2), 
      "to", round(max(t_grid), 2), "\n")
  cat("Monte Carlo simulations:", M, "each process\n")
  cat("Using", n_cores, "cores for parallelization\n\n")
  
  # Simulate limiting Gaussian process (once)
  cat("Simulating limit Gaussian process...\n")
  limit_values <- simulate_limit_gaussian_normal(
    omega_grid, t_grid, mu, sigma, M, n_cores, h0 = h0, unknown_param = unknown_param
  )

  # Simulate empirical process for each n
  empirical_data <- list()
  for (n in n_values) {
    cat("Simulating empirical process with n =", n, "...\n")
    empirical_values <- simulate_empirical_process_normal(
      omega_grid, t_grid, n, mu, sigma, M, n_cores,
      h0 = h0, unknown_param = unknown_param, empirical_ks_mode = empirical_ks_mode
    )
    empirical_data[[as.character(n)]] <- empirical_values
  }
  
  # Prepare data for plotting
  n_max <- max(n_values)
  
  # Density curves: All n values
  density_data <- data.frame()
  for (n in n_values) {
    n_char <- as.character(n)
    density_data <- rbind(
      density_data,
      data.frame(
        values = empirical_data[[n_char]],
        n = n,
        label = paste0("n=", n)
      )
    )
  }
  
  # Add limit process to density data
  limit_label <- "G"
  density_data_full <- rbind(
    density_data,
    data.frame(
      values = limit_values,
      n = Inf,
      label = limit_label
    )
  )
  
  # Create color palette (limit is black now)
  limit_color <- "#000000"
  n_colors <- length(n_values)
  empirical_colors <- scales::hue_pal()(n_colors)
  all_labels <- c(paste0("n=", n_values), limit_label)
  all_colors <- setNames(c(empirical_colors, limit_color), all_labels)
  linetype_vals <- setNames(c(rep("solid", n_colors), "dashed"), all_labels)
  density_data_full$label <- factor(density_data_full$label, levels = all_labels)
  
  # Create plot
  library(ggplot2)
  has_n50 <- any(n_values == 50)
  base_adjust <- 1.0
  adjust_n50 <- n50_adjust_multiplier * base_adjust
  if (has_n50) {
    cat(sprintf("[density] n=50 detected -> using adjust=%.1f for n=50 only (base adjust=%.1f, multiplier=%.1f)\n", adjust_n50, base_adjust, n50_adjust_multiplier))
  }
  p_convergence <- ggplot()
  if (has_n50) {
    p_convergence <- p_convergence +
      geom_density(
        data = subset(density_data_full, label != "n=50"),
        aes(x = values, color = label, linetype = label),
        bw = "bcv",
        fill = NA,
        linewidth = 1.0,
        adjust = base_adjust,
        show.legend = TRUE,
        key_glyph = draw_key_path
      ) +
      geom_density(
        data = subset(density_data_full, label == "n=50"),
        aes(x = values, color = label, linetype = label),
        bw = "bcv",
        fill = NA,
        linewidth = 1.0,
        adjust = adjust_n50,
        show.legend = TRUE,
        key_glyph = draw_key_path
      )
  } else {
    p_convergence <- p_convergence +
      geom_density(
        data = density_data_full,
        aes(x = values, color = label, linetype = label),
        bw = "bcv",
        fill = NA,
        linewidth = 1.0,
        adjust = base_adjust,
        show.legend = TRUE,
        key_glyph = draw_key_path
      )
  }
  p_convergence <- p_convergence +
        scale_color_manual(values = all_colors, breaks = all_labels, drop = FALSE) +
        scale_linetype_manual(values = linetype_vals, breaks = all_labels, drop = FALSE) +
    labs(
      x = "Supremum of the process",
      y = "Density",
      color = "Process",
      linetype = "Process"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.text = element_text(size = 19),
      legend.title = element_text(size = 19),
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19),
      axis.title.y = element_text(size = 19)
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = "solid", fill = NA, alpha = 1, linewidth = 1.5)),
      linetype = "none"
    )
  if (!is.null(xlim) && length(xlim) == 2 && !any(is.na(xlim))) {
    p_convergence <- p_convergence + coord_cartesian(xlim = xlim)
  }
  
  # No manual legend: using ggplot2 automatic legend
  
  # Statistical comparison
  cat("=== Kolmogorov-Smirnov Tests vs Limit ===\n")
  ks_results <- list()
  for (n in n_values) {
    n_char <- as.character(n)
    ks_test <- safe_ks_test(limit_values, empirical_data[[n_char]])
    ks_results[[n_char]] <- ks_test
    pval <- if (!is.null(ks_test$p.value)) ks_test$p.value else NA_real_
    if (is.na(pval)) {
      cat(sprintf("n=%5d: p-value = NA (test not run / invalid input)\n", n))
    } else {
      cat(sprintf("n=%5d: p-value = %.4f", n, pval))
      if (pval > 0.05) {
        cat(" -> Not significantly different ✓\n")
      } else {
        cat(" -> Significantly different\n")
      }
    }
  }
  cat("\n")
  # Optional QQ plot: use same simulated data
  qq_plot <- NULL
  if (isTRUE(qqplot)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) warning("ggplot2 required for QQ plot. Skipping QQ plot.")
    else {
      probs <- ppoints(M)
      limit_qs <- as.numeric(quantile(limit_values, probs = probs, type = 8, na.rm = TRUE))
      df_all <- data.frame()
      for (n in n_values) {
        empirical_qs <- as.numeric(quantile(empirical_data[[as.character(n)]], probs = probs, type = 8, na.rm = TRUE))
        df <- data.frame(sample_size = as.factor(n), p = probs, theoretical = limit_qs, empirical = empirical_qs)
        df_all <- rbind(df_all, df)
      }
      qq_plot <- ggplot2::ggplot(df_all, ggplot2::aes(x = theoretical, y = empirical, color = sample_size)) +
        ggplot2::geom_point(alpha = 0.7, size = 1.5) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        ggplot2::scale_color_manual(values = setNames(empirical_colors, as.character(n_values))) +
        ggplot2::labs(x = "Limit quantiles", y = "Empirical quantiles", color = "n") +
        ggplot2::theme_minimal()
      if (!is.null(qqplot_save)) ggplot2::ggsave(qqplot_save, plot = qq_plot)
    }
  }
  
  return(list(
    limit_values = limit_values,
    empirical_data = empirical_data,
    n_values = n_values,
    plot = p_convergence,
    qq_plot = qq_plot,
    ks_results = ks_results
  ))
}

cat("Normal distribution Gaussian process functions loaded successfully!\n")
