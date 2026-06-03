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

#' Build default sample size grid for Bahadur analysis
build_bahadur_sample_sizes <- function(min_n = 50, max_n = 100000) {
  full_grid <- c(
    seq(50, 1000, by = 100),
    seq(1250, 5000, by = 250),
    seq(5500, 10000, by = 500),
    seq(11000, 50000, by = 1000),
    seq(55000, 100000, by = 2500)
  )
  full_grid[full_grid >= min_n & full_grid <= max_n]
}

#' Default mu labels used in the AoS comparison figure
build_default_bahadur_mu_labels <- function() {
  c(
    expression("μ" * " = (1, 0, 0)"),
    expression("μ" * " = " * bgroup("(", list(frac(1, sqrt(3)), frac(1, sqrt(3)), frac(1, sqrt(3))), ")")),
    expression("μ" * " = " * bgroup("(", list(-frac(1, sqrt(3)), frac(1, sqrt(3)), frac(1, sqrt(3))), ")"))
  )
}

#' Encode one mu component as a readable filename token
build_bahadur_mu_component_slug <- function(x, tol = 1e-10) {
  inv_sqrt_3 <- 1 / sqrt(3)

  if (abs(x) < tol) {
    return("0")
  }
  if (abs(x - 1) < tol) {
    return("1")
  }
  if (abs(x + 1) < tol) {
    return("minus_1")
  }
  if (abs(x - inv_sqrt_3) < tol) {
    return("inv_sqrt_3")
  }
  if (abs(x + inv_sqrt_3) < tol) {
    return("minus_inv_sqrt_3")
  }

  component_slug <- format(round(x, 6), scientific = FALSE, trim = TRUE)
  component_slug <- gsub("-", "minus_", component_slug, fixed = TRUE)
  component_slug <- gsub("\\.", "p", component_slug)
  gsub("[^A-Za-z0-9_]", "", component_slug)
}

#' Build a readable filename suffix from mu
build_bahadur_mu_filename_suffix <- function(mu_true) {
  paste(vapply(mu_true, build_bahadur_mu_component_slug, character(1)), collapse = "_")
}

#' Identify the panel where the legend should be displayed
is_bahadur_legend_mu <- function(mu_true, tol = 1e-10) {
  reference_mu <- rep(1 / sqrt(3), 3)
  isTRUE(all.equal(as.numeric(mu_true), reference_mu, tolerance = tol))
}

#' Add default metadata when loading an old cache without filename information
enrich_bahadur_results_metadata <- function(combined_results) {
  if (all(c("mu_file_suffix", "show_legend") %in% names(combined_results))) {
    return(combined_results)
  }

  default_mu_list <- list(
    c(1, 0, 0),
    c(1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)),
    c(-1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3))
  )

  metadata <- data.frame(
    mu_index = seq_along(default_mu_list),
    mu_file_suffix = vapply(default_mu_list, build_bahadur_mu_filename_suffix, character(1)),
    show_legend = vapply(default_mu_list, is_bahadur_legend_mu, logical(1)),
    stringsAsFactors = FALSE
  )

  if (!("mu_file_suffix" %in% names(combined_results))) {
    combined_results <- combined_results %>%
      left_join(select(metadata, mu_index, mu_file_suffix), by = "mu_index")
  }

  if (!("show_legend" %in% names(combined_results))) {
    combined_results <- combined_results %>%
      left_join(select(metadata, mu_index, show_legend), by = "mu_index")
  }

  combined_results
}

#' Summarise pointwise Monte Carlo behaviour by mu and sample size
build_bahadur_pointwise_summary <- function(combined_results) {
  combined_results %>%
    group_by(n, mu_index, mu_label_factor, mu_file_suffix, show_legend) %>%
    summarise(
      mean_difference_norm = mean(difference_norm),
      sd_difference_norm = sd(difference_norm),
      n_mc = n(),
      .groups = "drop"
    ) %>%
    mutate(
      se_difference_norm = sd_difference_norm / sqrt(n_mc),
      ci_low = pmax(mean_difference_norm - 1.96 * se_difference_norm, 0),
      ci_high = mean_difference_norm + 1.96 * se_difference_norm
    )
}

#' Build the common theoretical benchmark 1/sqrt(n)
build_bahadur_theoretical_data <- function(combined_results) {
  x_range <- range(combined_results$n)
  theoretical_data <- data.frame(n = seq(x_range[1], x_range[2], length.out = 200))
  theoretical_data$inv_sqrt_n <- 1 / sqrt(theoretical_data$n)
  theoretical_data
}

#' Shared y-limits for log-log Bahadur plots
build_bahadur_loglog_limits <- function(pointwise_summary, theoretical_data, y_max) {
  means_df <- pointwise_summary %>%
    filter(mean_difference_norm > 0, n > 0) %>%
    mutate(log10_diff = log10(mean_difference_norm))

  ci_df <- pointwise_summary %>%
    filter(ci_low > 0, ci_high > 0, n > 0) %>%
    mutate(log10_ci_low = log10(ci_low))

  theo_df <- theoretical_data %>%
    filter(inv_sqrt_n > 0, n > 0) %>%
    mutate(log10_diff = log10(inv_sqrt_n))

  c(
    min(means_df$log10_diff, theo_df$log10_diff, ci_df$log10_ci_low, na.rm = TRUE),
    log10(y_max)
  )
}

#' Run heavy Bahadur simulation for one mu and cache results
#' @return Data frame with columns trajectory, n, difference_norm
run_bahadur_analysis <- function(mu_true, kappa_true = 1, n_trajectories = 500,
                                 n_cores = 10, min_n = 50, max_n = 100000,
                                 plot_suffix = "", seed = 123,
                                 output_dir = "output/bahadur/vmf",
                                 force_recompute = FALSE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(output_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  cache_file <- if (plot_suffix == "") {
    file.path(cache_dir, "vmf_bahadur_results_default.rds")
  } else {
    file.path(cache_dir, sprintf("vmf_bahadur_results_mu%s.rds", plot_suffix))
  }

  if (!force_recompute && file.exists(cache_file)) {
    cat("Loading cached results from", cache_file, "\n")
    return(readRDS(cache_file))
  }

  cat("=== Bahadur Representation Analysis for vMF Distribution ===\n\n")

  q <- length(mu_true) - 1
  xi_true <- kappa_true * mu_true
  sample_sizes <- build_bahadur_sample_sizes(min_n = min_n, max_n = max_n)

  cat("True parameters:\n")
  cat("Sphere dimension: S^", q, "\n", sep = "")
  cat("True mu:", paste(round(mu_true, 4), collapse = ", "), "\n")
  cat("True kappa:", kappa_true, "\n")
  cat("True xi:", paste(round(xi_true, 4), collapse = ", "), "\n\n")
  cat("Total sample points:", length(sample_sizes), "\n")
  cat("Number of trajectories:", n_trajectories, "\n")
  cat("Number of cores:", n_cores, "\n")

  set.seed(seed)
  cat("Setting up parallel cluster with", n_cores, "cores...\n")
  cl <- makeCluster(n_cores)

  trajectory_results <- tryCatch({
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(movMF)
        library(sphunif)
        library(pracma)
        library(rotasym)
      })
      source(file.path("utils.R"))
    })

    cat("Running", n_trajectories, "trajectories in parallel...\n")
    parLapply(cl, 1:n_trajectories, function(traj, mu, kappa, sizes, base_seed) {
      set.seed(base_seed + traj)
      analyze_single_trajectory(mu, kappa, sizes, traj)
    }, mu_true, kappa_true, sample_sizes, seed)
  }, finally = {
    if (exists("cl") && inherits(cl, "cluster")) {
      stopCluster(cl)
      cat("Cluster stopped successfully.\n")
    }
  })

  all_results <- do.call(rbind, trajectory_results)
  na_diff <- sum(is.na(all_results$difference_norm))
  inf_diff <- sum(is.infinite(all_results$difference_norm))

  if (na_diff > 0 || inf_diff > 0) {
    stop("Data contains NA or Inf values")
  }

  saveRDS(all_results, cache_file)
  cat("Cached results saved to", cache_file, "\n")
  all_results
}

#' Build one comparison plot from precomputed results
build_bahadur_comparison_plot <- function(combined_results, mu_labels, y_max = 0.10, loglog = FALSE) {
  pointwise_summary <- build_bahadur_pointwise_summary(combined_results)
  theoretical_data <- build_bahadur_theoretical_data(combined_results)

  if (!loglog) {
    return(
      ggplot(pointwise_summary, aes(x = n, y = mean_difference_norm)) +
        geom_ribbon(aes(x = n, ymin = ci_low, ymax = ci_high, fill = "95% CI"), alpha = 0.35, inherit.aes = FALSE) +
        geom_line(aes(color = "Mean"),
                  linetype = "dashed", linewidth = 0.5) +
        geom_line(data = theoretical_data, aes(x = n, y = inv_sqrt_n, color = "theoretical"),
                  linetype = "solid", linewidth = 0.5) +
        scale_color_manual(
          values = c(
            "Mean" = "red",
            "theoretical" = "black"
          ),
          labels = c(
            "Mean" = "Mean",
            "theoretical" = expression(frac(1, sqrt(n)))
          ),
          name = NULL
        ) +
        scale_fill_manual(
          values = c("95% CI" = "#b36565"),
          labels = c("95% CI" = "95% CI"),
          name = NULL
        ) +
        labs(x = expression("Sample size (" * italic(n) * ")"), y = NULL) +
        coord_cartesian(ylim = c(0, y_max)) +
        facet_wrap(~ mu_label_factor, nrow = 1, labeller = label_parsed) +
        theme_minimal() +
        theme(
          legend.position = c(0.65, 0.85),
          legend.justification = c("right", "top"),
          legend.text = element_text(size = 19),
          legend.key.size = unit(1.8, "lines"),
          strip.text = element_text(size = 17),
          axis.text.x = element_text(size = 19),
          axis.text.y = element_text(size = 19),
          axis.title.x = element_text(size = 19)
        )
    )
  }

  means_df <- pointwise_summary %>%
    filter(mean_difference_norm > 0, n > 0) %>%
    mutate(log10_n = log10(n), log10_diff = log10(mean_difference_norm))

  ci_df <- pointwise_summary %>%
    filter(ci_low > 0, ci_high > 0, n > 0) %>%
    mutate(
      log10_n = log10(n),
      log10_ci_low = log10(ci_low),
      log10_ci_high = log10(ci_high)
    )

  theo_df <- theoretical_data %>%
    filter(inv_sqrt_n > 0, n > 0) %>%
    mutate(log10_n = log10(n), log10_diff = log10(inv_sqrt_n))

  loglog_limits <- build_bahadur_loglog_limits(pointwise_summary, theoretical_data, y_max = y_max)

  ggplot(means_df, aes(x = log10_n, y = log10_diff)) +
    geom_ribbon(data = ci_df, aes(x = log10_n, ymin = log10_ci_low, ymax = log10_ci_high, fill = "95% CI"), alpha = 0.35, inherit.aes = FALSE) +
    geom_line(aes(color = "Mean"),
              linetype = "dashed", linewidth = 0.5) +
    geom_line(data = theo_df, aes(x = log10_n, y = log10_diff, color = "theoretical"),
              linetype = "solid", linewidth = 0.5) +
    scale_color_manual(
      values = c(
        "Mean" = "red",
        "theoretical" = "black"
      ),
      labels = c(
        "Mean" = expression(log[10]("Mean")),
        "theoretical" = expression(log[10] * bgroup("(", n^{-1/2}, ")"))
      ),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("95% CI" = "#b36565"),
      labels = c("95% CI" = "95% CI"),
      name = NULL
    ) +
    labs(x = expression(log[10](n)), y = NULL) +
    coord_cartesian(ylim = loglog_limits) +
    facet_wrap(~ mu_label_factor, nrow = 1, labeller = label_parsed) +
    theme_minimal() +
    theme(
      legend.position = c(0.05, 0.08),
      legend.justification = c("left", "bottom"),
      legend.text = element_text(size = 19),
      legend.key.size = unit(1.8, "lines"),
      strip.text = element_text(size = 17),
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19)
    )
}

#' Build one standalone linear-scale plot for a fixed mu, without title
build_bahadur_single_linear_plot <- function(combined_results, mu_index, y_max = 0.10) {
  pointwise_summary <- build_bahadur_pointwise_summary(combined_results)
  theoretical_data <- build_bahadur_theoretical_data(combined_results)

  pointwise_summary_mu <- pointwise_summary %>%
    filter(mu_index == !!mu_index)

  show_legend <- pointwise_summary_mu$show_legend[1]

  ggplot(pointwise_summary_mu, aes(x = n, y = mean_difference_norm)) +
    geom_ribbon(aes(x = n, ymin = ci_low, ymax = ci_high, fill = "95% CI"), alpha = 0.35, inherit.aes = FALSE) +
    geom_line(aes(color = "Mean"),
              linetype = "dashed", linewidth = 0.5) +
    geom_line(data = theoretical_data, aes(x = n, y = inv_sqrt_n, color = "theoretical"),
              linetype = "solid", linewidth = 0.5) +
    scale_color_manual(
      values = c(
        "Mean" = "red",
        "theoretical" = "black"
      ),
      labels = c(
        "Mean" = "Mean",
        "theoretical" = expression(frac(1, sqrt(n)))
      ),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("95% CI" = "#b36565"),
      labels = c("95% CI" = "95% CI"),
      name = NULL
    ) +
    labs(x = expression("Sample size (" * italic(n) * ")"), y = NULL) +
    coord_cartesian(ylim = c(0, y_max)) +
    theme_minimal() +
    theme(
      legend.position = if (isTRUE(show_legend)) c(0.65, 0.85) else "none",
      legend.justification = c("right", "top"),
      legend.text = element_text(size = 19),
      legend.key.size = unit(1.8, "lines"),
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19)
    )
}

#' Build one standalone log-log plot for a fixed mu, without title
build_bahadur_single_loglog_plot <- function(combined_results, mu_index, y_max = 0.10) {
  pointwise_summary <- build_bahadur_pointwise_summary(combined_results)
  theoretical_data <- build_bahadur_theoretical_data(combined_results)
  loglog_limits <- build_bahadur_loglog_limits(pointwise_summary, theoretical_data, y_max = y_max)

  pointwise_summary_mu <- pointwise_summary %>%
    filter(mu_index == !!mu_index)

  means_df <- pointwise_summary_mu %>%
    filter(mean_difference_norm > 0, n > 0) %>%
    mutate(log10_n = log10(n), log10_diff = log10(mean_difference_norm))

  ci_df <- pointwise_summary_mu %>%
    filter(ci_low > 0, ci_high > 0, n > 0) %>%
    mutate(
      log10_n = log10(n),
      log10_ci_low = log10(ci_low),
      log10_ci_high = log10(ci_high)
    )

  theo_df <- theoretical_data %>%
    filter(inv_sqrt_n > 0, n > 0) %>%
    mutate(log10_n = log10(n), log10_diff = log10(inv_sqrt_n))

  show_legend <- pointwise_summary_mu$show_legend[1]

  ggplot(means_df, aes(x = log10_n, y = log10_diff)) +
    geom_ribbon(data = ci_df, aes(x = log10_n, ymin = log10_ci_low, ymax = log10_ci_high, fill = "95% CI"), alpha = 0.35, inherit.aes = FALSE) +
    geom_line(aes(color = "Mean"),
              linetype = "dashed", linewidth = 0.5) +
    geom_line(data = theo_df, aes(x = log10_n, y = log10_diff, color = "theoretical"),
              linetype = "solid", linewidth = 0.5) +
    scale_color_manual(
      values = c(
        "Mean" = "red",
        "theoretical" = "black"
      ),
      labels = c(
        "Mean" = expression(log[10]("Mean")),
        "theoretical" = expression(log[10] * bgroup("(", n^{-1/2}, ")"))
      ),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("95% CI" = "#b36565"),
      labels = c("95% CI" = "95% CI"),
      name = NULL
    ) +
    labs(x = expression(log[10](n)), y = NULL) +
    coord_cartesian(ylim = loglog_limits) +
    theme_minimal() +
    theme(
      legend.position = if (isTRUE(show_legend)) c(0.05, 0.08) else "none",
      legend.justification = c("left", "bottom"),
      legend.text = element_text(size = 19),
      legend.key.size = unit(1.8, "lines"),
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      axis.title.x = element_text(size = 19)
    )
}

#' Generate all requested plot variants from precomputed combined results
generate_bahadur_plot_variants <- function(combined_results, output_dir = "output/bahadur/vmf") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  combined_results <- enrich_bahadur_results_metadata(combined_results)

  mu_labels <- build_default_bahadur_mu_labels()

  combined_results$mu_label_factor <- factor(combined_results$mu_index, levels = 1:3, labels = mu_labels)

  p_linear_010 <- build_bahadur_comparison_plot(combined_results, mu_labels, y_max = 0.10, loglog = FALSE)
  p_linear_005 <- build_bahadur_comparison_plot(combined_results, mu_labels, y_max = 0.05, loglog = FALSE)
  p_log_010 <- build_bahadur_comparison_plot(combined_results, mu_labels, y_max = 0.10, loglog = TRUE)
  p_log_005 <- build_bahadur_comparison_plot(combined_results, mu_labels, y_max = 0.05, loglog = TRUE)

  ggsave(file.path(output_dir, "vmf_bahadur_comparison_ylim_0p10.png"), p_linear_010, width = 18, height = 6, dpi = 600)
  ggsave(file.path(output_dir, "vmf_bahadur_comparison_ylim_0p05.png"), p_linear_005, width = 18, height = 6, dpi = 600)
  ggsave(file.path(output_dir, "vmf_bahadur_comparison_loglog_ymax_0p10.png"), p_log_010, width = 18, height = 6, dpi = 600)
  ggsave(file.path(output_dir, "vmf_bahadur_comparison_loglog_ymax_0p05.png"), p_log_005, width = 18, height = 6, dpi = 600)

  unique_mu_info <- combined_results %>%
    distinct(mu_index, mu_file_suffix)

  for (i in seq_len(nrow(unique_mu_info))) {
    mu_index_current <- unique_mu_info$mu_index[i]
    mu_file_suffix_current <- unique_mu_info$mu_file_suffix[i]
    p_single_linear <- build_bahadur_single_linear_plot(
      combined_results = combined_results,
      mu_index = mu_index_current,
      y_max = 0.10
    )
    p_single_loglog <- build_bahadur_single_loglog_plot(
      combined_results = combined_results,
      mu_index = mu_index_current,
      y_max = 0.10
    )
    ggsave(
      file.path(output_dir, sprintf("vmf_bahadur_linear_mu_%s.png", mu_file_suffix_current)),
      p_single_linear,
      width = 6,
      height = 6,
      dpi = 600
    )
    ggsave(
      file.path(output_dir, sprintf("vmf_bahadur_loglog_mu_%s.png", mu_file_suffix_current)),
      p_single_loglog,
      width = 6,
      height = 6,
      dpi = 600
    )
  }

  cat("Saved 4 comparison plot variants, 3 standalone linear plots, and 3 standalone log-log plots in", output_dir, "\n")
}

#' Regenerate plot variants from cached combined results (fast path, no simulation)
regenerate_bahadur_plots_from_cache <- function(output_dir = "output/bahadur/vmf") {
  combined_cache <- file.path(output_dir, "cache", "vmf_bahadur_comparison_results.rds")
  if (!file.exists(combined_cache)) {
    stop(sprintf("Cache file not found: %s", combined_cache))
  }
  combined_results <- readRDS(combined_cache)
  generate_bahadur_plot_variants(combined_results = combined_results, output_dir = output_dir)
  invisible(combined_results)
}

#' Run Bahadur analysis for 3 mu values, cache heavy results once, then plot
run_bahadur_comparison <- function(mu_list = list(
                                    c(1, 0, 0),
                                    c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)),
                                    c(-1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
                                  ),
                                  kappa_true = 1,
                                  n_trajectories = 500,
                                  n_cores = 10,
                                  min_n = 50,
                                  max_n = 100000,
                                  output_dir = "output/bahadur/vmf",
                                  force_recompute = FALSE) {
  cat("=== Bahadur Analysis Comparison for Multiple mu Values ===\n\n")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(output_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  all_comparison_results <- list()
  for (i in seq_along(mu_list)) {
    mu_current <- mu_list[[i]]
    mu_seed <- 123 + i * 1000
    cat("Processing mu", i, ":", paste(round(mu_current, 4), collapse = ", "), "\n")
    results <- run_bahadur_analysis(
      mu_true = mu_current,
      kappa_true = kappa_true,
      n_trajectories = n_trajectories,
      n_cores = n_cores,
      min_n = min_n,
      max_n = max_n,
      plot_suffix = i,
      seed = mu_seed,
      output_dir = output_dir,
      force_recompute = force_recompute
    )
    results$mu_index <- i
    results$mu_file_suffix <- build_bahadur_mu_filename_suffix(mu_current)
    results$show_legend <- is_bahadur_legend_mu(mu_current)
    all_comparison_results[[i]] <- results
  }

  combined_results <- do.call(rbind, all_comparison_results)
  combined_cache <- file.path(cache_dir, "vmf_bahadur_comparison_results.rds")
  saveRDS(combined_results, combined_cache)
  cat("Combined results cached to", combined_cache, "\n")

  generate_bahadur_plot_variants(combined_results = combined_results, output_dir = output_dir)
  combined_results
}

if (sys.nframe() == 0) {
  run_bahadur_comparison(force_recompute = TRUE)
}
