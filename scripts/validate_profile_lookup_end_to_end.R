#!/usr/bin/env Rscript

# Paired validation of experimental invariant profile lookup tables. The two
# paths use the same data, fit, centers, thresholds, information matrix, and
# multiplier draws. No production method selector is changed by this script.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("bootstrap/multiplier_bootstrap.R")
source("bootstrap/profile_lookup_interpolation.R")

parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- commandArgs(trailingOnly = TRUE)
  hit <- hit[startsWith(hit, prefix)]
  if (!length(hit)) default else substring(hit[[length(hit)]], nchar(prefix) + 1L)
}

n <- as.integer(parse_option("n", "200"))
B <- as.integer(parse_option("B", "4999"))
cores <- as.integer(parse_option("cores", "8"))
reference_cpp_kernel <- normalize_fast_multiplier_cpp_kernel(
  parse_option("reference-cpp-kernel", "legacy")
)
lookup_cpp_kernel <- normalize_fast_multiplier_cpp_kernel(
  parse_option("lookup-cpp-kernel", "contiguous_double")
)
data_seed <- as.integer(parse_option("data-seed", "2026080101"))
multiplier_seed <- as.integer(parse_option("multiplier-seed", "2026080102"))
output <- parse_option(
  "output",
  file.path("benchmarks", "profile_lookup_end_to_end_q10_n200_B4999.csv")
)
if (n < 20L || B < 99L || cores < 1L) stop("Invalid validation arguments.")

paired_metrics <- function(reference, comparison) {
  difference <- as.numeric(comparison) - as.numeric(reference)
  c(
    max_abs = max(abs(difference)),
    rmse = sqrt(mean(difference^2)),
    mean = mean(difference)
  )
}

metric_row <- function(model, quantity, statistic, reference, comparison,
                       comparison_path = "combined") {
  metrics <- paired_metrics(reference, comparison)
  data.frame(
    model = model,
    quantity = quantity,
    statistic = statistic,
    comparison_path = comparison_path,
    n_values = length(reference),
    max_abs_difference = metrics[["max_abs"]],
    rmse_difference = metrics[["rmse"]],
    mean_difference = metrics[["mean"]],
    reference_value = if (length(reference) == 1L) reference else NA_real_,
    lookup_value = if (length(comparison) == 1L) comparison else NA_real_,
    stringsAsFactors = FALSE
  )
}

fast_contribution_metrics <- function(centered_weights,
                                      score_matrix,
                                      correction_reference,
                                      correction_lookup,
                                      n,
                                      block_size = 100L) {
  difference_operator <- correction_lookup - correction_reference
  max_abs <- 0
  sum_squares <- 0
  sum_values <- 0
  count <- 0
  for (start in seq.int(1L, nrow(centered_weights), by = block_size)) {
    end <- min(start + block_size - 1L, nrow(centered_weights))
    score_block <- centered_weights[start:end, , drop = FALSE] %*% score_matrix
    difference <- -(score_block %*% t(difference_operator)) / sqrt(n)
    max_abs <- max(max_abs, max(abs(difference)))
    sum_squares <- sum_squares + sum(difference^2)
    sum_values <- sum_values + sum(difference)
    count <- count + length(difference)
  }
  c(
    max_abs = max_abs,
    rmse = sqrt(sum_squares / count),
    mean = sum_values / count,
    count = count
  )
}

build_table <- function(model) {
  if (identical(model, "vmf")) {
    return(profile_lookup_build(
      model = model,
      q = 10L,
      kappa_grid = seq(7, 14, length.out = 17L),
      geometry_grid = seq(-1, 1, length.out = 65L),
      t_grid = seq(0, pi, length.out = 129L),
      integration_grid_size = 16385L,
      cores = cores
    ))
  }
  profile_lookup_build(
    model = model,
    q = 10L,
    kappa_grid = seq(5, 12, length.out = 41L),
    geometry_grid = seq(0, 2.2, length.out = 65L),
    t_grid = seq(0, 3.6, length.out = 129L),
    integration_grid_size = 16385L,
    cores = cores
  )
}

run_case <- function(model, data, multiplier_matrix) {
  table <- build_table(model)
  spec <- if (identical(model, "vmf")) {
    make_vmf_spec(distance_type = "geodesic", unknown_param = "xi")
  } else {
    make_hvmf_spec(unknown_param = "both")
  }
  control <- if (identical(model, "vmf")) {
    list(
      derivative_method = "quadrature",
      vmf_profile_n_u = 32769L,
      vmf_derivative_n_u = 32769L
    )
  } else {
    list(
      derivative_method = "quadrature",
      hvmf_profile_n_y = 32769L,
      hvmf_derivative_n_y = 32769L
    )
  }
  theta_hat <- spec$fit_theta(
    data = data,
    weights = NULL,
    null = list(type = "composite"),
    control = control
  )
  data <- spec_normalize_data(spec, data, control)

  observed_start <- proc.time()[["elapsed"]]
  ks_prep <- prepare_ks_observed_data(
    data = data,
    spec = spec,
    theta_hat = theta_hat,
    ks_grid = make_sample_unique_distance_ks_grid(),
    control = control,
    light = TRUE,
    share_cvm_statistic = TRUE
  )
  observed_seconds <- proc.time()[["elapsed"]] - observed_start
  empirical_profile <- compute_sorted_empirical_profile_block(
    sorted_distance_matrix = ks_prep$sorted_distance_matrix,
    row_indices = seq_len(nrow(data))
  )
  direct_F_start <- proc.time()[["elapsed"]]
  direct_F <- compute_theoretical_sample_profile_sorted_block(
    spec = spec,
    normalized_data = data,
    sorted_distance_matrix = ks_prep$sorted_distance_matrix,
    theta = theta_hat,
    row_indices = seq_len(nrow(data)),
    control = control
  )
  direct_F_seconds <- proc.time()[["elapsed"]] - direct_F_start

  direct_start <- proc.time()[["elapsed"]]
  direct_fast <- spec_fast_multiplier_prepare(
    spec = spec,
    data = data,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    cvm_prep = prepare_cvm_observed_data_from_sample_ks(
      data = data,
      spec = spec,
      theta_hat = theta_hat,
      ks_prep = ks_prep,
      control = control
    ),
    control = control
  )
  direct_prep_seconds <- proc.time()[["elapsed"]] - direct_start

  lookup_start <- proc.time()[["elapsed"]]
  lookup <- profile_lookup_evaluate(
    table = table,
    xi = theta_hat$xi,
    centers = data,
    thresholds = ks_prep$sorted_distance_matrix
  )
  lookup_seconds <- proc.time()[["elapsed"]] - lookup_start

  direct_derivative <- as.matrix(direct_fast$D_ks$derivative_sorted)
  lookup_derivative <- as.matrix(lookup$derivative_sorted)
  direct_information <- as.matrix(direct_fast$Vhat)
  lookup_information <- direct_information
  information_inverse <- solve(direct_information)
  direct_correction <- direct_derivative %*% t(information_inverse)
  lookup_correction <- lookup_derivative %*% t(information_inverse)

  make_stream_prep <- function(derivative) {
    prepare_fast_ks_sample_stream_prep(
      S_obs = direct_fast$S_obs,
      Vhat = direct_information,
      Psi_aux = matrix(numeric(0), nrow = 0L, ncol = ncol(direct_fast$S_obs)),
      ks_prep = ks_prep,
      D_ks_info = list(
        mode = "sample_points_unique_distances",
        derivative_sorted = derivative
      ),
      control = list(fast_multiplier_cache_corrections = TRUE)
    )
  }
  direct_stream <- make_stream_prep(direct_derivative)
  lookup_stream <- make_stream_prep(lookup_derivative)
  centered_weights <- multiplier_matrix - 1

  run_kernel <- function(stream_prep, kernel) {
    started <- proc.time()[["elapsed"]]
    value <- compute_fast_sample_ks_cvm_stats_cpp(
      centered_weight_block = centered_weights,
      stream_prep = stream_prep,
      scale_factor = 1,
      compute_ks = TRUE,
      compute_cvm = TRUE,
      fuse_ks_cvm = TRUE,
      cpp_kernel = kernel
    )
    list(value = value, seconds = proc.time()[["elapsed"]] - started)
  }
  direct_legacy <- run_kernel(direct_stream, reference_cpp_kernel)
  direct_candidate <- run_kernel(direct_stream, lookup_cpp_kernel)
  lookup_legacy <- run_kernel(lookup_stream, reference_cpp_kernel)
  lookup_candidate <- run_kernel(lookup_stream, lookup_cpp_kernel)

  direct_process <- sqrt(nrow(data)) * (empirical_profile - direct_F)
  lookup_process <- sqrt(nrow(data)) * (empirical_profile - lookup$F)
  direct_observed <- c(ks = max(abs(direct_process)), cvm = mean(direct_process^2))
  lookup_observed <- c(ks = max(abs(lookup_process)), cvm = mean(lookup_process^2))
  direct_inference <- compute_inference_summary(
    as.list(direct_observed),
    direct_legacy$value,
    alpha = 0.05
  )
  lookup_inference <- compute_inference_summary(
    as.list(lookup_observed),
    lookup_candidate$value,
    alpha = 0.05
  )

  contribution <- fast_contribution_metrics(
    centered_weights = centered_weights,
    score_matrix = direct_fast$S_obs,
    correction_reference = direct_correction,
    correction_lookup = lookup_correction,
    n = nrow(data)
  )
  rows <- list(
    metric_row(model, "F_matrix", "all", direct_F, lookup$F, "table_only"),
    metric_row(model, "dot_F_matrix", "all", direct_derivative, lookup_derivative, "table_only"),
    metric_row(model, "information_matrix", "all", direct_information, lookup_information, "table_only"),
    metric_row(model, "dot_F_times_information_inverse", "all", direct_correction, lookup_correction, "table_only"),
    data.frame(
      model = model,
      quantity = "complete_fast_correction_contribution",
      statistic = "all",
      comparison_path = "table_only",
      n_values = contribution[["count"]],
      max_abs_difference = contribution[["max_abs"]],
      rmse_difference = contribution[["rmse"]],
      mean_difference = contribution[["mean"]],
      reference_value = NA_real_,
      lookup_value = NA_real_
    ),
    metric_row(model, "bootstrap_replicates", "ks",
               direct_legacy$value$ks, direct_candidate$value$ks, "kernel_only"),
    metric_row(model, "bootstrap_replicates", "cvm",
               direct_legacy$value$cvm, direct_candidate$value$cvm, "kernel_only"),
    metric_row(model, "bootstrap_replicates", "ks",
               direct_legacy$value$ks, lookup_legacy$value$ks, "table_only"),
    metric_row(model, "bootstrap_replicates", "cvm",
               direct_legacy$value$cvm, lookup_legacy$value$cvm, "table_only"),
    metric_row(model, "bootstrap_replicates", "ks",
               direct_legacy$value$ks, lookup_candidate$value$ks, "combined"),
    metric_row(model, "bootstrap_replicates", "cvm",
               direct_legacy$value$cvm, lookup_candidate$value$cvm, "combined")
  )
  for (statistic in c("ks", "cvm")) {
    rows <- c(rows, list(
      metric_row(
        model, "observed_statistic", statistic,
        direct_inference[[statistic]]$observed,
        lookup_inference[[statistic]]$observed,
        "combined"
      ),
      metric_row(
        model, "critical_value", statistic,
        direct_inference[[statistic]]$critical_value,
        lookup_inference[[statistic]]$critical_value,
        "combined"
      ),
      metric_row(
        model, "p_value", statistic,
        direct_inference[[statistic]]$p_value,
        lookup_inference[[statistic]]$p_value,
        "combined"
      )
    ))
  }
  summary <- do.call(rbind, rows)
  summary$n <- nrow(data)
  summary$B <- B
  summary$fitted_kappa <- theta_hat$kappa
  summary$reference_cpp_kernel <- reference_cpp_kernel
  summary$lookup_cpp_kernel <- lookup_cpp_kernel
  summary$table_build_seconds <- table$build_seconds
  summary$table_MiB <- table$bytes / 1024^2
  summary$observed_reference_seconds <- observed_seconds
  summary$direct_F_seconds <- direct_F_seconds
  summary$direct_derivative_prep_seconds <- direct_prep_seconds
  summary$lookup_F_and_derivative_seconds <- lookup_seconds
  summary$direct_legacy_kernel_seconds <- direct_legacy$seconds
  summary$direct_candidate_kernel_seconds <- direct_candidate$seconds
  summary$lookup_legacy_kernel_seconds <- lookup_legacy$seconds
  summary$lookup_candidate_kernel_seconds <- lookup_candidate$seconds

  list(
    summary = summary,
    theta_hat = theta_hat,
    table = table[c(
      "model", "q", "kappa_grid", "geometry_grid", "t_grid",
      "stencil_size", "integration_grid_size", "build_seconds", "bytes"
    )],
    inference = list(reference = direct_inference, lookup = lookup_inference),
    bootstrap = list(
      direct_legacy = direct_legacy$value,
      direct_candidate = direct_candidate$value,
      lookup_legacy = lookup_legacy$value,
      lookup_candidate = lookup_candidate$value
    )
  )
}

set.seed(data_seed)
vmf_data <- rotasym::r_vMF(
  n,
  mu = c(1, rep.int(0, 10L)),
  kappa = 10
)
hvmf_mu <- c(sqrt(2), 1, rep.int(0, 9L))
hvmf_data <- rhvmf_polar(n, mu = hvmf_mu, kappa = 10)
set.seed(multiplier_seed)
multiplier_matrix <- matrix(stats::rexp(B * n), nrow = B, ncol = n)

cases <- list(
  vmf = run_case("vmf", vmf_data, multiplier_matrix),
  hvmf = run_case("hvmf", hvmf_data, multiplier_matrix)
)
summary <- do.call(rbind, lapply(cases, `[[`, "summary"))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, output, row.names = FALSE)
details_output <- sub("\\.csv$", ".rds", output)
saveRDS(list(
  configuration = list(
    n = n,
    B = B,
    cores = cores,
    sample_seed = data_seed,
    multiplier_seed = multiplier_seed,
    reference_cpp_kernel = reference_cpp_kernel,
    lookup_cpp_kernel = lookup_cpp_kernel
  ),
  cases = cases,
  summary = summary
), details_output)
print(summary, row.names = FALSE)
cat(sprintf("Wrote %s and %s\n", output, details_output))
