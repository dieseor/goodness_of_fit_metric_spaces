#!/usr/bin/env Rscript

# Long, paired calibration/precision pilot for the vMF and HvMF scenarios
# appearing in the active AoS Section 6 table.  It is deliberately separate
# from the production runner and never writes under simulation_results/final*.
#
# Every comparison within a task uses the same simulated sample, the same MLE,
# and the same B exponential multiplier vectors.  The only changing quantity
# is the derivative calculation (grid resolution, historical score-MC size,
# or—in the small kernel audit—the accumulation kernel).

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

source("scripts/run_section6_new_scenarios.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

pilot_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  value <- commandArgs(trailingOnly = TRUE)
  value <- value[startsWith(value, prefix)]
  if (!length(value)) default else substring(value[[length(value)]], nchar(prefix) + 1L)
}

pilot_integer <- function(name, default, minimum = 0L) {
  value <- suppressWarnings(as.integer(pilot_option(name, default)))
  if (length(value) != 1L || !is.finite(value) || value < minimum) {
    stop(sprintf("`--%s` must be an integer of at least %d.", name, minimum), call. = FALSE)
  }
  value
}

pilot_logical <- function(name, default = FALSE) {
  value <- tolower(pilot_option(name, if (isTRUE(default)) "true" else "false"))
  if (!value %in% c("true", "false")) {
    stop(sprintf("`--%s` must be true or false.", name), call. = FALSE)
  }
  identical(value, "true")
}

fixed_exp_multipliers <- function(matrix_draws) {
  matrix_draws <- as.matrix(matrix_draws)
  cursor <- 0L
  list(
    name = "paired Exp(1) multipliers",
    mean = 1,
    sd = 1,
    generator = function(n) {
      cursor <<- cursor + 1L
      if (cursor > nrow(matrix_draws) || n != ncol(matrix_draws)) {
        stop("Paired multiplier matrix was consumed inconsistently.", call. = FALSE)
      }
      matrix_draws[cursor, ]
    }
  )
}

pilot_theta_signature <- function(theta) {
  values <- as.numeric(unlist(theta, recursive = TRUE, use.names = FALSE))
  paste(formatC(values, digits = 16L, format = "fg", flag = "#"), collapse = ";")
}

pilot_seed <- function(base_seed, design_id, rep, stream) {
  section6_seed(base_seed, design_id = design_id, rep = rep, stream = stream)
}

pilot_mode_specifications <- function(panel) {
  switch(panel,
    main = list(list(
      mode = "integral_4097", derivative_method = "quadrature", grid_size = 4097L,
      kernel = "contiguous_double", derivative_mc_size = NA_integer_
    )),
    refine = list(
      list(mode = "integral_4097", derivative_method = "quadrature", grid_size = 4097L,
           kernel = "contiguous_double", derivative_mc_size = NA_integer_),
      list(mode = "integral_16385", derivative_method = "quadrature", grid_size = 16385L,
           kernel = "contiguous_double", derivative_mc_size = NA_integer_)
    ),
    ultra = list(
      list(mode = "integral_4097", derivative_method = "quadrature", grid_size = 4097L,
           kernel = "contiguous_double", derivative_mc_size = NA_integer_),
      list(mode = "integral_16385", derivative_method = "quadrature", grid_size = 16385L,
           kernel = "contiguous_double", derivative_mc_size = NA_integer_),
      list(mode = "integral_32769", derivative_method = "quadrature", grid_size = 32769L,
           kernel = "contiguous_double", derivative_mc_size = NA_integer_)
    ),
    mc = list(
      list(mode = "integral_4097", derivative_method = "quadrature", grid_size = 4097L,
           kernel = "contiguous_double", derivative_mc_size = NA_integer_),
      list(mode = "score_mc_1000", derivative_method = "score_mc", grid_size = NA_integer_,
           kernel = "contiguous_double", derivative_mc_size = 1000L),
      list(mode = "score_mc_10000", derivative_method = "score_mc", grid_size = NA_integer_,
           kernel = "contiguous_double", derivative_mc_size = 10000L)
    ),
    kernel = list(
      list(mode = "integral_4097", derivative_method = "quadrature", grid_size = 4097L,
           kernel = "contiguous_double", derivative_mc_size = NA_integer_),
      list(mode = "integral_4097_legacy_kernel", derivative_method = "quadrature", grid_size = 4097L,
           kernel = "legacy", derivative_mc_size = NA_integer_)
    ),
    stop(sprintf("Unknown pilot panel '%s'.", panel), call. = FALSE)
  )
}

pilot_empty_results <- function() {
  data.frame(
    panel = character(), mode = character(), scenario = character(), family = character(),
    d = integer(), n = integer(), rep = integer(), status = character(), error_message = character(),
    data_seed = integer(), multiplier_seed = integer(), derivative_seed = integer(),
    grid_size = integer(), derivative_mc_size = integer(), kernel = character(),
    theta_signature = character(),
    ks_observed = numeric(), ks_critical = numeric(), ks_pvalue = numeric(), ks_reject = logical(),
    cvm_observed = numeric(), cvm_critical = numeric(), cvm_pvalue = numeric(), cvm_reject = logical(),
    ks_bootstrap_mean = numeric(), ks_bootstrap_sd = numeric(), ks_bootstrap_q95 = numeric(),
    cvm_bootstrap_mean = numeric(), cvm_bootstrap_sd = numeric(), cvm_bootstrap_q95 = numeric(),
    derivative_method_effective = character(), bootstrap_method_effective = character(),
    fast_multiplier_backend_effective = character(), fast_multiplier_cpp_kernel_effective = character(),
    fast_multiplier_fuse_ks_cvm_effective = logical(), elapsed_seconds = numeric(),
    stringsAsFactors = FALSE
  )
}

pilot_empty_comparisons <- function() {
  data.frame(
    panel = character(), scenario = character(), family = character(), d = integer(), n = integer(), rep = integer(),
    reference_mode = character(), comparison_mode = character(), statistic = character(),
    max_abs_bootstrap_difference = numeric(), rmse_bootstrap_difference = numeric(),
    observed_difference = numeric(), critical_difference = numeric(), pvalue_difference = numeric(),
    stringsAsFactors = FALSE
  )
}

pilot_control <- function(family, mode, derivative_seed, cvm_block_size) {
  control <- section6_control(
    derivative_mc_size = if (is.na(mode$derivative_mc_size)) 1000L else mode$derivative_mc_size,
    derivative_seed = derivative_seed,
    cvm_block_size = cvm_block_size,
    derivative_method = mode$derivative_method
  )
  control$fast_multiplier_cpp_kernel <- mode$kernel
  # Set both aliases because the derivative code intentionally prioritises
  # vmf/hvmf_derivative_n_* over vmf/hvmf_profile_n_*.
  if (identical(mode$derivative_method, "quadrature")) {
    if (identical(family, "vmf")) {
      control$vmf_profile_n_u <- mode$grid_size
      control$vmf_derivative_n_u <- mode$grid_size
    } else {
      control$hvmf_profile_n_y <- mode$grid_size
      control$hvmf_derivative_n_y <- mode$grid_size
    }
  }
  control
}

pilot_fit_one_mode <- function(data, spec, theta_hat, multipliers, B, mode,
                               derivative_seed, cvm_block_size) {
  family <- if (startsWith(spec$name, "vmf_")) "vmf" else "hvmf"
  control <- pilot_control(family, mode, derivative_seed, cvm_block_size)
  started <- proc.time()[["elapsed"]]
  result <- multiplier_bootstrap_gof(
    data = data,
    spec = spec,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    multipliers = fixed_exp_multipliers(multipliers),
    n_cores = 1L,
    observed_theta_hat = theta_hat,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = TRUE, bootstrap_thetas = FALSE),
    control = control
  )
  list(result = result, elapsed_seconds = proc.time()[["elapsed"]] - started)
}

pilot_result_row <- function(task, mode, seeds, theta_text, fitted) {
  result <- fitted$result
  diagnostics <- result$diagnostics
  inference <- result$inference
  boot_ks <- as.numeric(result$bootstrap$statistics$ks)
  boot_cvm <- as.numeric(result$bootstrap$statistics$cvm)
  effective_derivative <- diagnostics$derivative_method_effective %||%
    diagnostics$derivative_method %||% NA_character_
  effective_bootstrap <- diagnostics$effective_bootstrap_method %||% NA_character_
  effective_backend <- diagnostics$fast_multiplier_backend_effective %||% NA_character_
  effective_kernel <- diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_
  effective_fused <- isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective)
  conforming <- identical(effective_derivative, mode$derivative_method) &&
    identical(effective_bootstrap, "fast_multiplier") &&
    identical(effective_backend, "cpp") &&
    identical(effective_kernel, mode$kernel) && effective_fused
  data.frame(
    panel = task$panel, mode = mode$mode, scenario = task$scenario, family = task$family,
    d = task$d, n = task$n, rep = task$rep,
    status = if (conforming) "ok" else "nonconforming",
    error_message = if (conforming) NA_character_ else paste(
      sprintf("derivative=%s", effective_derivative), sprintf("bootstrap=%s", effective_bootstrap),
      sprintf("backend=%s", effective_backend), sprintf("kernel=%s", effective_kernel),
      sprintf("fused=%s", effective_fused), sep = "; "
    ),
    data_seed = seeds$data, multiplier_seed = seeds$multiplier, derivative_seed = seeds$derivative,
    grid_size = if (is.na(mode$grid_size)) NA_integer_ else as.integer(mode$grid_size),
    derivative_mc_size = if (is.na(mode$derivative_mc_size)) NA_integer_ else as.integer(mode$derivative_mc_size),
    kernel = mode$kernel, theta_signature = theta_text,
    ks_observed = inference$ks$observed, ks_critical = inference$ks$critical_value,
    ks_pvalue = inference$ks$p_value, ks_reject = inference$ks$reject,
    cvm_observed = inference$cvm$observed, cvm_critical = inference$cvm$critical_value,
    cvm_pvalue = inference$cvm$p_value, cvm_reject = inference$cvm$reject,
    ks_bootstrap_mean = mean(boot_ks), ks_bootstrap_sd = stats::sd(boot_ks),
    ks_bootstrap_q95 = as.numeric(stats::quantile(boot_ks, 0.95, names = FALSE, type = 7L)),
    cvm_bootstrap_mean = mean(boot_cvm), cvm_bootstrap_sd = stats::sd(boot_cvm),
    cvm_bootstrap_q95 = as.numeric(stats::quantile(boot_cvm, 0.95, names = FALSE, type = 7L)),
    derivative_method_effective = effective_derivative,
    bootstrap_method_effective = effective_bootstrap,
    fast_multiplier_backend_effective = effective_backend,
    fast_multiplier_cpp_kernel_effective = effective_kernel,
    fast_multiplier_fuse_ks_cvm_effective = effective_fused,
    elapsed_seconds = fitted$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

pilot_comparison_rows <- function(task, fitted, reference_mode = "integral_4097") {
  reference <- fitted[[reference_mode]]$result
  output <- list()
  index <- 1L
  for (mode_name in setdiff(names(fitted), reference_mode)) {
    comparison <- fitted[[mode_name]]$result
    for (statistic in c("ks", "cvm")) {
      reference_boot <- as.numeric(reference$bootstrap$statistics[[statistic]])
      comparison_boot <- as.numeric(comparison$bootstrap$statistics[[statistic]])
      difference <- comparison_boot - reference_boot
      output[[index]] <- data.frame(
        panel = task$panel, scenario = task$scenario, family = task$family,
        d = task$d, n = task$n, rep = task$rep,
        reference_mode = reference_mode, comparison_mode = mode_name, statistic = statistic,
        max_abs_bootstrap_difference = max(abs(difference)),
        rmse_bootstrap_difference = sqrt(mean(difference^2)),
        observed_difference = comparison$inference[[statistic]]$observed -
          reference$inference[[statistic]]$observed,
        critical_difference = comparison$inference[[statistic]]$critical_value -
          reference$inference[[statistic]]$critical_value,
        pvalue_difference = comparison$inference[[statistic]]$p_value -
          reference$inference[[statistic]]$p_value,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  if (!length(output)) pilot_empty_comparisons() else do.call(rbind, output)
}

pilot_run_task <- function(task, B, base_seed, cvm_block_size) {
  modes <- pilot_mode_specifications(task$panel)
  seeds <- list(
    data = pilot_seed(base_seed, task$design_id, task$rep, 0L),
    multiplier = pilot_seed(base_seed, task$design_id, task$rep, 1L),
    derivative = pilot_seed(base_seed, task$design_id, task$rep, 2L)
  )
  started <- proc.time()[["elapsed"]]
  blank_rows <- function(error_message) do.call(rbind, lapply(modes, function(mode) {
    row <- pilot_empty_results()[0L, , drop = FALSE]
    row[1L, ] <- NA
    row$panel <- task$panel; row$mode <- mode$mode; row$scenario <- task$scenario
    row$family <- task$family; row$d <- task$d; row$n <- task$n; row$rep <- task$rep
    row$status <- "error"; row$error_message <- error_message
    row$data_seed <- seeds$data; row$multiplier_seed <- seeds$multiplier; row$derivative_seed <- seeds$derivative
    row$grid_size <- if (is.na(mode$grid_size)) NA_integer_ else mode$grid_size
    row$derivative_mc_size <- if (is.na(mode$derivative_mc_size)) NA_integer_ else mode$derivative_mc_size
    row$kernel <- mode$kernel; row$elapsed_seconds <- NA_real_
    row
  }))
  tryCatch({
    set.seed(seeds$data)
    data <- generate_section6_sample(task)
    spec <- if (identical(task$family, "vmf")) {
      make_vmf_spec(distance_type = "geodesic", unknown_param = "xi")
    } else {
      make_hvmf_spec(unknown_param = "both")
    }
    normalized <- spec_normalize_data(spec, data, control = list())
    theta_hat <- spec$fit_theta(
      data = normalized, weights = NULL, null = list(type = "composite"), control = list()
    )
    theta_text <- pilot_theta_signature(theta_hat)
    set.seed(seeds$multiplier)
    multiplier_matrix <- matrix(stats::rexp(B * nrow(normalized)), nrow = B, ncol = nrow(normalized))
    fitted <- list()
    rows <- list()
    for (mode in modes) {
      one <- tryCatch(
        pilot_fit_one_mode(normalized, spec, theta_hat, multiplier_matrix, B, mode,
                           seeds$derivative, cvm_block_size),
        error = identity
      )
      if (inherits(one, "error")) {
        stop(sprintf("mode %s failed: %s", mode$mode, conditionMessage(one)), call. = FALSE)
      }
      fitted[[mode$mode]] <- one
      rows[[length(rows) + 1L]] <- pilot_result_row(task, mode, seeds, theta_text, one)
    }
    list(
      rows = do.call(rbind, rows),
      comparisons = pilot_comparison_rows(task, fitted),
      task_elapsed_seconds = proc.time()[["elapsed"]] - started
    )
  }, error = function(error) {
    list(
      rows = blank_rows(conditionMessage(error)),
      comparisons = pilot_empty_comparisons(),
      task_elapsed_seconds = proc.time()[["elapsed"]] - started
    )
  })
}

pilot_write_csv_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  utils::write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop(sprintf("Could not update '%s'.", path), call. = FALSE)
}

pilot_task_key <- function(x) {
  paste(x$panel, x$scenario, x$d, x$n, x$rep, sep = "|")
}

pilot_summarise <- function(results) {
  results <- results[results$status == "ok", , drop = FALSE]
  if (!nrow(results)) return(data.frame())
  groups <- split(results, interaction(results$panel, results$mode, results$scenario,
                                       results$d, results$n, drop = TRUE))
  do.call(rbind, lapply(groups, function(x) {
    ci <- function(reject) {
      if (!length(reject)) return(c(NA_real_, NA_real_))
      stats::binom.test(sum(reject), length(reject), conf.level = 0.95)$conf.int
    }
    ks_ci <- ci(x$ks_reject); cvm_ci <- ci(x$cvm_reject)
    data.frame(
      panel = x$panel[[1L]], mode = x$mode[[1L]], scenario = x$scenario[[1L]], family = x$family[[1L]],
      d = x$d[[1L]], n = x$n[[1L]], replications = nrow(x),
      rejection_ks = mean(x$ks_reject), rejection_ks_ci_lower = ks_ci[[1L]], rejection_ks_ci_upper = ks_ci[[2L]],
      rejection_cvm = mean(x$cvm_reject), rejection_cvm_ci_lower = cvm_ci[[1L]], rejection_cvm_ci_upper = cvm_ci[[2L]],
      mean_ks_pvalue = mean(x$ks_pvalue), median_ks_pvalue = stats::median(x$ks_pvalue),
      mean_cvm_pvalue = mean(x$cvm_pvalue), median_cvm_pvalue = stats::median(x$cvm_pvalue),
      mean_elapsed_seconds = mean(x$elapsed_seconds), stringsAsFactors = FALSE
    )
  }))
}

pilot_summarise_comparisons <- function(comparisons) {
  if (!nrow(comparisons)) return(data.frame())
  groups <- split(comparisons, interaction(comparisons$panel, comparisons$scenario,
                                            comparisons$d, comparisons$n,
                                            comparisons$reference_mode,
                                            comparisons$comparison_mode,
                                            comparisons$statistic, drop = TRUE))
  do.call(rbind, lapply(groups, function(x) data.frame(
    panel = x$panel[[1L]], scenario = x$scenario[[1L]], family = x$family[[1L]],
    d = x$d[[1L]], n = x$n[[1L]], reference_mode = x$reference_mode[[1L]],
    comparison_mode = x$comparison_mode[[1L]], statistic = x$statistic[[1L]],
    paired_replications = nrow(x),
    max_abs_bootstrap_difference = max(x$max_abs_bootstrap_difference),
    mean_rmse_bootstrap_difference = mean(x$rmse_bootstrap_difference),
    max_abs_observed_difference = max(abs(x$observed_difference)),
    max_abs_critical_difference = max(abs(x$critical_difference)),
    max_abs_pvalue_difference = max(abs(x$pvalue_difference)),
    mean_pvalue_difference = mean(x$pvalue_difference),
    stringsAsFactors = FALSE
  )))
}

pilot_write_status <- function(path, total, completed, results, started, active, cap_reached) {
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  task_count <- length(unique(pilot_task_key(results)))
  task_seconds <- if (task_count) elapsed / max(1L, completed) else NA_real_
  eta <- if (is.finite(task_seconds) && completed > 0L) task_seconds * (total - completed) / 10 else NA_real_
  writeLines(c(
    sprintf("updated: %s", format(Sys.time(), tz = "Europe/Madrid")),
    sprintf("completed_tasks: %d/%d", completed, total),
    sprintf("active_tasks: %d", active),
    sprintf("elapsed_seconds: %.1f", elapsed),
    sprintf("rough_eta_seconds_at_10_cores: %s", if (is.finite(eta)) format(round(eta), scientific = FALSE) else "NA"),
    sprintf("successful_mode_runs: %d", sum(results$status == "ok")),
    sprintf("nonconforming_mode_runs: %d", sum(results$status == "nonconforming")),
    sprintf("failed_mode_runs: %d", sum(results$status == "error")),
    sprintf("wall_time_cap_reached: %s", cap_reached)
  ), path)
}

run_aos_vmf_hvmf_integral_precision_pilot <- function(
    output_dir,
    main_M = 120L,
    refine_M = 20L,
    ultra_M = 8L,
    mc_M = 12L,
    kernel_M = 8L,
    B = 5000L,
    cores = 10L,
    base_seed = 20260801L,
    cvm_block_size = 50L,
    max_wall_minutes = 230L,
    checkpoint_tasks = 10L) {
  integer_values <- as.integer(c(main_M, refine_M, ultra_M, mc_M, kernel_M, B, cores,
                                 base_seed, cvm_block_size, max_wall_minutes, checkpoint_tasks))
  if (any(!is.finite(integer_values)) || any(integer_values[1:5] < 0L) ||
      any(integer_values[c(6:7, 9, 11)] < 1L) || max_wall_minutes < 0L) {
    stop("Invalid pilot controls.", call. = FALSE)
  }
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("This pilot requires Unix for outer parallelism when cores > 1.", call. = FALSE)
  }
  if (cores != 10L) warning("The requested pilot was calibrated for 10 outer workers.", call. = FALSE)

  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(lock), add = TRUE)

  raw_path <- file.path(output_dir, "raw_results.csv")
  comparison_path <- file.path(output_dir, "paired_bootstrap_comparisons.csv")
  summary_path <- file.path(output_dir, "rejection_summary.csv")
  comparison_summary_path <- file.path(output_dir, "paired_comparison_summary.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  manifest_path <- file.path(output_dir, "manifest.txt")
  log_path <- file.path(output_dir, "run.log")

  design <- rbind(
    make_section6_design("vmf", dimensions = c(2L, 10L), n_values = c(50L, 100L, 200L, 400L), beta_values = 0),
    make_section6_design("hvmf", dimensions = c(2L, 10L), n_values = c(50L, 100L, 200L, 400L), beta_values = 0)
  )
  # section6_seed assumes distinct design ids.  The individual family designs
  # each start at 1, so make the IDs globally distinct in this pilot only.
  design$design_id <- seq_len(nrow(design))

  make_panel <- function(name, M, subset) {
    rows <- design[subset, , drop = FALSE]
    if (!nrow(rows) || M == 0L) {
      out <- rows[0L, c("scenario", "family", "d", "n", "beta", "design_id"), drop = FALSE]
      out$rep <- integer()
      out$panel <- character()
      return(out)
    }
    out <- merge(rows, data.frame(rep = seq_len(M)), by = NULL)
    out$panel <- name
    out
  }
  panels <- list(
    make_panel("main", main_M, rep(TRUE, nrow(design))),
    make_panel("refine", refine_M, design$d == 10L),
    make_panel("ultra", ultra_M, design$d == 10L & design$n == 200L),
    make_panel("mc", mc_M, design$d == 10L & design$n == 200L),
    make_panel("kernel", kernel_M, design$d == 10L & design$n == 200L)
  )
  tasks <- do.call(rbind, panels)
  # Give every panel early coverage: within each replication the 32 default
  # cells come first, followed by the paired precision checks.  Thus a wall
  # cap still leaves both a broad null-calibration sample and completed grid
  # comparisons, rather than only the first part of the default panel.
  panel_order <- match(tasks$panel, c("main", "refine", "ultra", "mc", "kernel"))
  tasks <- tasks[order(tasks$rep, panel_order, tasks$family, tasks$scenario, tasks$d, tasks$n), , drop = FALSE]
  rownames(tasks) <- NULL

  manifest <- c(
    sprintf("created_at: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z")),
    "purpose: paired AoS vMF/HvMF null calibration and integral-grid precision pilot",
    "aos_table_design: Section 6 vMF/HvMF, d=2,10; n=50,100,200,400; beta=0",
    "production_result_directories_modified: FALSE",
    "tex_modified: FALSE",
    "integral_method_label_in_code: quadrature",
    "integral_default_grid_size: 4097",
    "strict_integral_grid_size: 16385",
    "ultra_integral_grid_size: 32769",
    "same_data_mle_and_multipliers_within_task: TRUE",
    "multiplier_law: Exp(1)",
    sprintf("B: %d", B), sprintf("cores: %d", cores), sprintf("base_seed: %d", base_seed),
    sprintf("main_M: %d", main_M), sprintf("refine_M: %d", refine_M),
    sprintf("ultra_M: %d", ultra_M), sprintf("mc_M: %d", mc_M), sprintf("kernel_M: %d", kernel_M),
    sprintf("max_wall_minutes: %d", max_wall_minutes),
    sprintf("git_head: %s", tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) "unavailable"))
  )
  if (file.exists(manifest_path)) {
    existing_manifest <- readLines(manifest_path, warn = FALSE)
    compare_prefixes <- c("aos_table_design:", "B:", "base_seed:", "main_M:", "refine_M:", "ultra_M:", "mc_M:", "kernel_M:")
    for (prefix in compare_prefixes) {
      old <- existing_manifest[startsWith(existing_manifest, prefix)]
      new <- manifest[startsWith(manifest, prefix)]
      if (length(old) != 1L || length(new) != 1L || !identical(old, new)) {
        stop("Existing output directory has an incompatible manifest; choose another --output-dir.", call. = FALSE)
      }
    }
  } else {
    writeLines(manifest, manifest_path)
  }

  results <- if (file.exists(raw_path)) utils::read.csv(raw_path, stringsAsFactors = FALSE) else pilot_empty_results()
  comparisons <- if (file.exists(comparison_path)) utils::read.csv(comparison_path, stringsAsFactors = FALSE) else pilot_empty_comparisons()
  missing_result_columns <- setdiff(names(pilot_empty_results()), names(results))
  missing_comparison_columns <- setdiff(names(pilot_empty_comparisons()), names(comparisons))
  if (length(missing_result_columns) || length(missing_comparison_columns)) {
    stop("Existing CSV schema is incompatible; use another --output-dir.", call. = FALSE)
  }
  results <- results[, names(pilot_empty_results()), drop = FALSE]
  comparisons <- comparisons[, names(pilot_empty_comparisons()), drop = FALSE]
  expected_mode_count <- vapply(tasks$panel, function(x) length(pilot_mode_specifications(x)), integer(1))
  existing_ok <- results[results$status == "ok", , drop = FALSE]
  counts <- table(pilot_task_key(existing_ok))
  task_keys <- pilot_task_key(tasks)
  done <- names(counts)[counts == expected_mode_count[match(names(counts), task_keys)]]
  # Count-based resumption alone is insufficient if a stale row has another mode.
  is_complete <- vapply(done, function(key) {
    task <- tasks[task_keys == key, , drop = FALSE]
    wanted <- vapply(pilot_mode_specifications(task$panel[[1L]]), `[[`, character(1), "mode")
    got <- existing_ok$mode[pilot_task_key(existing_ok) == key]
    setequal(wanted, got) && length(got) == length(wanted)
  }, logical(1))
  done <- names(is_complete)[is_complete]
  pending <- tasks[!task_keys %in% done, , drop = FALSE]

  write_outputs <- function() {
    results <<- results[order(results$panel, results$family, results$scenario, results$d, results$n,
                               results$rep, results$mode), , drop = FALSE]
    comparisons <<- comparisons[order(comparisons$panel, comparisons$family, comparisons$scenario,
                                      comparisons$d, comparisons$n, comparisons$rep,
                                      comparisons$comparison_mode, comparisons$statistic), , drop = FALSE]
    pilot_write_csv_atomic(results, raw_path)
    pilot_write_csv_atomic(comparisons, comparison_path)
    pilot_write_csv_atomic(pilot_summarise(results), summary_path)
    pilot_write_csv_atomic(pilot_summarise_comparisons(comparisons), comparison_summary_path)
  }

  ensure_distance_profile_cpp_loaded()
  total <- nrow(tasks)
  completed <- total - nrow(pending)
  started <- Sys.time()
  cap_seconds <- 60 * max_wall_minutes
  cap_reached <- FALSE
  last_checkpoint <- completed
  active <- list()
  active_index <- list()
  next_index <- 1L
  max_active <- min(cores, nrow(pending))
  last_progress_seconds <- -Inf
  report_progress <- function(force = FALSE) {
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (!force && elapsed - last_progress_seconds < 10) return(invisible(FALSE))
    last_progress_seconds <<- elapsed
    observed_elapsed <- results$elapsed_seconds[results$status == "ok"]
    mean_task_seconds <- if (length(observed_elapsed)) mean(observed_elapsed) else NA_real_
    # Mode-times are conservative because a paired task contains two or three modes.
    mean_modes <- if (completed > 0L) nrow(results) / completed else NA_real_
    eta <- if (is.finite(mean_task_seconds) && is.finite(mean_modes)) {
      mean_task_seconds * mean_modes * (total - completed) / max(1L, cores)
    } else NA_real_
    bar_width <- 30L
    filled <- floor(bar_width * completed / max(1L, total))
    bar <- paste0(strrep("#", filled), strrep("-", bar_width - filled))
    cat(sprintf("\r[%s] %d/%d tasks | %d active | %.1f min | ETA %s%s",
      bar, completed, total, length(active), elapsed / 60,
      if (is.finite(eta)) sprintf("%.1f min", eta / 60) else "estimating",
      if (cap_reached) " | wall cap reached" else ""))
    flush.console()
    invisible(TRUE)
  }
  cat(sprintf("%s AoS vMF/HvMF integral precision pilot: %d tasks pending, B=%d, cores=%d\n",
              format(started, tz = "Europe/Madrid"), nrow(pending), B, cores), file = log_path, append = TRUE)
  pilot_write_status(status_path, total, completed, results, started, length(active), cap_reached)
  report_progress(force = TRUE)

  collect_one <- function(pid, value, index) {
    if (inherits(value, "try-error")) {
      value <- list(rows = pilot_empty_results(), comparisons = pilot_empty_comparisons())
      cat(sprintf("\nERROR task=%s: %s\n", pilot_task_key(pending[index, , drop = FALSE]), as.character(value)), file = log_path, append = TRUE)
    }
    results <<- rbind(results, value$rows)
    comparisons <<- rbind(comparisons, value$comparisons)
    active[[pid]] <<- NULL
    active_index[[pid]] <<- NULL
    completed <<- completed + 1L
  }

  if (.Platform$OS.type != "unix" || cores <= 1L) {
    for (index in seq_len(nrow(pending))) {
      if (cap_seconds > 0 && as.numeric(difftime(Sys.time(), started, units = "secs")) >= cap_seconds) {
        cap_reached <- TRUE; break
      }
      value <- pilot_run_task(pending[index, , drop = FALSE], B, base_seed, cvm_block_size)
      collect_one("serial", value, index)
      if (completed - last_checkpoint >= checkpoint_tasks) {
        write_outputs(); last_checkpoint <- completed
        pilot_write_status(status_path, total, completed, results, started, 0L, cap_reached)
      }
      report_progress(force = TRUE)
    }
  } else {
    repeat {
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      if (cap_seconds > 0 && elapsed >= cap_seconds) cap_reached <- TRUE
      while (!cap_reached && next_index <= nrow(pending) && length(active) < max_active) {
        index <- next_index
        job <- parallel::mcparallel(
          pilot_run_task(pending[index, , drop = FALSE], B, base_seed, cvm_block_size),
          mc.set.seed = FALSE
        )
        pid <- as.character(job$pid)
        active[[pid]] <- job
        active_index[[pid]] <- index
        next_index <- next_index + 1L
      }
      if (!length(active)) {
        if (next_index > nrow(pending) || cap_reached) break
        next
      }
      collected <- parallel::mccollect(active, wait = FALSE)
      if (is.null(collected)) {
        Sys.sleep(0.5)
        report_progress()
        next
      }
      for (pid in names(collected)) collect_one(pid, collected[[pid]], active_index[[pid]])
      if (completed - last_checkpoint >= checkpoint_tasks) {
        write_outputs(); last_checkpoint <- completed
        pilot_write_status(status_path, total, completed, results, started, length(active), cap_reached)
        cat(sprintf("%s completed=%d/%d cap=%s\n", format(Sys.time(), tz = "Europe/Madrid"),
                    completed, total, cap_reached), file = log_path, append = TRUE)
      }
      report_progress(force = TRUE)
      if (next_index > nrow(pending) && !length(active)) break
    }
  }
  cat("\n")
  write_outputs()
  pilot_write_status(status_path, total, completed, results, started, length(active), cap_reached)
  cat(sprintf("%s finished completed=%d/%d cap=%s errors=%d nonconforming=%d\n",
              format(Sys.time(), tz = "Europe/Madrid"), completed, total, cap_reached,
              sum(results$status == "error"), sum(results$status == "nonconforming")),
      file = log_path, append = TRUE)
  if (cap_reached) {
    message("Wall-time cap reached: completed rows were checkpointed and the same command resumes safely.")
  }
  invisible(list(results = results, comparisons = comparisons, cap_reached = cap_reached))
}

if (sys.nframe() == 0L) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- pilot_option(
    "output-dir",
    file.path("benchmarks", sprintf("aos_vmf_hvmf_integral_precision_%s", timestamp))
  )
  run_aos_vmf_hvmf_integral_precision_pilot(
    output_dir = output_dir,
    main_M = pilot_integer("main-M", "120", 0L),
    refine_M = pilot_integer("refine-M", "20", 0L),
    ultra_M = pilot_integer("ultra-M", "8", 0L),
    mc_M = pilot_integer("mc-M", "12", 0L),
    kernel_M = pilot_integer("kernel-M", "8", 0L),
    B = pilot_integer("B", "5000", 2L),
    cores = pilot_integer("cores", "10", 1L),
    base_seed = pilot_integer("seed", "20260801", 1L),
    cvm_block_size = pilot_integer("cvm-block-size", "50", 1L),
    max_wall_minutes = pilot_integer("max-wall-minutes", "230", 0L),
    checkpoint_tasks = pilot_integer("checkpoint-tasks", "10", 1L)
  )
}
