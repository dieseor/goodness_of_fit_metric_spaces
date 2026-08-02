#!/usr/bin/env Rscript

# Paired null-calibration pilot for the paper vMF/HvMF scenarios.  For every
# data set it fits the null once, supplies exactly the same exponential
# multipliers to both derivative methods, and compares quadrature with the
# historical score-MC derivative.  It is a validation runner only: it never
# writes to a paper-results directory.

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

source("scripts/run_power_mixtures_paper.R")

option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  values <- commandArgs(trailingOnly = TRUE)
  values <- values[startsWith(values, prefix)]
  if (!length(values)) default else substring(values[[length(values)]], nchar(prefix) + 1L)
}

integer_csv <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(strsplit(option(name, default), ",", fixed = TRUE)[[1L]]))
  if (!length(value) || any(!is.finite(value)) || any(value < minimum)) {
    stop(sprintf("`--%s` must contain integers at least %d.", name, minimum), call. = FALSE)
  }
  sort(unique(value))
}

character_csv <- function(name, default) {
  value <- trimws(strsplit(option(name, default), ",", fixed = TRUE)[[1L]])
  value[nzchar(value)]
}

M <- as.integer(option("M", "800"))
B <- as.integer(option("B", "4999"))
cores <- as.integer(option("cores", "10"))
base_seed <- as.integer(option("seed", "20260830"))
derivative_mc_size <- as.integer(option("derivative-mc-size", "1000"))
checkpoint_every <- as.integer(option("checkpoint-every", "20"))
n_values <- integer_csv("n-values", "50,100,200")
scenarios <- character_csv("scenarios", "vmf_s2_antipodal,hvmf_h2_case1")
output_dir <- option(
  "output-dir",
  file.path("benchmarks", "paper_vmf_hvmf_null_derivative_pilot")
)

if (!is.finite(M) || M < 2L || !is.finite(B) || B < 2L ||
    !is.finite(cores) || cores < 1L || !is.finite(base_seed) ||
    !is.finite(derivative_mc_size) || derivative_mc_size < 2L ||
    !is.finite(checkpoint_every) || checkpoint_every < 1L) {
  stop("Invalid scalar pilot arguments.", call. = FALSE)
}

catalog <- paper_scenario_catalog()
allowed_scenarios <- c("vmf_s2_antipodal", "hvmf_h2_case1", "hvmf_h2_case2")
if (!length(scenarios) || !all(scenarios %in% allowed_scenarios)) {
  stop(
    sprintf(
      "`--scenarios` must be drawn from: %s.",
      paste(allowed_scenarios, collapse = ", ")
    ),
    call. = FALSE
  )
}
if (all(c("hvmf_h2_case1", "hvmf_h2_case2") %in% scenarios)) {
  warning(
    paste(
      "Both HvMF cases have the same null law at beta = 0. Pooling both is",
      "valid but doubles the cost; one case is sufficient for this pilot."
    ),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
raw_path <- file.path(output_dir, "raw_results.csv")
summary_path <- file.path(output_dir, "rejection_summary.csv")
paired_path <- file.path(output_dir, "paired_method_comparison.csv")
timing_path <- file.path(output_dir, "timing_summary.csv")
manifest_path <- file.path(output_dir, "manifest.txt")

manifest_lines <- c(
  sprintf("created_at: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z")),
  "purpose: paired paper-null calibration of quadrature versus score_mc",
  sprintf("scenarios: %s", paste(scenarios, collapse = ",")),
  sprintf("n_values: %s", paste(n_values, collapse = ",")),
  sprintf("M: %d", M),
  sprintf("B: %d", B),
  sprintf("outer_cores: %d", cores),
  sprintf("base_seed: %d", base_seed),
  sprintf("score_mc_derivative_mc_size: %d", derivative_mc_size),
  sprintf("checkpoint_every_tasks: %d", checkpoint_every),
  "data_law: beta = 0 exactly",
  "same_mle_within_pair: TRUE",
  "same_multiplier_matrix_within_pair: TRUE",
  "quadrature_fast_backend: cpp",
  "quadrature_fast_kernel: contiguous_double",
  "profile_lookup_table_used: FALSE"
)
if (file.exists(manifest_path)) {
  existing_manifest <- readLines(manifest_path, warn = FALSE)
  config_prefixes <- c("scenarios:", "n_values:", "M:", "B:", "base_seed:",
                       "score_mc_derivative_mc_size:", "same_mle_within_pair:",
                       "same_multiplier_matrix_within_pair:", "quadrature_fast_kernel:")
  for (prefix in config_prefixes) {
    old <- existing_manifest[startsWith(existing_manifest, prefix)]
    new <- manifest_lines[startsWith(manifest_lines, prefix)]
    if (length(old) != 1L || length(new) != 1L || !identical(old, new)) {
      stop(
        "Existing pilot directory has incompatible configuration; choose a new `--output-dir`.",
        call. = FALSE
      )
    }
  }
} else {
  writeLines(manifest_lines, manifest_path)
}

empty_results <- data.frame(
  scenario = character(), model = character(), n = integer(), rep = integer(),
  method = character(), status = character(), error_message = character(),
  data_seed = integer(), bootstrap_seed = integer(), derivative_seed = integer(),
  ks_pvalue = numeric(), cvm_pvalue = numeric(), ks_reject = logical(), cvm_reject = logical(),
  theta_signature = character(), multiplier_sum = numeric(), multiplier_sumsq = numeric(),
  derivative_method_effective = character(), fast_multiplier_backend_effective = character(),
  fast_multiplier_cpp_kernel_effective = character(),
  elapsed_seconds = numeric(), stringsAsFactors = FALSE
)

results <- if (file.exists(raw_path)) {
  existing <- utils::read.csv(raw_path, stringsAsFactors = FALSE)
  missing <- setdiff(names(empty_results), names(existing))
  if (length(missing)) {
    stop(
      sprintf("Existing pilot results lack columns: %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  existing[, names(empty_results), drop = FALSE]
} else {
  empty_results
}
if (nrow(results) && any(results$status != "ok")) {
  warning(
    "Discarding incomplete failed rows before resume; they will be recomputed.",
    call. = FALSE
  )
  results <- results[results$status == "ok", , drop = FALSE]
}

seed_for <- function(scenario_index, n_index, rep_id, stream) {
  modulus <- 2147483647
  value <- (as.numeric(base_seed) + 1000003 * scenario_index +
    10007 * n_index + 1009 * rep_id + 10000019 * stream) %% modulus
  as.integer(value) + 1L
}

make_paired_multipliers <- function(raw_matrix) {
  cursor <- 0L
  list(
    name = "paired_Exp(1)",
    mean = 1,
    sd = 1,
    generator = function(n) {
      cursor <<- cursor + 1L
      if (cursor > nrow(raw_matrix) || n != ncol(raw_matrix)) {
        stop("Paired multiplier generator was called inconsistently.")
      }
      raw_matrix[cursor, ]
    }
  )
}

theta_signature <- function(theta) {
  values <- unlist(theta, recursive = TRUE, use.names = FALSE)
  paste(formatC(as.numeric(values), digits = 16L, format = "fg", flag = "#"), collapse = ";")
}

draw_null_data <- function(scenario, n, data_seed) {
  spec <- catalog[[scenario]]
  set.seed(data_seed)
  if (identical(spec$model, "vmf")) {
    return(generate_vmf_antipodal_sample(
      n = n, beta = 0, mu = spec$mu0, kappa = spec$kappa
    ))
  }
  generate_hvmf_halfway_sample(
    n = n, beta = 0, mu0 = spec$mu0, mu1 = spec$mu1, kappa = spec$kappa
  )
}

run_paired_task <- function(task) {
  scenario <- as.character(task$scenario)
  scenario_spec <- catalog[[scenario]]
  data_seed <- seed_for(task$scenario_index, task$n_index, task$rep, 0L)
  bootstrap_seed <- seed_for(task$scenario_index, task$n_index, task$rep, 1L)
  derivative_seed <- seed_for(task$scenario_index, task$n_index, task$rep, 2L)
  x <- draw_null_data(scenario, as.integer(task$n), data_seed)
  model <- scenario_spec$model
  common_control <- list(
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fast_multiplier_fuse_ks_cvm = TRUE,
    fast_multiplier_cache_corrections = "auto",
    fast_multiplier_cvm_block_size = 50L,
    fast_bootstrap_chunk_size = 100L,
    vmf_profile_method = "tabulated",
    vmf_profile_n_u = 4097L,
    hvmf_profile_method = "tabulated",
    hvmf_profile_n_y = 4097L
  )
  model_spec <- if (identical(model, "vmf")) {
    make_vmf_spec(distance_type = "geodesic", unknown_param = "xi")
  } else {
    make_hvmf_spec(unknown_param = "both")
  }
  x_normalized <- spec_normalize_data(model_spec, x, common_control)
  theta_hat <- model_spec$fit_theta(
    data = x_normalized, weights = NULL, null = list(type = "composite"),
    control = common_control
  )
  theta_text <- theta_signature(theta_hat)
  set.seed(bootstrap_seed)
  raw_multipliers <- matrix(stats::rexp(B * nrow(x_normalized)), nrow = B)
  multiplier_sum <- sum(raw_multipliers)
  multiplier_sumsq <- sum(raw_multipliers^2)

  run_method <- function(method) {
    control <- common_control
    control$derivative_method <- method
    if (identical(method, "score_mc")) {
      control$derivative_mc_size <- derivative_mc_size
      control$derivative_mc_seed <- derivative_seed
    }
    started <- proc.time()[["elapsed"]]
    tryCatch({
      fit <- multiplier_bootstrap_gof(
        data = x_normalized,
        spec = model_spec,
        null = list(type = "composite"),
        statistics = c("ks", "cvm"),
        ks_grid = make_sample_unique_distance_ks_grid(),
        B = B,
        alpha = 0.05,
        multipliers = make_paired_multipliers(raw_multipliers),
        n_cores = 1L,
        observed_theta_hat = theta_hat,
        bootstrap_method = "fast_multiplier",
        keep = list(
          observed_process = FALSE,
          bootstrap_statistics = FALSE,
          bootstrap_thetas = FALSE
        ),
        control = control
      )
      data.frame(
        scenario = scenario, model = model, n = as.integer(task$n), rep = as.integer(task$rep),
        method = method, status = "ok", error_message = NA_character_,
        data_seed = data_seed, bootstrap_seed = bootstrap_seed, derivative_seed = derivative_seed,
        ks_pvalue = fit$inference$ks$p_value, cvm_pvalue = fit$inference$cvm$p_value,
        ks_reject = fit$inference$ks$reject, cvm_reject = fit$inference$cvm$reject,
        theta_signature = theta_text, multiplier_sum = multiplier_sum, multiplier_sumsq = multiplier_sumsq,
        derivative_method_effective = fit$diagnostics$derivative_method_effective %||%
          fit$diagnostics$derivative_method %||% NA_character_,
        fast_multiplier_backend_effective = fit$diagnostics$fast_multiplier_backend_effective %||% NA_character_,
        fast_multiplier_cpp_kernel_effective = fit$diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_,
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        stringsAsFactors = FALSE
      )
    }, error = function(error) {
      data.frame(
        scenario = scenario, model = model, n = as.integer(task$n), rep = as.integer(task$rep),
        method = method, status = "error", error_message = conditionMessage(error),
        data_seed = data_seed, bootstrap_seed = bootstrap_seed, derivative_seed = derivative_seed,
        ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA, cvm_reject = NA,
        theta_signature = theta_text, multiplier_sum = multiplier_sum, multiplier_sumsq = multiplier_sumsq,
        derivative_method_effective = NA_character_, fast_multiplier_backend_effective = NA_character_,
        fast_multiplier_cpp_kernel_effective = NA_character_,
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        stringsAsFactors = FALSE
      )
    })
  }
  rbind(run_method("quadrature"), run_method("score_mc"))
}

task_grid <- expand.grid(
  scenario = scenarios,
  n = n_values,
  rep = seq_len(M),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_grid$scenario_index <- match(task_grid$scenario, scenarios)
task_grid$n_index <- match(task_grid$n, n_values)
ensure_distance_profile_cpp_loaded()
task_key <- function(data) paste(data$scenario, data$n, data$rep, sep = "|")
existing_ok <- results[results$status == "ok", , drop = FALSE]
completed_methods <- table(task_key(existing_ok))
done_keys <- names(completed_methods)[completed_methods == 2L]
pending <- task_grid[!task_key(task_grid) %in% done_keys, , drop = FALSE]

write_outputs <- function() {
  results <<- results[order(results$scenario, results$n, results$rep, results$method), , drop = FALSE]
  utils::write.csv(results, raw_path, row.names = FALSE)
  valid <- results[results$status == "ok", , drop = FALSE]
  statistics <- c("ks", "cvm")
  summaries <- list()
  paired <- list()
  row <- 1L
  pair_row <- 1L
  groups <- split(valid, interaction(valid$scenario, valid$n, drop = TRUE))
  for (group in groups) {
    for (statistic in statistics) {
      reject_name <- paste0(statistic, "_reject")
      pvalue_name <- paste0(statistic, "_pvalue")
      for (method in c("quadrature", "score_mc")) {
        method_rows <- group[group$method == method, , drop = FALSE]
        count <- nrow(method_rows)
        rejected <- sum(method_rows[[reject_name]])
        ci <- stats::binom.test(rejected, count, conf.level = 0.95)$conf.int
        summaries[[row]] <- data.frame(
          scenario = method_rows$scenario[[1L]], model = method_rows$model[[1L]],
          n = method_rows$n[[1L]], statistic = statistic, method = method,
          replications = count, rejections = rejected, rejection_rate = rejected / count,
          rejection_ci_lower = ci[[1L]], rejection_ci_upper = ci[[2L]],
          mean_pvalue = mean(method_rows[[pvalue_name]]),
          median_pvalue = stats::median(method_rows[[pvalue_name]]),
          stringsAsFactors = FALSE
        )
        row <- row + 1L
      }
      wide <- reshape(
        group[, c("rep", "method", reject_name, pvalue_name), drop = FALSE],
        idvar = "rep", timevar = "method", direction = "wide"
      )
      wide <- wide[stats::complete.cases(wide), , drop = FALSE]
      reject_difference <- wide[[paste0(reject_name, ".quadrature")]] -
        wide[[paste0(reject_name, ".score_mc")]]
      pvalue_difference <- wide[[paste0(pvalue_name, ".quadrature")]] -
        wide[[paste0(pvalue_name, ".score_mc")]]
      difference_se <- if (nrow(wide) > 1L) stats::sd(reject_difference) / sqrt(nrow(wide)) else NA_real_
      paired[[pair_row]] <- data.frame(
        scenario = group$scenario[[1L]], model = group$model[[1L]], n = group$n[[1L]],
        statistic = statistic, paired_replications = nrow(wide),
        quadrature_minus_score_mc_rejection_rate = mean(reject_difference),
        paired_difference_ci_lower = mean(reject_difference) - 1.96 * difference_se,
        paired_difference_ci_upper = mean(reject_difference) + 1.96 * difference_se,
        quadrature_minus_score_mc_mean_pvalue = mean(pvalue_difference),
        quadrature_minus_score_mc_median_pvalue = stats::median(pvalue_difference),
        discordant_quadrature_only = sum(reject_difference == 1L),
        discordant_score_mc_only = sum(reject_difference == -1L),
        stringsAsFactors = FALSE
      )
      pair_row <- pair_row + 1L
    }
  }
  utils::write.csv(do.call(rbind, summaries), summary_path, row.names = FALSE)
  utils::write.csv(do.call(rbind, paired), paired_path, row.names = FALSE)
  timing <- stats::aggregate(
    elapsed_seconds ~ scenario + model + n + method,
    data = valid,
    FUN = function(value) c(mean = mean(value), median = stats::median(value))
  )
  timing <- do.call(data.frame, timing)
  names(timing)[names(timing) == "elapsed_seconds.mean"] <- "mean_elapsed_seconds"
  names(timing)[names(timing) == "elapsed_seconds.median"] <- "median_elapsed_seconds"
  utils::write.csv(timing, timing_path, row.names = FALSE)
}

total <- nrow(task_grid)
completed <- total - nrow(pending)
last_checkpoint <- completed
active_jobs <- list()
active_index <- list()
next_index <- 1L
max_active <- min(cores, nrow(pending))
started <- Sys.time()
report_progress <- function() {
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  rate <- if (elapsed > 0) completed / elapsed else 0
  remaining <- if (rate > 0) (total - completed) / rate else NA_real_
  cat(sprintf(
    "\r[null derivative pilot] %d/%d complete | %d active | elapsed %.1f min | ETA %s",
    completed, total, length(active_jobs), elapsed / 60,
    if (is.finite(remaining)) sprintf("%.1f min", remaining / 60) else "calculating"
  ))
  flush.console()
}

if (!nrow(pending)) {
  cat("All requested paired pilot tasks are already complete; refreshing summaries.\n")
} else if (.Platform$OS.type != "unix" || cores <= 1L) {
  for (index in seq_len(nrow(pending))) {
    results <- rbind(results, run_paired_task(pending[index, , drop = FALSE]))
    completed <- completed + 1L
    if (completed - last_checkpoint >= checkpoint_every) {
      write_outputs()
      last_checkpoint <- completed
    }
    report_progress()
  }
  cat("\n")
} else {
  while (completed < total) {
    while (next_index <= nrow(pending) && length(active_jobs) < max_active) {
      index <- next_index
      job <- parallel::mcparallel(
        run_paired_task(pending[index, , drop = FALSE]),
        mc.set.seed = FALSE
      )
      pid <- as.character(job$pid)
      active_jobs[[pid]] <- job
      active_index[[pid]] <- index
      next_index <- next_index + 1L
    }
    collected <- parallel::mccollect(active_jobs, wait = FALSE)
    if (is.null(collected)) {
      report_progress()
      Sys.sleep(0.5)
      next
    }
    for (pid in names(collected)) {
      index <- active_index[[pid]]
      value <- collected[[pid]]
      if (inherits(value, "try-error")) {
        value <- data.frame(
          scenario = pending$scenario[[index]], model = catalog[[pending$scenario[[index]]]]$model,
          n = pending$n[[index]], rep = pending$rep[[index]],
          method = c("quadrature", "score_mc"), status = "error",
          error_message = as.character(value), data_seed = NA_integer_, bootstrap_seed = NA_integer_,
          derivative_seed = NA_integer_, ks_pvalue = NA_real_, cvm_pvalue = NA_real_,
          ks_reject = NA, cvm_reject = NA, theta_signature = NA_character_,
          multiplier_sum = NA_real_, multiplier_sumsq = NA_real_,
          derivative_method_effective = NA_character_, fast_multiplier_backend_effective = NA_character_,
          fast_multiplier_cpp_kernel_effective = NA_character_, elapsed_seconds = NA_real_,
          stringsAsFactors = FALSE
        )
        cat(sprintf("\nERROR in %s n=%d rep=%d: %s\n", pending$scenario[[index]], pending$n[[index]], pending$rep[[index]], as.character(collected[[pid]])))
      }
      results <- rbind(results, value)
      active_jobs[[pid]] <- NULL
      active_index[[pid]] <- NULL
      completed <- completed + 1L
    }
    if (completed - last_checkpoint >= checkpoint_every) {
      write_outputs()
      last_checkpoint <- completed
    }
    report_progress()
  }
  cat("\n")
}

write_outputs()
failed <- results[results$status != "ok", , drop = FALSE]
if (nrow(failed)) {
  stop(sprintf("Pilot finished with %d failed method runs; inspect %s.", nrow(failed), raw_path), call. = FALSE)
}

cat(sprintf("Wrote %s\n", raw_path))
cat(sprintf("Wrote %s\n", summary_path))
cat(sprintf("Wrote %s\n", paired_path))
cat(sprintf("Wrote %s\n", timing_path))
