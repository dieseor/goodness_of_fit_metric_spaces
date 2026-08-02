#!/usr/bin/env Rscript

# Paired end-to-end validation of the production contiguous-double C++
# fast-multiplier kernel against the explicit historical `legacy` fallback.

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
source("bootstrap/multiplier_bootstrap.R")

option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  values <- commandArgs(trailingOnly = TRUE)
  match <- values[startsWith(values, prefix)]
  if (!length(match)) default else substring(match[[length(match)]], nchar(prefix) + 1L)
}

integer_vector <- function(name, default) {
  values <- suppressWarnings(as.integer(strsplit(option(name, default), ",", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!is.finite(values)) || any(values < 2L)) {
    stop(sprintf("`--%s` must be a nonempty vector of integers at least two.", name))
  }
  unique(values)
}

models <- strsplit(tolower(option("models", "vmf,hvmf")), ",", fixed = TRUE)[[1L]]
if (!all(models %in% c("vmf", "hvmf"))) stop("`--models` must use vMF and/or HvMF.")
dimensions <- integer_vector("dimensions", "2,10")
sizes <- integer_vector("sizes", "50,100,200,400")
B <- as.integer(option("B", "999"))
outer_cores <- as.integer(option("cores", "1"))
seed <- as.integer(option("seed", "20260821"))
derivative_method <- tolower(option("derivative-method", "quadrature"))
output_dir <- option(
  "output-dir",
  file.path("benchmarks", "fast_multiplier_cpp_kernel_end_to_end")
)
if (!is.finite(B) || B < 2L || !is.finite(outer_cores) || outer_cores < 1L ||
    !is.finite(seed) || !derivative_method %in% c("quadrature", "score_mc")) {
  stop("Invalid end-to-end validation arguments.")
}

draw_data <- function(model, q, n, seed) {
  set.seed(seed)
  kappa <- max(2, q)
  if (identical(model, "vmf")) {
    return(normalize_vmf_data(
      rotasym::r_vMF(n, mu = c(1, rep.int(0, q)), kappa = kappa)
    ))
  }
  mu <- c(cosh(0.25), sinh(0.25), rep.int(0, q - 1L))
  rhvmf_polar(n, mu = mu, kappa = kappa)
}

run_bootstrap <- function(model, q, x, B, bootstrap_seed, kernel,
                          derivative_method) {
  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    seed = bootstrap_seed,
    n_cores = 1L,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = derivative_method,
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = kernel,
      fast_bootstrap_chunk_size = 100L,
      derivative_mc_size = if (identical(derivative_method, "score_mc")) 5000L else NULL,
      derivative_mc_seed = if (identical(derivative_method, "score_mc")) {
        bootstrap_seed + 1L
      } else {
        NULL
      }
    )
  )
  if (identical(model, "vmf")) {
    return(do.call(multiplier_bootstrap_vmf, c(
      common,
      list(distance_type = "geodesic", unknown_param = "xi")
    )))
  }
  do.call(multiplier_bootstrap_hvmf, c(
    common,
    list(unknown_param = "both", fast_multiplier_backend = "cpp")
  ))
}

max_relative_difference <- function(first, second) {
  max(abs(first - second)) / max(1, max(abs(first)))
}

run_case <- function(item) {
  data_seed <- seed + 100000L * item$q + 1000L * item$n + item$case_index
  bootstrap_seed <- data_seed + 500L
  x <- draw_data(item$model, item$q, item$n, data_seed)
  legacy <- run_bootstrap(
    item$model, item$q, x, B, bootstrap_seed, "legacy", derivative_method
  )
  candidate <- run_bootstrap(
    item$model, item$q, x, B, bootstrap_seed, "contiguous_double",
    derivative_method
  )
  ks_abs <- max(abs(legacy$bootstrap$statistics$ks - candidate$bootstrap$statistics$ks))
  cvm_abs <- max(abs(legacy$bootstrap$statistics$cvm - candidate$bootstrap$statistics$cvm))
  row <- data.frame(
    model = item$model,
    q = item$q,
    n = item$n,
    B = B,
    derivative_method = derivative_method,
    data_seed = data_seed,
    bootstrap_seed = bootstrap_seed,
    max_abs_ks_difference = ks_abs,
    max_rel_ks_difference = max_relative_difference(
      legacy$bootstrap$statistics$ks, candidate$bootstrap$statistics$ks
    ),
    max_abs_cvm_difference = cvm_abs,
    max_rel_cvm_difference = max_relative_difference(
      legacy$bootstrap$statistics$cvm, candidate$bootstrap$statistics$cvm
    ),
    ks_critical_difference = abs(
      legacy$inference$ks$critical_value - candidate$inference$ks$critical_value
    ),
    cvm_critical_difference = abs(
      legacy$inference$cvm$critical_value - candidate$inference$cvm$critical_value
    ),
    ks_p_value_difference = abs(
      legacy$inference$ks$p_value - candidate$inference$ks$p_value
    ),
    cvm_p_value_difference = abs(
      legacy$inference$cvm$p_value - candidate$inference$cvm$p_value
    ),
    legacy_common_observed_seconds = legacy$diagnostics$common_observed_seconds,
    candidate_common_observed_seconds = candidate$diagnostics$common_observed_seconds,
    legacy_fast_prep_seconds = legacy$diagnostics$fast_prep_seconds,
    candidate_fast_prep_seconds = candidate$diagnostics$fast_prep_seconds,
    legacy_fast_loop_seconds = legacy$diagnostics$fast_loop_seconds,
    candidate_fast_loop_seconds = candidate$diagnostics$fast_loop_seconds,
    legacy_fast_total_seconds = legacy$diagnostics$fast_total_seconds,
    candidate_fast_total_seconds = candidate$diagnostics$fast_total_seconds,
    loop_speedup = if (candidate$diagnostics$fast_loop_seconds > 0) {
      legacy$diagnostics$fast_loop_seconds /
        candidate$diagnostics$fast_loop_seconds
    } else {
      NA_real_
    },
    legacy_kernel_effective = legacy$diagnostics$fast_multiplier_cpp_kernel_effective,
    candidate_kernel_effective = candidate$diagnostics$fast_multiplier_cpp_kernel_effective,
    stringsAsFactors = FALSE
  )
  if (!identical(row$legacy_kernel_effective, "legacy") ||
      !identical(row$candidate_kernel_effective, "contiguous_double")) {
    stop("The paired runs did not select the intended C++ kernels.")
  }
  row
}

case_grid <- expand.grid(
  model = models,
  q = dimensions,
  n = sizes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
case_grid$case_index <- seq_len(nrow(case_grid))
ensure_distance_profile_cpp_loaded()
items <- lapply(seq_len(nrow(case_grid)), function(index) case_grid[index, ])
run_cases_with_progress <- function(items, n_cores) {
  total <- length(items)
  rows <- vector("list", total)
  completed <- 0L

  report_progress <- function() {
    cat(sprintf(
      "\r[end-to-end kernel validation] %d/%d complete (%d active)",
      completed, total, length(active_jobs)
    ))
    flush.console()
  }

  if (.Platform$OS.type != "unix" || n_cores <= 1L) {
    active_jobs <- integer()
    for (index in seq_len(total)) {
      rows[[index]] <- try(run_case(items[[index]]), silent = TRUE)
      completed <- completed + 1L
      report_progress()
    }
    cat("\n")
    return(rows)
  }

  active_jobs <- list()
  active_index <- list()
  next_index <- 1L
  max_active <- min(as.integer(n_cores), total)
  while (completed < total) {
    while (next_index <= total && length(active_jobs) < max_active) {
      index <- next_index
      job <- parallel::mcparallel(
        try(run_case(items[[index]]), silent = TRUE),
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
      Sys.sleep(0.2)
      next
    }
    for (pid in names(collected)) {
      index <- active_index[[pid]]
      rows[[index]] <- collected[[pid]]
      if (inherits(rows[[index]], "try-error")) {
        cat(sprintf(
          "\nERROR in %s q=%s n=%s: %s\n",
          as.character(items[[index]]$model),
          as.character(items[[index]]$q),
          as.character(items[[index]]$n),
          as.character(rows[[index]])
        ))
      }
      active_jobs[[pid]] <- NULL
      active_index[[pid]] <- NULL
      completed <- completed + 1L
    }
    report_progress()
  }
  cat("\n")
  rows
}

rows <- run_cases_with_progress(items, outer_cores)
if (any(vapply(rows, inherits, logical(1), what = "try-error"))) {
  stop("At least one paired end-to-end kernel validation failed.")
}
summary <- do.call(rbind, rows)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
print(summary, row.names = FALSE)
cat(sprintf("Wrote %s\n", file.path(output_dir, "summary.csv")))

if (any(summary$max_rel_ks_difference > 1e-12) ||
    any(summary$max_rel_cvm_difference > 1e-12) ||
    any(summary$ks_critical_difference > 1e-12) ||
    any(summary$cvm_critical_difference > 1e-12) ||
    any(summary$ks_p_value_difference > 1e-12) ||
    any(summary$cvm_p_value_difference > 1e-12)) {
  stop("The candidate kernel exceeded the paired equivalence tolerance.")
}
