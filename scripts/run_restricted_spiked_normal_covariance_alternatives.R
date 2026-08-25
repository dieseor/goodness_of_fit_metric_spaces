#!/usr/bin/env Rscript

# Resumable paper-design runner for the restricted mean-aligned single-spiked
# Gaussian model. This file owns only this model's double-spike experiment.
#
# P0 = N_d(theta0, I_d + 2 u0 u0^T), u0 = theta0 / ||theta0||,
# Q  = N_d(theta0, I_d + 2 u0 u0^T + 2 q q^T), q perpendicular to u0,
# P_beta = (1-beta) P0 + beta Q.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

restricted_spiked_parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    fields <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[fields[[1L]]]] <- if (length(fields) == 1L) {
      "TRUE"
    } else {
      paste(fields[-1L], collapse = "=")
    }
  }
  out
}

restricted_spiked_parse_csv <- function(value, default, integer = FALSE) {
  if (is.null(value) || !nzchar(value)) return(default)
  out <- as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!is.finite(out))) stop("Invalid comma-separated argument.")
  if (isTRUE(integer)) {
    if (any(out != as.integer(out))) stop("Expected integer values.")
    out <- as.integer(out)
  }
  out
}

restricted_spiked_mean_catalog <- function() {
  list(
    axis_075 = list(
      label = "0.75 e1", direction = "axis", norm = 0.75, config_id = 1L,
      base_seed = 20260831L
    ),
    axis_100 = list(
      label = "e1", direction = "axis", norm = 1, config_id = 2L,
      base_seed = 20260833L
    ),
    diagonal_150 = list(
      label = "1.5 1/sqrt(d)", direction = "diagonal", norm = 1.5,
      config_id = 3L, base_seed = 20260831L
    ),
    diagonal_100 = list(
      label = "1/sqrt(d)", direction = "diagonal", norm = 1,
      config_id = 4L, base_seed = 20260833L
    )
  )
}

restricted_spiked_default_seed <- function(mean_config) {
  specification <- restricted_spiked_mean_catalog()[[mean_config]]
  if (is.null(specification)) stop("Unknown restricted-spiked mean configuration.")
  specification$base_seed
}

restricted_spiked_mean_vector <- function(d, mean_config) {
  specification <- restricted_spiked_mean_catalog()[[mean_config]]
  if (is.null(specification)) stop("Unknown restricted-spiked mean configuration.")
  if (identical(specification$direction, "axis")) {
    return(c(specification$norm, rep(0, d - 1L)))
  }
  rep(specification$norm / sqrt(d), d)
}

restricted_spiked_second_direction <- function(d, mean_config) {
  specification <- restricted_spiked_mean_catalog()[[mean_config]]
  if (identical(specification$direction, "axis")) {
    return(c(0, 1, rep(0, d - 2L)))
  }
  c(1 / sqrt(2), -1 / sqrt(2), rep(0, d - 2L))
}

restricted_spiked_design_key <- function(x, include_replication = FALSE) {
  pieces <- list(
    as.character(x$mean_config), as.integer(x$d), as.integer(x$n),
    formatC(as.numeric(x$beta), digits = 16L, format = "fg", flag = "#")
  )
  if (isTRUE(include_replication)) {
    pieces <- c(pieces, list(as.integer(x$replication)))
  }
  do.call(paste, c(pieces, sep = "|"))
}

make_restricted_spiked_design <- function(mean_config,
                                           dimensions = c(2L, 5L),
                                           n_values = c(50L, 100L, 200L, 400L),
                                           beta_values = c(0, 0.25, 0.5, 1),
                                           M = 1000L) {
  catalog <- restricted_spiked_mean_catalog()
  if (length(mean_config) != 1L || !mean_config %in% names(catalog)) {
    stop(sprintf(
      "`mean_config` must be one of: %s.", paste(names(catalog), collapse = ", ")
    ))
  }
  combinations <- expand.grid(
    mean_config = mean_config,
    d = sort(unique(as.integer(dimensions))),
    n = sort(unique(as.integer(n_values))),
    beta = sort(unique(as.numeric(beta_values))),
    stringsAsFactors = FALSE
  )
  combinations$mean_config_id <- catalog[[mean_config]]$config_id
  combinations$design_id <- seq_len(nrow(combinations)) +
    100L * (catalog[[mean_config]]$config_id - 1L)
  merge(combinations, data.frame(replication = seq_len(as.integer(M))), by = NULL)
}

restricted_spiked_seed <- function(base_seed, design_id, replication, stream = 0L) {
  modulus <- 2147483647
  value <- (as.numeric(base_seed) + 1000003 * as.numeric(design_id) +
    1009 * as.numeric(replication) + 10000019 * as.numeric(stream)) %% modulus
  as.integer(value) + 1L
}

restricted_spiked_manifest <- function(design, M, B, base_seed, lambda,
                                        derivative_mc_size, cvm_block_size) {
  combinations <- unique(design[c(
    "mean_config", "mean_config_id", "d", "n", "beta", "design_id"
  )])
  transform(
    combinations,
    M = as.integer(M), B = as.integer(B), base_seed = as.integer(base_seed),
    lambda = as.numeric(lambda), derivative_method = "score_mc",
    derivative_mc_size = as.integer(derivative_mc_size),
    vhat_method = "restricted_spiked_analytic_fisher",
    ks_grid = "sample_points_unique_distances",
    cvm_grid = "sample_points_unique_distances",
    ks_cvm_sample_grid_shared = TRUE, statistics = "ks,cvm",
    bootstrap_method = "fast_multiplier", fast_backend = "cpp",
    fast_kernel = "contiguous_double", fast_fused = TRUE,
    cvm_block_size = as.integer(cvm_block_size), fallback_allowed = FALSE
  )
}

restricted_spiked_write_atomic_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) {
    unlink(temporary, force = TRUE)
    stop(sprintf("Could not atomically update '%s'.", path))
  }
}

restricted_spiked_lock_path <- function(output_dir) {
  file.path(output_dir, ".restricted_spiked.lock")
}

restricted_spiked_acquire_lock <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock_path <- restricted_spiked_lock_path(output_dir)
  if (!dir.create(lock_path, showWarnings = FALSE)) {
    owner_path <- file.path(lock_path, "owner.txt")
    owner <- if (file.exists(owner_path)) {
      paste(readLines(owner_path, warn = FALSE), collapse = "; ")
    } else {
      "owner metadata unavailable"
    }
    stop(sprintf("Output directory is locked (%s).", owner), call. = FALSE)
  }
  token <- paste(Sys.info()[["nodename"]], Sys.getpid(), Sys.time(), sep = "|")
  writeLines(c(
    sprintf("slurm_job_id: %s", Sys.getenv("SLURM_JOB_ID", unset = "unavailable")),
    sprintf("hostname: %s", Sys.info()[["nodename"]]),
    sprintf("pid: %d", Sys.getpid()),
    sprintf("owner_token: %s", token)
  ), file.path(lock_path, "owner.txt"))
  list(path = lock_path, token = token)
}

restricted_spiked_release_lock <- function(lock) {
  if (is.null(lock$path) || !dir.exists(lock$path)) return(invisible(FALSE))
  owner <- readLines(file.path(lock$path, "owner.txt"), warn = FALSE)
  if (!sprintf("owner_token: %s", lock$token) %in% owner) {
    warning("Refusing to release a restricted-spiked lock owned by another process.")
    return(invisible(FALSE))
  }
  unlink(lock$path, recursive = TRUE, force = TRUE)
  invisible(!dir.exists(lock$path))
}

restricted_spiked_validate_manifest <- function(path, expected) {
  if (!file.exists(path)) return(invisible(TRUE))
  observed <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- names(expected)
  if (!all(required %in% names(observed))) {
    stop("Existing restricted-spiked manifest is incomplete; use a new output directory.")
  }
  row_key <- function(x) restricted_spiked_design_key(x)
  scalar_fields <- setdiff(required, c(
    "mean_config", "mean_config_id", "d", "n", "beta", "design_id"
  ))
  observed_order <- match(row_key(expected), row_key(observed))
  checks <- c(
    design = setequal(row_key(observed), row_key(expected)),
    design_id = !anyNA(observed_order) && identical(
      as.integer(observed$design_id[observed_order]), as.integer(expected$design_id)
    ),
    vapply(scalar_fields, function(field) {
      !anyNA(observed_order) && identical(
        as.character(observed[[field]][observed_order]),
        as.character(expected[[field]])
      )
    }, logical(1))
  )
  checks[is.na(checks)] <- FALSE
  if (!all(checks)) {
    stop(sprintf(
      "Existing restricted-spiked manifest is incompatible (%s); use a new output directory.",
      paste(names(checks)[!checks], collapse = ", ")
    ))
  }
  invisible(TRUE)
}

empty_restricted_spiked_results <- function() {
  data.frame(
    mean_config = character(), mean_label = character(), mean_config_id = integer(),
    d = integer(), n = integer(), beta = numeric(), design_id = integer(),
    replication = integer(), status = character(), error_message = character(),
    warning_message = character(), theta_true = character(), theta_norm_true = numeric(),
    u0 = character(), q = character(), lambda = numeric(), seed_data = integer(),
    seed_bootstrap = integer(), seed_derivative = integer(), theta_hat = character(),
    lambda_hat = numeric(), sigma_eigen_max_error = numeric(), score_mean_norm = numeric(),
    score_mean_max_abs = numeric(), ks_statistic = numeric(), cvm_statistic = numeric(),
    ks_pvalue = numeric(), cvm_pvalue = numeric(), ks_reject = logical(),
    cvm_reject = logical(), effective_bootstrap_method = character(),
    fallback_to_reestimated = logical(), fast_backend = character(),
    fast_kernel = character(), fast_fused = logical(), derivative_method = character(),
    vhat_method = character(), sample_grid_mode = character(),
    shared_sample_ks_cvm_cache = logical(), lightweight_ks_prep = logical(),
    lightweight_cvm_prep = logical(), elapsed_seconds = numeric(),
    stringsAsFactors = FALSE
  )
}

restricted_spiked_conforming <- function(results) {
  results$status == "ok" &
    results$effective_bootstrap_method == "fast_multiplier" &
    !results$fallback_to_reestimated & results$fast_backend == "cpp" &
    results$fast_kernel == "contiguous_double" & results$fast_fused &
    results$derivative_method == "score_mc" &
    results$vhat_method == "restricted_spiked_analytic_fisher" &
    results$sample_grid_mode == "sample_points_unique_distances" &
    results$shared_sample_ks_cvm_cache & results$lightweight_ks_prep &
    results$lightweight_cvm_prep
}

restricted_spiked_summarize <- function(results) {
  ok <- results[restricted_spiked_conforming(results), , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  split_key <- interaction(ok$mean_config, ok$d, ok$n, ok$beta, drop = TRUE)
  do.call(rbind, lapply(split(ok, split_key), function(group) data.frame(
    mean_config = group$mean_config[[1L]], mean_label = group$mean_label[[1L]],
    d = group$d[[1L]], n = group$n[[1L]], beta = group$beta[[1L]],
    completed = nrow(group), ks_rejection_percent = 100 * mean(group$ks_reject),
    cvm_rejection_percent = 100 * mean(group$cvm_reject),
    ks_pvalue_mean = mean(group$ks_pvalue), cvm_pvalue_mean = mean(group$cvm_pvalue),
    mean_elapsed_seconds = mean(group$elapsed_seconds), stringsAsFactors = FALSE
  )))
}

restricted_spiked_write_status <- function(path, total, completed, results, started, cores) {
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  conforming <- if (nrow(results)) restricted_spiked_conforming(results) else logical()
  successful <- results$elapsed_seconds[conforming]
  mean_job <- if (length(successful)) mean(successful) else NA_real_
  eta <- if (is.finite(mean_job)) mean_job * (total - completed) / cores else NA_real_
  writeLines(c(
    sprintf("updated: %s", format(Sys.time(), tz = "Europe/Madrid")),
    sprintf("completed: %d/%d", completed, total),
    sprintf("pending: %d", total - completed), sprintf("elapsed_seconds: %.1f", elapsed),
    sprintf("mean_seconds_per_conforming_job: %s", if (is.finite(mean_job)) mean_job else "NA"),
    sprintf("eta_seconds: %s", if (is.finite(eta)) round(eta) else "NA"),
    sprintf("ok: %d", sum(results$status == "ok")),
    sprintf("errors: %d", sum(results$status == "error")),
    sprintf("nonconforming: %d", sum(results$status == "nonconforming"))
  ), path)
}

run_restricted_spiked_covariance_alternatives <- function(
    mean_config, output_dir, M = 1000L, B = 5000L,
    dimensions = c(2L, 5L), n_values = c(50L, 100L, 200L, 400L),
    beta_values = c(0, 0.25, 0.5, 1), lambda = 2,
    derivative_mc_size = 10000L, cvm_block_size = 50L, cores = 4L,
    checkpoint_results = 2000L, base_seed = NULL, show_progress = TRUE) {
  if (is.null(base_seed)) base_seed <- restricted_spiked_default_seed(mean_config)
  numeric_settings <- c(M, B, dimensions, n_values, beta_values, lambda,
                        derivative_mc_size, cvm_block_size, cores,
                        checkpoint_results, base_seed)
  if (any(!is.finite(numeric_settings)) || M < 1L || B < 1L ||
      any(dimensions < 2L) || any(n_values < 2L) ||
      any(beta_values < 0 | beta_values > 1) || lambda <= 0 ||
      derivative_mc_size < 1L || cvm_block_size < 1L || cores < 1L ||
      checkpoint_results < 1L) {
    stop("Invalid restricted-spiked experiment settings.")
  }
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("Outer restricted-spiked parallelism requires a Unix platform.")
  }

  source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"),
         local = environment())
  design <- make_restricted_spiked_design(
    mean_config, dimensions, n_values, beta_values, M
  )
  expected_manifest <- restricted_spiked_manifest(
    design, M, B, base_seed, lambda, derivative_mc_size, cvm_block_size
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- restricted_spiked_acquire_lock(output_dir)
  on.exit(restricted_spiked_release_lock(lock), add = TRUE)

  manifest_path <- file.path(output_dir, "manifest.csv")
  results_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  log_path <- file.path(output_dir, "run.log")
  restricted_spiked_validate_manifest(manifest_path, expected_manifest)
  if (!file.exists(manifest_path)) {
    restricted_spiked_write_atomic_csv(expected_manifest, manifest_path)
  }

  existing <- if (file.exists(results_path)) {
    utils::read.csv(results_path, stringsAsFactors = FALSE)
  } else {
    empty_restricted_spiked_results()
  }
  required <- names(empty_restricted_spiked_results())
  if (!all(required %in% names(existing))) {
    stop("Existing restricted-spiked results have an incompatible schema.")
  }
  if (nrow(existing) && anyDuplicated(restricted_spiked_design_key(
      existing, include_replication = TRUE))) {
    stop("Existing restricted-spiked results contain duplicate keys.")
  }
  conforming <- if (nrow(existing)) restricted_spiked_conforming(existing) else logical()
  done <- restricted_spiked_design_key(
    existing[conforming, , drop = FALSE], include_replication = TRUE
  )
  pending <- design[!restricted_spiked_design_key(
    design, include_replication = TRUE
  ) %in% done, , drop = FALSE]
  started <- Sys.time()
  cat(sprintf(
    "%s mean_config=%s M=%d B=%d cores=%d pending=%d total=%d\n",
    format(started, tz = "Europe/Madrid"), mean_config, M, B, cores,
    nrow(pending), nrow(design)
  ), file = log_path, append = TRUE)
  restricted_spiked_write_status(
    status_path, nrow(design), nrow(design) - nrow(pending), existing, started, cores
  )
  if (!nrow(pending)) return(invisible(list(
    results = existing, summary = restricted_spiked_summarize(existing)
  )))

  catalog <- restricted_spiked_mean_catalog()
  call_fast <- function(x, bootstrap_seed, derivative_seed) {
    multiplier_bootstrap_restricted_spiked_normal(
      data = x, null = list(type = "composite"), statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = 0.05,
      n_cores = 1L, seed = bootstrap_seed, bootstrap_method = "fast_multiplier",
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                  bootstrap_thetas = FALSE),
      control = list(
        derivative_method = "score_mc", derivative_mc_size = derivative_mc_size,
        derivative_mc_seed = derivative_seed,
        fast_multiplier_cvm_block_size = cvm_block_size,
        fast_multiplier_backend = "cpp",
        fast_multiplier_cpp_kernel = "contiguous_double",
        fast_multiplier_fuse_ks_cvm = TRUE,
        fast_multiplier_cache_corrections = "auto", progress_bar = FALSE
      ),
      distance_profile_backend = "r", fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double", fuse_ks_cvm = TRUE,
      cache_block_corrections = "auto"
    )
  }
  run_job <- function(row) {
    elapsed_start <- proc.time()[["elapsed"]]
    d <- as.integer(row$d)
    theta0 <- restricted_spiked_mean_vector(d, mean_config)
    u0 <- theta0 / sqrt(sum(theta0^2))
    q <- restricted_spiked_second_direction(d, mean_config)
    data_seed <- restricted_spiked_seed(base_seed, row$design_id, row$replication, 0L)
    bootstrap_seed <- restricted_spiked_seed(base_seed, row$design_id, row$replication, 1L)
    derivative_seed <- restricted_spiked_seed(base_seed, row$design_id, row$replication, 2L)
    base <- data.frame(
      mean_config = mean_config, mean_label = catalog[[mean_config]]$label,
      mean_config_id = row$mean_config_id, d = d, n = as.integer(row$n),
      beta = as.numeric(row$beta), design_id = as.integer(row$design_id),
      replication = as.integer(row$replication), status = "ok",
      error_message = NA_character_, warning_message = NA_character_,
      theta_true = paste(theta0, collapse = ";"), theta_norm_true = sqrt(sum(theta0^2)),
      u0 = paste(u0, collapse = ";"), q = paste(q, collapse = ";"), lambda = lambda,
      seed_data = data_seed, seed_bootstrap = bootstrap_seed,
      seed_derivative = derivative_seed, stringsAsFactors = FALSE
    )
    tryCatch({
      set.seed(data_seed)
      sigma0 <- diag(d) + lambda * tcrossprod(u0)
      sigma_alt <- sigma0 + lambda * tcrossprod(q)
      labels <- stats::rbinom(base$n, 1L, base$beta) == 1L
      x <- mvtnorm::rmvnorm(base$n, mean = theta0, sigma = sigma0)
      if (any(labels)) {
        x[labels, ] <- mvtnorm::rmvnorm(sum(labels), mean = theta0, sigma = sigma_alt)
      }
      warnings <- character()
      result <- withCallingHandlers(
        call_fast(x, bootstrap_seed, derivative_seed),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      fitted <- result$observed$theta_hat
      score <- restricted_spiked_normal_score_matrix(x, fitted)
      eig <- eigen(fitted$Sigma, symmetric = TRUE, only.values = TRUE)$values
      diagnostics <- result$diagnostics
      output <- cbind(base, data.frame(
        theta_hat = paste(fitted$theta, collapse = ";"), lambda_hat = fitted$lambda,
        sigma_eigen_max_error = max(abs(sort(eig) - sort(c(1 + fitted$lambda, rep(1, d - 1L))))),
        score_mean_norm = sqrt(sum(colMeans(score)^2)),
        score_mean_max_abs = max(abs(colMeans(score))),
        ks_statistic = result$inference$ks$observed,
        cvm_statistic = result$inference$cvm$observed,
        ks_pvalue = result$inference$ks$p_value,
        cvm_pvalue = result$inference$cvm$p_value,
        ks_reject = result$inference$ks$reject, cvm_reject = result$inference$cvm$reject,
        effective_bootstrap_method = diagnostics$effective_bootstrap_method %||% NA_character_,
        fallback_to_reestimated = isTRUE(diagnostics$fallback_to_reestimated),
        fast_backend = diagnostics$fast_multiplier_backend_effective %||% NA_character_,
        fast_kernel = diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_,
        fast_fused = isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective),
        derivative_method = diagnostics$derivative_method_effective %||% NA_character_,
        vhat_method = diagnostics$vhat_method %||% NA_character_,
        sample_grid_mode = result$grid$mode %||% NA_character_,
        shared_sample_ks_cvm_cache = isTRUE(diagnostics$shared_sample_ks_cvm_cache),
        lightweight_ks_prep = isTRUE(diagnostics$lightweight_ks_prep),
        lightweight_cvm_prep = isTRUE(diagnostics$lightweight_cvm_prep),
        elapsed_seconds = proc.time()[["elapsed"]] - elapsed_start,
        stringsAsFactors = FALSE
      ))
      output$warning_message <- if (length(warnings)) {
        paste(unique(warnings), collapse = " | ")
      } else {
        NA_character_
      }
      if (!restricted_spiked_conforming(output)) {
        output$status <- "nonconforming"
        output$error_message <- "Requested restricted-spiked fast configuration was not effective."
      }
      output
    }, error = function(error) {
      failure <- cbind(base, data.frame(
        theta_hat = NA_character_, lambda_hat = NA_real_, sigma_eigen_max_error = NA_real_,
        score_mean_norm = NA_real_, score_mean_max_abs = NA_real_, ks_statistic = NA_real_,
        cvm_statistic = NA_real_, ks_pvalue = NA_real_, cvm_pvalue = NA_real_,
        ks_reject = NA, cvm_reject = NA, effective_bootstrap_method = NA_character_,
        fallback_to_reestimated = NA, fast_backend = NA_character_, fast_kernel = NA_character_,
        fast_fused = NA, derivative_method = NA_character_, vhat_method = NA_character_,
        sample_grid_mode = NA_character_, shared_sample_ks_cvm_cache = NA,
        lightweight_ks_prep = NA, lightweight_cvm_prep = NA,
        elapsed_seconds = proc.time()[["elapsed"]] - elapsed_start,
        stringsAsFactors = FALSE
      ))
      failure$status <- "error"
      failure$error_message <- conditionMessage(error)
      failure
    })
  }

  chunks <- split(seq_len(nrow(pending)),
                  ceiling(seq_len(nrow(pending)) / as.integer(checkpoint_results)))
  for (chunk_index in seq_along(chunks)) {
    indices <- chunks[[chunk_index]]
    rows <- parallel::mclapply(
      indices, function(i) run_job(pending[i, , drop = FALSE]),
      mc.cores = min(as.integer(cores), length(indices)),
      mc.preschedule = FALSE, mc.set.seed = FALSE
    )
    rows <- do.call(rbind, rows)
    new_keys <- restricted_spiked_design_key(rows, include_replication = TRUE)
    if (nrow(existing)) {
      existing <- existing[!restricted_spiked_design_key(
        existing, include_replication = TRUE
      ) %in% new_keys, , drop = FALSE]
    }
    existing <- rbind(existing, rows)
    existing <- existing[order(existing$design_id, existing$replication), , drop = FALSE]
    restricted_spiked_write_atomic_csv(existing, results_path)
    restricted_spiked_write_atomic_csv(restricted_spiked_summarize(existing), summary_path)
    completed <- nrow(design) - nrow(pending) + max(indices)
    restricted_spiked_write_status(
      status_path, nrow(design), completed, existing, started, cores
    )
    cat(sprintf("%s completed=%d/%d\n", format(Sys.time(), tz = "Europe/Madrid"),
                completed, nrow(design)), file = log_path, append = TRUE)
    if (isTRUE(show_progress)) {
      message(sprintf("completed %d/%d; pending %d", completed, nrow(design),
                      nrow(design) - completed))
    }
  }
  invisible(list(results = existing, summary = restricted_spiked_summarize(existing)))
}

if (sys.nframe() == 0L) {
  args <- restricted_spiked_parse_args(commandArgs(trailingOnly = TRUE))
  mean_config <- as.character(args$mean_config %||% "")
  if (!nzchar(mean_config)) {
    stop(sprintf(
      "`--mean_config` is required; choose one of: %s.",
      paste(names(restricted_spiked_mean_catalog()), collapse = ", ")
    ))
  }
  M <- as.integer(args$M %||% 1000L)
  B <- as.integer(args$B %||% 5000L)
  derivative_mc_size <- as.integer(args$derivative_mc_size %||% 10000L)
  output_dir <- as.character(args$output_dir %||% file.path(
    "simulation_results",
    sprintf(
      paste0(
        "restricted_spiked_normal_double_spike_lambda2_%s_",
        "d2_d5_n50_100_200_400_M%d_B%d_Nderiv%d"
      ),
      mean_config, M, B, derivative_mc_size
    )
  ))
  run_restricted_spiked_covariance_alternatives(
    mean_config = mean_config, output_dir = output_dir, M = M, B = B,
    dimensions = restricted_spiked_parse_csv(args$dimensions, c(2L, 5L), TRUE),
    n_values = restricted_spiked_parse_csv(args$n_values %||% args$n,
                                            c(50L, 100L, 200L, 400L), TRUE),
    beta_values = restricted_spiked_parse_csv(args$beta_values %||% args$betas,
                                               c(0, 0.25, 0.5, 1)),
    lambda = as.numeric(args$lambda %||% 2),
    derivative_mc_size = derivative_mc_size,
    cvm_block_size = as.integer(args$cvm_block_size %||% 50L),
    cores = as.integer(args$cores %||% 4L),
    checkpoint_results = as.integer(args$checkpoint_results %||% 2000L),
    base_seed = as.integer(args$seed %||% restricted_spiked_default_seed(mean_config)),
    show_progress = !identical(tolower(as.character(args$show_progress %||% "true")),
                               "false")
  )
}
