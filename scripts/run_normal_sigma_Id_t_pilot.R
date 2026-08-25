#!/usr/bin/env Rscript

# Resumable GOF pilot for the normal-sigma-Id versus multivariate-t scenario:
#
#   P_beta = (1-beta) N_d(0, I_d) + beta T_{nu,d}^{std}.
#
# The fitted null is N_d(mu, sigma^2 I_d), with the normal-sigma-Id adapter.
# It does not modify the unrestricted multivariate-normal experiment runner.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    fields <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[fields[[1L]]]] <- if (length(fields) == 1L) "TRUE" else paste(fields[-1L], collapse = "=")
  }
  out
}

parse_csv <- function(value, default, integer = FALSE) {
  if (is.null(value) || !nzchar(value)) return(default)
  out <- as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!is.finite(out))) stop("Invalid comma-separated option.")
  if (isTRUE(integer)) {
    if (any(out != as.integer(out))) stop("Expected integer values.")
    out <- as.integer(out)
  }
  out
}

standardized_multivariate_t <- function(n, d, nu) {
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  w <- stats::rchisq(n, df = nu)
  z * sqrt((nu - 2) / w)
}

normal_sigma_Id_t_design <- function(dimensions, n_values, beta_values, M) {
  design <- expand.grid(
    d = sort(unique(as.integer(dimensions))),
    n = sort(unique(as.integer(n_values))),
    beta = sort(unique(as.numeric(beta_values))),
    stringsAsFactors = FALSE
  )
  design$design_id <- seq_len(nrow(design))
  merge(design, data.frame(replication = seq_len(as.integer(M))), by = NULL)
}

normal_sigma_Id_t_key <- function(x) {
  paste(x$d, x$n, formatC(x$beta, digits = 16L, format = "fg", flag = "#"),
        x$replication, sep = "|")
}

normal_sigma_Id_t_seed <- function(base_seed, design_id, replication, stream) {
  as.integer((as.numeric(base_seed) + 1000003 * design_id +
    1009 * replication + 10000019 * stream) %% 2147483647) + 1L
}

empty_results <- function() {
  data.frame(
    d = integer(), n = integer(), beta = numeric(), design_id = integer(), replication = integer(),
    status = character(), error_message = character(), warning_message = character(),
    seed_data = integer(), seed_bootstrap = integer(), seed_derivative = integer(),
    mu_hat = character(), sigma_hat = numeric(), score_mean_norm = numeric(),
    ks_statistic = numeric(), cvm_statistic = numeric(), ks_pvalue = numeric(), cvm_pvalue = numeric(),
    ks_reject = logical(), cvm_reject = logical(), effective_bootstrap_method = character(),
    fallback_to_reestimated = logical(), fast_backend = character(), fast_kernel = character(),
    fast_fused = logical(), derivative_method = character(), derivative_mc_size = integer(),
    vhat_method = character(), elapsed_seconds = numeric(), stringsAsFactors = FALSE
  )
}

write_atomic_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop(sprintf("Could not update '%s'.", path))
}

summarize_results <- function(x) {
  x <- x[x$status == "ok", , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  groups <- split(x, interaction(x$d, x$n, x$beta, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) data.frame(
    d = z$d[[1L]], n = z$n[[1L]], beta = z$beta[[1L]], M = nrow(z),
    ks_rejection_percent = 100 * mean(z$ks_reject),
    cvm_rejection_percent = 100 * mean(z$cvm_reject),
    ks_pvalue_mean = mean(z$ks_pvalue), cvm_pvalue_mean = mean(z$cvm_pvalue),
    mean_elapsed_seconds = mean(z$elapsed_seconds), stringsAsFactors = FALSE
  )))
}

conforming <- function(x, derivative_mc_size) {
  x$status == "ok" & x$effective_bootstrap_method == "fast_multiplier" &
    !x$fallback_to_reestimated & x$fast_backend == "cpp" &
    x$fast_kernel == "contiguous_double" & x$fast_fused &
    x$derivative_method == "score_mc" &
    x$derivative_mc_size == as.integer(derivative_mc_size) &
    x$vhat_method == "normal_sigma_Id_analytic_fisher"
}

write_status <- function(path, total, results, started, cores, derivative_mc_size) {
  ok <- if (nrow(results)) conforming(results, derivative_mc_size) else logical()
  per_job <- if (any(ok)) mean(results$elapsed_seconds[ok]) else NA_real_
  remaining <- total - sum(ok)
  eta <- if (is.finite(per_job)) per_job * remaining / cores else NA_real_
  writeLines(c(
    sprintf("updated: %s", format(Sys.time(), tz = "Europe/Madrid")),
    sprintf("completed: %d/%d", sum(ok), total), sprintf("pending: %d", remaining),
    sprintf("elapsed_seconds: %.1f", as.numeric(difftime(Sys.time(), started, units = "secs"))),
    sprintf("mean_seconds_per_conforming_job: %s", if (is.finite(per_job)) per_job else "NA"),
    sprintf("eta_seconds: %s", if (is.finite(eta)) round(eta) else "NA"),
    sprintf("ok: %d", sum(results$status == "ok")),
    sprintf("errors: %d", sum(results$status == "error")),
    sprintf("nonconforming: %d", sum(results$status == "nonconforming"))
  ), path)
}

run_normal_sigma_Id_t_pilot <- function(output_dir, M = 100L, B = 499L,
                                           dimensions = 5L,
                                           n_values = c(50L, 100L, 200L, 400L),
                                           beta_values = c(0, 0.5, 1), nu = 3,
                                           derivative_mc_size = 10000L,
                                           cvm_block_size = 50L, cores = 10L,
                                           checkpoint_results = 100L,
                                           base_seed = 20260825L,
                                           show_progress = TRUE) {
  if (nu <= 2 || any(dimensions < 2L) || any(n_values < 2L) ||
      any(beta_values < 0 | beta_values > 1) || M < 1L || B < 1L ||
      derivative_mc_size < 1L || cvm_block_size < 1L || cores < 1L || checkpoint_results < 1L) {
    stop("Invalid normal-sigma-Id t pilot settings.")
  }
  if (.Platform$OS.type != "unix" && cores > 1L) stop("Outer parallelism requires a Unix platform.")
  source(file.path("bootstrap", "normal_sigma_Id_bootstrap.R"), local = environment())
  design <- normal_sigma_Id_t_design(dimensions, n_values, beta_values, M)
  manifest <- unique(design[c("d", "n", "beta", "design_id")])
  manifest <- transform(
    manifest, M = as.integer(M), B = as.integer(B), nu = as.numeric(nu),
    base_seed = as.integer(base_seed), derivative_mc_size = as.integer(derivative_mc_size),
    cvm_block_size = as.integer(cvm_block_size), null_model = "N_d(mu,sigma^2 I_d)",
    alternative = "standardized_multivariate_t", derivative_method = "score_mc",
    statistics = "ks,cvm", ks_grid = "sample_points_unique_distances",
    bootstrap_method = "fast_multiplier", fast_backend = "cpp",
    fast_kernel = "contiguous_double", fast_fused = TRUE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(output_dir, ".normal_sigma_Id_t.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Output directory is locked by another run.")
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)
  manifest_path <- file.path(output_dir, "manifest.csv")
  results_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  if (file.exists(manifest_path)) {
    old <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(old, manifest, check.attributes = FALSE))) {
      stop("Existing manifest is incompatible; use a new output directory.")
    }
  } else write_atomic_csv(manifest, manifest_path)
  results <- if (file.exists(results_path)) utils::read.csv(results_path, stringsAsFactors = FALSE) else empty_results()
  if (!all(names(empty_results()) %in% names(results))) stop("Existing raw results have an incompatible schema.")
  done <- if (nrow(results)) normal_sigma_Id_t_key(results[conforming(results, derivative_mc_size), , drop = FALSE]) else character()
  pending <- design[!normal_sigma_Id_t_key(design) %in% done, , drop = FALSE]
  started <- Sys.time()
  run_one <- function(row) {
    began <- proc.time()[["elapsed"]]
    data_seed <- normal_sigma_Id_t_seed(base_seed, row$design_id, row$replication, 0L)
    bootstrap_seed <- normal_sigma_Id_t_seed(base_seed, row$design_id, row$replication, 1L)
    derivative_seed <- normal_sigma_Id_t_seed(base_seed, row$design_id, row$replication, 2L)
    base <- data.frame(
      d = as.integer(row$d), n = as.integer(row$n), beta = as.numeric(row$beta),
      design_id = as.integer(row$design_id), replication = as.integer(row$replication),
      status = "ok", error_message = NA_character_, warning_message = NA_character_,
      seed_data = data_seed, seed_bootstrap = bootstrap_seed, seed_derivative = derivative_seed,
      stringsAsFactors = FALSE
    )
    tryCatch({
      set.seed(data_seed)
      x <- matrix(stats::rnorm(base$n * base$d), nrow = base$n, ncol = base$d)
      is_t <- stats::runif(base$n) < base$beta
      if (any(is_t)) x[is_t, ] <- standardized_multivariate_t(sum(is_t), base$d, nu)
      warnings <- character()
      fit <- withCallingHandlers(
        multiplier_bootstrap_normal_sigma_Id(
          data = x, null = list(type = "composite"), statistics = c("ks", "cvm"),
          ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = 0.05,
          n_cores = 1L, seed = bootstrap_seed, bootstrap_method = "fast_multiplier",
          keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE),
          control = list(derivative_method = "score_mc", derivative_mc_size = derivative_mc_size,
            derivative_mc_seed = derivative_seed, fast_multiplier_cvm_block_size = cvm_block_size,
            fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
            fast_multiplier_fuse_ks_cvm = TRUE, fast_multiplier_cache_corrections = "auto",
            fast_multiplier_stream_chunk_size = 100L),
          distance_profile_backend = "r", fast_multiplier_backend = "cpp",
          fast_multiplier_cpp_kernel = "contiguous_double", fuse_ks_cvm = TRUE,
          cache_block_corrections = "auto"
        ), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
        }
      )
      theta_hat <- fit$observed$theta_hat
      score <- normal_sigma_Id_score_matrix(x, theta_hat)
      diagnostics <- fit$diagnostics
      out <- cbind(base, data.frame(
        mu_hat = paste(theta_hat$mu, collapse = ";"), sigma_hat = theta_hat$sigma,
        score_mean_norm = sqrt(sum(colMeans(score)^2)),
        ks_statistic = fit$inference$ks$observed, cvm_statistic = fit$inference$cvm$observed,
        ks_pvalue = fit$inference$ks$p_value, cvm_pvalue = fit$inference$cvm$p_value,
        ks_reject = fit$inference$ks$reject, cvm_reject = fit$inference$cvm$reject,
        effective_bootstrap_method = diagnostics$effective_bootstrap_method %||% NA_character_,
        fallback_to_reestimated = isTRUE(diagnostics$fallback_to_reestimated),
        fast_backend = diagnostics$fast_multiplier_backend_effective %||% NA_character_,
        fast_kernel = diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_,
        fast_fused = isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective),
        derivative_method = diagnostics$derivative_method_effective %||% NA_character_,
        derivative_mc_size = diagnostics$derivative_mc_size %||% NA_integer_,
        vhat_method = diagnostics$vhat_method %||% NA_character_,
        elapsed_seconds = proc.time()[["elapsed"]] - began, stringsAsFactors = FALSE
      ))
      out$warning_message <- if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_
      if (!conforming(out, derivative_mc_size)) {
        out$status <- "nonconforming"
        out$error_message <- "Requested normal-sigma-Id fast configuration was not effective."
      }
      out
    }, error = function(e) {
      cbind(base, data.frame(
        mu_hat = NA_character_, sigma_hat = NA_real_, score_mean_norm = NA_real_,
        ks_statistic = NA_real_, cvm_statistic = NA_real_, ks_pvalue = NA_real_, cvm_pvalue = NA_real_,
        ks_reject = NA, cvm_reject = NA, effective_bootstrap_method = NA_character_,
        fallback_to_reestimated = NA, fast_backend = NA_character_, fast_kernel = NA_character_,
        fast_fused = NA, derivative_method = NA_character_, derivative_mc_size = NA_integer_,
        vhat_method = NA_character_, elapsed_seconds = proc.time()[["elapsed"]] - began,
        stringsAsFactors = FALSE
      ), status = "error", error_message = conditionMessage(e))
    })
  }
  total <- nrow(design)
  write_status(status_path, total, results, started, cores, derivative_mc_size)
  while (nrow(pending)) {
    batch_n <- min(as.integer(checkpoint_results), nrow(pending))
    batch <- pending[seq_len(batch_n), , drop = FALSE]
    pending <- pending[-seq_len(batch_n), , drop = FALSE]
    rows <- parallel::mclapply(seq_len(nrow(batch)), function(i) run_one(batch[i, , drop = FALSE]),
      mc.cores = min(as.integer(cores), nrow(batch)), mc.preschedule = FALSE, mc.set.seed = FALSE)
    rows <- do.call(rbind, rows)
    rows <- rows[order(rows$design_id, rows$replication), , drop = FALSE]
    keys <- normal_sigma_Id_t_key(rows)
    if (nrow(results)) results <- results[!normal_sigma_Id_t_key(results) %in% keys, , drop = FALSE]
    results <- rbind(results, rows)
    results <- results[order(results$design_id, results$replication), , drop = FALSE]
    write_atomic_csv(results, results_path)
    write_atomic_csv(summarize_results(results), summary_path)
    write_status(status_path, total, results, started, cores, derivative_mc_size)
    if (isTRUE(show_progress)) {
      cat(sprintf("\rcompleted %d/%d", sum(conforming(results, derivative_mc_size)), total)); flush.console()
    }
  }
  if (isTRUE(show_progress)) cat("\n")
  invisible(list(results = results, summary = summarize_results(results)))
}

if (sys.nframe() == 0L) {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  M <- as.integer(args$M %||% 100L)
  B <- as.integer(args$B %||% 499L)
  output_dir <- as.character(args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios",
    sprintf("pilot_normal_sigma_Id_t_d5_nu3_M%d_B%d_scoremc_Nderiv10000", M, B)
  ))
  run_normal_sigma_Id_t_pilot(
    output_dir = output_dir, M = M, B = B,
    dimensions = parse_csv(args$dimensions, 5L, integer = TRUE),
    n_values = parse_csv(args$n_values %||% args$n, c(50L, 100L, 200L, 400L), integer = TRUE),
    beta_values = parse_csv(args$beta_values %||% args$betas, c(0, 0.5, 1)),
    nu = as.numeric(args$nu %||% 3),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 10000L),
    cvm_block_size = as.integer(args$cvm_block_size %||% 50L),
    cores = as.integer(args$cores %||% 10L),
    checkpoint_results = as.integer(args$checkpoint_results %||% 100L),
    base_seed = as.integer(args$seed %||% 20260825L),
    show_progress = !identical(tolower(as.character(args$show_progress %||% "true")), "false")
  )
}
