#!/usr/bin/env Rscript

# Resumable Section 6 pilot for the orthogonal vMF mixture when the null is
# {vMF(mu, kappa): mu in S^d}, with kappa fixed and only mu estimated.
# This runner and bootstrap/vmf_fixed_kappa_model_spec.R are isolated from the
# existing vMF model, which still estimates xi = kappa * mu.

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

source("scripts/run_section6_new_scenarios.R", encoding = "UTF-8")
source("bootstrap/vmf_fixed_kappa_model_spec.R", encoding = "UTF-8")

`%||%` <- function(x, y) if (is.null(x)) y else x

vmf_mu_only_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  values <- commandArgs(trailingOnly = TRUE)
  values <- values[startsWith(values, prefix)]
  if (!length(values)) return(default)
  substring(values[[length(values)]], nchar(prefix) + 1L)
}

vmf_mu_only_csv <- function(name, default, storage.mode = c("integer", "numeric")) {
  storage.mode <- match.arg(storage.mode)
  raw <- vmf_mu_only_arg(name, NULL)
  if (is.null(raw) || !nzchar(raw)) return(default)
  values <- trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])
  switch(storage.mode, integer = as.integer(values), numeric = as.numeric(values))
}

vmf_mu_only_fixed_kappa_design <- function(kappa_values = c(1, 2),
                                            dimensions = c(2L, 5L),
                                            n_values = c(50L, 100L, 200L, 400L),
                                            beta_values = c(0, 0.5, 1)) {
  grid <- expand.grid(
    kappa = sort(unique(as.numeric(kappa_values))),
    d = sort(unique(as.integer(dimensions))),
    n = sort(unique(as.integer(n_values))),
    beta = sort(unique(as.numeric(beta_values))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(grid$kappa)) || any(grid$kappa <= 0) ||
      any(!is.finite(grid$d)) || any(grid$d < 2L) ||
      any(!is.finite(grid$n)) || any(grid$n < 1L) ||
      any(!is.finite(grid$beta)) || any(grid$beta < 0 | grid$beta > 1)) {
    stop("Invalid fixed-kappa vMF design grid.", call. = FALSE)
  }
  grid$scenario <- paste0(
    "vmf_31_orthogonal_mu_only_kappa_",
    gsub("\\.", "p", formatC(grid$kappa, digits = 12L, format = "fg", flag = "#"))
  )
  grid$family <- "vmf_fixed_kappa"
  grid$alternative <- "orthogonal_vmf_mixture"
  grid$description <- sprintf(
    "(1-beta/2) vMF(e1,%.12g) + (beta/2) vMF(e2,%.12g); null fits mu only",
    grid$kappa, grid$kappa
  )
  grid$design_id <- seq_len(nrow(grid))
  grid[, c("scenario", "family", "alternative", "description", "kappa", "d", "n", "beta", "design_id")]
}

generate_vmf_mu_only_fixed_kappa_sample <- function(job) {
  ambient_dim <- as.integer(job$d) + 1L
  mu0 <- section6_e(ambient_dim, index = 1L)
  mu1 <- section6_e(ambient_dim, index = 2L)
  n <- as.integer(job$n)
  choose_alt <- stats::runif(n) < as.numeric(job$beta) / 2
  x <- rotasym::r_vMF(n, mu = mu0, kappa = as.numeric(job$kappa))
  if (any(choose_alt)) {
    x[choose_alt, ] <- rotasym::r_vMF(
      sum(choose_alt), mu = mu1, kappa = as.numeric(job$kappa)
    )
  }
  x
}

vmf_mu_only_empty_results <- function() {
  output <- empty_section6_results()
  output$kappa <- numeric()
  output$unknown_param <- character()
  output
}

vmf_mu_only_fixed_kappa_job <- function(job, B, base_seed, cvm_block_size) {
  started <- proc.time()[["elapsed"]]
  data_seed <- section6_seed(base_seed, job$design_id, job$rep, 0L)
  bootstrap_seed <- section6_seed(base_seed, job$design_id, job$rep, 1L)
  out <- data.frame(
    scenario = as.character(job$scenario), family = "vmf_fixed_kappa",
    alternative = as.character(job$alternative), d = as.integer(job$d), n = as.integer(job$n),
    beta = as.numeric(job$beta), design_id = as.integer(job$design_id), rep = as.integer(job$rep),
    status = "ok", error_message = NA_character_, ks_pvalue = NA_real_, cvm_pvalue = NA_real_,
    ks_reject = NA, cvm_reject = NA, bootstrap_method_requested = "fast_multiplier",
    bootstrap_method_effective = NA_character_, fallback_to_reestimated = NA,
    derivative_method_requested = "quadrature", derivative_method_effective = NA_character_,
    derivative_method_selection_source = "fixed_kappa_quadrature", quadrature_algorithm = NA_character_,
    quadrature_abs_tol = NA_real_, quadrature_max_terms = NA_integer_,
    quadrature_max_terms_used = NA_integer_, quadrature_max_residual_error_estimate = NA_real_,
    ks_grid = "sample_unique_distances", fast_multiplier_backend_requested = "cpp",
    fast_multiplier_backend_effective = NA_character_,
    fast_multiplier_cpp_kernel_requested = "contiguous_double",
    fast_multiplier_cpp_kernel_effective = NA_character_,
    fast_multiplier_fuse_ks_cvm_requested = TRUE,
    fast_multiplier_fuse_ks_cvm_effective = NA,
    seed_data = data_seed, seed_bootstrap = bootstrap_seed, seed_derivative = NA_integer_,
    elapsed_seconds = NA_real_, kappa = as.numeric(job$kappa),
    unknown_param = "mu_fixed_kappa", stringsAsFactors = FALSE
  )
  out <- tryCatch({
    set.seed(data_seed)
    x <- generate_vmf_mu_only_fixed_kappa_sample(job)
    fit <- multiplier_bootstrap_vmf_fixed_kappa(
      data = x, kappa = job$kappa, null = list(type = "composite"),
      statistics = c("ks", "cvm"), ks_grid = make_sample_unique_distance_ks_grid(),
      B = as.integer(B), alpha = 0.05, n_cores = 1L, seed = bootstrap_seed,
      bootstrap_method = "fast_multiplier",
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE),
      control = list(
        derivative_method = "quadrature", fast_multiplier_cvm_block_size = as.integer(cvm_block_size),
        fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
        fast_multiplier_fuse_ks_cvm = TRUE, fast_multiplier_cache_corrections = "auto",
        fast_multiplier_stream_chunk_size = 100L, vmf_profile_method = "tabulated",
        vmf_profile_n_u = 4097L
      ),
      distance_type = "geodesic", distance_profile_backend = "r",
      fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
      fuse_ks_cvm = TRUE, cache_block_corrections = "auto"
    )
    diagnostics <- fit$diagnostics
    out$ks_pvalue <- fit$inference$ks$p_value
    out$cvm_pvalue <- fit$inference$cvm$p_value
    out$ks_reject <- fit$inference$ks$reject
    out$cvm_reject <- fit$inference$cvm$reject
    out$bootstrap_method_effective <- diagnostics$effective_bootstrap_method %||% NA_character_
    out$derivative_method_effective <- diagnostics$derivative_method_effective %||%
      diagnostics$derivative_method %||% NA_character_
    out$quadrature_algorithm <- diagnostics$quadrature_algorithm %||% NA_character_
    out$quadrature_abs_tol <- diagnostics$quadrature_abs_tol %||% NA_real_
    out$quadrature_max_terms <- diagnostics$quadrature_max_terms %||% NA_integer_
    out$quadrature_max_terms_used <- diagnostics$quadrature_max_terms_used %||% NA_integer_
    out$quadrature_max_residual_error_estimate <-
      diagnostics$quadrature_max_residual_error_estimate %||% NA_real_
    out$fallback_to_reestimated <- isTRUE(diagnostics$fallback_to_reestimated)
    out$fast_multiplier_backend_effective <- diagnostics$fast_multiplier_backend_effective %||% NA_character_
    out$fast_multiplier_cpp_kernel_effective <-
      diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_
    out$fast_multiplier_fuse_ks_cvm_effective <-
      isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective)
    conforming <- identical(out$bootstrap_method_effective, "fast_multiplier") &&
      identical(out$derivative_method_effective, "quadrature") &&
      identical(out$fast_multiplier_backend_effective, "cpp") &&
      identical(out$fast_multiplier_cpp_kernel_effective, "contiguous_double") &&
      isTRUE(out$fast_multiplier_fuse_ks_cvm_effective) && !isTRUE(out$fallback_to_reestimated)
    if (!conforming) {
      out$status <- "nonconforming"
      out$error_message <- "Requested fixed-kappa quadrature, fused C++ fast bootstrap was not effective."
    }
    out
  }, error = function(error) {
    out$status <- "error"
    out$error_message <- conditionMessage(error)
    out
  })
  out$elapsed_seconds <- proc.time()[["elapsed"]] - started
  out
}

summarize_vmf_mu_only_results <- function(results) {
  ok <- results[results$status == "ok", , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  keys <- interaction(ok$scenario, ok$kappa, ok$d, ok$n, ok$beta, drop = TRUE)
  do.call(rbind, lapply(split(ok, keys), function(x) data.frame(
    scenario = x$scenario[[1L]], family = x$family[[1L]], alternative = x$alternative[[1L]],
    kappa = x$kappa[[1L]], d = x$d[[1L]], n = x$n[[1L]], beta = x$beta[[1L]],
    M = nrow(x), rejection_ks = mean(x$ks_reject), rejection_cvm = mean(x$cvm_reject),
    mean_elapsed_seconds = mean(x$elapsed_seconds), stringsAsFactors = FALSE
  )))
}

vmf_mu_only_validate_manifest <- function(path, design, M, B, seed, cvm_block_size) {
  if (!file.exists(path)) return(invisible(TRUE))
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("scenario", "kappa", "d", "n", "beta", "M", "B", "base_seed", "cvm_block_size",
                "unknown_param", "derivative_method", "ks_grid", "bootstrap_method",
                "fast_multiplier_backend", "fast_multiplier_cpp_kernel", "fused_ks_cvm_kernel")
  if (!all(required %in% names(manifest))) {
    stop("Existing fixed-kappa vMF manifest is incomplete; use a new output directory.", call. = FALSE)
  }
  current <- design[order(design$design_id), , drop = FALSE]
  saved <- manifest[order(manifest$design_id), , drop = FALSE]
  design_match <- nrow(saved) == nrow(current) &&
    identical(as.character(saved$scenario), as.character(current$scenario)) &&
    isTRUE(all.equal(as.numeric(saved$kappa), as.numeric(current$kappa), tolerance = 0)) &&
    identical(as.integer(saved$d), as.integer(current$d)) &&
    identical(as.integer(saved$n), as.integer(current$n)) &&
    isTRUE(all.equal(as.numeric(saved$beta), as.numeric(current$beta), tolerance = 0))
  controls_match <- all(as.integer(saved$M) == as.integer(M)) &&
    all(as.integer(saved$B) == as.integer(B)) &&
    all(as.integer(saved$base_seed) == as.integer(seed)) &&
    all(as.integer(saved$cvm_block_size) == as.integer(cvm_block_size)) &&
    all(saved$unknown_param == "mu_fixed_kappa") &&
    all(saved$derivative_method == "quadrature") &&
    all(saved$ks_grid == "sample_unique_distances") &&
    all(saved$bootstrap_method == "fast_multiplier") &&
    all(saved$fast_multiplier_backend == "cpp") &&
    all(saved$fast_multiplier_cpp_kernel == "contiguous_double") &&
    all(as.logical(saved$fused_ks_cvm_kernel))
  if (!isTRUE(design_match) || !isTRUE(controls_match)) {
    stop("Existing fixed-kappa vMF output directory has an incompatible manifest; do not resume into it.", call. = FALSE)
  }
  invisible(TRUE)
}

run_vmf_mu_only_fixed_kappa_pilot <- function(
    output_dir,
    kappa_values = c(1, 2),
    dimensions = c(2L, 5L),
    n_values = c(50L, 100L, 200L, 400L),
    beta_values = c(0, 0.5, 1),
    M = 100L,
    B = 1000L,
    cores = 6L,
    seed = 20260907L,
    cvm_block_size = 50L,
    checkpoint_results = 12L,
    show_progress = TRUE) {
  cores <- as.integer(cores)
  M <- as.integer(M)
  B <- as.integer(B)
  checkpoint_results <- as.integer(checkpoint_results)
  if (!is.finite(cores) || cores < 1L || !is.finite(M) || M < 1L ||
      !is.finite(B) || B < 1L || !is.finite(checkpoint_results) || checkpoint_results < 1L) {
    stop("`M`, `B`, `cores`, and `checkpoint_results` must be positive integers.", call. = FALSE)
  }
  ensure_distance_profile_cpp_loaded()
  design <- vmf_mu_only_fixed_kappa_design(kappa_values, dimensions, n_values, beta_values)
  jobs <- merge(design, data.frame(rep = seq_len(M)), by = NULL)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(lock), add = TRUE)
  manifest_path <- file.path(output_dir, "manifest.csv")
  results_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  vmf_mu_only_validate_manifest(manifest_path, design, M, B, seed, cvm_block_size)
  if (!file.exists(manifest_path)) {
    manifest <- transform(
      design, M = M, B = B, cores = cores, base_seed = as.integer(seed),
      unknown_param = "mu_fixed_kappa", derivative_method = "quadrature",
      ks_grid = "sample_unique_distances", bootstrap_method = "fast_multiplier",
      fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
      fused_ks_cvm_kernel = TRUE, cvm_block_size = cvm_block_size
    )
    section6_write_atomic_csv(manifest, manifest_path)
  }
  existing <- if (file.exists(results_path)) {
    utils::read.csv(results_path, stringsAsFactors = FALSE)
  } else {
    vmf_mu_only_empty_results()
  }
  completed_keys <- if (nrow(existing)) {
    section6_design_key(existing[existing$status == "ok", , drop = FALSE], include_rep = TRUE)
  } else {
    character()
  }
  pending <- jobs[!section6_design_key(jobs, include_rep = TRUE) %in% completed_keys, , drop = FALSE]
  started <- Sys.time()
  completed <- nrow(jobs) - nrow(pending)
  section6_write_status(status_path, "vmf_fixed_kappa", nrow(jobs), completed, existing, started, cores)
  if (isTRUE(show_progress)) section6_progress(completed, nrow(jobs), started, existing, cores)
  if (!nrow(pending)) return(invisible(list(results = existing, summary = summarize_vmf_mu_only_results(existing))))

  # A dynamic worker queue keeps all requested outer workers busy despite
  # heterogeneous dimensions and sample sizes; bootstrap work remains serial.
  workers <- min(cores, nrow(pending))
  next_job <- 1L
  finished <- 0L
  buffer <- vmf_mu_only_empty_results()
  active <- list()
  on.exit({
    if (length(active)) try(parallel:::mckill(unname(active), signal = 2L), silent = TRUE)
  }, add = TRUE)
  submit <- function(index) {
    process <- parallel::mcparallel(
      vmf_mu_only_fixed_kappa_job(pending[index, , drop = FALSE], B, seed, cvm_block_size),
      mc.set.seed = FALSE, silent = TRUE
    )
    active[[as.character(process$pid)]] <<- process
  }
  for (i in seq_len(workers)) {
    submit(next_job)
    next_job <- next_job + 1L
  }
  while (finished < nrow(pending)) {
    received <- parallel::mccollect(active, wait = FALSE, timeout = -1)
    if (is.null(received)) next
    for (pid in names(received)) {
      active[[pid]] <- NULL
      row <- received[[pid]]
      if (!is.data.frame(row) || nrow(row) != 1L) {
        stop("A fixed-kappa vMF worker returned an invalid result.", call. = FALSE)
      }
      buffer <- rbind(buffer, row)
      finished <- finished + 1L
      if (next_job <= nrow(pending)) {
        submit(next_job)
        next_job <- next_job + 1L
      }
      if (nrow(buffer) >= checkpoint_results || finished == nrow(pending)) {
        existing <- rbind(existing, buffer)
        existing <- existing[order(existing$design_id, existing$rep), , drop = FALSE]
        section6_write_atomic_csv(existing, results_path)
        section6_write_atomic_csv(summarize_vmf_mu_only_results(existing), summary_path)
        completed <- nrow(jobs) - nrow(pending) + finished
        section6_write_status(status_path, "vmf_fixed_kappa", nrow(jobs), completed, existing, started, cores)
        buffer <- vmf_mu_only_empty_results()
      }
      if (isTRUE(show_progress)) {
        section6_progress(nrow(jobs) - nrow(pending) + finished, nrow(jobs), started, existing, cores)
      }
    }
  }
  invisible(list(results = existing, summary = summarize_vmf_mu_only_results(existing)))
}

if (sys.nframe() == 0L) {
  run_vmf_mu_only_fixed_kappa_pilot(
    output_dir = vmf_mu_only_arg(
      "output_dir",
      "simulation_results/section6_new_scenarios/pilot_vmf_31_orthogonal_mu_only_kappa1_2_M100_B1000"
    ),
    kappa_values = vmf_mu_only_csv("kappa_values", c(1, 2), "numeric"),
    dimensions = vmf_mu_only_csv("dimensions", c(2L, 5L), "integer"),
    n_values = vmf_mu_only_csv("n_values", c(50L, 100L, 200L, 400L), "integer"),
    beta_values = vmf_mu_only_csv("beta_values", c(0, 0.5, 1), "numeric"),
    M = as.integer(vmf_mu_only_arg("M", "100")),
    B = as.integer(vmf_mu_only_arg("B", "1000")),
    cores = as.integer(vmf_mu_only_arg("cores", "6")),
    seed = as.integer(vmf_mu_only_arg("seed", "20260907")),
    cvm_block_size = as.integer(vmf_mu_only_arg("cvm_block_size", "50")),
    checkpoint_results = as.integer(vmf_mu_only_arg("checkpoint_results", "12")),
    show_progress = !identical(tolower(vmf_mu_only_arg("show_progress", "true")), "false")
  )
}
