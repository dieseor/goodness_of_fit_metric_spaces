#!/usr/bin/env Rscript

# Resumable Section 6 pilot for the dimension-indexed antipodal vMF mixture
# with a concentration fixed across dimensions.  This is deliberately kept
# outside run_section6_new_scenarios.R: it does not alter the paper catalogue
# (whose current vMF 3.1 parametrisation has kappa = d).

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

source("scripts/run_section6_new_scenarios.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

pilot_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  values <- commandArgs(trailingOnly = TRUE)
  values <- values[startsWith(values, prefix)]
  if (!length(values)) return(default)
  substring(values[[length(values)]], nchar(prefix) + 1L)
}

pilot_csv <- function(name, default, type = c("integer", "numeric")) {
  type <- match.arg(type)
  raw <- pilot_arg(name, NULL)
  if (is.null(raw) || !nzchar(raw)) return(default)
  values <- strsplit(raw, ",", fixed = TRUE)[[1L]]
  values <- trimws(values)
  switch(type, integer = as.integer(values), numeric = as.numeric(values))
}

fixed_kappa_design <- function(dimensions, n_values, beta_values, kappa,
                               kappa_rule = c("fixed", "sqrt_d"),
                               scenario_type = c(
                                 "antipodal",
                                 "projected_normal_half_concentration",
                                 "orthogonal_vmf_mixture",
                                 "projected_normal_mean_d",
                                 "projected_normal_sqrt_d",
                                 "projected_normal_sqrt_d_kappa_2d",
                                 "projected_normal_2sqrt_d_kappa_2d"
)) {
  kappa_rule <- match.arg(kappa_rule)
  scenario_type <- match.arg(scenario_type)
  design <- expand.grid(
    d = sort(unique(as.integer(dimensions))),
    n = sort(unique(as.integer(n_values))),
    beta = sort(unique(as.numeric(beta_values))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(design$d)) || any(design$d < 1L) ||
      any(!is.finite(design$n)) || any(design$n < 1L) ||
      any(!is.finite(design$beta)) || any(design$beta < 0 | design$beta > 1)) {
    stop("Invalid dimension, sample-size, or beta grid.", call. = FALSE)
  }
  if (identical(scenario_type, "projected_normal_half_concentration")) {
    design$kappa <- as.numeric(design$d)
    design$projected_normal_mean_norm <- sqrt(design$d / 2)
    design$alternative_mu_index <- NA_integer_
    scenario <- "vmf_2_projected_normal_half_concentration"
    alternative <- "projected_normal_mixture"
    description <- sprintf(
      "(1-beta) vMF(e1,%.12g) + beta Law(Z/||Z||), Z~N(sqrt(%.12g/2)e1,I)",
      design$kappa, design$d
    )
  } else if (identical(scenario_type, "projected_normal_mean_d")) {
    design$kappa <- 2 * as.numeric(design$d)
    design$projected_normal_mean_norm <- as.numeric(design$d)
    design$alternative_mu_index <- NA_integer_
    scenario <- "vmf_2_projected_normal_mean_d"
    alternative <- "projected_normal_mixture"
    description <- sprintf(
      "(1-beta) vMF(e1,2*%d) + beta Law(Z/||Z||), Z~N(%d e1,I)",
      design$d, design$d
    )
  } else if (identical(scenario_type, "projected_normal_2sqrt_d_kappa_2d")) {
    design$kappa <- 2 * as.numeric(design$d)
    design$projected_normal_mean_norm <- 2 * sqrt(as.numeric(design$d))
    design$alternative_mu_index <- NA_integer_
    scenario <- "vmf_2_projected_normal_2sqrt_d_kappa_2d"
    alternative <- "projected_normal_mixture"
    description <- sprintf(
      "(1-beta) vMF(e1,2*%d) + beta Law(Z/||Z||), Z~N(2sqrt(%d)e1,I)",
      design$d, design$d
    )
  } else if (identical(scenario_type, "projected_normal_sqrt_d_kappa_2d")) {
    design$kappa <- 2 * as.numeric(design$d)
    design$projected_normal_mean_norm <- sqrt(as.numeric(design$d))
    design$alternative_mu_index <- NA_integer_
    scenario <- "vmf_2_projected_normal_sqrt_d_kappa_2d"
    alternative <- "projected_normal_mixture"
    description <- sprintf(
      "(1-beta) vMF(e1,2*%d) + beta Law(Z/||Z||), Z~N(sqrt(%d)e1,I)",
      design$d, design$d
    )
  } else if (identical(scenario_type, "projected_normal_sqrt_d")) {
    design$kappa <- as.numeric(design$d)
    design$projected_normal_mean_norm <- sqrt(as.numeric(design$d))
    design$alternative_mu_index <- NA_integer_
    scenario <- "vmf_2_projected_normal_sqrt_d"
    alternative <- "projected_normal_mixture"
    description <- sprintf(
      "(1-beta) vMF(e1,%d) + beta Law(Z/||Z||), Z~N(sqrt(%d)e1,I)",
      design$d, design$d
    )
  } else if (identical(scenario_type, "orthogonal_vmf_mixture")) {
    design$kappa <- as.numeric(design$d)
    design$projected_normal_mean_norm <- NA_real_
    design$alternative_mu_index <- 2L
    scenario <- "vmf_1_orthogonal_mixture"
    alternative <- "orthogonal_vmf_mixture"
    description <- sprintf(
      "(1-beta/2) vMF(e1,%d) + (beta/2) vMF(e2,%d)",
      design$d, design$d
    )
  } else {
    design$kappa <- if (identical(kappa_rule, "sqrt_d")) sqrt(design$d) else as.numeric(kappa)
    design$projected_normal_mean_norm <- NA_real_
    design$alternative_mu_index <- NA_integer_
    scenario <- if (identical(kappa_rule, "sqrt_d")) {
      "vmf_1_antipodal_sqrt_d"
    } else {
      "vmf_1_antipodal_fixed_kappa"
    }
    alternative <- "antipodal_vmf_mixture"
    description <- sprintf(
      "(1-beta/2) vMF(e1,%.12g) + (beta/2) vMF(-e1,%.12g)",
      design$kappa, design$kappa
    )
  }
  transform(
    design,
    scenario = scenario,
    family = "vmf",
    alternative = alternative,
    description = description,
    design_id = seq_len(nrow(design))
  )[, c("scenario", "family", "alternative", "description", "kappa", "projected_normal_mean_norm", "alternative_mu_index", "d", "n", "beta", "design_id")]
}

generate_fixed_kappa_antipodal <- function(job) {
  mu <- section6_e(as.integer(job$d) + 1L)
  if (identical(as.character(job$alternative), "projected_normal_mixture")) {
    choose_alt <- stats::runif(as.integer(job$n)) < as.numeric(job$beta)
    x <- rotasym::r_vMF(as.integer(job$n), mu = mu, kappa = as.numeric(job$kappa))
    if (any(choose_alt)) {
      ambient_dim <- as.integer(job$d) + 1L
      z <- matrix(stats::rnorm(sum(choose_alt) * ambient_dim), ncol = ambient_dim)
      z[, 1L] <- z[, 1L] + as.numeric(job$projected_normal_mean_norm)
      x[choose_alt, ] <- z / sqrt(rowSums(z^2))
    }
    return(x)
  }
  choose_alt <- stats::runif(as.integer(job$n)) < as.numeric(job$beta) / 2
  x <- rotasym::r_vMF(as.integer(job$n), mu = mu, kappa = as.numeric(job$kappa))
  if (any(choose_alt)) {
    alternative_mu <- if (identical(as.character(job$alternative), "orthogonal_vmf_mixture")) {
      section6_e(as.integer(job$d) + 1L, index = 2L)
    } else {
      -mu
    }
    x[choose_alt, ] <- rotasym::r_vMF(sum(choose_alt), mu = alternative_mu, kappa = as.numeric(job$kappa))
  }
  x
}

fixed_kappa_job <- function(job, B, base_seed, cvm_block_size) {
  started <- proc.time()[["elapsed"]]
  data_seed <- section6_seed(base_seed, job$design_id, job$rep, 0L)
  bootstrap_seed <- section6_seed(base_seed, job$design_id, job$rep, 1L)
  out <- data.frame(
    scenario = as.character(job$scenario), family = "vmf", alternative = as.character(job$alternative),
    d = as.integer(job$d), n = as.integer(job$n), beta = as.numeric(job$beta),
    design_id = as.integer(job$design_id), rep = as.integer(job$rep), status = "ok", error_message = NA_character_,
    ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA, cvm_reject = NA,
    bootstrap_method_requested = "fast_multiplier", bootstrap_method_effective = NA_character_, fallback_to_reestimated = NA,
    derivative_method_requested = "quadrature", derivative_method_effective = NA_character_,
    derivative_method_selection_source = "explicit", quadrature_algorithm = NA_character_, quadrature_abs_tol = NA_real_,
    quadrature_max_terms = NA_integer_, quadrature_max_terms_used = NA_integer_, quadrature_max_residual_error_estimate = NA_real_,
    ks_grid = "sample_unique_distances", fast_multiplier_backend_requested = "cpp", fast_multiplier_backend_effective = NA_character_,
    fast_multiplier_cpp_kernel_requested = "contiguous_double", fast_multiplier_cpp_kernel_effective = NA_character_,
    fast_multiplier_fuse_ks_cvm_requested = TRUE, fast_multiplier_fuse_ks_cvm_effective = NA,
    seed_data = data_seed, seed_bootstrap = bootstrap_seed, seed_derivative = NA_integer_,
    elapsed_seconds = NA_real_, stringsAsFactors = FALSE
  )
  out <- tryCatch({
    set.seed(data_seed)
    x <- generate_fixed_kappa_antipodal(job)
    fit <- run_section6_bootstrap(
      design_row = job, x = x, B = B, seed = bootstrap_seed,
      derivative_seed = NA_integer_, derivative_mc_size = 1000L,
      cvm_block_size = cvm_block_size, derivative_method = "quadrature", n_cores = 1L
    )
    diagnostics <- fit$diagnostics
    out$ks_pvalue <- fit$inference$ks$p_value
    out$cvm_pvalue <- fit$inference$cvm$p_value
    out$ks_reject <- fit$inference$ks$reject
    out$cvm_reject <- fit$inference$cvm$reject
    out$bootstrap_method_effective <- diagnostics$effective_bootstrap_method %||% NA_character_
    out$derivative_method_effective <- diagnostics$derivative_method_effective %||% diagnostics$derivative_method %||% NA_character_
    out$quadrature_algorithm <- diagnostics$quadrature_algorithm %||% NA_character_
    out$quadrature_abs_tol <- diagnostics$quadrature_abs_tol %||% NA_real_
    out$quadrature_max_terms <- diagnostics$quadrature_max_terms %||% NA_integer_
    out$quadrature_max_terms_used <- diagnostics$quadrature_max_terms_used %||% NA_integer_
    out$quadrature_max_residual_error_estimate <- diagnostics$quadrature_max_residual_error_estimate %||% NA_real_
    out$fallback_to_reestimated <- isTRUE(diagnostics$fallback_to_reestimated)
    out$fast_multiplier_backend_effective <- diagnostics$fast_multiplier_backend_effective %||% NA_character_
    out$fast_multiplier_cpp_kernel_effective <- diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_
    out$fast_multiplier_fuse_ks_cvm_effective <- isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective)
    conforming <- identical(out$bootstrap_method_effective, "fast_multiplier") &&
      identical(out$derivative_method_effective, "quadrature") &&
      identical(out$fast_multiplier_backend_effective, "cpp") &&
      identical(out$fast_multiplier_cpp_kernel_effective, "contiguous_double") &&
      isTRUE(out$fast_multiplier_fuse_ks_cvm_effective) && !isTRUE(out$fallback_to_reestimated)
    if (!conforming) {
      out$status <- "nonconforming"
      out$error_message <- "Requested quadrature, fused C++ fast bootstrap was not effective."
    }
    out
  }, error = function(e) {
    out$status <- "error"
    out$error_message <- conditionMessage(e)
    out
  })
  out$elapsed_seconds <- proc.time()[["elapsed"]] - started
  out
}

run_fixed_kappa_pilot <- function(output_dir, kappa = 1, dimensions = c(2L, 5L),
                                  n_values = c(50L, 100L, 200L, 400L),
                                  beta_values = c(0, 0.5, 1), M = 100L, B = 1000L,
                                  cores = 6L, seed = 20260831L, cvm_block_size = 50L,
                                  checkpoint_results = 12L,
                                  kappa_rule = c("fixed", "sqrt_d"),
                                  scenario_type = c(
                                    "antipodal",
                                    "projected_normal_half_concentration",
                                    "orthogonal_vmf_mixture",
                                    "projected_normal_mean_d",
                                    "projected_normal_sqrt_d",
                                 "projected_normal_sqrt_d_kappa_2d",
                                 "projected_normal_2sqrt_d_kappa_2d"
)) {
  kappa_rule <- match.arg(kappa_rule)
  scenario_type <- match.arg(scenario_type)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) stop("`kappa` must be strictly positive.", call. = FALSE)
  design <- fixed_kappa_design(
    dimensions, n_values, beta_values, kappa,
    kappa_rule = kappa_rule, scenario_type = scenario_type
  )
  jobs <- merge(design, data.frame(rep = seq_len(as.integer(M))), by = NULL)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(lock), add = TRUE)
  manifest_path <- file.path(output_dir, "manifest.csv")
  results_path <- file.path(output_dir, "raw_results.csv")
  if (!file.exists(manifest_path)) {
    manifest <- transform(design, M = as.integer(M), B = as.integer(B), cores = as.integer(cores),
      base_seed = as.integer(seed), derivative_method = "quadrature", ks_grid = "sample_unique_distances",
      bootstrap_method = "fast_multiplier", fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double", fused_ks_cvm_kernel = TRUE,
      cvm_block_size = as.integer(cvm_block_size))
    section6_write_atomic_csv(manifest, manifest_path)
  }
  existing <- if (file.exists(results_path)) utils::read.csv(results_path, stringsAsFactors = FALSE) else empty_section6_results()
  completed_keys <- if (nrow(existing)) section6_design_key(existing[existing$status == "ok", ], include_rep = TRUE) else character()
  pending <- jobs[!section6_design_key(jobs, include_rep = TRUE) %in% completed_keys, , drop = FALSE]
  started <- Sys.time()
  total <- nrow(jobs)
  completed <- total - nrow(pending)
  section6_write_status(file.path(output_dir, "progress_status.txt"), "vmf", total, completed, existing, started, cores)
  section6_progress(completed, total, started, existing, cores)
  if (nrow(pending) == 0L) {
    cat("\n")
    return(invisible(existing))
  }
  for (start in seq.int(1L, nrow(pending), by = as.integer(checkpoint_results))) {
    index <- start:min(start + as.integer(checkpoint_results) - 1L, nrow(pending))
    rows <- parallel::mclapply(index, function(i) fixed_kappa_job(pending[i, , drop = FALSE], B, seed, cvm_block_size),
      mc.cores = min(as.integer(cores), length(index)), mc.preschedule = FALSE)
    existing <- rbind(existing, do.call(rbind, rows))
    existing <- existing[order(existing$design_id, existing$rep), , drop = FALSE]
    completed <- nrow(existing[existing$status == "ok", , drop = FALSE])
    section6_write_atomic_csv(existing, results_path)
    section6_write_atomic_csv(summarize_section6_results(existing), file.path(output_dir, "summary.csv"))
    section6_write_status(file.path(output_dir, "progress_status.txt"), "vmf", total, completed, existing, started, cores)
    section6_progress(completed, total, started, existing, cores)
  }
  cat("\n")
}

if (sys.nframe() == 0L) {
  run_fixed_kappa_pilot(
    output_dir = pilot_arg("output_dir", "simulation_results/section6_new_scenarios/pilot_vmf_31_fixed_kappa1_M100_B1000"),
    kappa = as.numeric(pilot_arg("kappa", "1")),
    kappa_rule = tolower(as.character(pilot_arg("kappa_rule", "fixed"))),
    scenario_type = tolower(as.character(pilot_arg("scenario_type", "antipodal"))),
    dimensions = pilot_csv("dimensions", c(2L, 5L), "integer"),
    n_values = pilot_csv("n_values", c(50L, 100L, 200L, 400L), "integer"),
    beta_values = pilot_csv("beta_values", c(0, 0.5, 1), "numeric"),
    M = as.integer(pilot_arg("M", "100")), B = as.integer(pilot_arg("B", "1000")),
    cores = as.integer(pilot_arg("cores", "6")), seed = as.integer(pilot_arg("seed", "20260831")),
    cvm_block_size = as.integer(pilot_arg("cvm_block_size", "50")),
    checkpoint_results = as.integer(pilot_arg("checkpoint_results", "12"))
  )
}
