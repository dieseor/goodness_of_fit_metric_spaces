#!/usr/bin/env Rscript

# Paired diagnostic only.  It compares the production fast bootstrap with a
# hybrid that retains the observed MLE in the statistic, but replaces the
# derivative correction by D_{theta0} and V_{theta0}.  The latter derivative is
# approximated by an independent, large auxiliary sample.  It is not a valid
# bootstrap calibration to be used in the paper.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/run_section6_new_scenarios.R")

normal_score_both <- function(x, theta) {
  x <- as.matrix(x)
  d <- ncol(x)
  Sinv <- solve(theta$Sigma)
  centered <- sweep(x, 2L, theta$mu, "-")
  score_mu <- centered %*% t(Sinv)
  score_sigma <- t(vapply(seq_len(nrow(x)), function(i) {
    rr <- centered[i, , drop = FALSE]
    fast_multiplier_sym_score_to_vech(
      0.5 * (Sinv %*% crossprod(rr) %*% Sinv - Sinv)
    )
  }, numeric(d * (d + 1L) / 2L)))
  cbind(score_mu, score_sigma)
}

normal_theta0_section6 <- function(job) {
  d <- as.integer(job$d)
  if (identical(as.character(job$scenario), "normal_1_mixture")) {
    return(normalize_mvnormal_theta(list(
      mu = 0.5 * section6_e(d), Sigma = section6_sigma(d, "plus")
    )))
  }
  if (identical(as.character(job$scenario), "normal_2_t3")) {
    return(normalize_mvnormal_theta(list(mu = rep(0, d), Sigma = diag(d))))
  }
  stop("Unsupported Normal scenario in the theta0-hybrid diagnostic.")
}

make_theta0_hybrid_normal_spec <- function(theta0) {
  spec <- make_mvnormal_spec(unknown_param = "both")
  spec$fast_multiplier_prepare <- function(data, theta_hat, ks_prep = NULL,
                                            cvm_prep = NULL, control = list()) {
    n_aux <- as.integer(control$hybrid_reference_mc_size)
    if (!is.finite(n_aux) || n_aux < 1000L) {
      stop("`hybrid_reference_mc_size` must be at least 1000.")
    }
    if (!is.null(control$hybrid_reference_seed)) set.seed(as.integer(control$hybrid_reference_seed))

    # This is the requested hybrid: psi remains evaluated at the observed MLE,
    # whereas D and V are evaluated at the known generating parameter theta0.
    score_obs_hat <- normal_score_both(data, theta_hat)
    V0 <- fast_multiplier_gaussian_paper_vhat(theta0$Sigma, "both")
    S_obs <- -score_obs_hat %*% t(solve(V0))
    aux <- rmvnormal_euclidean(n_aux, theta0$mu, theta0$Sigma)
    Psi_aux <- normal_score_both(aux, theta0)
    D_ks <- fast_multiplier_compute_D_ks(spec, aux, Psi_aux, ks_prep, control)
    D_cvm <- fast_multiplier_reuse_sample_ks_derivative_for_cvm(D_ks, cvm_prep)
    if (is.null(D_cvm)) {
      D_cvm <- fast_multiplier_compute_D_cvm(
        spec, aux, Psi_aux, data, cvm_prep, control
      )
    }
    list(
      S_obs = S_obs,
      Vhat = diag(ncol(S_obs)),
      Psi_aux = Psi_aux,
      correction_representation = "hybrid_Dtheta0_Vtheta0_score_thetahat",
      paper_Vhat = V0,
      paper_Vhat_inverse = solve(V0),
      paper_Vhat_method = "analytic_at_theta0",
      vhat_method = "hybrid_analytic_at_theta0",
      vhat_diagnostics = fast_multiplier_vhat_diagnostics(
        S_obs = S_obs, Psi_aux = Psi_aux, Vhat = diag(ncol(S_obs)), par0 = theta0$mu
      ),
      observed_cvm_distance_matrix = NULL,
      derivative_method = "score_mc_reference_at_theta0",
      derivative_mc_size = n_aux,
      derivative_mc_seed = control$hybrid_reference_seed %||% NA_integer_,
      D_ks = D_ks,
      D_cvm = D_cvm
    )
  }
  spec
}

run_hybrid_bootstrap <- function(x, theta_hat, theta0, B, bootstrap_seed,
                                 reference_seed, reference_mc_size,
                                 cvm_block_size) {
  spec <- make_theta0_hybrid_normal_spec(theta0)
  control <- section6_control(
    derivative_mc_size = reference_mc_size,
    derivative_seed = reference_seed,
    cvm_block_size = cvm_block_size
  )
  control$hybrid_reference_mc_size <- as.integer(reference_mc_size)
  control$hybrid_reference_seed <- as.integer(reference_seed)
  multiplier_bootstrap_gof(
    data = x, spec = spec, null = list(type = "composite"),
    statistics = c("ks", "cvm"), ks_grid = make_sample_unique_distance_ks_grid(),
    B = as.integer(B), alpha = .05, n_cores = 1L, seed = as.integer(bootstrap_seed),
    observed_theta_hat = theta_hat, bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                bootstrap_thetas = FALSE), control = control
  )
}

diagnostic_row <- function(job, rep, procedure, fit, elapsed, seeds) {
  data.frame(
    scenario = as.character(job$scenario), d = as.integer(job$d), n = as.integer(job$n),
    rep = as.integer(rep), procedure = procedure, status = "ok", error_message = NA_character_,
    ks_pvalue = fit$inference$ks$p_value, cvm_pvalue = fit$inference$cvm$p_value,
    ks_reject = fit$inference$ks$reject, cvm_reject = fit$inference$cvm$reject,
    method = fit$diagnostics$effective_bootstrap_method,
    backend = fit$diagnostics$fast_multiplier_backend_effective,
    fused = fit$diagnostics$fast_multiplier_fuse_ks_cvm_effective,
    seed_data = seeds$data, seed_bootstrap = seeds$bootstrap,
    seed_reference = seeds$reference, elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )
}

run_one_pair <- function(job, rep, B, base_seed, reference_mc_size, cvm_block_size) {
  started <- proc.time()[["elapsed"]]
  seeds <- list(
    data = section6_seed(base_seed, job$design_id, rep, 0L),
    bootstrap = section6_seed(base_seed, job$design_id, rep, 1L),
    derivative = section6_seed(base_seed, job$design_id, rep, 2L),
    reference = section6_seed(base_seed, job$design_id, rep, 3L)
  )
  tryCatch({
    set.seed(seeds$data)
    x <- generate_section6_sample(job)
    standard <- run_section6_bootstrap(
      job, x, B, seeds$bootstrap, seeds$derivative,
      derivative_mc_size = 1000L, cvm_block_size = cvm_block_size
    )
    theta_hat <- standard$observed$theta_hat
    hybrid <- run_hybrid_bootstrap(
      x, theta_hat, normal_theta0_section6(job), B, seeds$bootstrap,
      seeds$reference, reference_mc_size, cvm_block_size
    )
    rbind(
      diagnostic_row(job, rep, "production_Nderiv1000", standard,
                     proc.time()[["elapsed"]] - started, seeds),
      diagnostic_row(job, rep, "hybrid_Dtheta0_Vtheta0", hybrid,
                     proc.time()[["elapsed"]] - started, seeds)
    )
  }, error = function(e) {
    data.frame(
      scenario = as.character(job$scenario), d = as.integer(job$d), n = as.integer(job$n),
      rep = as.integer(rep), procedure = "pair", status = "error",
      error_message = conditionMessage(e), ks_pvalue = NA_real_, cvm_pvalue = NA_real_,
      ks_reject = NA, cvm_reject = NA, method = NA_character_, backend = NA_character_,
      fused = NA, seed_data = seeds$data, seed_bootstrap = seeds$bootstrap,
      seed_reference = seeds$reference, elapsed_seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    )
  })
}

summarize_hybrid <- function(x) {
  ok <- x[x$status == "ok", , drop = FALSE]
  do.call(rbind, lapply(split(ok, interaction(ok$scenario, ok$procedure, drop = TRUE)), function(z) {
    data.frame(
      scenario = z$scenario[[1L]], procedure = z$procedure[[1L]], M = nrow(z),
      ks_size = mean(z$ks_reject), cvm_size = mean(z$cvm_reject),
      median_ks_pvalue = stats::median(z$ks_pvalue),
      median_cvm_pvalue = stats::median(z$cvm_pvalue),
      all_fast = all(z$method == "fast_multiplier"),
      all_cpp = all(z$backend == "cpp"), all_fused = all(z$fused),
      stringsAsFactors = FALSE
    )
  }))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  M <- as.integer(args$M %||% 50L)
  B <- as.integer(args$B %||% 299L)
  cores <- as.integer(args$cores %||% 10L)
  reference_mc_size <- as.integer(args$reference_mc_size %||% 100000L)
  cvm_block_size <- as.integer(args$cvm_block_size %||% 50L)
  base_seed <- as.integer(args$seed %||% 20260808L)
  output_dir <- as.character(args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios",
    "diagnostic_normal_d10_n100_theta0_hybrid_M50_B299_Nref100k"
  ))
  if (M < 2L || B < 99L || cores < 1L || reference_mc_size < 10000L) stop("Invalid controls.")
  ensure_distance_profile_cpp_loaded()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- section6_acquire_output_lock(output_dir)
  on.exit(section6_release_output_lock(lock), add = TRUE)
  out_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  design <- make_section6_design(
    "normal", dimensions = 10L, n_values = 100L, beta_values = 0
  )
  existing <- if (file.exists(out_path)) utils::read.csv(out_path, stringsAsFactors = FALSE) else data.frame()
  tasks <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    cbind(design[rep(i, M), , drop = FALSE], rep = seq_len(M))
  }))
  done_key <- if (nrow(existing)) with(existing[existing$status == "ok" & existing$procedure == "hybrid_Dtheta0_Vtheta0", ],
    paste(scenario, rep, sep = "|")) else character()
  task_key <- paste(tasks$scenario, tasks$rep, sep = "|")
  pending <- tasks[!task_key %in% done_key, , drop = FALSE]
  writeLines(c(
    "purpose: diagnostic only; not a valid calibration for the paper",
    "comparison: production fast bootstrap Nderiv=1000 versus hybrid D(theta0), V(theta0)",
    "observed statistic and its MLE: retained in both procedures",
    sprintf("reference auxiliary size: %d", reference_mc_size),
    sprintf("M per scenario: %d; B: %d; outer cores: %d", M, B, cores)
  ), file.path(output_dir, "README.txt"))
  total <- nrow(tasks)
  completed <- total - nrow(pending)
  message(sprintf("Paired theta0-hybrid diagnostic: %d/%d completed; %d pending.", completed, total, nrow(pending)))
  while (nrow(pending)) {
    batch <- pending[seq_len(min(cores, nrow(pending))), , drop = FALSE]
    rows <- parallel::mclapply(seq_len(nrow(batch)), function(i) {
      run_one_pair(batch[i, , drop = FALSE], batch$rep[[i]], B, base_seed,
                   reference_mc_size, cvm_block_size)
    }, mc.cores = min(cores, nrow(batch)), mc.set.seed = FALSE)
    existing <- rbind(existing, do.call(rbind, rows))
    section6_write_atomic_csv(existing, out_path)
    section6_write_atomic_csv(summarize_hybrid(existing), summary_path)
    pending <- pending[-seq_len(nrow(batch)), , drop = FALSE]
    completed <- completed + nrow(batch)
    message(sprintf("completed pairs: %d/%d", completed, total))
  }
}
