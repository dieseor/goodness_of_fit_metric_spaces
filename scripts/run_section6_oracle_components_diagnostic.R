#!/usr/bin/env Rscript

# Paired, stepwise diagnosis of the Normal composite fast correction.  These
# hybrids are diagnostic devices only; none except `production` and `simple`
# is a valid calibration procedure for reporting.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/run_section6_theta0_hybrid_diagnostic.R")
source("scripts/run_section6_simple_null_diagnostic.R")

make_component_normal_spec <- function(theta0, score_at = c("hat", "theta0"),
                                       V_at = c("hat", "theta0"),
                                       D_at = c("hat", "theta0")) {
  score_at <- match.arg(score_at)
  V_at <- match.arg(V_at)
  D_at <- match.arg(D_at)
  spec <- make_mvnormal_spec("both")
  spec$fast_multiplier_prepare <- function(data, theta_hat, ks_prep = NULL,
                                            cvm_prep = NULL, control = list()) {
    theta_score <- if (score_at == "hat") theta_hat else theta0
    theta_V <- if (V_at == "hat") theta_hat else theta0
    theta_D <- if (D_at == "hat") theta_hat else theta0
    n_aux <- if (D_at == "hat") as.integer(control$derivative_mc_size) else
      as.integer(control$reference_mc_size)
    seed <- if (D_at == "hat") control$derivative_mc_seed else control$reference_mc_seed
    if (!is.null(seed)) set.seed(as.integer(seed))
    score_obs <- normal_score_both(data, theta_score)
    V <- fast_multiplier_gaussian_paper_vhat(theta_V$Sigma, "both")
    S_obs <- -score_obs %*% t(solve(V))
    aux <- rmvnormal_euclidean(n_aux, theta_D$mu, theta_D$Sigma)
    Psi_aux <- normal_score_both(aux, theta_D)
    D_ks <- fast_multiplier_compute_D_ks(spec, aux, Psi_aux, ks_prep, control)
    D_cvm <- fast_multiplier_reuse_sample_ks_derivative_for_cvm(D_ks, cvm_prep)
    if (is.null(D_cvm)) D_cvm <- fast_multiplier_compute_D_cvm(
      spec, aux, Psi_aux, data, cvm_prep, control
    )
    list(
      S_obs = S_obs, Vhat = diag(ncol(S_obs)), Psi_aux = Psi_aux,
      correction_representation = sprintf("score_%s_V_%s_D_%s", score_at, V_at, D_at),
      paper_Vhat = V, paper_Vhat_inverse = solve(V),
      paper_Vhat_method = paste0("analytic_at_", V_at),
      vhat_method = paste0("diagnostic_", V_at),
      vhat_diagnostics = fast_multiplier_vhat_diagnostics(
        S_obs, Psi_aux, diag(ncol(S_obs)), par0 = theta_D$mu
      ),
      observed_cvm_distance_matrix = NULL,
      derivative_method = paste0("score_mc_at_", D_at), derivative_mc_size = n_aux,
      derivative_mc_seed = seed %||% NA_integer_, D_ks = D_ks, D_cvm = D_cvm
    )
  }
  spec
}

# Diagnostic only: retain the fitted parameter in the fast correction, but
# evaluate the observed distance-profile discrepancy against F^{theta0}.
# This is not a valid composite-null procedure; it isolates the fitted-profile
# term P_n y_{omega,t} - F^{theta_hat}_omega(t).
make_profile0_normal_spec <- function(theta0, theta_correction) {
  base <- make_mvnormal_spec("both")
  spec <- base
  base_profile_eval <- base$profile_eval
  base_profile_matrix_eval <- base$extras$profile_matrix_eval
  base_sample_profile_matrix_eval <- base$extras$sample_profile_matrix_eval

  spec$profile_eval <- function(omega, t, theta, control = list()) {
    base_profile_eval(omega, t, theta0, control)
  }
  spec$extras$profile_matrix_eval <- function(omega_grid, t_grid, theta, control = list()) {
    base_profile_matrix_eval(omega_grid, t_grid, theta0, control)
  }
  spec$extras$sample_profile_matrix_eval <- function(data, distance_matrix, theta, control = list()) {
    base_sample_profile_matrix_eval(data, distance_matrix, theta0, control)
  }
  spec$fast_multiplier_prepare <- function(data, theta_hat, ks_prep = NULL,
                                            cvm_prep = NULL, control = list()) {
    prepare_mvnormal_fast_multiplier(
      spec = base, data = data, theta_hat = theta_correction,
      ks_prep = ks_prep, cvm_prep = cvm_prep, control = control,
      unknown_param = "both"
    )
  }
  spec
}

run_component_bootstrap <- function(x, theta_hat, theta0, component, B,
                                    bootstrap_seed, derivative_seed, reference_seed,
                                    reference_mc_size, cvm_block_size) {
  if (component == "production") {
    job <- data.frame(family = "normal")
    return(run_section6_bootstrap(
      job, x, B, bootstrap_seed, derivative_seed, 1000L, cvm_block_size
    ))
  }
  if (component == "simple") {
    return(multiplier_bootstrap_gof(
      data = x, spec = make_mvnormal_spec("both"),
      null = list(type = "simple", theta = theta0), statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = .05,
      n_cores = 1L, seed = bootstrap_seed, bootstrap_method = "reestimated",
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                  bootstrap_thetas = FALSE),
      control = list(cvm_block_size = cvm_block_size)
    ))
  }
  if (component == "F0") {
    return(multiplier_bootstrap_gof(
      data = x, spec = make_profile0_normal_spec(theta0, theta_hat),
      null = list(type = "composite"), statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = .05,
      n_cores = 1L, seed = bootstrap_seed, observed_theta_hat = theta0,
      bootstrap_method = "fast_multiplier",
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                  bootstrap_thetas = FALSE),
      control = section6_control(1000L, derivative_seed, cvm_block_size)
    ))
  }
  lookup <- list(
    score0 = c("theta0", "hat", "hat"),
    V0 = c("hat", "theta0", "hat"),
    D0 = c("hat", "hat", "theta0"),
    DV0 = c("hat", "theta0", "theta0")
  )
  choice <- lookup[[component]]
  if (is.null(choice)) stop("Unknown diagnostic component.")
  control <- section6_control(1000L, derivative_seed, cvm_block_size)
  control$reference_mc_size <- reference_mc_size
  control$reference_mc_seed <- reference_seed
  multiplier_bootstrap_gof(
    data = x, spec = make_component_normal_spec(theta0, choice[[1]], choice[[2]], choice[[3]]),
    null = list(type = "composite"), statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(), B = B, alpha = .05,
    n_cores = 1L, seed = bootstrap_seed, observed_theta_hat = theta_hat,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                bootstrap_thetas = FALSE), control = control
  )
}

component_row <- function(job, rep, component, fit, seeds, started) {
  data.frame(
    scenario = as.character(job$scenario), d = as.integer(job$d), n = as.integer(job$n),
    rep = rep, component = component, status = "ok", error_message = NA_character_,
    observed_ks = fit$inference$ks$observed, observed_cvm = fit$inference$cvm$observed,
    critical_ks = fit$inference$ks$critical_value, critical_cvm = fit$inference$cvm$critical_value,
    ks_pvalue = fit$inference$ks$p_value, cvm_pvalue = fit$inference$cvm$p_value,
    ks_reject = fit$inference$ks$reject, cvm_reject = fit$inference$cvm$reject,
    effective_method = fit$diagnostics$effective_bootstrap_method,
    seed_data = seeds$data, seed_bootstrap = seeds$bootstrap,
    elapsed_seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE
  )
}

run_component_pair <- function(job, rep, components, B, seed, reference_mc_size, cvm_block_size) {
  started <- proc.time()[["elapsed"]]
  seeds <- list(data = section6_seed(seed, job$design_id, rep, 0L),
                bootstrap = section6_seed(seed, job$design_id, rep, 1L),
                derivative = section6_seed(seed, job$design_id, rep, 2L),
                reference = section6_seed(seed, job$design_id, rep, 3L))
  tryCatch({
    set.seed(seeds$data); x <- generate_section6_sample(job)
    theta_hat <- make_mvnormal_spec("both")$fit_theta(x, NULL, list(type = "composite"), list())
    theta0 <- normal_theta0_section6(job)
    do.call(rbind, lapply(components, function(component) component_row(
      job, rep, component,
      run_component_bootstrap(x, theta_hat, theta0, component, B, seeds$bootstrap,
        seeds$derivative, seeds$reference, reference_mc_size, cvm_block_size),
      seeds, started
    )))
  }, error = function(e) data.frame(
    scenario = as.character(job$scenario), d = as.integer(job$d), n = as.integer(job$n),
    rep = rep, component = "pair", status = "error", error_message = conditionMessage(e),
    observed_ks = NA_real_, observed_cvm = NA_real_, critical_ks = NA_real_, critical_cvm = NA_real_,
    ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA, cvm_reject = NA,
    effective_method = NA_character_, seed_data = seeds$data, seed_bootstrap = seeds$bootstrap,
    elapsed_seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE
  ))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  M <- as.integer(args$M %||% 50L); B <- as.integer(args$B %||% 299L)
  cores <- as.integer(args$cores %||% 12L)
  reference_mc_size <- as.integer(args$reference_mc_size %||% 100000L)
  cvm_block_size <- as.integer(args$cvm_block_size %||% 50L)
  seed <- as.integer(args$seed %||% 20260810L)
  output_dir <- args$output_dir %||% file.path("simulation_results", "section6_new_scenarios", "diagnostic_normal_d10_components_M50_B299")
  available_components <- c("production", "F0", "score0", "V0", "D0", "DV0", "simple")
  components_arg <- args$components %||% paste(available_components, collapse = ",")
  components <- trimws(strsplit(as.character(components_arg), ",", fixed = TRUE)[[1L]])
  if (!length(components) || any(!nzchar(components)) || any(!components %in% available_components)) {
    stop("`--components` must be a non-empty comma-separated subset of: ",
      paste(available_components, collapse = ", "), call. = FALSE)
  }
  components <- unique(components)
  ensure_distance_profile_cpp_loaded(); dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- section6_acquire_output_lock(output_dir); on.exit(section6_release_output_lock(lock), add = TRUE)
  path <- file.path(output_dir, "raw_results.csv")
  existing <- if (file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else data.frame()
  design <- make_section6_design("normal", dimensions = 10L, n_values = 100L, beta_values = 0)
  tasks <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) cbind(design[rep(i, M), , drop = FALSE], rep = seq_len(M))))
  complete <- if (nrow(existing)) {
    ok <- existing[existing$status == "ok" & existing$component %in% components, , drop = FALSE]
    by_task <- split(ok$component, paste(ok$scenario, ok$rep, sep = "|"))
    names(by_task)[vapply(by_task, function(x) all(components %in% x), logical(1L))]
  } else character()
  pending <- tasks[!paste(tasks$scenario, tasks$rep, sep = "|") %in% complete, , drop = FALSE]
  writeLines(c("diagnostic only; components are not reportable calibrations", paste("components:", paste(components, collapse = ", ")), sprintf("M=%d, B=%d, reference_mc_size=%d", M, B, reference_mc_size)), file.path(output_dir, "README.txt"))
  total <- nrow(tasks); message(sprintf("Component diagnostic: %d/%d pairs complete.", total - nrow(pending), total))
  while (nrow(pending)) {
    batch <- pending[seq_len(min(cores, nrow(pending))), , drop = FALSE]
    rows <- parallel::mclapply(seq_len(nrow(batch)), function(i) run_component_pair(batch[i, , drop = FALSE], batch$rep[[i]], components, B, seed, reference_mc_size, cvm_block_size), mc.cores = min(cores, nrow(batch)), mc.set.seed = FALSE)
    existing <- rbind(existing, do.call(rbind, rows)); section6_write_atomic_csv(existing, path)
    pending <- pending[-seq_len(nrow(batch)), , drop = FALSE]
    message(sprintf("completed pairs: %d/%d", total - nrow(pending), total))
  }
  ok <- existing[existing$status == "ok", ]
  summary <- do.call(rbind, lapply(split(ok, interaction(ok$scenario, ok$component, drop = TRUE)), function(z) data.frame(scenario=z$scenario[[1]], component=z$component[[1]], M=nrow(z), ks_size=mean(z$ks_reject), cvm_size=mean(z$cvm_reject), median_observed_ks=median(z$observed_ks), median_critical_ks=median(z$critical_ks), median_observed_cvm=median(z$observed_cvm), median_critical_cvm=median(z$critical_cvm), median_ks_pvalue=median(z$ks_pvalue), median_cvm_pvalue=median(z$cvm_pvalue))))
  section6_write_atomic_csv(summary, file.path(output_dir, "summary.csv"))
}
