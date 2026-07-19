#!/usr/bin/env Rscript

# Reproducible power experiment for the first HvMF alternative, generated with
# the polar sampler.  No GIG sampler is used here.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source("scripts/run_second_scenarios_power.R")

mu0_hvmf_alt1 <- c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2))
mu1_hvmf_alt1 <- c(sqrt(2), -1 / sqrt(2), 1 / sqrt(2))

as_numeric_vector <- function(x, name) {
  value <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (!length(value) || any(!is.finite(value))) stop(sprintf("`%s` must be a non-empty numeric vector.", name))
  value
}

hvmf_alt1_label <- function(kappa) paste0("kappa_", format(kappa, trim = TRUE))

generate_hvmf_alt1_polar_sample <- function(n, beta, kappa) {
  # P_beta = (1 - beta / 2) HvMF(mu0, kappa) + (beta / 2) HvMF(mu1, kappa).
  take_alternative <- stats::runif(n) < beta / 2
  x <- rhvmf_h2_polar(n, mu = mu0_hvmf_alt1, kappa = kappa)
  if (any(take_alternative)) {
    x[take_alternative, ] <- rhvmf_h2_polar(sum(take_alternative), mu = mu1_hvmf_alt1, kappa = kappa)
  }
  x
}

make_hvmf_alt1_design <- function(kappas, n_values, beta_values) {
  rows <- list(); i <- 1L
  for (kappa in kappas) for (n in n_values) for (beta in beta_values) {
    rows[[i]] <- data.frame(
      family = "hvmf_alt1_polar", candidate = hvmf_alt1_label(kappa),
      kappa = kappa, n = as.integer(n), beta = beta, stringsAsFactors = FALSE
    )
    i <- i + 1L
  }
  out <- do.call(rbind, rows)
  out$design_id <- seq_len(nrow(out))
  out
}

empty_hvmf_alt1_results <- function() {
  data.frame(
    family = character(), candidate = character(), kappa = numeric(), n = integer(), beta = numeric(),
    design_id = integer(), rep = integer(), method = character(), ks_pvalue = numeric(), cvm_pvalue = numeric(),
    ks_reject = logical(), cvm_reject = logical(), elapsed_seconds = numeric(), status = character(),
    error_message = character(), stringsAsFactors = FALSE
  )
}

run_one_hvmf_alt1_job <- function(job, B, base_seed) {
  started <- proc.time()[["elapsed"]]
  out <- data.frame(
    family = as.character(job$family), candidate = as.character(job$candidate), kappa = as.numeric(job$kappa),
    n = as.integer(job$n), beta = as.numeric(job$beta), design_id = as.integer(job$design_id), rep = as.integer(job$rep),
    method = "fast_multiplier", ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA,
    cvm_reject = NA, elapsed_seconds = NA_real_, status = "ok", error_message = NA_character_, stringsAsFactors = FALSE
  )
  out <- tryCatch({
    set.seed(power_seed(base_seed, out$design_id, out$rep, 0L))
    x <- generate_hvmf_alt1_polar_sample(out$n, out$beta, out$kappa)
    fit <- run_power_bootstrap(
      "hvmf", x, B = B, seed = power_seed(base_seed, out$design_id, out$rep, 1L),
      bootstrap_method = "fast_multiplier", hvmf_small_grid = FALSE
    )
    out$ks_pvalue <- fit$inference$ks$p_value; out$cvm_pvalue <- fit$inference$cvm$p_value
    out$ks_reject <- fit$inference$ks$reject; out$cvm_reject <- fit$inference$cvm$reject
    out
  }, error = function(error) {
    out$status <- "error"; out$error_message <- conditionMessage(error); out
  })
  out$elapsed_seconds <- proc.time()[["elapsed"]] - started
  out
}

run_hvmf_alt1_jobs <- function(design, M, B, cores, output_dir, stage, base_seed = 20260716L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "raw_results.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  manifest_path <- file.path(output_dir, "manifest.csv")
  log_path <- file.path(output_dir, "run.log")
  if (!file.exists(manifest_path)) {
    manifest <- transform(
      design, M = as.integer(M), B = as.integer(B), cores = as.integer(cores), stage = stage,
      bootstrap_method = "fast_multiplier", base_seed = as.integer(base_seed), hvmf_small_grid = FALSE,
      null_mu = paste(format(mu0_hvmf_alt1, digits = 17L), collapse = ";"),
      alternative_mu = paste(format(mu1_hvmf_alt1, digits = 17L), collapse = ";"),
      mixture_weight = "beta/2", sampler = "rhvmf_h2_polar"
    )
    utils::write.csv(manifest, manifest_path, row.names = FALSE)
  }
  cat(sprintf("%s stage=%s M=%d B=%d cores=%d sampler=polar\n", format(Sys.time(), tz = "Europe/Madrid"), stage, M, B, cores), file = log_path, append = TRUE)
  existing <- if (file.exists(result_path)) utils::read.csv(result_path, stringsAsFactors = FALSE) else empty_hvmf_alt1_results()
  jobs <- merge(design, data.frame(rep = seq_len(M)), by = NULL)
  key <- function(x) paste(x$design_id, x$rep, "fast_multiplier", sep = "_")
  done <- if (nrow(existing)) key(existing[existing$status == "ok" & existing$method == "fast_multiplier", , drop = FALSE]) else character()
  pending <- jobs[!key(jobs) %in% done, , drop = FALSE]
  started <- Sys.time(); total <- nrow(jobs); completed <- total - nrow(pending)
  write_power_status(status_path, completed, total, started, existing, stage, cores, design)
  print_power_progress(completed, total, started, existing, cores)
  if (!nrow(pending)) return(existing)
  chunks <- split(pending, ceiling(seq_len(nrow(pending)) / max(1L, 2L * as.integer(cores))))
  for (chunk in chunks) {
    rows <- if (.Platform$OS.type == "unix" && cores > 1L) {
      parallel::mclapply(seq_len(nrow(chunk)), function(i) run_one_hvmf_alt1_job(chunk[i, , drop = FALSE], B, base_seed), mc.cores = as.integer(cores), mc.preschedule = FALSE)
    } else lapply(seq_len(nrow(chunk)), function(i) run_one_hvmf_alt1_job(chunk[i, , drop = FALSE], B, base_seed))
    existing <- rbind(existing, do.call(rbind, rows))
    write_atomic_csv(existing, result_path)
    completed <- completed + nrow(chunk)
    write_power_status(status_path, completed, total, started, existing, stage, cores, design)
    print_power_progress(completed, total, started, existing, cores)
    cat(sprintf("%s completed=%d/%d\n", format(Sys.time(), tz = "Europe/Madrid"), completed, total), file = log_path, append = TRUE)
  }
  existing
}

summarize_hvmf_alt1 <- function(results) {
  good <- results[results$status == "ok" & results$method == "fast_multiplier", , drop = FALSE]
  do.call(rbind, lapply(split(good, interaction(good$candidate, good$n, good$beta, drop = TRUE)), function(x) {
    data.frame(candidate = x$candidate[[1L]], kappa = x$kappa[[1L]], n = x$n[[1L]], beta = x$beta[[1L]], M = nrow(x),
               power_ks = mean(x$ks_reject), power_cvm = mean(x$cvm_reject), stringsAsFactors = FALSE)
  }))
}

if (identical(Sys.getenv("RUN_HVMF_ALT1_POLAR_POWER"), "1")) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  M <- as.integer(cli$M %||% 100L); B <- as.integer(cli$B %||% 299L); cores <- as.integer(cli$cores %||% 10L)
  if (M < 1L || B < 1L || cores < 1L) stop("`M`, `B`, and `cores` must be positive integers.")
  kappas <- as_numeric_vector(cli$kappas %||% "1,2", "kappas")
  n_values <- as_numeric_vector(cli$n %||% "200", "n")
  beta_values <- as_numeric_vector(cli$beta %||% "1", "beta")
  if (any(n_values < 1 | n_values != as.integer(n_values)) || any(beta_values < 0 | beta_values > 1)) stop("Invalid `n` or `beta`.")
  output_dir <- cli$output_dir %||% file.path("simulation_results", "second_scenarios_power", "hvmf_alt1_polar_screen_M100_B299")
  design <- make_hvmf_alt1_design(kappas, as.integer(n_values), beta_values)
  results <- run_hvmf_alt1_jobs(design, M = M, B = B, cores = cores, output_dir = output_dir, stage = cli$stage %||% "screen", base_seed = as.integer(cli$seed %||% 20260716L))
  utils::write.csv(summarize_hvmf_alt1(results), file.path(output_dir, "summary.csv"), row.names = FALSE)
}
