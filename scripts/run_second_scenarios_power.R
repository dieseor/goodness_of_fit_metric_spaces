#!/usr/bin/env Rscript

# Reproducible pilot and final power runs for the four second scenarios.
# Outer parallelism is the only parallelism used here; every bootstrap call is
# deliberately single-core to respect the global core budget.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source("bootstrap/calibration_study.R")
source("bootstrap/cardioid_model_spec.R")

power_root <- file.path("simulation_results", "second_scenarios_power")
power_mu_hvmf <- c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2))

parse_cli <- function(args) {
  parsed <- list()
  for (arg in args[startsWith(args, "--")]) {
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    parsed[[pieces[[1L]]]] <- if (length(pieces) == 1L) "TRUE" else paste(pieces[-1L], collapse = "=")
  }
  parsed
}

power_seed <- function(base_seed, design_id, rep_id, stream = 0L) {
  as.integer((as.numeric(base_seed) + 1000003 * design_id + 1009 * rep_id + 10000019 * stream) %% 2147483647) + 1L
}

r_standardized_t <- function(n, nu) {
  stats::rt(n, df = nu) * sqrt((nu - 2) / nu)
}

r_standardized_t_ilr <- function(n, nu) {
  g <- matrix(stats::rnorm(n * 2L), ncol = 2L)
  g * sqrt((nu - 2) / stats::rchisq(n, df = nu))
}

r_lg_t <- function(n, nu) {
  logistic_gaussian_ilr_to_simplex(r_standardized_t_ilr(n, nu), ambient_dim = 3L)
}

generate_power_sample <- function(family, candidate, n, beta) {
  if (identical(family, "normal")) {
    return(if (beta == 0) stats::rnorm(n) else {
      take_h1 <- stats::runif(n) < beta
      x <- stats::rnorm(n)
      x[take_h1] <- r_standardized_t(sum(take_h1), candidate)
      x
    })
  }
  if (identical(family, "lg")) {
    take_h1 <- stats::runif(n) < beta
    z <- matrix(stats::rnorm(n * 2L), ncol = 2L)
    if (any(take_h1)) z[take_h1, ] <- r_standardized_t_ilr(sum(take_h1), candidate)
    return(logistic_gaussian_ilr_to_simplex(z, ambient_dim = 3L))
  }
  if (identical(family, "vmf")) {
    take_h1 <- stats::runif(n) < beta
    x <- rotasym::r_vMF(n, mu = c(1, 0, 0), kappa = 0.5)
    if (any(take_h1)) x[take_h1, ] <- r_sph_car(sum(take_h1), mu = c(1, 0, 0), rho = candidate, k = 1L)
    return(x)
  }
  if (identical(family, "hvmf")) {
    take_h1 <- stats::runif(n) < beta
    x <- rhvmf_h2_polar(n, mu = power_mu_hvmf, kappa = candidate$kappa)
    if (any(take_h1)) {
      x[take_h1, ] <- rhvmf_h2_angular_mixture(
        sum(take_h1), mu = power_mu_hvmf, kappa = candidate$kappa, delta = candidate$delta
      )
    }
    return(x)
  }
  stop(sprintf("Unsupported family '%s'.", family))
}

run_power_bootstrap <- function(family, x, B, seed, bootstrap_method, hvmf_small_grid = FALSE) {
  ks_grid <- if (identical(family, "hvmf") && isTRUE(hvmf_small_grid)) {
    make_hvmf_ks_grid(x, mu = power_mu_hvmf, n_omega = 4L, n_t = 4L)
  } else {
    make_sample_unique_distance_ks_grid()
  }
  common <- list(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = as.integer(B), alpha = 0.05, n_cores = 1L, seed = seed,
    bootstrap_method = bootstrap_method,
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE)
  )
  if (identical(family, "normal")) {
    return(do.call(multiplier_bootstrap_normal, c(common, list(
      unknown_param = "both", control = list(derivative_mc_size = 1000L, derivative_mc_seed = seed + 7L)
    ))))
  }
  if (identical(family, "lg")) {
    return(do.call(multiplier_bootstrap_logistic_gaussian, c(common, list(
      unknown_param = "both", control = list(derivative_method = "quadrature")
    ))))
  }
  if (identical(family, "vmf")) {
    return(do.call(multiplier_bootstrap_vmf, c(common, list(
      unknown_param = "xi", distance_type = "geodesic",
      control = list(
        derivative_method = "quadrature",
        vmf_profile_method = "tabulated",
        vmf_profile_n_u = 4097L
      )
    ))))
  }
  if (identical(family, "hvmf")) {
    return(do.call(multiplier_bootstrap_hvmf, c(common, list(
      unknown_param = "both", control = list(
        derivative_method = "quadrature",
        hvmf_profile_method = "tabulated",
        hvmf_profile_n_y = if (isTRUE(hvmf_small_grid)) 257L else 4097L
      )
    ))))
  }
  stop(sprintf("Unsupported family '%s'.", family))
}

candidate_label <- function(family, candidate) {
  if (identical(family, "hvmf")) sprintf("kappa_%s_delta_%s", candidate$kappa, format(candidate$delta, trim = TRUE))
  else paste0(if (identical(family, "lg")) "nu_" else if (identical(family, "vmf")) "rho_" else "nu_", format(candidate, trim = TRUE))
}

candidate_from_label <- function(family, label) {
  values <- as.numeric(unlist(regmatches(label, gregexpr("[0-9]+(?:\\.[0-9]+)?", label, perl = TRUE))))
  if (identical(family, "hvmf")) list(kappa = values[[1L]], delta = values[[2L]]) else values[[1L]]
}

candidate_catalog <- function() {
  hvmf_candidates <- unlist(
    lapply(c(5, 25, 50), function(kappa) {
      lapply(c(0.05, 0.10, 0.20), function(delta) list(kappa = kappa, delta = delta))
    }),
    recursive = FALSE
  )
  list(
    normal = as.list(c(3, 4, 6, 10, 20)),
    lg = as.list(c(3, 4, 6, 10, 20)),
    vmf = as.list(c(0.1, 0.3, 0.5)),
    hvmf = unname(hvmf_candidates)
  )
}

make_design <- function(catalog, n_values, beta_values) {
  rows <- list(); index <- 1L
  for (family in names(catalog)) for (candidate in catalog[[family]]) for (n in n_values) for (beta in beta_values) {
    rows[[index]] <- data.frame(
      family = family, candidate = candidate_label(family, candidate), n = as.integer(n), beta = as.numeric(beta), stringsAsFactors = FALSE
    ); index <- index + 1L
  }
  design <- do.call(rbind, rows); design$design_id <- seq_len(nrow(design)); design
}

empty_power_results <- function() {
  data.frame(family = character(), candidate = character(), n = integer(), beta = numeric(), design_id = integer(), rep = integer(),
             method = character(), ks_pvalue = numeric(), cvm_pvalue = numeric(), ks_reject = logical(), cvm_reject = logical(),
             derivative_method_requested = character(), derivative_method_effective = character(),
             elapsed_seconds = numeric(), status = character(), error_message = character(), stringsAsFactors = FALSE)
}

run_one_power_job <- function(job, B, base_seed, bootstrap_method, hvmf_small_grid = FALSE) {
  family <- as.character(job$family); candidate <- candidate_from_label(family, as.character(job$candidate))
  started <- proc.time()[["elapsed"]]
  out <- data.frame(family = family, candidate = as.character(job$candidate), n = as.integer(job$n), beta = as.numeric(job$beta),
                    design_id = as.integer(job$design_id), rep = as.integer(job$rep), method = bootstrap_method,
                    ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA, cvm_reject = NA,
                    derivative_method_requested = if (family %in% c("vmf", "hvmf")) "quadrature" else "score_mc",
                    derivative_method_effective = NA_character_,
                    elapsed_seconds = NA_real_, status = "ok", error_message = NA_character_, stringsAsFactors = FALSE)
  out <- tryCatch({
    set.seed(power_seed(base_seed, out$design_id, out$rep, 0L))
    x <- generate_power_sample(family, candidate, out$n, out$beta)
    fit <- run_power_bootstrap(family, x, B = B, seed = power_seed(base_seed, out$design_id, out$rep, 1L), bootstrap_method = bootstrap_method, hvmf_small_grid = hvmf_small_grid)
    out$ks_pvalue <- fit$inference$ks$p_value; out$cvm_pvalue <- fit$inference$cvm$p_value
    out$ks_reject <- fit$inference$ks$reject; out$cvm_reject <- fit$inference$cvm$reject
    out$derivative_method_effective <- fit$diagnostics$derivative_method_effective %||%
      fit$diagnostics$derivative_method %||% NA_character_
    out
  }, error = function(error) {
    out$status <- "error"; out$error_message <- conditionMessage(error); out
  })
  out$elapsed_seconds <- proc.time()[["elapsed"]] - started
  out
}

write_atomic_csv <- function(data, path) {
  temporary <- paste0(path, ".tmp")
  utils::write.csv(data, temporary, row.names = FALSE)
  file.rename(temporary, path)
}

write_power_status <- function(path, completed, total, started, results, stage, cores, design) {
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  mean_seconds <- if (nrow(results) == 0L) NA_real_ else mean(results$elapsed_seconds[results$status == "ok"], na.rm = TRUE)
  eta <- if (is.finite(mean_seconds)) mean_seconds * (total - completed) / max(1L, as.integer(cores)) else NA_real_
  writeLines(c(
    paste("stage:", stage), paste("updated:", format(Sys.time(), tz = "Europe/Madrid")),
    paste("families:", paste(unique(design$family), collapse = ", ")),
    paste("candidates:", paste(unique(design$candidate), collapse = ", ")),
    sprintf("completed: %d/%d", completed, total), sprintf("elapsed_seconds: %.1f", elapsed),
    sprintf("mean_seconds_per_job: %.3f", mean_seconds), sprintf("eta_seconds: %.1f", eta)
  ), path)
}

print_power_progress <- function(completed, total, started, results, cores, width = 30L) {
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  proportion <- if (total > 0L) completed / total else 1
  filled <- as.integer(floor(width * proportion))
  bar <- paste0(strrep("=", filled), strrep("-", max(0L, width - filled)))
  mean_seconds <- if (nrow(results) == 0L) NA_real_ else {
    mean(results$elapsed_seconds[results$status == "ok"], na.rm = TRUE)
  }
  eta <- if (is.finite(mean_seconds)) {
    mean_seconds * (total - completed) / max(1L, as.integer(cores))
  } else {
    NA_real_
  }
  eta_text <- if (is.finite(eta)) sprintf("ETA %s", format(round(eta), units = "secs")) else "ETA --"
  cat(sprintf("\r[%s] %6.2f%%  %d/%d  elapsed %s  %s",
              bar, 100 * proportion, completed, total,
              format(round(elapsed), units = "secs"), eta_text))
  if (completed >= total) cat("\n")
  flush.console()
}

run_power_jobs <- function(design, M, B, cores, output_dir, stage, bootstrap_method = "fast_multiplier", base_seed = 20260714L, hvmf_small_grid = FALSE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "raw_results.csv")
  status_path <- file.path(output_dir, "progress_status.txt")
  manifest_path <- file.path(output_dir, "manifest.csv")
  log_path <- file.path(output_dir, "run.log")
  expected_derivative_method <- ifelse(
    design$family %in% c("vmf", "hvmf"),
    "quadrature",
    "score_mc"
  )
  if (file.exists(manifest_path)) {
    manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    method_matches <- "derivative_method" %in% names(manifest) &&
      nrow(manifest) == nrow(design) &&
      all(as.character(manifest$derivative_method) == expected_derivative_method)
    if (!isTRUE(method_matches)) {
      stop(
        paste(
          "The existing second-scenarios manifest does not explicitly match",
          "the requested derivative methods. Use a new output directory;",
          "score-MC and quadrature replications must not be mixed silently."
        ),
        call. = FALSE
      )
    }
  }
  if (!file.exists(manifest_path)) {
    manifest <- transform(design, M = as.integer(M), B = as.integer(B), cores = as.integer(cores),
                          stage = stage, bootstrap_method = bootstrap_method, base_seed = as.integer(base_seed),
                          derivative_method = expected_derivative_method,
                          hvmf_small_grid = isTRUE(hvmf_small_grid))
    utils::write.csv(manifest, manifest_path, row.names = FALSE)
  }
  cat(sprintf("%s stage=%s M=%d B=%d cores=%d\\n", format(Sys.time(), tz = "Europe/Madrid"), stage, M, B, cores),
      file = log_path, append = TRUE)
  existing <- if (file.exists(result_path)) utils::read.csv(result_path, stringsAsFactors = FALSE) else empty_power_results()
  jobs <- merge(design, data.frame(rep = seq_len(M)), by = NULL)
  key <- function(x) paste(x$design_id, x$rep, bootstrap_method, sep = "_")
  done <- if (nrow(existing)) key(existing[existing$status == "ok" & existing$method == bootstrap_method, , drop = FALSE]) else character()
  pending <- jobs[!key(jobs) %in% done, , drop = FALSE]
  started <- Sys.time(); total <- nrow(jobs); completed <- total - nrow(pending)
  write_power_status(status_path, completed, total, started, existing, stage, cores, design)
  print_power_progress(completed, total, started, existing, cores)
  if (nrow(pending) == 0L) return(existing)
  chunk_size <- max(1L, 2L * as.integer(cores))
  chunks <- split(pending, ceiling(seq_len(nrow(pending)) / chunk_size))
  for (chunk in chunks) {
    rows <- if (.Platform$OS.type == "unix" && cores > 1L) {
      parallel::mclapply(
        seq_len(nrow(chunk)),
        function(i) run_one_power_job(chunk[i, , drop = FALSE], B, base_seed, bootstrap_method, hvmf_small_grid),
        mc.cores = as.integer(cores), mc.preschedule = FALSE
      )
    } else lapply(seq_len(nrow(chunk)), function(i) run_one_power_job(chunk[i, , drop = FALSE], B, base_seed, bootstrap_method, hvmf_small_grid))
    existing <- rbind(existing, do.call(rbind, rows))
    write_atomic_csv(existing, result_path)
    completed <- completed + nrow(chunk)
    write_power_status(status_path, completed, total, started, existing, stage, cores, design)
    print_power_progress(completed, total, started, existing, cores)
    cat(sprintf("%s completed=%d/%d\\n", format(Sys.time(), tz = "Europe/Madrid"), completed, total),
        file = log_path, append = TRUE)
  }
  existing
}

summarize_power <- function(results) {
  good <- results[results$status == "ok" & results$method == "fast_multiplier", , drop = FALSE]
  do.call(rbind, lapply(split(good, interaction(good$family, good$candidate, good$n, good$beta, drop = TRUE)), function(x) {
    data.frame(family = x$family[[1L]], candidate = x$candidate[[1L]], n = x$n[[1L]], beta = x$beta[[1L]], M = nrow(x),
               power_ks = mean(x$ks_reject), power_cvm = mean(x$cvm_reject), stringsAsFactors = FALSE)
  }))
}

select_candidates <- function(summary_df) {
  screening <- summary_df[summary_df$n == 100L & abs(summary_df$beta - 0.5) < 1e-12, , drop = FALSE]
  selected <- lapply(split(screening, screening$family), function(x) {
    x$mean_power <- rowMeans(x[, c("power_ks", "power_cvm")])
    x$score <- abs(x$mean_power - 0.5) + ifelse(x$mean_power < 0.35 | x$mean_power > 0.65, 1, 0)
    x[which.min(x$score), c("family", "candidate"), drop = FALSE]
  })
  do.call(rbind, selected)
}

run_fast_vs_reestimated_validation <- function(cores = 10L) {
  output_dir <- file.path(power_root, "hvmf_fast_validation_M100_B299")
  catalog <- list(hvmf = unname(lapply(c(5, 25, 50), function(kappa) list(kappa = kappa, delta = 0))))
  design <- make_design(catalog, n_values = c(50L, 100L), beta_values = 0)
  fast <- run_power_jobs(design, M = 100L, B = 299L, cores = cores, output_dir = output_dir, stage = "fast_validation", bootstrap_method = "fast_multiplier", hvmf_small_grid = TRUE)
  slow <- run_power_jobs(design, M = 100L, B = 299L, cores = cores, output_dir = output_dir, stage = "reestimated_validation", bootstrap_method = "reestimated", hvmf_small_grid = TRUE)
  all_results <- rbind(fast, slow[slow$method == "reestimated", , drop = FALSE])
  split_key <- interaction(all_results$candidate, all_results$n, all_results$rep, drop = TRUE)
  paired <- lapply(split(all_results[all_results$status == "ok", , drop = FALSE], split_key), function(x) {
    if (!all(c("fast_multiplier", "reestimated") %in% x$method)) return(NULL)
    fast_row <- x[x$method == "fast_multiplier", , drop = FALSE][1L, ]; slow_row <- x[x$method == "reestimated", , drop = FALSE][1L, ]
    data.frame(candidate = fast_row$candidate, n = fast_row$n, ks_diff = as.numeric(fast_row$ks_reject) - as.numeric(slow_row$ks_reject), cvm_diff = as.numeric(fast_row$cvm_reject) - as.numeric(slow_row$cvm_reject))
  })
  paired <- do.call(rbind, Filter(Negate(is.null), paired))
  validation <- do.call(rbind, lapply(split(paired, interaction(paired$candidate, paired$n, drop = TRUE)), function(x) {
    data.frame(candidate = x$candidate[[1L]], n = x$n[[1L]], M = nrow(x),
               ks_mean_difference = mean(x$ks_diff), ks_mc_se = stats::sd(x$ks_diff) / sqrt(nrow(x)),
               cvm_mean_difference = mean(x$cvm_diff), cvm_mc_se = stats::sd(x$cvm_diff) / sqrt(nrow(x)))
  }))
  validation$passes <- abs(validation$ks_mean_difference) <= pmax(0.02, 3 * validation$ks_mc_se) &
    abs(validation$cvm_mean_difference) <= pmax(0.02, 3 * validation$cvm_mc_se)
  utils::write.csv(validation, file.path(output_dir, "fast_vs_reestimated_summary.csv"), row.names = FALSE)
  if (!all(validation$passes)) stop("HvMF fast-versus-reestimated validation did not agree within Monte Carlo uncertainty.")
  invisible(validation)
}

run_pilot_and_finals <- function() {
  run_fast_vs_reestimated_validation(cores = 10L)
  screening_dir <- file.path(power_root, "pilot_screen_M100_B299")
  catalog <- candidate_catalog()
  screening_design <- make_design(catalog, n_values = 100L, beta_values = 0.5)
  screening_results <- run_power_jobs(screening_design, M = 100L, B = 299L, cores = 10L, output_dir = screening_dir, stage = "pilot_screen")
  screening_summary <- summarize_power(screening_results)
  utils::write.csv(screening_summary, file.path(screening_dir, "summary.csv"), row.names = FALSE)
  selected <- select_candidates(screening_summary)
  utils::write.csv(selected, file.path(screening_dir, "selected_candidates.csv"), row.names = FALSE)

  selected_catalog <- lapply(seq_len(nrow(selected)), function(i) {
    family <- selected$family[[i]]; candidate <- candidate_from_label(family, selected$candidate[[i]])
    setNames(list(list(candidate)), family)
  })
  selected_catalog <- do.call(c, selected_catalog)
  curve_dir <- file.path(power_root, "pilot_curves_M100_B299")
  curve_design <- make_design(selected_catalog, n_values = c(50L, 100L, 200L), beta_values = c(0.25, 0.5, 1))
  curve_results <- run_power_jobs(curve_design, M = 100L, B = 299L, cores = 10L, output_dir = curve_dir, stage = "pilot_curves")
  utils::write.csv(summarize_power(curve_results), file.path(curve_dir, "summary.csv"), row.names = FALSE)

  for (family in c("normal", "lg", "vmf", "hvmf")) {
    final_catalog <- setNames(list(selected_catalog[[family]]), family)
    final_dir <- file.path(power_root, "final_M1000_B5000", family)
    final_design <- make_design(final_catalog, n_values = c(50L, 100L, 200L), beta_values = c(0.25, 0.5, 1))
    final_results <- run_power_jobs(final_design, M = 1000L, B = 5000L, cores = 5L, output_dir = final_dir, stage = paste0("final_", family))
    utils::write.csv(summarize_power(final_results), file.path(final_dir, "summary.csv"), row.names = FALSE)
  }
}

if (identical(Sys.getenv("RUN_SECOND_SCENARIOS_POWER"), "1")) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  mode <- cli$mode %||% "all"
  if (!mode %in% c("validation", "all")) stop("`--mode` must be `validation` or `all`.")
  if (identical(mode, "validation")) run_fast_vs_reestimated_validation(cores = 10L) else run_pilot_and_finals()
}
