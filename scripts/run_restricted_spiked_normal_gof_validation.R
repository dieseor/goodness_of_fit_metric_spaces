#!/usr/bin/env Rscript

# Resumable exploratory GOF validation for the restricted mean-aligned
# single-spiked Gaussian model.  It deliberately studies the null only: no
# alternative is defined here, so this runner cannot silently impose one.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  value <- args[startsWith(args, key)]
  if (!length(value)) return(default)
  sub(key, "", value[[1L]], fixed = TRUE)
}
parse_csv_numeric <- function(value, name, integer = FALSE) {
  out <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (!length(out) || any(!is.finite(out))) stop("Invalid `--", name, "`.")
  if (integer && any(out != as.integer(out))) stop("`--", name, "` must contain integers.")
  if (integer) as.integer(out) else out
}

M <- as.integer(arg_value("M", "100"))
B <- as.integer(arg_value("B", "499"))
cores <- as.integer(arg_value("cores", "4"))
n_values <- parse_csv_numeric(arg_value("n", "50,100,200"), "n", integer = TRUE)
lambda_values <- parse_csv_numeric(arg_value("lambdas", "0.5,2"), "lambdas")
beta_values <- parse_csv_numeric(arg_value("betas", "0,0.5,1"), "betas")
comparison_M <- as.integer(arg_value("comparison_M", "10"))
comparison_B <- as.integer(arg_value("comparison_B", "99"))
derivative_mc_size <- as.integer(arg_value("derivative_mc_size", "1000"))
cvm_block_size <- as.integer(arg_value("cvm_block_size", "50"))
axis_norm <- as.numeric(arg_value("axis_norm", "0.75"))
diagonal_norm <- as.numeric(arg_value("diagonal_norm", "1.50"))
seed <- as.integer(arg_value("seed", "20260830"))
output_dir <- arg_value(
  "output_dir",
  "simulation_results/restricted_spiked_normal_gof_validation"
)
if (any(!is.finite(c(M, B, cores, n_values, lambda_values, comparison_M,
                     comparison_B, derivative_mc_size, cvm_block_size, axis_norm,
                     diagonal_norm, seed))) ||
    M < 1L || B < 1L || cores < 1L || comparison_M < 1L || comparison_B < 1L ||
    derivative_mc_size < 1L || cvm_block_size < 1L || axis_norm <= 0 || diagonal_norm <= 0 ||
    any(n_values < 2L) || any(lambda_values <= 0) || any(beta_values < 0 | beta_values > 1)) {
  stop("M, B, cores, and sample sizes must be positive; lambda values must be strictly positive and beta values must lie in [0,1].")
}
if (any(abs(lambda_values - 0.1) <= 100 * .Machine$double.eps)) {
  stop("lambda = 0.1 is excluded from this GOF validation because its MLE reaches the forbidden boundary too frequently under the null.")
}

source(file.path("bootstrap", "restricted_spiked_normal_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theta_design <- function(d) {
  list(
    axis_moderate = c(axis_norm, rep(0, d - 1L)),
    diagonal_large = rep(diagonal_norm / sqrt(d), d)
  )
}

make_design <- function(M) {
  pieces <- lapply(c(2L, 5L), function(d) {
    theta_names <- names(theta_design(d))
    expand.grid(
      d = d,
      theta_name = theta_names,
      n = n_values,
      lambda = lambda_values,
      beta = beta_values,
      replication = seq_len(M),
      stringsAsFactors = FALSE
    )
  })
  design <- do.call(rbind, pieces)
  design$theta <- vapply(seq_len(nrow(design)), function(i) {
    paste(format(theta_design(design$d[[i]])[[design$theta_name[[i]]]],
                 scientific = FALSE, trim = TRUE), collapse = ";")
  }, character(1L))
  design$theta_norm <- vapply(seq_len(nrow(design)), function(i) {
    sqrt(sum(theta_design(design$d[[i]])[[design$theta_name[[i]]]]^2))
  }, numeric(1L))
  design$job_id <- seq_len(nrow(design))
  design
}

design <- make_design(M)
utils::write.csv(design, file.path(output_dir, "design.csv"), row.names = FALSE)

call_fast <- function(x, job_seed, B) {
  multiplier_bootstrap_restricted_spiked_normal(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = B,
    alpha = 0.05,
    n_cores = 1L,
    seed = job_seed + 1L,
    bootstrap_method = "fast_multiplier",
    keep = list(observed_process = FALSE, bootstrap_statistics = FALSE,
                bootstrap_thetas = FALSE),
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = derivative_mc_size,
      derivative_mc_seed = job_seed + 2L,
      fast_multiplier_cvm_block_size = cvm_block_size,
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_multiplier_cache_corrections = "auto",
      progress_bar = FALSE
    ),
    distance_profile_backend = "r",
    fast_multiplier_backend = "cpp",
    fast_multiplier_cpp_kernel = "contiguous_double",
    fuse_ks_cvm = TRUE,
    cache_block_corrections = "auto"
  )
}

run_job <- function(row) {
  started <- proc.time()[["elapsed"]]
  job_seed <- seed + 100003L * as.integer(row$job_id)
  set.seed(job_seed)
  d <- as.integer(row$d)
  theta_true <- theta_design(d)[[as.character(row$theta_name)]]
  lambda_true <- as.numeric(row$lambda)
  base <- data.frame(
    job_id = as.integer(row$job_id), d = d, theta_name = as.character(row$theta_name),
    theta_true = paste(theta_true, collapse = ";"), theta_norm_true = sqrt(sum(theta_true^2)),
    lambda_true = lambda_true, n = as.integer(row$n),
    beta = as.numeric(row$beta),
    data_generating_distribution = if (as.numeric(row$beta) == 0) {
      "restricted_spiked_normal_null"
    } else {
      "(1-beta/2)N(theta,Sigma)+beta/2 N(-theta,Sigma)"
    },
    replication = as.integer(row$replication), seed = job_seed,
    stringsAsFactors = FALSE
  )
  tryCatch({
    n_minus <- stats::rbinom(1L, size = base$n, prob = base$beta / 2)
    x <- rrestricted_spiked_normal(base$n, theta_true, lambda_true)
    if (n_minus > 0L) {
      x[seq_len(n_minus), ] <- rrestricted_spiked_normal(n_minus, -theta_true, lambda_true)
    }
    warnings <- character()
    result <- withCallingHandlers(
      call_fast(x, job_seed, B),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
    fitted <- result$observed$theta_hat
    score <- restricted_spiked_normal_score_matrix(x, fitted)
    eig <- eigen(fitted$Sigma, symmetric = TRUE, only.values = TRUE)$values
    expected_eig <- c(1 + fitted$lambda, rep(1, d - 1L))
    cbind(base, data.frame(
      status = "ok", error_message = NA_character_,
      warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      theta_hat = paste(fitted$theta, collapse = ";"),
      theta_error = sqrt(sum((fitted$theta - theta_true)^2)),
      lambda_hat = fitted$lambda, lambda_error = fitted$lambda - lambda_true,
      sigma_eigen_max_error = max(abs(sort(eig) - sort(expected_eig))),
      score_mean_norm = sqrt(sum(colMeans(score)^2)),
      score_mean_max_abs = max(abs(colMeans(score))),
      profile_tau = fitted$fit_diagnostics$profile_tau,
      profile_radius = fitted$fit_diagnostics$profiled_radius,
      ks_statistic = result$inference$ks$observed,
      cvm_statistic = result$inference$cvm$observed,
      ks_pvalue = result$inference$ks$p_value,
      cvm_pvalue = result$inference$cvm$p_value,
      effective_bootstrap_method = result$diagnostics$effective_bootstrap_method,
      distance_profile_backend = result$diagnostics$distance_profile_backend_effective,
      fast_backend = result$diagnostics$fast_multiplier_backend_effective,
      fast_kernel = result$diagnostics$fast_multiplier_cpp_kernel_effective,
      fast_fused = result$diagnostics$fast_multiplier_fuse_ks_cvm_effective,
      derivative_method = result$diagnostics$derivative_method_effective,
      vhat_method = result$diagnostics$vhat_method,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    ))
  }, error = function(error) {
    cbind(base, data.frame(
      status = "error", error_message = conditionMessage(error), warning_message = NA_character_,
      theta_hat = NA_character_, theta_error = NA_real_, lambda_hat = NA_real_,
      lambda_error = NA_real_, sigma_eigen_max_error = NA_real_,
      score_mean_norm = NA_real_, score_mean_max_abs = NA_real_,
      profile_tau = NA_real_, profile_radius = NA_real_, ks_statistic = NA_real_,
      cvm_statistic = NA_real_, ks_pvalue = NA_real_, cvm_pvalue = NA_real_,
      effective_bootstrap_method = NA_character_, fast_backend = NA_character_,
      distance_profile_backend = NA_character_, fast_kernel = NA_character_, fast_fused = NA, derivative_method = NA_character_,
      vhat_method = NA_character_, elapsed_seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    ))
  })
}

results_path <- file.path(output_dir, "gof_results.csv")
existing <- if (file.exists(results_path)) utils::read.csv(results_path, stringsAsFactors = FALSE) else data.frame()
done <- if (nrow(existing)) as.integer(existing$job_id) else integer()
pending <- design[!design$job_id %in% done, , drop = FALSE]
message(sprintf("Restricted-spiked GOF validation: %d pending / %d total jobs; M=%d, B=%d, outer cores=%d.",
                nrow(pending), nrow(design), M, B, cores))

started_all <- proc.time()[["elapsed"]]
if (nrow(pending)) {
  for (begin in seq.int(1L, nrow(pending), by = cores)) {
    end <- min(nrow(pending), begin + cores - 1L)
    batch <- pending[begin:end, , drop = FALSE]
    batch_results <- parallel::mclapply(
      seq_len(nrow(batch)), function(i) run_job(batch[i, , drop = FALSE]),
      mc.cores = min(cores, nrow(batch)), mc.preschedule = FALSE
    )
    existing <- rbind(existing, do.call(rbind, batch_results))
    utils::write.csv(existing, results_path, row.names = FALSE)
    elapsed <- proc.time()[["elapsed"]] - started_all
    completed <- nrow(existing)
    rate <- elapsed / max(1L, end)
    eta <- rate * (nrow(pending) - end)
    message(sprintf("completed %d/%d GOF jobs (this run %d/%d); elapsed %.1fs; ETA %.1fs",
                    completed, nrow(design), end, nrow(pending), elapsed, eta))
  }
}

summarize_null <- function(group) {
  successful <- group[group$status == "ok", , drop = FALSE]
  data.frame(
    d = group$d[[1L]], theta_name = group$theta_name[[1L]],
    theta_norm_true = group$theta_norm_true[[1L]], lambda_true = group$lambda_true[[1L]],
    beta = group$beta[[1L]], n = group$n[[1L]], attempted = nrow(group), successful = nrow(successful),
    failures = sum(group$status != "ok"),
    ks_rejection_rate = if (nrow(successful)) mean(successful$ks_pvalue <= 0.05) else NA_real_,
    cvm_rejection_rate = if (nrow(successful)) mean(successful$cvm_pvalue <= 0.05) else NA_real_,
    ks_mean_pvalue = if (nrow(successful)) mean(successful$ks_pvalue) else NA_real_,
    cvm_mean_pvalue = if (nrow(successful)) mean(successful$cvm_pvalue) else NA_real_,
    theta_rmse = if (nrow(successful)) sqrt(mean(successful$theta_error^2)) else NA_real_,
    lambda_bias = if (nrow(successful)) mean(successful$lambda_error) else NA_real_,
    lambda_rmse = if (nrow(successful)) sqrt(mean(successful$lambda_error^2)) else NA_real_,
    max_score_mean_norm = if (nrow(successful)) max(successful$score_mean_norm) else NA_real_,
    max_sigma_eigen_error = if (nrow(successful)) max(successful$sigma_eigen_max_error) else NA_real_,
    mean_elapsed_seconds = if (nrow(successful)) mean(successful$elapsed_seconds) else NA_real_
  )
}
groups <- split(existing, interaction(existing$d, existing$theta_name, existing$lambda_true,
                                      existing$beta, existing$n, drop = TRUE))
summary_table <- do.call(rbind, lapply(groups, summarize_null))
utils::write.csv(summary_table, file.path(output_dir, "gof_summary.csv"), row.names = FALSE)

# A deliberately small paired fast/reestimated diagnostic.  It uses only an
# interior, moderate spike and is separate from the null calibration grid.
comparison_design <- expand.grid(
  d = c(2L, 5L), n = c(50L, 200L), replication = seq_len(comparison_M),
  stringsAsFactors = FALSE
)
comparison_design$job_id <- seq_len(nrow(comparison_design))
comparison_path <- file.path(output_dir, "fast_vs_reestimated.csv")
comparison_existing <- if (file.exists(comparison_path)) {
  utils::read.csv(comparison_path, stringsAsFactors = FALSE)
} else data.frame()
comparison_done <- if (nrow(comparison_existing)) as.integer(comparison_existing$job_id) else integer()

run_comparison_job <- function(row) {
  started <- proc.time()[["elapsed"]]
  job_seed <- seed + 7000003L + 100003L * as.integer(row$job_id)
  set.seed(job_seed)
  d <- as.integer(row$d)
  theta <- theta_design(d)$axis_moderate
  x <- rrestricted_spiked_normal(as.integer(row$n), theta, 0.5)
  base <- data.frame(job_id = as.integer(row$job_id), d = d, n = as.integer(row$n),
                     replication = as.integer(row$replication), lambda_true = 0.5,
                     seed = job_seed, stringsAsFactors = FALSE)
  warnings <- character()
  tryCatch({
    common <- list(
      data = x, null = list(type = "composite"), statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(), B = comparison_B, alpha = 0.05,
      n_cores = 1L, seed = job_seed + 1L,
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE),
      control = list(derivative_method = "score_mc", derivative_mc_size = derivative_mc_size,
                     derivative_mc_seed = job_seed + 2L, fast_multiplier_cvm_block_size = cvm_block_size,
                     fast_multiplier_backend = "cpp", fast_multiplier_cpp_kernel = "contiguous_double",
                     fast_multiplier_fuse_ks_cvm = TRUE, fast_multiplier_cache_corrections = "auto"),
      distance_profile_backend = "r"
    )
    run_with_warning_capture <- function(expression) {
      withCallingHandlers(expression, warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      })
    }
    fast <- run_with_warning_capture(do.call(
      multiplier_bootstrap_restricted_spiked_normal,
      c(common, list(bootstrap_method = "fast_multiplier",
                     fast_multiplier_backend = "cpp",
                     fast_multiplier_cpp_kernel = "contiguous_double",
                     fuse_ks_cvm = TRUE, cache_block_corrections = "auto"))
    ))
    slow <- run_with_warning_capture(do.call(
      multiplier_bootstrap_restricted_spiked_normal,
      c(common, list(bootstrap_method = "reestimated"))
    ))
    cbind(base, data.frame(
      status = "ok", error_message = NA_character_,
      warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      fast_ks_pvalue = fast$inference$ks$p_value, fast_cvm_pvalue = fast$inference$cvm$p_value,
      slow_ks_pvalue = slow$inference$ks$p_value, slow_cvm_pvalue = slow$inference$cvm$p_value,
      pvalue_difference_ks = fast$inference$ks$p_value - slow$inference$ks$p_value,
      pvalue_difference_cvm = fast$inference$cvm$p_value - slow$inference$cvm$p_value,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    ))
  }, error = function(error) {
    cbind(base, data.frame(status = "error", error_message = conditionMessage(error),
      warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      fast_ks_pvalue = NA_real_, fast_cvm_pvalue = NA_real_, slow_ks_pvalue = NA_real_,
      slow_cvm_pvalue = NA_real_, pvalue_difference_ks = NA_real_, pvalue_difference_cvm = NA_real_,
      elapsed_seconds = proc.time()[["elapsed"]] - started, stringsAsFactors = FALSE))
  })
}

comparison_pending <- comparison_design[!comparison_design$job_id %in% comparison_done, , drop = FALSE]
message(sprintf("Restricted-spiked paired fast/reestimated diagnostic: %d pending / %d total jobs (M=%d, B=%d).",
                nrow(comparison_pending), nrow(comparison_design), comparison_M, comparison_B))
if (nrow(comparison_pending)) {
  for (begin in seq.int(1L, nrow(comparison_pending), by = cores)) {
    end <- min(nrow(comparison_pending), begin + cores - 1L)
    batch <- comparison_pending[begin:end, , drop = FALSE]
    rows <- parallel::mclapply(seq_len(nrow(batch)), function(i) run_comparison_job(batch[i, , drop = FALSE]),
                               mc.cores = min(cores, nrow(batch)), mc.preschedule = FALSE)
    comparison_existing <- rbind(comparison_existing, do.call(rbind, rows))
    utils::write.csv(comparison_existing, comparison_path, row.names = FALSE)
    message(sprintf("completed %d/%d paired diagnostic jobs", nrow(comparison_existing), nrow(comparison_design)))
  }
}

comparison_ok <- comparison_existing[comparison_existing$status == "ok", , drop = FALSE]
comparison_summary <- if (nrow(comparison_ok)) {
  stats::aggregate(cbind(pvalue_difference_ks, pvalue_difference_cvm, elapsed_seconds) ~ d + n,
                   comparison_ok, mean)
} else comparison_ok
utils::write.csv(comparison_summary, file.path(output_dir, "fast_vs_reestimated_summary.csv"), row.names = FALSE)

writeLines(c(
  "Restricted mean-aligned single-spiked normal GOF validation",
  sprintf("M: %d", M), sprintf("B: %d", B), sprintf("outer cores: %d; inner bootstrap cores: 1", cores),
  sprintf("dimensions: 2,5; n: %s", paste(n_values, collapse = ",")),
  sprintf("strictly positive lambda grid: %s", paste(lambda_values, collapse = ",")),
  sprintf("beta grid: %s", paste(beta_values, collapse = ",")),
  sprintf("theta designs: %.8g e1 and %.8g (1,...,1)/sqrt(d)", axis_norm, diagonal_norm),
  "beta=0: restricted-spiked null; beta>0: (1-beta/2)N(theta,Sigma)+beta/2 N(-theta,Sigma)",
  "KS: sample unique-distance grid; CvM: fast multiplier; backend: C++ contiguous_double fused KS-CvM",
  sprintf("distance profiles: current pipeline default requested as R; derivative: score_mc with auxiliary size %d", derivative_mc_size),
  sprintf("CvM fast block size: %d", cvm_block_size),
  sprintf("paired fast/reestimated diagnostic: lambda=0.5, M=%d, B=%d, n=50,200", comparison_M, comparison_B),
  sprintf("seed: %d", seed)
), file.path(output_dir, "manifest.txt"))

message(sprintf("Restricted-spiked GOF validation complete. Results: %s", output_dir))
