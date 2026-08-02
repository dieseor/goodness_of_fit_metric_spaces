#!/usr/bin/env Rscript

# One-point null calibration for the composite vMF GOF tests.
# The data-generating point (mu, kappa) is explicit, while the fitted null
# remains composite in (mu, kappa).

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source("bootstrap/calibration_study.R")

parse_cli <- function(args) {
  parsed <- list()
  for (arg in args[startsWith(args, "--")]) {
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    parsed[[pieces[[1L]]]] <- if (length(pieces) == 1L) "TRUE" else paste(pieces[-1L], collapse = "=")
  }
  parsed
}

require_arg <- function(x, name) {
  if (is.null(x) || !nzchar(x)) stop(sprintf("Supply `--%s=...`.", name))
  x
}

parse_numeric_csv <- function(x, name) {
  values <- as.numeric(strsplit(require_arg(x, name), ",", fixed = TRUE)[[1L]])
  if (any(!is.finite(values))) stop(sprintf("`--%s` must contain finite numbers.", name))
  values
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
mu <- parse_numeric_csv(args$mu, "mu")
if (length(mu) != 3L) stop("`--mu` must contain exactly three coordinates.")
if (sum(mu^2) <= 0) stop("`--mu` must be nonzero.")
mu <- mu / sqrt(sum(mu^2))

kappa <- as.numeric(require_arg(args$kappa, "kappa"))
n <- as.integer(require_arg(args$n, "n"))
M <- as.integer(require_arg(args$M, "M"))
B <- as.integer(require_arg(args$B, "B"))
cores <- as.integer(require_arg(args$cores, "cores"))
output_dir <- require_arg(args$output_dir, "output_dir")
base_seed <- as.integer(args$seed %||% 20260725L)

if (!is.finite(kappa) || kappa <= 0) stop("`--kappa` must be strictly positive.")
if (!is.finite(n) || n <= 0L || !is.finite(M) || M <= 0L || !is.finite(B) || B <= 0L || !is.finite(cores) || cores <= 0L) {
  stop("`n`, `M`, `B`, and `cores` must be strictly positive integers.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(output_dir, "raw_results.csv")
summary_path <- file.path(output_dir, "summary.csv")
scenario_path <- file.path(output_dir, "scenario.txt")

if (file.exists(raw_path)) {
  existing_names <- names(utils::read.csv(
    raw_path,
    nrows = 1L,
    stringsAsFactors = FALSE
  ))
  required_existing <- c(
    "derivative_method_effective",
    "fast_multiplier_backend_effective",
    "fast_multiplier_cpp_kernel_effective",
    "fast_multiplier_fuse_ks_cvm_effective"
  )
  missing_existing <- setdiff(required_existing, existing_names)
  if (length(missing_existing)) {
    stop(
      paste(
        "The existing calibration results predate explicit fast-kernel tracking.",
        "Use a new output directory; bootstrap implementations must not be",
        "mixed silently."
      ),
      call. = FALSE
    )
  }
}

writeLines(c(
  sprintf("null: vMF((%.8f, %.8f, %.8f), kappa = %.8f)", mu[[1L]], mu[[2L]], mu[[3L]], kappa),
  sprintf("n: %d", n),
  "beta: 0",
  sprintf("M: %d", M),
  sprintf("B: %d", B),
  sprintf("cores: %d", cores),
  "statistics: KS sample, CvM",
  "bootstrap: fast_multiplier",
  "derivative_method_requested: quadrature",
  "fast_multiplier_backend_requested: cpp",
  "fast_multiplier_cpp_kernel_requested: contiguous_double",
  "fast_multiplier_fuse_ks_cvm_requested: TRUE",
  "fast_multiplier_cvm_block_size: 50"
), scenario_path)

empty_results <- data.frame(
  rep = integer(), ks_pvalue = numeric(), cvm_pvalue = numeric(),
  ks_reject = logical(), cvm_reject = logical(), status = character(),
  derivative_method_requested = character(), derivative_method_effective = character(),
  fast_multiplier_backend_requested = character(), fast_multiplier_backend_effective = character(),
  fast_multiplier_cpp_kernel_requested = character(), fast_multiplier_cpp_kernel_effective = character(),
  fast_multiplier_fuse_ks_cvm_requested = logical(), fast_multiplier_fuse_ks_cvm_effective = logical(),
  error_message = character(), elapsed_seconds = numeric(), stringsAsFactors = FALSE
)
results <- if (file.exists(raw_path)) utils::read.csv(raw_path, stringsAsFactors = FALSE) else empty_results

seed_for <- function(rep_id, stream) {
  as.integer((as.numeric(base_seed) + 1000003 * as.numeric(rep_id) + 10000019 * as.numeric(stream)) %% 2147483646) + 1L
}

one_replication <- function(rep_id) {
  started <- proc.time()[["elapsed"]]
  tryCatch({
    data_seed <- seed_for(rep_id, 0L)
    bootstrap_seed <- seed_for(rep_id, 1L)
    set.seed(data_seed)
    x <- rotasym::r_vMF(n = n, mu = mu, kappa = kappa)
    fit <- multiplier_bootstrap_vmf(
      data = x,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = B,
      alpha = 0.05,
      n_cores = 1L,
      seed = bootstrap_seed,
      bootstrap_method = "fast_multiplier",
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE, bootstrap_thetas = FALSE),
      control = list(
        derivative_method = "quadrature",
        fast_multiplier_cvm_block_size = 50L,
        fast_multiplier_backend = "cpp",
        fast_multiplier_cpp_kernel = "contiguous_double",
        fast_multiplier_fuse_ks_cvm = TRUE,
        fast_multiplier_cache_corrections = "auto",
        vmf_profile_method = "tabulated",
        vmf_profile_n_u = 4097L
      ),
      distance_type = "geodesic",
      unknown_param = "xi"
    )
    data.frame(
      rep = rep_id,
      ks_pvalue = fit$inference$ks$p_value,
      cvm_pvalue = fit$inference$cvm$p_value,
      ks_reject = fit$inference$ks$reject,
      cvm_reject = fit$inference$cvm$reject,
      derivative_method_requested = "quadrature",
      derivative_method_effective = fit$diagnostics$derivative_method_effective %||%
        fit$diagnostics$derivative_method %||% NA_character_,
      fast_multiplier_backend_requested = "cpp",
      fast_multiplier_backend_effective = fit$diagnostics$fast_multiplier_backend_effective %||% NA_character_,
      fast_multiplier_cpp_kernel_requested = "contiguous_double",
      fast_multiplier_cpp_kernel_effective = fit$diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_,
      fast_multiplier_fuse_ks_cvm_requested = TRUE,
      fast_multiplier_fuse_ks_cvm_effective = isTRUE(fit$diagnostics$fast_multiplier_fuse_ks_cvm_effective),
      status = "ok",
      error_message = NA_character_,
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )
  }, error = function(e) {
    data.frame(
      rep = rep_id,
      ks_pvalue = NA_real_, cvm_pvalue = NA_real_, ks_reject = NA, cvm_reject = NA,
      derivative_method_requested = "quadrature", derivative_method_effective = NA_character_,
      fast_multiplier_backend_requested = "cpp", fast_multiplier_backend_effective = NA_character_,
      fast_multiplier_cpp_kernel_requested = "contiguous_double", fast_multiplier_cpp_kernel_effective = NA_character_,
      fast_multiplier_fuse_ks_cvm_requested = TRUE, fast_multiplier_fuse_ks_cvm_effective = NA,
      status = "error", error_message = conditionMessage(e),
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )
  })
}

done <- results$rep[results$status == "ok"]
pending <- setdiff(seq_len(M), done)
batches <- split(pending, ceiling(seq_along(pending) / max(1L, cores)))

for (batch in batches) {
  if (length(batch) == 0L) next
  rows <- parallel::mclapply(batch, one_replication, mc.cores = min(cores, length(batch)), mc.preschedule = FALSE)
  results <- rbind(results, do.call(rbind, rows))
  results <- results[order(results$rep), , drop = FALSE]
  utils::write.csv(results, raw_path, row.names = FALSE)
  completed <- sum(results$status == "ok")
  cat(sprintf("\r[%6.2f%%] %d/%d completed", 100 * completed / M, completed, M))
  flush.console()
}
if (length(pending) > 0L) cat("\n")

valid <- results[results$status == "ok", , drop = FALSE]
summary <- data.frame(
  null_mu_x = mu[[1L]], null_mu_y = mu[[2L]], null_mu_z = mu[[3L]],
  kappa = kappa, n = n, beta = 0, M = nrow(valid), B = B,
  power_ks = mean(valid$ks_reject), power_cvm = mean(valid$cvm_reject)
)
utils::write.csv(summary, summary_path, row.names = FALSE)
print(summary)
