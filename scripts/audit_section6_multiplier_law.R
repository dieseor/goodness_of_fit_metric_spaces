#!/usr/bin/env Rscript

# Sensitivity of the *fast* multiplier bootstrap to the positive multiplier
# law.  All laws below have mean and standard deviation one.  The two-point
# law W in {0,2} is exactly the non-negative version of centred Rademacher
# multipliers (W - 1 in {-1,1}); it is admissible for the weighted empirical
# measure, unlike raw signed Gaussian/Rademacher weights.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    bits <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[bits[[1L]]]] <- if (length(bits) == 1L) "TRUE" else paste(bits[-1L], collapse = "=")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
d <- as.integer(args$d %||% 10L)
n <- as.integer(args$n %||% 100L)
M <- as.integer(args$M %||% 100L)
B <- as.integer(args$B %||% 299L)
cores <- as.integer(args$cores %||% 10L)
derivative_mc_size <- as.integer(args$derivative_mc_size %||% 1000L)
seed <- as.integer(args$seed %||% 20260807L)
replicate_offset <- as.integer(args$replicate_offset %||% 0L)
output_dir <- args$output_dir %||% file.path(
  "simulation_results", "section6_new_scenarios",
  sprintf("audit_multiplier_law_normal_d%d_n%d_M%d_B%d", d, n, M, B)
)

if (d < 2L || n < 10L || M < 2L || B < 99L || cores < 1L || derivative_mc_size < 100L) {
  stop("Invalid diagnostic sizes.")
}

source(file.path("bootstrap", "multiplier_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

laws <- list(
  exponential = list(
    name = "Exp(1)", generator = function(n) stats::rexp(n), mean = 1, sd = 1
  ),
  two_point = list(
    name = "two-point {0,2}", generator = function(n) 2 * stats::rbinom(n, 1, .5),
    mean = 1, sd = 1
  ),
  lognormal = list(
    name = "Lognormal(mean=1, sd=1)",
    generator = function(n) stats::rlnorm(n, meanlog = -.5 * log(2), sdlog = sqrt(log(2))),
    mean = 1, sd = 1
  )
)

one_replication <- function(i) {
  replication_id <- replicate_offset + i
  set.seed(seed + 1009L * replication_id)
  x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = diag(d))
  spec <- make_mvnormal_spec("both")
  theta <- spec$fit_theta(x, NULL, list(type = "composite"), list())
  rows <- lapply(seq_along(laws), function(j) {
    result <- multiplier_bootstrap_gof(
      data = x, spec = spec, null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = B, alpha = .05, multipliers = laws[[j]], n_cores = 1L,
      seed = seed + 100000L * replication_id + j,
      observed_theta_hat = theta,
      bootstrap_method = "fast_multiplier",
      keep = list(observed_process = FALSE, bootstrap_statistics = FALSE),
      control = list(
        derivative_mc_size = derivative_mc_size,
        derivative_mc_seed = seed + 1000000L * replication_id,
        # Keep this diagnostic independent of a concurrent sourceCpp build.
        # The multiplier law changes no numerical formula, so the R backend is
        # the appropriate stable route for this comparison.
        fast_multiplier_backend = "r",
        fast_multiplier_fuse_ks_cvm = TRUE,
        fast_multiplier_cvm_block_size = 50L
      )
    )
    data.frame(
      replication = replication_id, multiplier = names(laws)[[j]],
      ks_p_value = result$inference$ks$p_value,
      cvm_p_value = result$inference$cvm$p_value,
      ks_reject = result$inference$ks$reject,
      cvm_reject = result$inference$cvm$reject,
      backend = result$diagnostics$fast_multiplier_backend_effective,
      effective_method = result$diagnostics$effective_bootstrap_method,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

message(sprintf("Multiplier-law audit: %d Normal d=%d samples, B=%d, using %d outer cores ...", M, d, B, min(M, cores)))
if (.Platform$OS.type == "unix" && cores > 1L) {
  values <- parallel::mclapply(seq_len(M), one_replication,
    mc.cores = min(M, cores), mc.set.seed = FALSE
  )
} else values <- lapply(seq_len(M), one_replication)
results <- do.call(rbind, values)
summary <- do.call(rbind, lapply(split(results, results$multiplier), function(x) {
  data.frame(
    multiplier = x$multiplier[[1L]], M = nrow(x), B = B,
    ks_size = mean(x$ks_reject), cvm_size = mean(x$cvm_reject),
    ks_mcse = sqrt(mean(x$ks_reject) * (1 - mean(x$ks_reject)) / nrow(x)),
    cvm_mcse = sqrt(mean(x$cvm_reject) * (1 - mean(x$cvm_reject)) / nrow(x)),
    median_ks_p_value = stats::median(x$ks_p_value),
    median_cvm_p_value = stats::median(x$cvm_p_value),
    backend_all_r = all(x$backend == "r"),
    all_fast = all(x$effective_method == "fast_multiplier"),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(results, file.path(output_dir, "replications.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(d = d, n = n, M = M, B = B, cores = cores,
  derivative_mc_size = derivative_mc_size, seed = seed,
  replicate_offset = replicate_offset
), file.path(output_dir, "config.csv"), row.names = FALSE)
message("Audit written to: ", normalizePath(output_dir))
