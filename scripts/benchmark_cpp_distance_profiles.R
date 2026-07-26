#!/usr/bin/env Rscript

# Confirmatory paired benchmark for the retained optional C++ backends.
# Example:
#   Rscript scripts/benchmark_cpp_distance_profiles.R --reps=20 --sizes=50,100,200

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  matches <- commandArgs(trailingOnly = TRUE)
  matches <- matches[startsWith(matches, prefix)]
  if (length(matches) == 0L) return(default)
  substring(matches[[length(matches)]], nchar(prefix) + 1L)
}

repetitions <- as.integer(parse_option("reps", "20"))
sizes <- as.integer(strsplit(parse_option("sizes", "50,100,200"), ",", fixed = TRUE)[[1L]])
normal_B <- as.integer(parse_option("normal-b", "200"))
weighted_B <- as.integer(parse_option("weighted-b", "10"))
output_dir <- parse_option("output", file.path("output", "cpp_distance_profile_benchmark"))

if (!is.finite(repetitions) || repetitions < 2L ||
    any(!is.finite(sizes)) || any(sizes < 2L) ||
    !is.finite(normal_B) || normal_B < 1L ||
    !is.finite(weighted_B) || weighted_B < 1L) {
  stop("Invalid benchmark arguments.")
}

source(file.path("bootstrap", "multiplier_bootstrap.R"))

strip_backend_timings <- function(result) {
  timing_names <- grep("_seconds$", names(result$diagnostics), value = TRUE)
  result$diagnostics[c(
    timing_names,
    "distance_profile_backend_requested",
    "distance_profile_backend_effective"
  )] <- NULL
  result
}

time_one <- function(run, backend) {
  answer <- NULL
  elapsed <- system.time(answer <- run(backend))[["elapsed"]]
  list(elapsed = unname(elapsed), answer = answer)
}

paired_benchmark <- function(family, n, run, repetitions) {
  invisible(run("cpp"))
  rows <- vector("list", repetitions)
  for (pair_index in seq_len(repetitions)) {
    backend_order <- if (pair_index %% 2L == 1L) c("r", "cpp") else c("cpp", "r")
    timed <- setNames(vector("list", 2L), c("r", "cpp"))
    for (backend in backend_order) timed[[backend]] <- time_one(run, backend)
    exact <- identical(
      strip_backend_timings(timed$r$answer),
      strip_backend_timings(timed$cpp$answer)
    )
    rows[[pair_index]] <- data.frame(
      family = family,
      n = n,
      pair = pair_index,
      r_seconds = timed$r$elapsed,
      cpp_seconds = timed$cpp$elapsed,
      gain = 1 - timed$cpp$elapsed / timed$r$elapsed,
      speedup = timed$r$elapsed / timed$cpp$elapsed,
      exact = exact
    )
    message(sprintf(
      "%s n=%d pair=%02d/%02d: R=%.3fs C++=%.3fs gain=%+.1f%% exact=%s",
      family, n, pair_index, repetitions,
      timed$r$elapsed, timed$cpp$elapsed,
      100 * rows[[pair_index]]$gain, exact
    ))
  }
  do.call(rbind, rows)
}

make_normal_run <- function(n) {
  set.seed(10000L + n)
  x <- rnorm(n, mean = 1.25, sd = 2.1)
  grid <- list(
    omega_grid = seq(-4, 6, length.out = 17L),
    t_grid = seq(0, 7, length.out = 19L)
  )
  function(backend) {
    multiplier_bootstrap_normal(
      data = x,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = grid,
      B = normal_B,
      seed = 901L,
      n_cores = 1L,
      unknown_param = "both",
      keep = list(
        observed_process = TRUE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      distance_profile_backend = backend
    )
  }
}

make_weighted_run <- function(n) {
  set.seed(20000L + n)
  mu <- c(0, 0, 1)
  x <- r_sph_small_circle_weighted_mixture2(n, mu, 0.6, 10, 0.45, 8, 0.25)
  grid <- list(
    omega_grid = generate_canonical_lattice(5L, dim = 3L),
    t_grid = seq(1e-8, pi - 1e-8, length.out = 7L)
  )
  function(backend) {
    multiplier_bootstrap_small_circle_weighted_mixture2(
      data = x,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = grid,
      B = weighted_B,
      seed = 902L,
      n_cores = 1L,
      distance_type = "geodesic",
      keep = list(
        observed_process = TRUE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      control = list(
        small_circle_weighted_mixture2_L_max = 120L,
        small_circle_weighted_mixture2_quad_n = 300L,
        small_circle_weighted_mixture2_n_starts = 1L,
        small_circle_weighted_mixture2_optim_control = list(maxit = 80L, reltol = 1e-6)
      ),
      distance_profile_backend = backend
    )
  }
}

source_file <- resolve_distance_profile_backend_path("cpp", "distance_profile_backend.cpp")
compile_environment <- new.env(parent = globalenv())
compile_cache <- tempfile("distance-profile-cold-compile-")
dir.create(compile_cache)
cold_compile_seconds <- unname(system.time(Rcpp::sourceCpp(
  file = source_file,
  env = compile_environment,
  cacheDir = compile_cache,
  rebuild = TRUE,
  showOutput = FALSE,
  verbose = FALSE
))[["elapsed"]])
hot_load_seconds <- unname(system.time(ensure_distance_profile_cpp_loaded())[["elapsed"]])

results <- list()
result_index <- 1L
for (n in sizes) {
  results[[result_index]] <- paired_benchmark("normal", n, make_normal_run(n), repetitions)
  result_index <- result_index + 1L
  results[[result_index]] <- paired_benchmark(
    "small_circle_weighted_mixture2", n, make_weighted_run(n), repetitions
  )
  result_index <- result_index + 1L
}
pair_results <- do.call(rbind, results)

summarize_family <- function(rows) {
  set.seed(20260722L)
  bootstrap_medians <- replicate(
    10000L,
    median(sample(rows$speedup, nrow(rows), replace = TRUE))
  )
  median_delta <- median(rows$r_seconds - rows$cpp_seconds)
  data.frame(
    family = rows$family[[1L]],
    pairs = nrow(rows),
    all_exact = all(rows$exact),
    median_r_seconds = median(rows$r_seconds),
    median_cpp_seconds = median(rows$cpp_seconds),
    median_gain = median(rows$gain),
    median_speedup = median(rows$speedup),
    speedup_ci_low = unname(quantile(bootstrap_medians, 0.025)),
    speedup_ci_high = unname(quantile(bootstrap_medians, 0.975)),
    worst_gain = min(rows$gain),
    cold_compile_seconds = cold_compile_seconds,
    hot_load_seconds = hot_load_seconds,
    approximate_runs_to_amortize = if (median_delta > 0) cold_compile_seconds / median_delta else Inf
  )
}

summary_results <- do.call(rbind, lapply(
  split(pair_results, pair_results$family),
  summarize_family
))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(pair_results, file.path(output_dir, "paired_results.csv"), row.names = FALSE)
write.csv(summary_results, file.path(output_dir, "summary.csv"), row.names = FALSE)
saveRDS(
  list(
    pairs = pair_results,
    summary = summary_results,
    session_info = sessionInfo(),
    arguments = list(
      repetitions = repetitions,
      sizes = sizes,
      normal_B = normal_B,
      weighted_B = weighted_B
    )
  ),
  file.path(output_dir, "benchmark.rds")
)

print(summary_results, row.names = FALSE)
if (!all(pair_results$exact)) stop("At least one benchmark pair was not bitwise identical.")
