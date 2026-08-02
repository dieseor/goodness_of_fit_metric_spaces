#!/usr/bin/env Rscript

# Prototype benchmark for the rank-two derivative structure of deterministic
# vMF/HvMF profile derivatives. No production backend is changed here.

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
source("bootstrap/multiplier_bootstrap.R")

parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- commandArgs(trailingOnly = TRUE)
  hit <- hit[startsWith(hit, prefix)]
  if (!length(hit)) default else substring(hit[[length(hit)]], nchar(prefix) + 1L)
}

sizes <- as.integer(strsplit(parse_option("sizes", "50,100,200,400"), ",", fixed = TRUE)[[1L]])
dimensions <- as.integer(strsplit(parse_option("dimensions", "2,10"), ",", fixed = TRUE)[[1L]])
B_chunk <- as.integer(parse_option("B-chunk", "100"))
repetitions <- as.integer(parse_option("reps", "10"))
tile_size <- as.integer(parse_option("tile-size", "8"))
seed <- as.integer(parse_option("seed", "20260805"))
output <- parse_option(
  "output",
  file.path("benchmarks", "structured_fast_kernel.csv")
)
if (any(sizes < 2L) || any(dimensions < 2L) || B_chunk < 1L ||
    repetitions < 2L || tile_size < 1L || tile_size > 64L) {
  stop("Invalid benchmark arguments.")
}

if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
prototype_environment <- new.env(parent = globalenv())
Rcpp::sourceCpp(
  code = '
    #include <Rcpp.h>
    #include <R_ext/BLAS.h>
    #include <R_ext/RS.h>
    using namespace Rcpp;

    // [[Rcpp::export]]
    Rcpp::List cpp_structured_fast_sample_ks_cvm(
        Rcpp::NumericMatrix centered_weights,
        Rcpp::NumericMatrix center_projection,
        Rcpp::NumericVector mu_projection,
        Rcpp::IntegerMatrix obs_order_matrix,
        Rcpp::IntegerMatrix tie_end_matrix,
        Rcpp::NumericMatrix coefficient_omega,
        Rcpp::NumericMatrix coefficient_mu,
        const double scale_factor = 1.0) {
      const int n_reps = centered_weights.nrow();
      const int n = centered_weights.ncol();
      const int n_centers = obs_order_matrix.nrow();
      if (center_projection.nrow() != n_reps ||
          center_projection.ncol() != n_centers ||
          mu_projection.size() != n_reps ||
          obs_order_matrix.ncol() != n ||
          tie_end_matrix.nrow() != n_centers ||
          tie_end_matrix.ncol() != n ||
          coefficient_omega.nrow() != n_centers ||
          coefficient_omega.ncol() != n ||
          coefficient_mu.nrow() != n_centers ||
          coefficient_mu.ncol() != n) {
        Rcpp::stop("Structured kernel inputs have incompatible dimensions.");
      }

      Rcpp::NumericVector ks(n_reps, 0.0);
      Rcpp::NumericVector cvm_sum(n_reps, 0.0);
      std::vector<double> running(static_cast<std::size_t>(n_reps));
      std::vector<double> center_cvm_sum(static_cast<std::size_t>(n_reps));
      const double root_n = std::sqrt(static_cast<double>(n));
      const double cvm_denominator =
        static_cast<double>(n) * static_cast<double>(n);

      for (int center = 0; center < n_centers; ++center) {
        std::fill(running.begin(), running.end(), 0.0);
        std::fill(center_cvm_sum.begin(), center_cvm_sum.end(), 0.0);
        int threshold = 0;
        const double* projected_center = center_projection.begin() +
          static_cast<R_xlen_t>(center) * n_reps;
        for (int sorted_index = 0; sorted_index < n; ++sorted_index) {
          const int observation = obs_order_matrix(center, sorted_index) - 1;
          const double* weight_column = centered_weights.begin() +
            static_cast<R_xlen_t>(observation) * n_reps;
          for (int rep = 0; rep < n_reps; ++rep) {
            running[static_cast<std::size_t>(rep)] += weight_column[rep];
          }
          while (threshold < n &&
                 tie_end_matrix(center, threshold) - 1 == sorted_index) {
            const double coefficient_center = coefficient_omega(center, threshold);
            const double coefficient_location = coefficient_mu(center, threshold);
            for (int rep = 0; rep < n_reps; ++rep) {
              const double correction =
                coefficient_center * projected_center[rep] +
                coefficient_location * mu_projection[rep];
              const double process = scale_factor *
                (running[static_cast<std::size_t>(rep)] -
                 correction) / root_n;
              const double absolute_process = std::abs(process);
              if (absolute_process > ks[rep]) ks[rep] = absolute_process;
              center_cvm_sum[static_cast<std::size_t>(rep)] += process * process;
            }
            ++threshold;
          }
        }
        if (threshold != n) Rcpp::stop("Tie endpoints must be nondecreasing.");
        for (int rep = 0; rep < n_reps; ++rep) {
          cvm_sum[rep] += center_cvm_sum[static_cast<std::size_t>(rep)];
        }
      }
      for (int rep = 0; rep < n_reps; ++rep) {
        cvm_sum[rep] /= cvm_denominator;
      }
      return Rcpp::List::create(
        Rcpp::Named("ks") = ks,
        Rcpp::Named("cvm") = cvm_sum
      );
    }

    // [[Rcpp::export]]
    Rcpp::List cpp_cache_friendly_fast_sample_ks_cvm(
        Rcpp::NumericMatrix centered_weights,
        Rcpp::NumericMatrix score_block,
        Rcpp::IntegerMatrix obs_order_matrix,
        Rcpp::IntegerMatrix tie_end_matrix,
        Rcpp::NumericMatrix correction_matrix,
        const double scale_factor = 1.0) {
      const int n_reps = centered_weights.nrow();
      const int n = centered_weights.ncol();
      const int n_centers = obs_order_matrix.nrow();
      const int n_parameters = score_block.ncol();
      if (score_block.nrow() != n_reps || obs_order_matrix.ncol() != n ||
          tie_end_matrix.nrow() != n_centers || tie_end_matrix.ncol() != n ||
          correction_matrix.nrow() != n_centers * n ||
          correction_matrix.ncol() != n_parameters) {
        Rcpp::stop("Cache-friendly kernel inputs have incompatible dimensions.");
      }

      Rcpp::NumericVector ks(n_reps, 0.0);
      Rcpp::NumericVector cvm(n_reps, 0.0);
      Rcpp::NumericMatrix correction(n_reps, n);
      std::vector<double> running(static_cast<std::size_t>(n_reps));
      std::vector<double> center_cvm_sum(static_cast<std::size_t>(n_reps));
      const double root_n = std::sqrt(static_cast<double>(n));
      const double cvm_denominator =
        static_cast<double>(n) * static_cast<double>(n);
      const char trans_n = \'N\';
      const char trans_t = \'T\';
      const double alpha = 1.0;
      const double beta = 0.0;
      const int lda = n_reps;
      const int ldb = correction_matrix.nrow();
      const int ldc = n_reps;

      for (int center = 0; center < n_centers; ++center) {
        const double* correction_center = correction_matrix.begin() +
          static_cast<R_xlen_t>(center) * n;
        F77_CALL(dgemm)(
          &trans_n, &trans_t, &n_reps, &n, &n_parameters,
          &alpha, score_block.begin(), &lda,
          correction_center, &ldb, &beta,
          correction.begin(), &ldc FCONE FCONE
        );
        std::fill(running.begin(), running.end(), 0.0);
        std::fill(center_cvm_sum.begin(), center_cvm_sum.end(), 0.0);
        int threshold = 0;
        for (int sorted_index = 0; sorted_index < n; ++sorted_index) {
          const int observation = obs_order_matrix(center, sorted_index) - 1;
          const double* weight_column = centered_weights.begin() +
            static_cast<R_xlen_t>(observation) * n_reps;
          for (int rep = 0; rep < n_reps; ++rep) {
            running[static_cast<std::size_t>(rep)] += weight_column[rep];
          }
          while (threshold < n &&
                 tie_end_matrix(center, threshold) - 1 == sorted_index) {
            const double* correction_column = correction.begin() +
              static_cast<R_xlen_t>(threshold) * n_reps;
            for (int rep = 0; rep < n_reps; ++rep) {
              const double process = scale_factor *
                (running[static_cast<std::size_t>(rep)] -
                 correction_column[rep]) / root_n;
              const double absolute_process = std::abs(process);
              if (absolute_process > ks[rep]) ks[rep] = absolute_process;
              center_cvm_sum[static_cast<std::size_t>(rep)] += process * process;
            }
            ++threshold;
          }
        }
        if (threshold != n) Rcpp::stop("Tie endpoints must be nondecreasing.");
        for (int rep = 0; rep < n_reps; ++rep) {
          cvm[rep] += center_cvm_sum[static_cast<std::size_t>(rep)];
        }
      }
      for (int rep = 0; rep < n_reps; ++rep) cvm[rep] /= cvm_denominator;
      return Rcpp::List::create(
        Rcpp::Named("ks") = ks,
        Rcpp::Named("cvm") = cvm
      );
    }

    // [[Rcpp::export]]
    Rcpp::List cpp_tiled_fast_sample_ks_cvm(
        Rcpp::NumericMatrix centered_weights,
        Rcpp::NumericMatrix score_block,
        Rcpp::IntegerMatrix obs_order_matrix,
        Rcpp::IntegerMatrix tie_end_matrix,
        Rcpp::NumericMatrix correction_matrix,
        const double scale_factor = 1.0,
        const int tile_size = 8) {
      const int n_reps = centered_weights.nrow();
      const int n = centered_weights.ncol();
      const int n_centers = obs_order_matrix.nrow();
      const int n_parameters = score_block.ncol();
      if (tile_size < 1 || tile_size > 64 || score_block.nrow() != n_reps ||
          obs_order_matrix.ncol() != n ||
          tie_end_matrix.nrow() != n_centers || tie_end_matrix.ncol() != n ||
          correction_matrix.nrow() != n_centers * n ||
          correction_matrix.ncol() != n_parameters) {
        Rcpp::stop("Tiled kernel inputs have incompatible dimensions.");
      }

      Rcpp::NumericVector ks(n_reps, 0.0);
      Rcpp::NumericVector cvm(n_reps, 0.0);
      Rcpp::NumericMatrix correction(n_reps, n);
      const double root_n = std::sqrt(static_cast<double>(n));
      const double cvm_denominator =
        static_cast<double>(n) * static_cast<double>(n);
      const char trans_n = \'N\';
      const char trans_t = \'T\';
      const double alpha = 1.0;
      const double beta = 0.0;
      const int lda = n_reps;
      const int ldb = correction_matrix.nrow();
      const int ldc = n_reps;

      for (int center = 0; center < n_centers; ++center) {
        const double* correction_center = correction_matrix.begin() +
          static_cast<R_xlen_t>(center) * n;
        F77_CALL(dgemm)(
          &trans_n, &trans_t, &n_reps, &n, &n_parameters,
          &alpha, score_block.begin(), &lda,
          correction_center, &ldb, &beta,
          correction.begin(), &ldc FCONE FCONE
        );
        for (int rep_start = 0; rep_start < n_reps; rep_start += tile_size) {
          const int tile_n = std::min(tile_size, n_reps - rep_start);
          long double running[64] = {0.0L};
          long double center_cvm_sum[64] = {0.0L};
          int threshold = 0;
          for (int sorted_index = 0; sorted_index < n; ++sorted_index) {
            const int observation = obs_order_matrix(center, sorted_index) - 1;
            const double* weight_column = centered_weights.begin() +
              static_cast<R_xlen_t>(observation) * n_reps + rep_start;
            for (int offset = 0; offset < tile_n; ++offset) {
              running[offset] += static_cast<long double>(weight_column[offset]);
            }
            while (threshold < n &&
                   tie_end_matrix(center, threshold) - 1 == sorted_index) {
              const double* correction_column = correction.begin() +
                static_cast<R_xlen_t>(threshold) * n_reps + rep_start;
              for (int offset = 0; offset < tile_n; ++offset) {
                const int rep = rep_start + offset;
                const double process = scale_factor *
                  (static_cast<double>(running[offset]) -
                   correction_column[offset]) / root_n;
                const double absolute_process = std::abs(process);
                if (absolute_process > ks[rep]) ks[rep] = absolute_process;
                center_cvm_sum[offset] +=
                  static_cast<long double>(process * process);
              }
              ++threshold;
            }
          }
          if (threshold != n) Rcpp::stop("Tie endpoints must be nondecreasing.");
          for (int offset = 0; offset < tile_n; ++offset) {
            cvm[rep_start + offset] +=
              static_cast<double>(center_cvm_sum[offset]);
          }
        }
      }
      for (int rep = 0; rep < n_reps; ++rep) cvm[rep] /= cvm_denominator;
      return Rcpp::List::create(
        Rcpp::Named("ks") = ks,
        Rcpp::Named("cvm") = cvm
      );
    }
  ',
  env = prototype_environment,
  rebuild = TRUE,
  showOutput = FALSE,
  verbose = FALSE
)
ensure_distance_profile_cpp_loaded()

median_elapsed <- function(fun, repetitions) {
  timings <- numeric(repetitions)
  value <- NULL
  for (index in seq_len(repetitions)) {
    timings[[index]] <- system.time(value <- fun())[["elapsed"]]
  }
  list(value = value, median = median(timings), minimum = min(timings))
}

set.seed(seed)
rows <- list()
row_index <- 0L
for (q in dimensions) {
  p <- q + 1L
  for (n in sizes) {
    centers <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    centers <- centers / sqrt(rowSums(centers^2))
    mu <- c(1, rep.int(0, q))
    coefficient_omega <- matrix(stats::rnorm(n * n, sd = 0.15), n, n)
    coefficient_mu <- matrix(stats::rnorm(n * n, sd = 0.15), n, n)
    center_rows <- rep(seq_len(n), each = n)
    derivative <- as.vector(t(coefficient_omega)) * centers[center_rows, , drop = FALSE] +
      as.vector(t(coefficient_mu)) * matrix(mu, n * n, p, byrow = TRUE)

    information_factor <- matrix(stats::rnorm(p * p), p, p)
    information <- crossprod(information_factor) + diag(p)
    information_inverse <- solve(information)
    correction_matrix <- derivative %*% t(information_inverse)
    score <- matrix(stats::rnorm(n * p), n, p)
    centered_weights <- matrix(stats::rexp(B_chunk * n) - 1, B_chunk, n)
    order_matrix <- t(replicate(n, sample.int(n)))
    tie_end_matrix <- matrix(rep(seq_len(n), n), nrow = n, byrow = TRUE)

    current_call <- function() {
      score_block <- centered_weights %*% score
      distance_profile_cpp_call(
        "cpp_fast_sample_ks_cvm_stats",
        centered_weights = centered_weights,
        score_block = score_block,
        obs_order_matrix = order_matrix,
        tie_end_matrix = tie_end_matrix,
        correction_matrix = correction_matrix,
        scale_factor = 1,
        compute_ks = TRUE,
        compute_cvm = TRUE
      )
    }
    cache_friendly_call <- function() {
      score_block <- centered_weights %*% score
      prototype_environment$cpp_cache_friendly_fast_sample_ks_cvm(
        centered_weights = centered_weights,
        score_block = score_block,
        obs_order_matrix = order_matrix,
        tie_end_matrix = tie_end_matrix,
        correction_matrix = correction_matrix,
        scale_factor = 1
      )
    }
    tiled_call <- function() {
      score_block <- centered_weights %*% score
      prototype_environment$cpp_tiled_fast_sample_ks_cvm(
        centered_weights = centered_weights,
        score_block = score_block,
        obs_order_matrix = order_matrix,
        tie_end_matrix = tie_end_matrix,
        correction_matrix = correction_matrix,
        scale_factor = 1,
        tile_size = tile_size
      )
    }
    structured_call <- function() {
      score_block <- centered_weights %*% score
      solved_score <- score_block %*% information_inverse
      prototype_environment$cpp_structured_fast_sample_ks_cvm(
        centered_weights = centered_weights,
        center_projection = solved_score %*% t(centers),
        mu_projection = drop(solved_score %*% mu),
        obs_order_matrix = order_matrix,
        tie_end_matrix = tie_end_matrix,
        coefficient_omega = coefficient_omega,
        coefficient_mu = coefficient_mu,
        scale_factor = 1
      )
    }

    invisible(current_call())
    invisible(cache_friendly_call())
    invisible(tiled_call())
    invisible(structured_call())
    current <- median_elapsed(current_call, repetitions)
    cache_friendly <- median_elapsed(cache_friendly_call, repetitions)
    tiled <- median_elapsed(tiled_call, repetitions)
    structured <- median_elapsed(structured_call, repetitions)
    max_error <- max(
      abs(current$value$ks - cache_friendly$value$ks),
      abs(current$value$cvm - cache_friendly$value$cvm),
      abs(current$value$ks - tiled$value$ks),
      abs(current$value$cvm - tiled$value$cvm),
      abs(current$value$ks - structured$value$ks),
      abs(current$value$cvm - structured$value$cvm)
    )
    row_index <- row_index + 1L
    rows[[row_index]] <- data.frame(
      q = q,
      p = p,
      n = n,
      B_chunk = B_chunk,
      repetitions = repetitions,
      tile_size = tile_size,
      current_median_seconds = current$median,
      cache_friendly_double_median_seconds = cache_friendly$median,
      tiled_median_seconds = tiled$median,
      structured_double_median_seconds = structured$median,
      cache_friendly_double_speedup = current$median / cache_friendly$median,
      tiled_speedup = current$median / tiled$median,
      structured_double_speedup = current$median / structured$median,
      current_minimum_seconds = current$minimum,
      cache_friendly_double_minimum_seconds = cache_friendly$minimum,
      tiled_minimum_seconds = tiled$minimum,
      structured_double_minimum_seconds = structured$minimum,
      max_abs_statistic_error = max_error,
      current_derivative_MiB = as.numeric(object.size(correction_matrix)) / 1024^2,
      structured_derivative_MiB = as.numeric(
        object.size(coefficient_omega) + object.size(coefficient_mu) +
          object.size(centers) + object.size(mu)
      ) / 1024^2,
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output, row.names = FALSE)
print(results, row.names = FALSE)
cat(sprintf("Wrote %s\n", output))
