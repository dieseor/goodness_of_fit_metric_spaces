#!/usr/bin/env Rscript

# Paired microbenchmark for the sample-to-sample vMF/HvMF distance matrices.
# It compares the current full BLAS product with a symmetric C++ loop that
# evaluates only the upper triangle. No production backend is changed here.

Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  args <- commandArgs(trailingOnly = TRUE)
  matches <- args[startsWith(args, prefix)]
  if (!length(matches)) return(default)
  substring(matches[[length(matches)]], nchar(prefix) + 1L)
}

repetitions <- as.integer(parse_option("reps", "50"))
inner_iterations <- as.integer(parse_option("inner", "20"))
sizes <- as.integer(strsplit(parse_option("sizes", "50,100,200,400"), ",", fixed = TRUE)[[1L]])
dimensions <- as.integer(strsplit(parse_option("dimensions", "2,10"), ",", fixed = TRUE)[[1L]])
output_path <- parse_option(
  "output",
  file.path("benchmarks", "vmf_hvmf_pairwise_distance_symmetry.csv")
)

if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
source("utils.R")

cpp_environment <- new.env(parent = globalenv())
Rcpp::sourceCpp(
  code = '
    #include <Rcpp.h>
    #include <algorithm>
    #include <cmath>
    using namespace Rcpp;

    // [[Rcpp::export]]
    NumericMatrix symmetric_vmf_geodesic(NumericMatrix x) {
      const int n = x.nrow();
      const int p = x.ncol();
      NumericMatrix out(n, n);
      for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
          double inner = 0.0;
          for (int k = 0; k < p; ++k) inner += x(i, k) * x(j, k);
          inner = std::min(1.0, std::max(-1.0, inner));
          const double value = std::acos(inner);
          out(i, j) = value;
          out(j, i) = value;
        }
      }
      return out;
    }

    // [[Rcpp::export]]
    NumericMatrix symmetric_hvmf_geodesic(NumericMatrix x) {
      const int n = x.nrow();
      const int p = x.ncol();
      NumericMatrix out(n, n);
      for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
          double minkowski = -x(i, 0) * x(j, 0);
          for (int k = 1; k < p; ++k) minkowski += x(i, k) * x(j, k);
          const double value = std::acosh(std::max(-minkowski, 1.0));
          out(i, j) = value;
          out(j, i) = value;
        }
      }
      return out;
    }
  ',
  env = cpp_environment,
  showOutput = FALSE,
  verbose = FALSE
)

time_repeated <- function(fun, repetitions, inner_iterations) {
  values <- numeric(repetitions)
  answer <- NULL
  for (i in seq_len(repetitions)) {
    values[[i]] <- system.time({
      for (j in seq_len(inner_iterations)) answer <- fun()
    })[["elapsed"]] / inner_iterations
  }
  list(times = values, answer = answer)
}

results <- list()
result_index <- 1L
for (model in c("vmf", "hvmf")) {
  for (d in dimensions) {
    for (n in sizes) {
      set.seed(910000L + 1000L * d + n + if (model == "hvmf") 500000L else 0L)
      if (model == "vmf") {
        x <- rotasym::r_vMF(n, mu = c(1, rep.int(0, d)), kappa = d)
        x <- x / sqrt(rowSums(x^2))
        full <- function() {
          inner <- pmin(pmax(x %*% t(x), -1), 1)
          acos(inner)
        }
        symmetric <- function() cpp_environment$symmetric_vmf_geodesic(x)
      } else {
        x <- rhvmf_polar(n, mu = c(sqrt(2), 1, rep.int(0, d - 1L)), kappa = d)
        x <- normalize_hvmf_hq_data(x)
        full <- function() hvmf_distance_matrix_hq(x, x)
        symmetric <- function() cpp_environment$symmetric_hvmf_geodesic(x)
      }

      invisible(full())
      invisible(symmetric())
      full_timing <- time_repeated(full, repetitions, inner_iterations)
      symmetric_timing <- time_repeated(symmetric, repetitions, inner_iterations)
      max_error <- max(abs(full_timing$answer - symmetric_timing$answer))
      median_full <- median(full_timing$times)
      median_symmetric <- median(symmetric_timing$times)
      results[[result_index]] <- data.frame(
        model = model,
        d = d,
        n = n,
        repetitions = repetitions,
        inner_iterations = inner_iterations,
        full_blas_median_seconds = median_full,
        symmetric_cpp_median_seconds = median_symmetric,
        symmetric_speedup = median_full / median_symmetric,
        symmetric_gain = 1 - median_symmetric / median_full,
        max_abs_error = max_error
      )
      result_index <- result_index + 1L
    }
  }
}

output <- do.call(rbind, results)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(output, output_path, row.names = FALSE)
print(output, row.names = FALSE)
