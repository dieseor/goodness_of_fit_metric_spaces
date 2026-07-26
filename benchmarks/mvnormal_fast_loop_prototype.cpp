#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::plugins(cpp17)]]

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

// Prototype for benchmarking only.  It deliberately leaves the score matrix
// product in R and uses the same BLAS dgemm operation as `%*%` for each centre.
// The empirical weighted profile, tie selection, and KS/CvM reductions are
// fused so that KS and CvM do not construct the same intermediate matrices
// independently.

// [[Rcpp::export]]
Rcpp::List cpp_mvnormal_fast_sample_ks_cvm(
    Rcpp::NumericMatrix centered_weights,
    Rcpp::NumericMatrix score_block,
    Rcpp::IntegerMatrix obs_order_matrix,
    Rcpp::IntegerMatrix tie_end_matrix,
    Rcpp::NumericMatrix correction_cache,
    const double scale_factor) {
  const int n_reps = centered_weights.nrow();
  const int n = centered_weights.ncol();
  const int n_centers = obs_order_matrix.nrow();
  const int n_parameters = score_block.ncol();

  if (score_block.nrow() != n_reps) {
    Rcpp::stop("`score_block` has an incompatible number of rows.");
  }
  if (obs_order_matrix.ncol() != n ||
      tie_end_matrix.nrow() != n_centers ||
      tie_end_matrix.ncol() != n) {
    Rcpp::stop("Ordering and tie matrices have incompatible dimensions.");
  }
  if (correction_cache.nrow() != n_centers * n ||
      correction_cache.ncol() != n_parameters) {
    Rcpp::stop("The correction cache has incompatible dimensions.");
  }

  Rcpp::NumericVector ks(n_reps, 0.0);
  Rcpp::NumericVector cvm_sum(n_reps, 0.0);
  Rcpp::NumericMatrix correction(n_reps, n);
  std::vector<double> cumulative(static_cast<std::size_t>(n));

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta = 0.0;
  const int lda = n_reps;
  const int ldb = correction_cache.nrow();
  const int ldc = n_reps;
  const double root_n = std::sqrt(static_cast<double>(n));
  const double cvm_denominator = static_cast<double>(n) * static_cast<double>(n);

  for (int center = 0; center < n_centers; ++center) {
    const double* correction_center = correction_cache.begin() +
      static_cast<R_xlen_t>(center) * n;

    F77_CALL(dgemm)(
      &trans_n,
      &trans_t,
      &n_reps,
      &n,
      &n_parameters,
      &alpha,
      score_block.begin(),
      &lda,
      correction_center,
      &ldb,
      &beta,
      correction.begin(),
      &ldc FCONE FCONE
    );

    for (int rep = 0; rep < n_reps; ++rep) {
      long double running = 0.0L;
      for (int sorted_index = 0; sorted_index < n; ++sorted_index) {
        const int observation = obs_order_matrix(center, sorted_index) - 1;
        if (observation < 0 || observation >= n) {
          Rcpp::stop("`obs_order_matrix` contains an invalid index.");
        }
        running = running + static_cast<long double>(centered_weights(rep, observation));
        cumulative[static_cast<std::size_t>(sorted_index)] = static_cast<double>(running);
      }

      long double center_cvm_sum = 0.0L;
      for (int threshold_index = 0; threshold_index < n; ++threshold_index) {
        const int tie_end = tie_end_matrix(center, threshold_index) - 1;
        if (tie_end < 0 || tie_end >= n) {
          Rcpp::stop("`tie_end_matrix` contains an invalid index.");
        }
        const double empirical = cumulative[static_cast<std::size_t>(tie_end)];
        const double difference = empirical - correction(rep, threshold_index);
        const double scaled = scale_factor * difference;
        const double process = scaled / root_n;
        const double absolute_process = std::abs(process);
        if (absolute_process > ks[rep]) {
          ks[rep] = absolute_process;
        }
        center_cvm_sum = center_cvm_sum +
          static_cast<long double>(process * process);
      }
      cvm_sum[rep] = cvm_sum[rep] + static_cast<double>(center_cvm_sum);
    }
  }

  Rcpp::NumericVector cvm(n_reps);
  for (int rep = 0; rep < n_reps; ++rep) {
    cvm[rep] = cvm_sum[rep] / cvm_denominator;
  }

  return Rcpp::List::create(
    Rcpp::Named("ks") = ks,
    Rcpp::Named("cvm") = cvm
  );
}
