#include <Rcpp.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::plugins(cpp17)]]

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

namespace {

inline double clip_unit(const double x) {
  return std::min(1.0, std::max(-1.0, x));
}

inline double clip_probability(const double x) {
  return std::min(1.0, std::max(0.0, x));
}

Rcpp::NumericMatrix legendre_matrix_impl(const Rcpp::NumericVector& x,
                                         const int l_max) {
  const int n = x.size();
  Rcpp::NumericMatrix out(n, l_max + 1);
  for (int i = 0; i < n; ++i) {
    out(i, 0) = 1.0;
  }
  if (l_max == 0) {
    return out;
  }
  for (int i = 0; i < n; ++i) {
    out(i, 1) = x[i];
  }
  for (int ell = 1; ell < l_max; ++ell) {
    const double factor_one = static_cast<double>(2 * ell + 1);
    const double factor_two = static_cast<double>(ell);
    const double denominator = static_cast<double>(ell + 1);
    for (int i = 0; i < n; ++i) {
      const double x_i = x[i];
      const double first_product = factor_one * x_i;
      const double first_term = first_product * out(i, ell);
      const double second_term = factor_two * out(i, ell - 1);
      const double numerator = first_term - second_term;
      out(i, ell + 1) = numerator / denominator;
    }
  }
  return out;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dp_legendre_matrix(Rcpp::NumericVector x,
                                           const int l_max) {
  if (l_max < 0) {
    Rcpp::stop("`l_max` must be a nonnegative integer.");
  }
  Rcpp::NumericVector x_clipped = Rcpp::clone(x);
  for (R_xlen_t i = 0; i < x_clipped.size(); ++i) {
    x_clipped[i] = clip_unit(x_clipped[i]);
  }
  return legendre_matrix_impl(x_clipped, l_max);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dp_projection_cdf_legendre_matrix(
    Rcpp::NumericMatrix x_matrix,
    Rcpp::NumericVector r,
    Rcpp::NumericVector coefficients,
    const bool enforce_bounds = true) {
  const int n_rows = x_matrix.nrow();
  const int n_cols = x_matrix.ncol();
  const int l_max = coefficients.size() - 1;
  if (r.size() != n_rows) {
    Rcpp::stop("`r` must have length nrow(`x_matrix`).");
  }
  if (l_max < 0) {
    Rcpp::stop("`coefficients` cannot be empty.");
  }

  Rcpp::NumericVector r_clipped = Rcpp::clone(r);
  for (int i = 0; i < n_rows; ++i) {
    r_clipped[i] = clip_unit(r_clipped[i]);
  }
  Rcpp::NumericVector x_flat = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(x_matrix));
  for (R_xlen_t i = 0; i < x_flat.size(); ++i) {
    x_flat[i] = clip_unit(x_flat[i]);
  }

  Rcpp::NumericMatrix out(n_rows, n_cols);
  for (int col = 0; col < n_cols; ++col) {
    for (int row = 0; row < n_rows; ++row) {
      const int flat = row + n_rows * col;
      out(row, col) = (x_flat[flat] + 1.0) / 2.0;
    }
  }

  Rcpp::NumericVector p_r_previous(n_rows, 1.0);
  Rcpp::NumericVector p_r_current = Rcpp::clone(r_clipped);
  Rcpp::NumericVector p_r_next(n_rows);
  Rcpp::NumericVector p_x_previous(x_flat.size(), 1.0);
  Rcpp::NumericVector p_x_current = Rcpp::clone(x_flat);
  Rcpp::NumericVector p_x_next(x_flat.size());

  for (int ell = 1; ell <= l_max; ++ell) {
    const double recurrence_factor_one = static_cast<double>(2 * ell + 1);
    const double recurrence_factor_two = static_cast<double>(ell);
    const double recurrence_denominator = static_cast<double>(ell + 1);
    for (int row = 0; row < n_rows; ++row) {
      const double first_product = recurrence_factor_one * r_clipped[row];
      const double first_term = first_product * p_r_current[row];
      const double second_term = recurrence_factor_two * p_r_previous[row];
      const double numerator = first_term - second_term;
      p_r_next[row] = numerator / recurrence_denominator;
    }
    for (R_xlen_t flat = 0; flat < x_flat.size(); ++flat) {
      const double first_product = recurrence_factor_one * x_flat[flat];
      const double first_term = first_product * p_x_current[flat];
      const double second_term = recurrence_factor_two * p_x_previous[flat];
      const double numerator = first_term - second_term;
      p_x_next[flat] = numerator / recurrence_denominator;
    }

    const double denominator = 2.0 * static_cast<double>(2 * ell + 1);
    for (int col = 0; col < n_cols; ++col) {
      for (int row = 0; row < n_rows; ++row) {
        const double weight = coefficients[ell] * p_r_current[row];
        const int flat = row + n_rows * col;
        const double difference = p_x_next[flat] - p_x_previous[flat];
        const double basis = difference / denominator;
        const double increment = basis * weight;
        out(row, col) = out(row, col) + increment;
      }
    }

    std::swap(p_r_previous, p_r_current);
    std::swap(p_r_current, p_r_next);
    std::swap(p_x_previous, p_x_current);
    std::swap(p_x_current, p_x_next);
  }

  if (enforce_bounds) {
    for (R_xlen_t i = 0; i < out.size(); ++i) {
      out[i] = clip_probability(out[i]);
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_dp_normal_profile(Rcpp::NumericVector omega,
                                          Rcpp::NumericVector t_values,
                                          const double mu,
                                          const double sigma) {
  const R_xlen_t n = std::max(omega.size(), t_values.size());
  Rcpp::NumericVector out(n);
  if (omega.size() == 0 || t_values.size() == 0) {
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    const double t = t_values[i % t_values.size()];
    if (t > 0.0) {
      const double omega_i = omega[i % omega.size()];
      const double upper_numerator = omega_i + t - mu;
      const double lower_numerator = omega_i - t - mu;
      const double upper = upper_numerator / sigma;
      const double lower = lower_numerator / sigma;
      const double upper_prob = R::pnorm5(upper, 0.0, 1.0, 1, 0);
      const double lower_prob = R::pnorm5(lower, 0.0, 1.0, 1, 0);
      out[i] = upper_prob - lower_prob;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dp_normal_profile_matrix(
    Rcpp::NumericVector omega,
    Rcpp::NumericMatrix t_matrix,
    const double mu,
    const double sigma) {
  const int n_rows = t_matrix.nrow();
  const int n_cols = t_matrix.ncol();
  if (omega.size() != n_rows) {
    Rcpp::stop("`omega` must have length nrow(`t_matrix`).");
  }
  Rcpp::NumericMatrix out(n_rows, n_cols);
  for (int col = 0; col < n_cols; ++col) {
    for (int row = 0; row < n_rows; ++row) {
      const double t = t_matrix(row, col);
      if (t > 0.0) {
        const double upper_numerator = omega[row] + t - mu;
        const double lower_numerator = omega[row] - t - mu;
        const double upper = upper_numerator / sigma;
        const double lower = lower_numerator / sigma;
        const double upper_prob = R::pnorm5(upper, 0.0, 1.0, 1, 0);
        const double lower_prob = R::pnorm5(lower, 0.0, 1.0, 1, 0);
        out(row, col) = upper_prob - lower_prob;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dp_weighted_sample_profile_rows(
    Rcpp::IntegerMatrix order_matrix,
    Rcpp::IntegerMatrix rank_matrix,
    Rcpp::NumericVector normalized_weights) {
  const int n_rows = order_matrix.nrow();
  const int n_total = order_matrix.ncol();
  if (rank_matrix.nrow() != n_rows || rank_matrix.ncol() != n_total) {
    Rcpp::stop("Order and rank matrices have incompatible dimensions.");
  }
  Rcpp::NumericMatrix cumulative(n_rows, n_total);
  Rcpp::NumericMatrix out(n_rows, n_total);
  for (int row = 0; row < n_rows; ++row) {
    double running = 0.0;
    for (int col = 0; col < n_total; ++col) {
      const int index = order_matrix(row, col) - 1;
      running = running + normalized_weights[index];
      cumulative(row, col) = running;
    }
    for (int col = 0; col < n_total; ++col) {
      const int rank = rank_matrix(row, col) - 1;
      out(row, col) = cumulative(row, rank) / static_cast<double>(n_total);
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dp_weighted_sample_profile_linear(
    Rcpp::IntegerMatrix order_matrix,
    Rcpp::IntegerVector rank_linear_index,
    Rcpp::NumericVector normalized_weights) {
  const int n_rows = order_matrix.nrow();
  const int n_total = order_matrix.ncol();
  if (rank_linear_index.size() != static_cast<R_xlen_t>(n_rows) * n_total) {
    Rcpp::stop("The linear rank index has incompatible dimensions.");
  }
  Rcpp::NumericMatrix cumulative(n_rows, n_total);
  Rcpp::NumericMatrix out(n_rows, n_total);
  for (int row = 0; row < n_rows; ++row) {
    double running = 0.0;
    for (int col = 0; col < n_total; ++col) {
      const int index = order_matrix(row, col) - 1;
      running = running + normalized_weights[index];
      cumulative(row, col) = running;
    }
  }
  for (R_xlen_t i = 0; i < rank_linear_index.size(); ++i) {
    out[i] = cumulative[rank_linear_index[i] - 1] / static_cast<double>(n_total);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_fast_sample_ks_cvm_stats(
    Rcpp::NumericMatrix centered_weights,
    Rcpp::NumericMatrix score_block,
    Rcpp::IntegerMatrix obs_order_matrix,
    Rcpp::IntegerMatrix tie_end_matrix,
    Rcpp::NumericMatrix correction_matrix,
    const double scale_factor,
    const bool compute_ks = true,
    const bool compute_cvm = true) {
  const int n_reps = centered_weights.nrow();
  const int n = centered_weights.ncol();
  const int n_centers = obs_order_matrix.nrow();
  const int n_parameters = score_block.ncol();

  if (!compute_ks && !compute_cvm) {
    Rcpp::stop("At least one of `compute_ks` and `compute_cvm` must be true.");
  }
  if (score_block.nrow() != n_reps) {
    Rcpp::stop("`score_block` has an incompatible number of rows.");
  }
  if (obs_order_matrix.ncol() != n ||
      tie_end_matrix.nrow() != n_centers ||
      tie_end_matrix.ncol() != n) {
    Rcpp::stop("Ordering and tie matrices have incompatible dimensions.");
  }
  if (correction_matrix.nrow() != n_centers * n ||
      correction_matrix.ncol() != n_parameters) {
    Rcpp::stop("The correction matrix has incompatible dimensions.");
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
  const int ldb = correction_matrix.nrow();
  const int ldc = n_reps;
  const double root_n = std::sqrt(static_cast<double>(n));
  const double cvm_denominator =
    static_cast<double>(n) * static_cast<double>(n);

  for (int center = 0; center < n_centers; ++center) {
    const double* correction_center = correction_matrix.begin() +
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
        running += static_cast<long double>(
          centered_weights(rep, observation)
        );
        cumulative[static_cast<std::size_t>(sorted_index)] =
          static_cast<double>(running);
      }

      if (compute_ks && compute_cvm) {
        long double center_cvm_sum = 0.0L;
        for (int threshold_index = 0; threshold_index < n;
             ++threshold_index) {
          const int tie_end =
            tie_end_matrix(center, threshold_index) - 1;
          if (tie_end < 0 || tie_end >= n) {
            Rcpp::stop("`tie_end_matrix` contains an invalid index.");
          }
          const double empirical =
            cumulative[static_cast<std::size_t>(tie_end)];
          const double difference =
            empirical - correction(rep, threshold_index);
          const double scaled = scale_factor * difference;
          const double process = scaled / root_n;
          const double absolute_process = std::abs(process);
          if (absolute_process > ks[rep]) {
            ks[rep] = absolute_process;
          }
          center_cvm_sum +=
            static_cast<long double>(process * process);
        }
        cvm_sum[rep] += static_cast<double>(center_cvm_sum);
      } else if (compute_ks) {
        for (int threshold_index = 0; threshold_index < n;
             ++threshold_index) {
          const int tie_end =
            tie_end_matrix(center, threshold_index) - 1;
          if (tie_end < 0 || tie_end >= n) {
            Rcpp::stop("`tie_end_matrix` contains an invalid index.");
          }
          const double empirical =
            cumulative[static_cast<std::size_t>(tie_end)];
          const double difference =
            empirical - correction(rep, threshold_index);
          const double scaled = scale_factor * difference;
          const double process = scaled / root_n;
          const double absolute_process = std::abs(process);
          if (absolute_process > ks[rep]) {
            ks[rep] = absolute_process;
          }
        }
      } else {
        long double center_cvm_sum = 0.0L;
        for (int threshold_index = 0; threshold_index < n;
             ++threshold_index) {
          const int tie_end =
            tie_end_matrix(center, threshold_index) - 1;
          if (tie_end < 0 || tie_end >= n) {
            Rcpp::stop("`tie_end_matrix` contains an invalid index.");
          }
          const double empirical =
            cumulative[static_cast<std::size_t>(tie_end)];
          const double difference =
            empirical - correction(rep, threshold_index);
          const double scaled = scale_factor * difference;
          const double process = scaled / root_n;
          center_cvm_sum +=
            static_cast<long double>(process * process);
        }
        cvm_sum[rep] += static_cast<double>(center_cvm_sum);
      }
    }
  }

  Rcpp::NumericVector cvm(n_reps, 0.0);
  if (compute_cvm) {
    for (int rep = 0; rep < n_reps; ++rep) {
      cvm[rep] = cvm_sum[rep] / cvm_denominator;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("ks") = ks,
    Rcpp::Named("cvm") = cvm,
    Rcpp::Named("cvm_sum") = cvm_sum
  );
}
