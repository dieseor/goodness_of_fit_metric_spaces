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
Rcpp::NumericMatrix cpp_profile_lookup_tensor_local_polynomial(
    Rcpp::NumericMatrix values,
    Rcpp::NumericVector t_grid,
    Rcpp::NumericVector geometry_grid,
    Rcpp::NumericVector kappa_grid,
    Rcpp::NumericVector t_query,
    Rcpp::NumericVector geometry_query,
    Rcpp::NumericVector kappa_query,
    const int stencil_size = 6) {
  const int n_t = t_grid.size();
  const int n_geometry = geometry_grid.size();
  const int n_kappa = kappa_grid.size();
  const R_xlen_t n_grid =
    static_cast<R_xlen_t>(n_t) * n_geometry * n_kappa;
  const R_xlen_t n_query = t_query.size();
  const int n_fields = values.ncol();
  if (stencil_size != 4 && stencil_size != 6) {
    Rcpp::stop("`stencil_size` must be either 4 or 6.");
  }
  if (n_t < stencil_size || n_geometry < stencil_size ||
      n_kappa < stencil_size) {
    Rcpp::stop("Every lookup grid must contain at least `stencil_size` nodes.");
  }
  if (values.nrow() != n_grid || n_fields < 1) {
    Rcpp::stop("`values` has incompatible lookup-table dimensions.");
  }
  if (geometry_query.size() != n_query || kappa_query.size() != n_query) {
    Rcpp::stop("Lookup query vectors must have equal lengths.");
  }

  auto validate_grid = [](const Rcpp::NumericVector& grid,
                          const char* label) {
    for (R_xlen_t i = 1; i < grid.size(); ++i) {
      if (!R_finite(grid[i]) || !(grid[i] > grid[i - 1])) {
        Rcpp::stop("%s must be finite and strictly increasing.", label);
      }
    }
    if (!R_finite(grid[0])) {
      Rcpp::stop("%s must be finite and strictly increasing.", label);
    }
  };
  validate_grid(t_grid, "`t_grid`");
  validate_grid(geometry_grid, "`geometry_grid`");
  validate_grid(kappa_grid, "`kappa_grid`");

  auto stencil = [stencil_size](const Rcpp::NumericVector& grid,
                                const double value,
                                std::vector<int>& indices,
                                std::vector<double>& weights) {
    const int n = grid.size();
    if (!R_finite(value) || value < grid[0] || value > grid[n - 1]) {
      Rcpp::stop("Lookup query lies outside its certified grid.");
    }
    const auto upper = std::upper_bound(grid.begin(), grid.end(), value);
    int left = static_cast<int>(upper - grid.begin()) - 1;
    left = std::max(0, std::min(left, n - 2));
    int start = left - (stencil_size / 2 - 1);
    start = std::max(0, std::min(start, n - stencil_size));
    for (int i = 0; i < stencil_size; ++i) {
      indices[i] = start + i;
      const double node = grid[start + i];
      double weight = 1.0;
      for (int j = 0; j < stencil_size; ++j) {
        if (j == i) continue;
        weight *= (value - grid[start + j]) /
          (node - grid[start + j]);
      }
      weights[i] = weight;
    }
  };

  Rcpp::NumericMatrix out(n_query, n_fields);
  std::vector<int> t_indices(stencil_size);
  std::vector<int> geometry_indices(stencil_size);
  std::vector<int> kappa_indices(stencil_size);
  std::vector<double> t_weights(stencil_size);
  std::vector<double> geometry_weights(stencil_size);
  std::vector<double> kappa_weights(stencil_size);
  for (R_xlen_t query = 0; query < n_query; ++query) {
    stencil(t_grid, t_query[query], t_indices, t_weights);
    stencil(
      geometry_grid,
      geometry_query[query],
      geometry_indices,
      geometry_weights
    );
    stencil(kappa_grid, kappa_query[query], kappa_indices, kappa_weights);
    for (int field = 0; field < n_fields; ++field) {
      long double answer = 0.0L;
      for (int kappa_offset = 0; kappa_offset < stencil_size;
           ++kappa_offset) {
        const int kappa_index = kappa_indices[kappa_offset];
        for (int geometry_offset = 0; geometry_offset < stencil_size;
             ++geometry_offset) {
          const int geometry_index = geometry_indices[geometry_offset];
          const double outer_weight =
            kappa_weights[kappa_offset] *
            geometry_weights[geometry_offset];
          for (int t_offset = 0; t_offset < stencil_size; ++t_offset) {
            const int t_index = t_indices[t_offset];
            const R_xlen_t flat = t_index +
              static_cast<R_xlen_t>(n_t) *
              (geometry_index +
               static_cast<R_xlen_t>(n_geometry) * kappa_index);
            answer += static_cast<long double>(
              values(flat, field) * outer_weight * t_weights[t_offset]
            );
          }
        }
      }
      out(query, field) = static_cast<double>(answer);
    }
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

// Cache-friendly implementation of the same fast sample-based KS/CvM kernel
// as cpp_fast_sample_ks_cvm_stats().  It preserves the order of the sums
// within every replicate, but traverses replicates together after reading each
// contiguous column of centered_weights and correction.  This is the
// production default; cpp_fast_sample_ks_cvm_stats() remains the explicit
// historical fallback for paired validation.
// [[Rcpp::export]]
Rcpp::List cpp_fast_sample_ks_cvm_stats_contiguous_double(
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
  std::vector<double> running(static_cast<std::size_t>(n_reps));
  std::vector<double> center_cvm_sum(static_cast<std::size_t>(n_reps));

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

    std::fill(running.begin(), running.end(), 0.0);
    if (compute_cvm) {
      std::fill(center_cvm_sum.begin(), center_cvm_sum.end(), 0.0);
    }
    int threshold = 0;

    for (int sorted_index = 0; sorted_index < n; ++sorted_index) {
      const int observation = obs_order_matrix(center, sorted_index) - 1;
      if (observation < 0 || observation >= n) {
        Rcpp::stop("`obs_order_matrix` contains an invalid index.");
      }
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
            (running[static_cast<std::size_t>(rep)] - correction_column[rep]) /
            root_n;
          if (compute_ks) {
            const double absolute_process = std::abs(process);
            if (absolute_process > ks[rep]) {
              ks[rep] = absolute_process;
            }
          }
          if (compute_cvm) {
            center_cvm_sum[static_cast<std::size_t>(rep)] += process * process;
          }
        }
        ++threshold;
      }
    }
    if (threshold != n) {
      Rcpp::stop("`tie_end_matrix` must contain valid nondecreasing endpoints.");
    }
    if (compute_cvm) {
      for (int rep = 0; rep < n_reps; ++rep) {
        cvm_sum[rep] += center_cvm_sum[static_cast<std::size_t>(rep)];
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

// Joint profile for d((x, s), (omega, t)) =
// 0.5 * { arccos(x'omega) / pi + |s - t| }.  The conditional spherical
// profiles are represented by a row of Legendre coefficients at each time
// quadrature node.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_sunspots_joint_profile_block(
    Rcpp::NumericMatrix radii,
    Rcpp::NumericVector rho,
    Rcpp::NumericVector center_s,
    Rcpp::NumericVector time_nodes,
    Rcpp::NumericVector time_weights,
    Rcpp::NumericMatrix coefficients) {
  const int n_rows = radii.nrow();
  const int n_cols = radii.ncol();
  const int n_time = time_nodes.size();
  const int n_coefficients = coefficients.ncol();
  const int l_max = n_coefficients - 1;
  const double pi = 3.141592653589793238462643383279502884;

  if (rho.size() != n_rows || center_s.size() != n_rows) {
    Rcpp::stop("`rho` and `center_s` must have length nrow(`radii`).");
  }
  if (time_weights.size() != n_time || coefficients.nrow() != n_time) {
    Rcpp::stop("Time nodes, time weights, and coefficients have incompatible dimensions.");
  }
  if (n_coefficients < 1) {
    Rcpp::stop("`coefficients` must have at least one column.");
  }

  Rcpp::NumericMatrix out(n_rows, n_cols);
  Rcpp::NumericVector rho_clipped = Rcpp::clone(rho);
  for (int row = 0; row < n_rows; ++row) {
    if (!R_FINITE(center_s[row]) || !R_FINITE(rho_clipped[row])) {
      Rcpp::stop("`rho` and `center_s` must be finite.");
    }
    rho_clipped[row] = clip_unit(rho_clipped[row]);
  }

  // P_ell(rho) is shared by all temporal quadrature nodes.
  Rcpp::NumericMatrix rho_legendre(n_rows, l_max + 1);
  for (int row = 0; row < n_rows; ++row) {
    rho_legendre(row, 0) = 1.0;
    if (l_max >= 1) rho_legendre(row, 1) = rho_clipped[row];
  }
  for (int ell = 1; ell < l_max; ++ell) {
    const double factor_one = static_cast<double>(2 * ell + 1);
    const double factor_two = static_cast<double>(ell);
    const double denominator = static_cast<double>(ell + 1);
    for (int row = 0; row < n_rows; ++row) {
      rho_legendre(row, ell + 1) =
        (factor_one * rho_clipped[row] * rho_legendre(row, ell) -
          factor_two * rho_legendre(row, ell - 1)) / denominator;
    }
  }

  for (int node = 0; node < n_time; ++node) {
    if (!R_FINITE(time_nodes[node]) || !R_FINITE(time_weights[node])) {
      Rcpp::stop("Time nodes and weights must be finite.");
    }
    const double weight = time_weights[node];
    if (weight == 0.0) continue;
    for (int col = 0; col < n_cols; ++col) {
      for (int row = 0; row < n_rows; ++row) {
        const double radius = radii(row, col);
        if (!R_FINITE(radius)) Rcpp::stop("`radii` must be finite.");
        const double spherical_fraction = std::min(
          1.0,
          std::max(0.0, 2.0 * radius - std::abs(center_s[row] - time_nodes[node]))
        );
        const double x = std::cos(pi * spherical_fraction);
        double cdf = (x + 1.0) / 2.0;
        if (l_max >= 1) {
          double p_x_previous = 1.0;
          double p_x_current = x;
          for (int ell = 1; ell <= l_max; ++ell) {
            const double factor_one = static_cast<double>(2 * ell + 1);
            const double factor_two = static_cast<double>(ell);
            const double denominator = static_cast<double>(ell + 1);
            const double p_x_next =
              (factor_one * x * p_x_current - factor_two * p_x_previous) / denominator;
            cdf += coefficients(node, ell) * rho_legendre(row, ell) *
              (p_x_next - p_x_previous) / (2.0 * factor_one);
            p_x_previous = p_x_current;
            p_x_current = p_x_next;
          }
        }
        out(row, col) += weight * (1.0 - clip_probability(cdf));
      }
    }
  }

  for (int col = 0; col < n_cols; ++col) {
    for (int row = 0; row < n_rows; ++row) {
      out(row, col) = clip_probability(out(row, col));
    }
  }
  return out;
}
