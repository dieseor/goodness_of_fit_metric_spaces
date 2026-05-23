# S2 Validation Outputs

This directory contains **validation artifacts** for the S^2 exact-integral benchmarks in the simple-null vMF case.

## Files

- `s2_joint_probability_benchmark_geodesic.csv`
  Scalar benchmark for generic S^2 joint probabilities using a high-accuracy integral reference.
  Purpose: compare the exact-integral route and Monte Carlo on representative geometric configurations.

- `s2_joint_probability_error_geodesic.png`
- `s2_joint_probability_time_geodesic.png`
  Visual summaries of scalar error and timing for the generic S^2 benchmark.

- `s2_joint_probability_closed_form_benchmark_geodesic.csv`
  Scalar benchmark on special S^2 cases with genuinely closed-form truth.
  Purpose: compare Monte Carlo and the exact-integral route against probabilities known without numerical integration.

- `s2_joint_probability_closed_form_error_geodesic.png`
- `s2_joint_probability_closed_form_time_geodesic.png`
  Visual summaries for the closed-form S^2 benchmark.

- `s2_covariance_benchmark_geodesic_kappa0p50.csv`
  Matrix-level covariance benchmark on S^2 using an exact-integral reference matrix.
  Purpose: compare the computational cost and matrix accuracy of Monte Carlo versus the integral route.

- `s2_covariance_error_geodesic_kappa0p50.png`
- `s2_covariance_time_geodesic_kappa0p50.png`
  Visual summaries for the S^2 covariance benchmark.
