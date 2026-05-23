# S1 Validation Outputs

This directory contains **validation artifacts** for the S^1 vMF work.

## Files

- `s1_joint_probability_benchmark.csv`
  Scalar benchmark for joint probabilities on S^1.
  Purpose: compare Monte Carlo against the deterministic S^1 route on individual probabilities.

- `s1_joint_probability_error.png`
- `s1_joint_probability_time.png`
  Visual summaries of scalar error and timing for the previous benchmark.

- `s1_covariance_benchmark.csv`
  Matrix-level benchmark for the baseline S^1 covariance comparison.
  Purpose: compare the Monte Carlo covariance matrix against the deterministic S^1 covariance matrix on the same grids.

- `s1_covariance_error.png`
- `s1_covariance_time.png`
  Visual summaries of matrix error and timing for the baseline covariance benchmark.

- `s1_covariance_benchmark_realistic_mc.csv`
  Matrix-level benchmark with Monte Carlo sizes closer to the main experiments.
  Purpose: check whether large Monte Carlo sizes recover positive semidefiniteness and how fast the matrix error decreases.

- `s1_covariance_error_realistic_mc.png`
- `s1_covariance_time_realistic_mc.png`
  Visual summaries for the realistic-MC covariance benchmark.

- `s1_covariance_reference_benchmark_angles6_t10.csv`
  Matrix-level benchmark against an external reference built entrywise by adaptive quadrature.
  Purpose: validate the deterministic S^1 covariance route against a more independent numerical reference.

- `s1_covariance_reference_error_angles6_t10.png`
- `s1_covariance_reference_time_angles6_t10.png`
  Visual summaries for the external-reference covariance benchmark.
