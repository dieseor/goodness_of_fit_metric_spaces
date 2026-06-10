# Benchmark And Validation Outputs

This directory stores artifacts whose purpose is **testing, validation, and numerical benchmarking**.

These files are not primary scientific outputs. They are diagnostic products used to answer questions such as:

- whether a deterministic formula matches a trusted reference;
- whether Monte Carlo remains positive semidefinite on a given grid;
- how numerical error changes with tuning parameters;
- how computational cost compares across methods.

## Subdirectories

- `s1_validation`
  S^1 validation artifacts for the vMF Gaussian-process work.
  These files check the deterministic S^1 route, Monte Carlo behavior, and covariance-matrix diagnostics.

- `s2_validation`
  S^2 validation artifacts for the exact-integral benchmarks in the simple-null vMF case.
  These include:
  - scalar joint-probability benchmarks;
  - closed-form special-case benchmarks;
  - covariance-matrix benchmarks.

- `logistic_gaussian_quadform`
  Benchmarks for the weighted noncentral chi-squared evaluations used by the logistic Gaussian distance-profile code.
  These files compare exact and approximate methods, identify failure regimes, and justify the chosen fallback strategy.

## What Should Go Here

Store here:

- benchmark CSV files;
- diagnostic plots comparing methods;
- numerical-validation outputs used to justify implementation choices.

Do not store here:

- final plots intended for the paper;
- exploratory figures that are not tied to a test or validation question;
- routine production outputs from the main simulation pipelines.
