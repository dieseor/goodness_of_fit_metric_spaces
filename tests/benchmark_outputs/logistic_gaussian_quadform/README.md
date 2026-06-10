# Logistic Gaussian Quadratic-Form Benchmark

This directory stores benchmark artifacts for the scalar distribution-function evaluations that appear in the logistic Gaussian distance-profile calculations.

The benchmark compares:

- `CompQuadForm::farebrother`
- `CompQuadForm::davies`
- `CompQuadForm::imhof`
- `sphunif::p_wschisq(method = "SW")`
- `sphunif::p_wschisq(method = "HBE")`

The benchmark mixes:

- structured stress cases with controlled spectra and noncentralities;
- realistic cases induced by the logistic Gaussian calibration scenarios already used in the repo;
- random stress cases intended to expose unstable regions.

Saved artifacts should include:

- `quadform_cases.csv`
- `quadform_results.csv`
- `quadform_method_summary.csv`
- `dispatcher_rule_summary.csv`
- `farebrother_regime_summary.csv`
- `farebrother_slow_cases.csv`
- `recommendation.md`
- `sessionInfo.txt`

For the subset of extreme cases where the benchmark could not certify an exact
reference from `farebrother`/`imhof`, a separate Monte Carlo reference study is
available. It estimates the target probability directly from simulated weighted
noncentral chi-squared samples and compares all methods against that Monte Carlo
reference.

Run with:

```bash
Rscript scripts/run_logistic_gaussian_quadform_benchmark.R --n_cores=12
Rscript scripts/run_logistic_gaussian_quadform_mc_reference.R --n_cores=12
```
