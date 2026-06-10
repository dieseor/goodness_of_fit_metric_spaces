# Logistic Gaussian Quadratic-Form Benchmark

This benchmark compares candidate methods for evaluating weighted sums of noncentral chi-squared variables arising in the logistic Gaussian distance-profile calculations.

Reference rule:
- use `imhof` as the reference when both `imhof` and `farebrother` are valid;
- otherwise fall back to the available exact method among `imhof` and `farebrother`.

Farebrother slow-case threshold: 3.679000 seconds (empirical 95th percentile).

Headline findings:
- `farebrother`: median 0.000000 s, 95th percentile 3.679000 s, max abs. error 1.000e+00, 389 non-ok cases.
- `imhof`: median 0.002000 s, 95th percentile 0.087850 s, used as exact reference.
- `davies`: 419 non-ok cases out of 1084, so it is not suitable as the default exact route without explicit failure handling.
- `sphunif_sw`: median abs. error 8.886e-04; `sphunif_hbe`: median abs. error 1.114e-04.
- Best dispatcher among the benchmarked simple rules: `cond_gt_1e4_or_dim5_delta10_qhigh`, with total benchmark time 78.60 s, 95th percentile 0.09570 s, and maximum absolute deviation 1.015e-04 versus the benchmark reference.

Interpretation:
- `farebrother` is the fastest exact method whenever it succeeds and remains numerically stable.
- `imhof` is slower but provides a robust exact fallback and a reliable benchmark reference.
- `davies` can fail with non-zero `ifault` even in moderate cases, so it should not be the default route.
- `SW` and `HBE` are useful only as fast approximations, not as exact replacements.
- The benchmark favours the dispatcher `cond_gt_1e4_or_dim5_delta10_qhigh`: send clearly ill-conditioned or strongly noncentral medium/high-dimensional cases directly to `imhof`, and otherwise try `farebrother` first.
