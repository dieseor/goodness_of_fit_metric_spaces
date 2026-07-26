# Logistic Gaussian Quadratic-Form Benchmark

This benchmark compares candidate methods for evaluating weighted sums of noncentral chi-squared variables arising in the logistic Gaussian distance-profile calculations.

Reference rule:
- use `imhof` as the reference when both `imhof` and `farebrother` are valid;
- otherwise fall back to the available exact method among `imhof` and `farebrother`.

Farebrother slow-case threshold: 5.048200 seconds (empirical 95th percentile).

Headline findings:
- `farebrother`: median 0.000000 s, 95th percentile 5.048200 s, max abs. error 1.000e+00, 40 non-ok cases.
- `imhof`: median 0.038500 s, 95th percentile 0.130050 s, used as exact reference.
- `davies`: 67 non-ok cases out of 100, so it is not suitable as the default exact route without explicit failure handling.
- `sphunif_sw`: median abs. error 7.773e-04; `sphunif_hbe`: median abs. error 3.448e-05.
- Best dispatcher among the benchmarked simple rules: `cond_gt_1e4`, with total benchmark time 1.77 s, 95th percentile 0.09910 s, and maximum absolute deviation 1.328e-06 versus the benchmark reference.

Interpretation:
- `farebrother` is the fastest exact method whenever it succeeds and remains numerically stable.
- `imhof` is slower but provides a robust exact fallback and a reliable benchmark reference.
- `davies` can fail with non-zero `ifault` even in moderate cases, so it should not be the default route.
- `SW` and `HBE` are useful only as fast approximations, not as exact replacements.
- The benchmark favours the dispatcher `cond_gt_1e4`: send clearly ill-conditioned or strongly noncentral medium/high-dimensional cases directly to `imhof`, and otherwise try `farebrother` first.
