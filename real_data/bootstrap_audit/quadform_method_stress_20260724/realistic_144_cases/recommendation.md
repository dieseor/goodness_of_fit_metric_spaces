# Logistic Gaussian Quadratic-Form Benchmark

This benchmark compares candidate methods for evaluating weighted sums of noncentral chi-squared variables arising in the logistic Gaussian distance-profile calculations.

Reference rule:
- use `imhof` as the reference when both `imhof` and `farebrother` are valid;
- otherwise fall back to the available exact method among `imhof` and `farebrother`.

Farebrother slow-case threshold: 0.001000 seconds (empirical 95th percentile).

Headline findings:
- `farebrother`: median 0.000000 s, 95th percentile 0.001000 s, max abs. error 4.624e-05, 0 non-ok cases.
- `imhof`: median 0.063000 s, 95th percentile 0.085850 s, used as exact reference.
- `davies`: 60 non-ok cases out of 144, so it is not suitable as the default exact route without explicit failure handling.
- `sphunif_sw`: median abs. error 2.273e-05; `sphunif_hbe`: median abs. error 2.433e-05.
- Best dispatcher among the benchmarked simple rules: `fallback_after_farebrother`, with total benchmark time 0.04 s, 95th percentile 0.00100 s, and maximum absolute deviation 4.624e-05 versus the benchmark reference.

Interpretation:
- `farebrother` is the fastest exact method whenever it succeeds and remains numerically stable.
- `imhof` is slower but provides a robust exact fallback and a reliable benchmark reference.
- `davies` can fail with non-zero `ifault` even in moderate cases, so it should not be the default route.
- `SW` and `HBE` are useful only as fast approximations, not as exact replacements.
- The benchmark favours the dispatcher `fallback_after_farebrother`: send clearly ill-conditioned or strongly noncentral medium/high-dimensional cases directly to `imhof`, and otherwise try `farebrother` first.
