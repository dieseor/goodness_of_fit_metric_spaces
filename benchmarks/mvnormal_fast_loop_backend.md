# Multivariate-normal fast loop: R versus C++

Date: 2026-07-23

The comparison uses the same cached correction matrices and the same fused
KS/CvM algorithm in both backends.  Timings below are warm medians; compilation
is excluded.  Backend order was alternated and all compared KS and CvM vectors
were bitwise identical.

## Isolated fused loop

| dimension | n | B | R (s) | C++ (s) | speedup | gain |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 50 | 300 | 0.085 | 0.060 | 1.42 | 29.4% |
| 10 | 50 | 300 | 0.100 | 0.074 | 1.35 | 26.0% |
| 2 | 100 | 300 | 0.288 | 0.239 | 1.21 | 17.0% |
| 10 | 100 | 300 | 0.349 | 0.294 | 1.19 | 15.8% |
| 2 | 200 | 300 | 1.056 | 0.972 | 1.09 | 8.0% |
| 10 | 200 | 200 | 0.820 | 0.753 | 1.09 | 8.2% |

## Complete fast-multiplier contrast

For \(n=100\), \(B=500\), both requested statistics, automatic caching, and
500 auxiliary derivative draws:

| dimension | R (s) | C++ (s) | speedup | gain |
|---:|---:|---:|---:|---:|
| 2 | 1.160 | 1.061 | 1.09 | 8.5% |
| 10 | 1.441 | 1.326 | 1.09 | 8.0% |

At \(B=200\), the corresponding complete-test gains were 4.2% and 3.2%:
fixed preparation costs dilute the loop improvement when the bootstrap load is
small.  The first source compilation took about three seconds on the test
machine and is reported separately from warm execution.

## Exactness checks

- 500 reproducible randomized differential cases, with dimensions 2 and 10,
  random block sizes, and deliberately introduced ties: exact KS and CvM.
- Complete contrasts with common data, multipliers, derivative samples, and
  seeds: exact observed results, bootstrap distributions, inference, and
  decisions.
- Cached, uncached, fused, unfused, sequential, and parallel routes: exact.

