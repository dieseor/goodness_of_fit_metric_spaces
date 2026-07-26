# Optional C++ distance-profile backend: validation and benchmark

Date: 2026-07-22. Platform: R 4.4.2, Apple clang, x86_64 macOS. BLAS-related thread counts were fixed at one. Jones–Pewsey was excluded from implementation and benchmarking.

## Retention decision

The requested list contains 15 non-Jones–Pewsey families, not 14. Candidate kernels were tested for all 15. Only the following two implementations passed both the bitwise-equality gate and the end-to-end performance gate:

| Family | Confirmatory pairs | Median gain | Median speedup | Bootstrap 95% CI for speedup | Worst paired gain | Decision |
|---|---:|---:|---:|---:|---:|---|
| Normal | 20 | 74.56% | 3.9314 | [3.8770, 3.9587] | 73.63% | Retain |
| Weighted small-circle mixture | 20 | 27.24% | 1.3745 | [1.3486, 1.3960] | 16.38% | Retain |

Every confirmatory pair gave bitwise-identical observed profiles, statistics, bootstrap distributions, critical values, p-values, decisions, grids, and fitted parameters. Only timing diagnostics and the new backend metadata differed.

The weighted-mixture confirmation pools paired scenarios with sample sizes 50, 100, and 200 and bootstrap loads of 5 and 10. Its slowest improvement was positive; therefore no relevant scenario approached the allowed 5% regression limit. The normal confirmation used a composite null, KS and CvM, sample size 100, and 200 bootstrap replicates per run.

## Rejected candidates

The preliminary screen used paired alternating-order evaluations and checked `identical()` before considering timing. Approximate warm gains below reflect the screen's resolution; all corresponding production C++ routes were removed or made unavailable.

| Family | Main observed bottleneck | Warm gain in profile screen | Decision |
|---|---|---:|---|
| Multivariate normal | numerical quadratic-form CDF | -1.9% | Reject |
| Logistic Gaussian | numerical quadratic-form CDF | approximately 0% | Reject |
| vMF | Bessel evaluations and CDF tables | 0–12% | Reject |
| HvMF | Bessel evaluations and spline construction | -0.7% | Reject |
| Spherical Cauchy | existing vectorized/compiled work | approximately 0% | Reject |
| Small circle | Legendre recurrence plus matrix work | slower by about 25% | Reject |
| Watson | Legendre recurrence plus matrix work | slower by about 31% | Reject |
| Cardioid | closed vectorized R expressions | approximately 0% | Reject |
| Beta mixture | dense eigendecomposition already compiled | +0.5% | Reject |
| Uniform–beta mixture | Legendre/integration work | +0.9% | Reject |
| Logitnormal mixture | numerical integration | +0.4% | Reject |
| Symmetric small-circle mixture | Legendre recurrence and MLE | slower by about 30% | Reject |
| Axial truncated-normal mixture | stable Rmath CDF calculations | approximately +1% | Reject |

For beta and the rotational mixtures, the weighted MLE was as costly as or costlier than profile evaluation. It was deliberately left unchanged in this iteration.

## Compilation and amortization

A forced clean compilation of the final retained translation unit took 3.038 seconds; loading an already cached build took 0.019 seconds. At the confirmatory workloads, the clean compilation cost is recovered after approximately 4 normal tests or 5 weighted-mixture tests. Normal use remains on the R backend and pays no compilation or loading cost.

## Exactness coverage

The CI test contains 1,000 seeded differential scalar cases for the normal family, 50 random matrix/CvM cases, 40 weighted-small-circle scalar/grid cases, and 25 weighted empirical-profile cases with ties. It includes zero and extreme scales/concentrations, aligned and antipodal directions, boundary radii, mixed vector lengths, empty and invalid inputs, `NA` values, value/dimension/attribute checks, and error/warning comparisons.

Complete simple and composite tests use common seeds for KS and CvM. Complete weighted-mixture composite tests and normal sequential/parallel tests also require bitwise identity. Thirteen rejected families must fail explicitly when C++ is requested. Jones–Pewsey functions expose no backend selector and their existing regression tests remain authoritative.

The reproducible confirmatory runner is `scripts/benchmark_cpp_distance_profiles.R`; it writes pair-level timings, summaries, session information, cold compilation time, cached-load time, and estimated amortization counts under `output/cpp_distance_profile_benchmark/` by default.

All 44 existing `testthat` files were exercised. The full run exposed one legacy PSOCK worker that exported the normal profile without the backend helpers; the wrapper now detects that environment and preserves the R path. That failing file and the complete new differential suite passed after the fix; all other files had passed in the full run. The only remaining messages were five pre-existing package/data warnings.
