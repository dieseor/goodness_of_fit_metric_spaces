# Empirical audit of the comet Uniform--Beta fit

This directory is diagnostic only. It does not replace any official result and no production or manuscript file was changed.

## Configuration and reproduction

The data are the `sphunif::comets` orbital normals after the repository's current complete-case, period, class, and fragment filters: short-period n = 784 and long-period n = 610.

The production fit maximizes the mean projected log-density with BFGS, using the first 12 of 16 deterministic candidate starts, `maxit = 350`, and `reltol = 1e-9`. This differs from the full spherical log-likelihood only by multiplication by n and addition of the constant `-n log(4 pi)`. The winning initial candidate was 2 for short-period and 6 for long-period. There were no warnings.

| Dataset | Saved full log-likelihood | Reproduced full log-likelihood | AIC | BIC | convergence | maximum parameter discrepancy |
|---|---:|---:|---:|---:|---:|---:|
| Short | -168.854071945 | -168.854071944 | 347.708144 | 371.030189 | 0 | 1.23e-9 |
| Long | -1520.620732542 | -1520.620700817 | 3051.241402 | 3073.308696 | 0 | 1.21e-7 |

The implemented spherical density is exactly `1/(4 pi)` times the stated uniform--Beta mixture for interior projected coordinates. For numerical evaluation the code clips y to `[1e-12, 1 - 1e-12]`, clips w to `[0.01, 0.99]`, and clips each shape to `[0.05, 1000]`. Consequently, the mathematical likelihood is unbounded when a shape is below one, but the implemented objective has a finite endpoint plateau. Neither saved fit is itself on that plateau.

## Local analysis

The local chart uses two exponential-map coordinates for the direction and `(logit(w), log(alpha), log(beta))` for the remaining parameters. Gradients are for the full spherical log-likelihood. Central finite-difference steps were kept at least 50 times smaller than the closest endpoint distance. Independent analytic score calculations agreed with the numerical gradients to relative error below 7e-5.

| Dataset | gradient norm | Hessian eigenvalues (descending) | endpoint distance used to set step |
|---|---:|---|---:|
| Short | 5577.139 | 7.018e7, -64.09, -145.93, -1223.33, -7.017e7 | 8.124e-5 |
| Long | 16399.65 | 4.290e8, -1.50, -140.96, -464.95, -4.288e8 | 3.832e-5 |

The positive Hessian eigenvalue persisted when the finite-difference step varied from endpoint-distance/20 to endpoint-distance/500. Thus the reported points are not stationary local maxima, despite `optim` returning convergence code 0. Of 50 small perturbations, 16 improved the short fit and 45 improved the long fit. All 50 optimizations returned code 0.

## Multistart audit

There were 240 starts per dataset, including production starts, random directions, data and antipodal directions, weights from 0.03 to 0.8, and shapes on both sides of one. Basin counts below use a numerical equivalence rule and identify the exact symmetry `(mu, alpha, beta) = (-mu, beta, alpha)`; they are descriptive clusters, not a claim about the topology of the likelihood.

| Dataset | code 0 | code 1 | solutions better than reported | endpoint-plateau solutions | numerical basins | frequency of reported structural basin |
|---|---:|---:|---:|---:|---:|---:|
| Short | 224 | 16 | 132 | 16 | 46 | 166/240 |
| Long | 212 | 28 | 15 | 27 | 136 | 3/240 |

The best short run had log-likelihood -167.0284, AIC 344.0568, BIC 367.3788, and nearest-point distance 1.998e-6 radians. It belongs to the same broad north-pole/concentrated-Beta structural basin as the reported fit, but lies on the numerical endpoint plateau. Hence AIC 347.7 is reproducibly associated with a broad structural basin, yet the precise reported point is not its local maximizer.

The best long run had log-likelihood -1518.621, AIC 3047.242, and BIC 3069.309. The multistart landscape was much more fragmented and contained many observation-specific endpoint solutions.

## Observation contributions

For short-period, the maximum individual log-density is 4.660, versus a median of 0.457 and a 99th percentile of 1.959. Relative to the uniform model, the total log-likelihood gain is 1815.47; the largest observation contributes 0.40% of that gain, the largest five 1.46%, and the largest ten 2.70%. The exceptionally small AIC is therefore a distributed consequence of strong concentration, not a one-point spike.

For long-period, the maximum individual log-density is 3.307, versus a median of -2.599 and a 99th percentile of -1.487. The total gain over the uniform model is only 23.304; the largest observation contributes 25.1% of that gain, the largest five 53.6%, and the largest ten 74.3%. There is no astronomically large single contribution, but the AIC advantage over uniformity has strong endpoint leverage.

## Stable profiles toward singularities

Profiles use analytic log-coordinates for the target observation, avoiding floating-point rounding of `1 - y` or `y` to zero. Conditional profiles retain the production bounds on w, alpha, and beta.

- Short, beta endpoint: the fixed-parameter path first exceeds the reported fit at 5.84e-5 radians; gains of 1, 2, and 5 log-likelihood units require approximately 8.03e-6, 7.93e-7, and 1.48e-9 radians. Conditional reoptimization improves already at the reported distance because the reported point is not stationary.
- Long, alpha endpoint: the fixed path exceeds the report at 2.78e-5 radians; gains of 1, 2, and 5 require 7.74e-6, 1.56e-6, and 1.29e-8 radians.
- Long, beta endpoint: from the reported direction the target is 0.04165 radians away. The fixed path does not exceed the report until 6.01e-6 radians and needs 1.12e-6 radians for a one-unit gain. Conditional reoptimization improves immediately for the separate reason that the reported shape/weight coordinates are not stationary.

Thus the singular regions producing material fixed-parameter gains are narrow, but the reported fits are already close enough to one relevant endpoint that a direct ascent exists.

## Re-estimated multiplier bootstrap

The detailed refit audit uses the exact production reestimation controls: normalized Exp(1) weights, observed-fit warm start only, one initial start, BFGS `maxit = 80`, and `reltol = 1e-6`. It uses the paper's KS seeds and B = 200. The full GOF audit separately runs KS and CvM with the paper seeds, sample-points/unique-distances KS grid, and the production `reestimated` engine.

| Dataset | detailed refits | code 0 | code 1 | unchanged at start | endpoint refits |
|---|---:|---:|---:|---:|---:|
| Short | 200 | 200 | 0 | 16 | 0 |
| Long | 200 | 197 | 3 | 51 | 16 |

The full B = 200 results are:

| Dataset | statistic | p-value | 95% critical value | endpoint refits |
|---|---|---:|---:|---:|
| Short | KS | 0.5423 | 2.2576 | 0 |
| Short | CvM | 0.5423 | 0.3639 | 0 |
| Long | KS | 0.0896 | 1.8277 | 16 |
| Long | CvM | 0.1592 | 0.3216 | 15 |

The long-period endpoint runs did not create the upper tail: their KS statistics ranged from 1.027 to 1.718, and their CvM statistics from 0.052 to 0.238. No singular KS run exceeded the observed statistic; two singular CvM runs did. Removing endpoint runs diagnostically changes the long p-values only from 0.0896 to 0.0973 (KS) and from 0.1592 to 0.1613 (CvM). This deletion is a sensitivity calculation, not a proposed procedure.

The official B = 1000 results are 0.5445/0.4925 for short KS/CvM and 0.0929/0.1868 for long KS/CvM. The B = 200 diagnostics are compatible with them at ordinary Monte Carlo resolution.

## Answers A--G

**A. Global likelihood.** The mathematical UB likelihood is unbounded for these fitted shape regimes. The production evaluator implements the same interior density but makes its numerical optimization problem bounded through endpoint and parameter clipping.

**B. Local maxima.** Neither reported fit is a genuine stationary local maximum of the implemented or mathematical interior likelihood. Both have a large nonzero tangent score, an indefinite Hessian, and improving arbitrarily small perturbations. This is distinct from merely failing to be a global maximum.

**C. Short-period AIC 347.7.** It represents a broad, repeatedly reached structural basin, not a fit whose reported likelihood is dominated by one observation. However, the exact reported point is an optimizer stopping point rather than a local maximum, and nearby endpoint-directed fits improve it modestly before reaching the numerical plateau.

**D. Long-period local regularity.** The visual fit can be reasonable, but the fitted point is not numerically regular as a local maximum. It is closer to an alpha-endpoint singularity, has a strongly nonzero score, and multistart behavior is fragmented.

**E. Bootstrap stability.** Short-period reestimation remained away from endpoints in 200/200 runs. Long-period did not: 16/200 KS-seed refits and 15/200 CvM-seed refits reached the numerical endpoint region. Some long fits are also sensitive at machine-level branching near the plateau.

**F. Reliability of reported UB p-values.** There is no empirical evidence here that the reported short-period p-values are numerically unreliable. Long-period reestimation is demonstrably unstable in parameter space, but in B = 200 the endpoint refits were not responsible for the relevant upper tails, trimmed p-values barely changed, and the pilot p-values agree with B = 1000. Thus there is no observed material p-value distortion in these runs, while long-period deserves a numerical-fragility qualification rather than a clean stability claim.

**G. AIC/BIC influence.** Short-period AIC/BIC are not driven by a few observations. Long-period has substantial endpoint leverage: a few points account for most of its modest gain over the uniform model, although no individual log-density contribution is itself enormous.

## Reproduction commands

Run the likelihood/local/multistart/refit audit with:

`Rscript --vanilla scripts/diagnose_comets_uniform_beta_mle.R --output_root=real_data/comets/diagnostics/uniform_beta_mle_audit_final_exact_bootstrap_20260915 --n_multistart=240 --n_local_perturb=50 --n_bootstrap=200 --n_cores=3`

Run the full bootstrap GOF audit with:

`Rscript --vanilla scripts/diagnose_comets_uniform_beta_bootstrap_gof.R --B=200 --n_cores=3 --output_root=real_data/comets/diagnostics/uniform_beta_bootstrap_gof_B200_20260915`

