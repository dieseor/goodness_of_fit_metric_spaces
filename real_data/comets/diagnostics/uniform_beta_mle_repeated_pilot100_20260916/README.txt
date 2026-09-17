Repeated Uniform-Beta MLE diagnostic (separate from production outputs).
n_fits per dataset: 100; cores: 8; seed: 20260916
Every fit calls uniform_beta_mixture_mle_s2_weighted() directly.
The only added controls provide one external warm start and request warm_start_only=TRUE.
If that direct production routine returns non-convergence, its own documented fallback to its ordinary candidate starts remains active.
Density components use a 512-point sphere grid and Hellinger distance threshold 0.025; parameter distances are not used for grouping.
