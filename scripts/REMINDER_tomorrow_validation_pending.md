# Pending For Tomorrow

These AoS experiments were intentionally left out of the rerun launcher for tonight:

1. KS convergence/validation figures for Normal and vMF.
The current convergence code in `convergence_empirical_process/` still works on fixed `omega × t` grids, whereas the paper text now states that the empirical KS approximation should use sample points and distinct pairwise distances. This needs a careful redesign rather than a blind rerun.

2. Bahadur validation for vMF.
The paper currently reports `M = 500`, and updating that to `M = 1000` may be expensive because the default run spans sample sizes up to `n = 100000`. Feasibility should be checked before launching it.

3. Giant H0 calibration table.
Explicitly excluded for this rerun by request.

4. Jones--Pewsey reruns.
Explicitly excluded for this rerun by request.
