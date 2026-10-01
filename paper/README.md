# Reproducing the paper

Run everything from the root of the repository. First install the package and the package versions we used:

```sh
Rscript -e 'renv::restore()'
R CMD INSTALL .
```

The small inputs needed for the figures and tables are in `paper/artifacts/` (hashes and origin of each file in `paper/artifacts/manifest.tsv`). To check the saved results against the manuscript tables:

```sh
python3 paper/check_tables.py path/to/manuscript.tex
```

To redraw all the figures from the saved inputs (a few minutes):

```sh
sh paper/rebuild_figures.sh paper/rendered
```

The sections below say where each result comes from. The Monte Carlo experiments were run on a cluster; the figures and tables can be rebuilt without rerunning them.

## Figure 1: limit process on the circle

`scripts/gaussian_process_s1_visualization.R` simulates and draws the process. The simulated process used in the paper is saved in `artifacts/s1/figure1_plot_input.rds` and can be redrawn with `replot_limit_gaussian_s1_vmf_from_result()`.

## Simulation study: rejection rates

Each scenario has its own runner. All campaigns used `M=1000` and `B=1000`; the sample sizes, values of β, seeds and derivative method (`score_mc` for scenarios 1–4, `quadrature` for 5–8) of each campaign are in `artifacts/tables/simulation_audit.json`.

| Scenario | Runner |
| --- | --- |
| 1 | `scripts/run_restricted_spiked_normal_covariance_alternatives.R --mean_config=diagonal_100` |
| 2 | `scripts/run_normal_sigma_Id_t_pilot.R` |
| 3 | `scripts/run_logistic_gaussian_dirichlet15_mu_only.R` |
| 4 | `scripts/run_logistic_gaussian_sigma_shape_scenarios.R --scenario=t4` |
| 5 | `scripts/run_vmf_mu_only_fixed_kappa_pilot.R --kappa_values=2` |
| 6 | `scripts/run_vmf_antipodal_fixed_kappa_pilot.R` with `--scenario_type=projected_normal_mean_d` (β=0) and `--scenario_type=projected_normal_2sqrt_d_kappa_2d_beta_half` (β>0) |
| 7, 8 | `scripts/run_section6_new_scenarios.R --family=hvmf` |

## Fast versus re-estimated bootstrap

The replications were run with `scripts/run_paper_fast_matched_batch.R` and `scripts/run_paper_slow_reestimated_batch.R` (arguments `SCENARIO D REP_START REP_END`; `scripts/run_paper_fast_matched_batch.sbatch` submits the fast batches with Slurm). The raw batches are in `artifacts/table3/`, and the two summaries in the table are rebuilt with

```sh
python3 paper/summarize_table3.py fast paper/artifacts/table3/fast_batches fast.csv
python3 paper/summarize_table3.py slow paper/artifacts/table3/slow_batches slow.csv
```

## Real data: sunspots, comets, wind and compositions

```sh
bash paper/run_real_data.sh
```

runs the goodness-of-fit tests of the four real-data tables with the seeds of the paper. The BHEP column of the compositional table comes from `scripts/run_logistic_gaussian_hz_pvalues.R`.

The figures of the real-data section are drawn by `rebuild_figures.sh` from the saved samples: sunspots with `scripts/plot_spherical_parametric_hdr_bands.R`, wind with `real_data/wind/plot_risoe_extended_r_density_panels.R` and the simplex contours with `scripts/run_sediments_simplex_contours.R`. The wind samples were selected from the open DTU Risø archive cited in the paper.

## Supplementary material: convergence of the vMF process

`scripts/regenerate_paper_vmf_convergence_grid.R` runs the simulations, and `rebuild_figures.sh` draws the figures from the saved results in `artifacts/vmf_convergence/` with `scripts/plot_paper_vmf_convergence_from_rds.R`.
