# Reproducing the paper results

Run these commands from a clone of this repository. The manuscript lives in a separate repository; set `PAPER_TEX` to its main `.tex` file when comparing numbers. `paper/manuscript_inventory.tsv` indexes the 7 active tables and 28 active image references in the manuscript snapshot with SHA-256 `627029f39220747dbe53a2713446cfadd06ad8679b486666550cf22640a7a67a`.

```sh
Rscript -e 'renv::restore()'
R CMD INSTALL .
Rscript --vanilla paper/check_package.R
sh paper/rebuild_figures.sh /tmp/dpgof-figures
python3 paper/check_tables.py "$PAPER_TEX"
```

The figure command reads only the small saved inputs in `paper/artifacts/` and writes to the output directory supplied. It does not rerun Monte Carlo experiments or write to the manuscript. The artifact hashes and origins are in `paper/artifacts/manifest.tsv`.

## Tables

| Manuscript label | Saved result and effective source |
| --- | --- |
| `tab:simulation_scenarios` | Scenario definitions in the manuscript; effective per-campaign settings, result paths and hashes are in `artifacts/tables/simulation_audit.json`. The simulation runners vary by model; `scripts/run_section6_new_scenarios.R` is the main catalogue. |
| `tab:empirical_null_calibration` | 468 cells in `simulation_audit.json` plus 12 scenario 1, `n=800` cells in `scenario1_n800.csv`. The audit records campaign manifests, raw-result paths, bootstrap methods and 43 file hashes. |
| `tab:bootstrap_procedure_comparison` | `table3_fast_rates.csv` and `table3_slow_summary.csv`. The fast summary is rebuilt by `paper/summarize_table3_fast.py` from `simulation_results/paper_fast_matched_c3_20260929_024311/fast_batches/`; the slow summary uses `simulation_results/paper_fast_vs_reestimated_n100_beta0/slow_batches/`. Each of the 16 groups has 1,000 distinct, seed-matched replications. All 80 displayed values, including both timing rows and speed-ups, match the manuscript. The raw batches remain local research outputs; the processed summaries are included here. |
| `tab:sunspots_cycle23_temporal_models` | `sunspots_gof.csv`, step 1 of `paper/run_real_data.sh`. |
| `tab:comets-gof` | `comets_c2_sc.csv` and `comets_ub.csv`, steps 2-4, 7 and 8 of `paper/run_real_data.sh`. |
| `tab:risoe-wind-gof` | `wind_gof.csv`, step 5 of `paper/run_real_data.sh`. |
| `tab:logistic-gaussian-real-data` | `compositions_ks_cvm.csv`, step 6 of `paper/run_real_data.sh`. The BHEP column comes from `scripts/run_logistic_gaussian_hz_pvalues.R` and is not recomputed by `check_tables.py`. |

`python3 paper/check_tables.py "$PAPER_TEX"` compares 480 simulation cells, all 80 values in the bootstrap comparison table, and the stored GOF values in the four real-data tables with the current manuscript. It does not rerun the statistical experiments.

## Figures

| Manuscript figures | Saved input | Rendering code | Verification against the manuscript snapshot |
| --- | --- | --- | --- |
| `fig:s1_limit_process_visualization` (2 PNGs) | `artifacts/s1/figure1_plot_input.rds`, reduced from the saved simulation result, and copies of both current manuscript PNGs | `scripts/gaussian_process_s1_visualization.R` | The saved input is the simulated process shown in the figure. |
| `fig:sunspots_cycle23_spatial_windows` (5 PDFs) | `artifacts/sunspots/` | `scripts/plot_spherical_parametric_hdr_bands.R` | All 5 rendered PDFs matched pixel for pixel. |
| `fig:risoe-wind-densities` (4 PDFs) | `artifacts/wind/`, cleaned 77 m subsamples from DTU `Risoe_m_all.nc` | `real_data/wind/plot_risoe_extended_r_density_panels.R` | All 4 rendered PDFs matched pixel for pixel. The large original NetCDF is not required to redraw them. |
| `fig:logistic-gaussian-simplex-d3` (5 PDFs) | `artifacts/simplex/`, closed subsets of datasets from the R package `compositions` | `scripts/run_sediments_simplex_contours.R` | All 5 rendered PDFs matched pixel for pixel. |
| `fig:convergence_process_vmf_simple_mu_100` and `fig:convergence_process_vmf_comp_mu_100` (12 PNGs) | `artifacts/vmf_convergence/`, six saved fixed-grid Monte Carlo results | `scripts/plot_paper_vmf_convergence_from_rds.R` | All 12 PNGs matched byte for byte using the bandwidth adjustments in `rebuild_figures.sh`. |

The comparison was against the manuscript image files on 2026-09-30. PDF matching means matching rasterized pages at 100 dpi; PDF metadata may differ.

## Source data and scope

The sunspot retained observations and fitted parameters are saved in `artifacts/sunspots/`; the original catalogue and processing code are under `real_data/sunspots/`. Wind samples are derived from the open DTU Risø archive cited in the article. The five simplex inputs are copied from the `compositions` datasets after the same closure used for the article. The vMF files are saved simulation results, not fresh simulations. The complete simulations, cluster logs and large raw data are intentionally separate from the small publication artifacts.
