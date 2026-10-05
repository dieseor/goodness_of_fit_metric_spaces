# Goodness-of-fit for distributions on metric spaces

Reproducibility material accompanying
[*Goodness-of-fit for distributions on metric spaces*](https://arxiv.org/abs/2609.38158).
This is the repository linked in the article:
[goodness_of_fit_metric_spaces](https://github.com/dieseor/goodness_of_fit_metric_spaces).

The R package is maintained separately at [gofmetric](https://github.com/dieseor/gofmetric).
No package source is duplicated here. The final dependency will be the immutable
`v0.1.0` tag and its exact commit SHA, created only after CRAN accepts 0.1.0.
The current lockfile intentionally has no provisional gofmetric pin.

## Environment

The saved `renv.lock` records R 4.4.2. From this root, restore it with
`Rscript -e 'renv::restore()'`. Until the final release is available, install the candidate gofmetric package
into this project library separately with
`remotes::install_github("dieseor/gofmetric")`. This temporary installation is
not an immutable release pin. After CRAN acceptance, record the
final tag and SHA in the lockfile and validate a clean restoration.

`utils.R` loads the installed package namespace for historical scripts, then
loads research-only helpers. It never sources or compiles package code locally.
NetCDF preprocessing requires NetCDF development libraries. Long campaigns
require the computational resources specified in their runners/Slurm scripts.

## Saved artifacts and reruns

- Rebuild saved figures: `sh paper/rebuild_figures.sh paper/rendered`.
- Run simulations: `bash paper/run_simulations.sh` (long campaign).
- Run real-data analyses: `bash paper/run_real_data.sh`.
- Summarize timing batches with `paper/summarize_table3.py`.
- Compare saved tables with the manuscript using
  `python3 paper/check_tables.py /path/to/AoS/gof_metric_spaces_aos.tex`.

The manuscript is maintained separately. `paper/manuscript_inventory.tsv` maps
its figures/tables to candidate producers; `paper/artifacts/manifest.tsv` records
saved inputs, provenance, and SHA-256 hashes. The figure rebuild script does not
yet cover the S1/Figure 1 panels; their saved inputs and reference images are kept.

Selected saved inputs and results are published under `paper/artifacts/`,
with provenance and hashes in `paper/artifacts/manifest.tsv`. Complete historical
outputs, pilots, duplicate campaigns, and raw wind NetCDF files remain preserved
in local backups and are not included in this repository. Rerunning the wind
preprocessing requires obtaining those raw data separately; acquisition and
redistribution instructions remain to be documented.

## Historical environment evidence

Preserve seeds, numerical controls, manifests, and each run's `sessionInfo.txt`.
The lockfile pins sphunif 1.4.2 whereas several saved runs report 1.4.3; this
historical discrepancy is not resolved by moving the files. Some optional research
workflows may require packages not captured by the existing lockfile. Saved-output
preservation does not establish successful end-to-end recomputation. The final
release pin and validation of a clean environment restoration are still pending.

Local diagnostics are excluded. MIT licensing notices are retained.
