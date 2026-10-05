# Goodness-of-fit for distributions on metric spaces

Reproducibility material for
[*Goodness-of-fit for distributions on metric spaces*](https://arxiv.org/abs/2609.38158).

The R package is maintained separately at [gofmetric](https://github.com/dieseor/gofmetric).

## Environment

The saved `renv.lock` records R 4.4.2. From this root, restore it with
`Rscript -e 'renv::restore()'`. Install gofmetric with
`remotes::install_github("dieseor/gofmetric")`.

## Saved artifacts and reruns

- Rebuild saved figures: `sh paper/rebuild_figures.sh paper/rendered`.
- Run simulations: `bash paper/run_simulations.sh` (long campaign).
- Run real-data analyses: `bash paper/run_real_data.sh`.
- Summarize timing batches with `paper/summarize_table3.py`.
- Compare saved tables with the manuscript using
  `python3 paper/check_tables.py /path/to/AoS/gof_metric_spaces_aos.tex`.
