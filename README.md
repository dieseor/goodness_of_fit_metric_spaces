
# gofmetric

[![R CI](https://github.com/dieseor/gofmetric/actions/workflows/ci.yml/badge.svg)](https://github.com/dieseor/gofmetric/actions/workflows/ci.yml) [![Coverage](https://codecov.io/gh/dieseor/gofmetric/branch/main/graph/badge.svg)](https://app.codecov.io/gh/dieseor/gofmetric)

`gofmetric` implements distance-profile Kolmogorov–Smirnov and Cramér–von Mises goodness-of-fit tests for distributions on metric spaces. It accompanies [*Goodness-of-fit for distributions on metric spaces*](https://arxiv.org/abs/2609.38158). The badges become informative when this branch and its CI workflow are published under the new repository name.

## Install

The declared minimum is R 4.4.0 with compatible older dependencies; CI checks R 4.4.0 and 4.4.2 using [`renv.lock`](renv.lock). For a fresh install with current CRAN dependencies, use R 4.5 or newer: the current `gsl` version, required indirectly through `sphunif`, needs R 4.5. An ordinary install does not restore the historical dependencies. Install the development version from GitHub:

``` r
install.packages("remotes")
remotes::install_github("dieseor/gofmetric")
```

## Example

``` r
library(gofmetric)
set.seed(1)
x <- rnorm(100)
fit <- gof_test(x, h0 = "normal", B = 999)
fit$inference$ks$p_value
fit$inference$cvm$p_value
```

The `h0` argument chooses among twelve distribution families and the joint sunspots time-and-location model. See `?gof_test` for their names and each wrapper’s help page for its data format and parameters. For `mvnormal`, pass a matrix even for one-dimensional data: `matrix(x, ncol = 1)`. The joint sunspots model accepts a list with `x` (unit-sphere positions) and `s` (times in `(0, 1)`), or a four-column matrix. It currently supports a fitted composite null with the fast bootstrap.

## Reproduce the paper

The historical project environment is recorded in [`renv.lock`](renv.lock) and is kept separate from the current-dependency installation checks. It uses R 4.4.2. Individual saved runs can have different package versions: for example, the [sunspots run’s `sessionInfo.txt`](real_data/reruns/paper_main_realdata_B1000_3cores_20260831_113532/sunspots/sessionInfo.txt) records `sphunif` 1.4.3, while the lockfile records 1.4.2. Consult each run’s session information when checking its exact provenance. After restoring the lockfile environment, run the scripts from the repository root:

``` sh
Rscript -e 'renv::restore()'
R CMD INSTALL .
bash paper/run_simulations.sh   # long cluster campaign
bash paper/run_real_data.sh
sh paper/rebuild_figures.sh paper/rendered
```

The repository also contains data processing and figure generation code. The simulation campaign is computationally expensive.

## Citation

Run `citation("gofmetric")` for the versioned software citation and the three-author article citation. The package is distributed under the MIT license.
