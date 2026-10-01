# dpgof

This repository contains the code needed to reproduce the results of the paper [*Goodness-of-fit for distributions on metric spaces*](https://arxiv.org/abs/2609.38158).

The tests are built from distance profiles and calibrated with a multiplier bootstrap. There are models for Euclidean, directional, hyperbolic and compositional data.

## Install

You need R 4.1 or newer (required by ggplot2). We have run the code with R 4.2.1 and 4.4.2.

```r
install.packages("remotes")
remotes::install_github("dieseor/goodness_of_fit_metric_spaces")
```

The package versions used to obtain the results of the paper are in `renv.lock`.

## Example

Test whether a sample is normal, with the mean and variance estimated from the data:

```r
library(dpgof)
set.seed(1)
x <- rnorm(100)

fast <- gof_test(x, h0 = "normal")                           # fast bootstrap
slow <- gof_test(x, h0 = "normal", bootstrap = "reestimated")  # re-estimates the parameters in each replicate

fast$inference$ks$p_value
fast$inference$cvm$p_value
```

In the simulations of the paper, the fast bootstrap was on average 438 times faster than the re-estimated one, with almost the same rejection rates. Other null distributions are chosen with `h0`, for example `"vmf"`, `"hvmf"` or `"logistic_gaussian"` (see `?gof_test`).

## Reproducing the paper

To reproduce the experiments from the paper, install the package versions we used and run the scripts in `paper/` from the root of the repository:

```sh
Rscript -e 'renv::restore()'
R CMD INSTALL .

bash paper/run_simulations.sh            # simulation study (long, we ran it on a cluster)
bash paper/run_real_data.sh              # real data
sh paper/rebuild_figures.sh paper/rendered   # figures
```

## Citation

If you use this code, please cite the paper: D. Serrano, E. García-Portugués and I. Van Keilegom (2026), *Goodness-of-fit for distributions on metric spaces*, [arXiv:2609.38158](https://arxiv.org/abs/2609.38158).
