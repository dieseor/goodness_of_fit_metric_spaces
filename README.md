# dpgof

This repository contains the code needed to reproduce the results of the paper *Goodness-of-fit for distributions on metric spaces*.

The tests are built from distance profiles and calibrated with a multiplier bootstrap. There are models for Euclidean, directional, hyperbolic and compositional data.

## Install

You need R 4.4 or newer.

```r
install.packages("remotes")
remotes::install_github("dieseor/dpgof")
```

The package versions we used are in `renv.lock`.

## Example

Test whether a sample is normal, with the mean and variance estimated from the data:

```r
library(dpgof)
set.seed(1)
x <- rnorm(100)

# Fast bootstrap
fast <- multiplier_bootstrap_normal(x, null = list(type = "composite"),
                                    unknown_param = "both",
                                    bootstrap_method = "fast_multiplier",
                                    keep = list(observed_process = FALSE))
fast$inference$ks$p_value
fast$inference$cvm$p_value

# Bootstrap that re-estimates the parameters in each replicate
slow <- multiplier_bootstrap_normal(x, null = list(type = "composite"),
                                    unknown_param = "both",
                                    bootstrap_method = "reestimated")
slow$inference$ks$p_value
slow$inference$cvm$p_value
```

Both give similar p-values, but the fast bootstrap is much quicker: on this example it takes well under a second, against about half a minute for the re-estimated one. Other models work the same way through the `multiplier_bootstrap_*` functions (`multiplier_bootstrap_vmf()`, `multiplier_bootstrap_logistic_gaussian()`, ...).

## Reproducing the paper

[paper/README.md](paper/README.md) lists the tables and figures of the paper together with the scripts and saved inputs that produce them.
