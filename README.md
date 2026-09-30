# dpgof

`dpgof` computes goodness-of-fit tests from distance profiles using a multiplier bootstrap. It includes Euclidean, directional, hyperbolic, and compositional models.

## Install

Install R 4.4 or newer, then run:

```r
install.packages("remotes")
remotes::install_github("dieseor/dpgof")
```

For the package versions used in the research project, see `renv.lock`.

## A first test

```r
library(dpgof)
set.seed(1)
x <- rnorm(30)
fit <- multiplier_bootstrap_normal(
  x,
  null = list(type = "simple", theta = list(mu = 0, sigma = 1)),
  statistics = c("ks", "cvm"),
  B = 199
)
fit$inference$ks$p_value
fit$inference$cvm$p_value
```

`B = 199` keeps this example short; use more replications for an analysis. The help pages for `multiplier_bootstrap_*` show the available model wrappers and their arguments. Models outside the associated article may have less empirical validation.

## Article reproduction

[Paper files and their provenance](paper/README.md) lists the active tables and figures, their generating scripts, and the saved inputs. Long simulations are separate from figure and table rendering.

The package code is in `R/` and `src/`. The old `utils.R` and bootstrap engine load the same R modules for existing research scripts.
