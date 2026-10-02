\name{multiplier_bootstrap_normal}
\alias{multiplier_bootstrap_normal}
\title{Multiplier bootstrap: normal}
\description{Run distance-profile goodness-of-fit tests for the normal model.}
\usage{
multiplier_bootstrap_normal(data, null, statistics = c("ks", "cvm"), ks_grid = NULL,
  B = 5000, alpha = 0.05, multipliers = NULL, n_cores = 1, seed = NULL,
  bootstrap_method = c("reestimated", "fast_multiplier"),
  keep = list(observed_process = TRUE, bootstrap_statistics = TRUE,
  bootstrap_thetas = FALSE), control = list(), unknown_param = NULL,
  distance_profile_backend = c("r", "cpp"))
}
\arguments{
  \item{data}{
Observed sample: a numeric vector for the univariate normal wrapper, or a numeric matrix with observations in rows for multivariate models.
}
  \item{null}{
Use \code{list(type = "simple", theta = ...)} for fixed model parameters or \code{list(type = "composite")} to estimate parameters from the sample. The required fields of \code{theta} depend on the model.
}
  \item{statistics}{
One or both of \code{"ks"} and \code{"cvm"} (Cramer--von Mises).
}
  \item{ks_grid}{
Optional KS evaluation grid. \code{NULL} selects the model default.
}
  \item{B}{
Positive number of multiplier bootstrap replications. Monte Carlo p-values have resolution \code{1 / (B + 1)}.
}
  \item{alpha}{
Nominal level used for the rejection decision.
}
  \item{multipliers}{
Optional multiplier weights; \code{NULL} generates them internally.
}
  \item{n_cores}{
Number of parallel workers. \code{1} runs serially.
}
  \item{seed}{
Optional seed for reproducible bootstrap weights and model simulation.
}
  \item{bootstrap_method}{
\code{"reestimated"} refits the model in each bootstrap replicate; \code{"fast_multiplier"} uses the linear correction when supported.
}
  \item{keep}{
A list with logical fields \code{observed_process}, \code{bootstrap_statistics}, and \code{bootstrap_thetas}.
}
  \item{control}{
List of numerical controls passed to the model specification, including derivative and integration settings where applicable.
}
  \item{unknown_param}{
Parameters estimated under the null.
}
  \item{distance_profile_backend}{
\code{"r"} or \code{"cpp"} for distance-profile evaluation where supported. This is separate from the fast-multiplier kernel backend.
}
}
\value{A list with \code{observed}, \code{bootstrap}, \code{inference}, \code{grid}, and \code{diagnostics}. Each requested statistic has a p-value, critical value, and rejection flag under \code{inference}.}
\details{Univariate normal data. For a simple null, use \code{theta = list(mu = ..., sigma = ...)}. For the shared bootstrap arguments and result structure, see \code{\link{multiplier_bootstrap_gof}}.}
\examples{
set.seed(1)
x <- rnorm(12)
fit <- multiplier_bootstrap_normal(x, null = list(type = "simple",
  theta = list(mu = 0, sigma = 1)), statistics = "ks", B = 19)
fit$inference$ks$p_value
}
