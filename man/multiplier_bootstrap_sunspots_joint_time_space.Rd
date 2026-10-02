\name{multiplier_bootstrap_sunspots_joint_time_space}
\alias{multiplier_bootstrap_sunspots_joint_time_space}
\title{Joint time and location goodness-of-fit test for sunspots}
\description{Test the joint model with a two-component beta mixture for time and a time-varying symmetric small-circle mixture for location.}
\usage{
multiplier_bootstrap_sunspots_joint_time_space(data, null,
  statistics = c("ks", "cvm"), ks_grid = NULL, B = 5000, alpha = 0.05,
  multipliers = NULL, n_cores = 1, seed = NULL,
  bootstrap_method = "fast_multiplier",
  keep = list(observed_process = FALSE, bootstrap_statistics = TRUE,
              bootstrap_thetas = FALSE), control = list(),
  distance_profile_backend = "cpp")
}
\arguments{
  \item{data}{A list with \code{x}, an \eqn{n \times 3} matrix of unit sphere positions, and \code{s}, \eqn{n} times in \eqn{(0,1)}; or a four-column matrix with positions followed by times.}
  \item{null}{\code{list(type = "composite")} to fit both components. The current fast implementation does not support a simple null.}
  \item{statistics}{\code{"ks"}, \code{"cvm"}, or both.}
  \item{ks_grid}{Optional evaluation grid; \code{NULL} uses unique sample distances.}
  \item{B}{Positive number of bootstrap replicates.}
  \item{alpha}{Nominal test level.}
  \item{multipliers}{Optional supplied multiplier weights.}
  \item{n_cores}{Number of workers.}
  \item{seed}{Optional random seed.}
  \item{bootstrap_method}{Only \code{"fast_multiplier"} is implemented. Re-estimation raises an error.}
  \item{keep}{Logical retention flags for the observed process, bootstrap statistics and bootstrap parameter estimates.}
  \item{control}{Numerical controls. \code{hemisphere_regression} must be \code{"shared"}.}
  \item{distance_profile_backend}{\code{"cpp"} (default) or \code{"r"}.}
}
\details{The temporal and conditional spatial components are fitted separately, but the test concerns their joint distribution. The spatial model uses a common linear time regression for the two hemispheres.}
\value{A goodness-of-fit result with observed statistics, bootstrap statistics, p-values under \code{inference}, and diagnostics.}
\seealso{\code{\link{gof_test}}}
