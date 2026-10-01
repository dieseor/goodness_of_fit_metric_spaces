\name{gof_test}
\alias{gof_test}
\title{Goodness-of-fit test}
\description{Test whether a sample comes from a given family of distributions, using distance profiles and a multiplier bootstrap.}
\usage{
gof_test(data, h0, theta = NULL, bootstrap = c("fast", "reestimated"), B = 1000,
         statistics = c("ks", "cvm"), seed = NULL, n_cores = 1, ...)
}
\arguments{
  \item{data}{The sample: a numeric vector or a matrix with one observation per row.}
  \item{h0}{The distribution under the null hypothesis, for example \code{"normal"}, \code{"mvnormal"}, \code{"vmf"}, \code{"hvmf"}, \code{"cardioid"} or \code{"logistic_gaussian"}. Each value corresponds to a \code{multiplier_bootstrap_*} function.}
  \item{theta}{\code{NULL} (default) estimates the parameters from the sample. Otherwise, a list with the parameters of a fully specified null.}
  \item{bootstrap}{\code{"fast"} (default) or \code{"reestimated"}, which refits the model in each bootstrap replicate and is much slower.}
  \item{B}{Number of bootstrap replicates.}
  \item{statistics}{\code{"ks"}, \code{"cvm"} or both.}
  \item{seed}{Optional seed.}
  \item{n_cores}{Number of cores.}
  \item{\dots}{Further arguments for the model, passed to the corresponding \code{multiplier_bootstrap_*} function (for example \code{k} for the cardioid or \code{distance_type}).}
}
\value{The result of the corresponding \code{multiplier_bootstrap_*} function. The p-values are in \code{$inference$ks$p_value} and \code{$inference$cvm$p_value}. The element \code{$bootstrap} gives the bootstrap requested and used; if the fast bootstrap is not valid for the fitted model (for example at a boundary of the parameter space), the re-estimated one is used, \code{$bootstrap$fallback} is \code{TRUE} and a warning is given.}
\examples{
set.seed(1)
x <- rnorm(50)
fit <- gof_test(x, h0 = "normal", B = 99)
fit$inference$ks$p_value
}
