\name{gof_test}
\alias{gof_test}
\title{Goodness-of-fit test}
\description{Test whether a sample comes from a given family of distributions, using distance profiles and a multiplier bootstrap.}
\usage{
gof_test(data, h0, theta = NULL, bootstrap = c("fast", "reestimated"), B = 1000,
         statistics = c("ks", "cvm"), seed = NULL, n_cores = 1, ...)
}
\arguments{
  \item{data}{The sample. Use a vector for \code{normal}, a matrix with observations in rows for \code{mvnormal} (including dimension one), and a list \code{list(x = positions, s = times)} or a four-column matrix for \code{sunspots_joint_time_space}. See each model's help page.}
  \item{h0}{One of \code{"normal"}, \code{"mvnormal"}, \code{"restricted_spiked_normal"}, \code{"normal_sigma_Id"}, \code{"vmf"}, \code{"vmf_fixed_kappa"}, \code{"hvmf"}, \code{"logistic_gaussian"}, \code{"logistic_gaussian_ar1"}, \code{"cardioid"}, \code{"small_circle"}, \code{"uniform_beta_mixture"}, or \code{"sunspots_joint_time_space"}.}
  \item{theta}{\code{NULL} (default) estimates the parameters from the sample. Otherwise, a list with the parameters of a fully specified null. The joint sunspots model currently requires \code{theta = NULL}.}
  \item{bootstrap}{\code{"fast"} (default) or \code{"reestimated"}, which refits the model in each bootstrap replicate and is much slower. Only \code{"fast"} is available for the joint sunspots model.}
  \item{B}{Number of bootstrap replicates.}
  \item{statistics}{\code{"ks"}, \code{"cvm"} or both.}
  \item{seed}{Optional seed.}
  \item{n_cores}{Number of cores.}
  \item{\dots}{Further arguments for the model, passed to the corresponding \code{multiplier_bootstrap_*} function (for example \code{k} for the cardioid or \code{distance_type}).}
}
\value{The result of the corresponding \code{multiplier_bootstrap_*} function. The p-values are in \code{$inference$ks$p_value} and \code{$inference$cvm$p_value}. The element \code{$bootstrap} gives the bootstrap requested and used. For supported models, a fast bootstrap that is invalid at a fitted boundary may fall back to re-estimation with a warning. The joint sunspots model raises an error if re-estimation is requested.}
\examples{
set.seed(1)
x <- rnorm(50)
fit <- gof_test(x, h0 = "normal", B = 99)
fit$inference$ks$p_value
}
