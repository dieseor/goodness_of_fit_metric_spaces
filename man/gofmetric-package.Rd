\name{gofmetric-package}
\alias{gofmetric-package}
\alias{gofmetric}
\title{Goodness-of-fit tests on metric spaces}
\description{Kolmogorov--Smirnov and Cramer--von Mises tests based on distance profiles and calibrated by multiplier bootstrap.}
\details{
The public families are \code{normal}, \code{mvnormal},
\code{restricted_spiked_normal}, \code{normal_sigma_Id}, \code{vmf},
\code{vmf_fixed_kappa}, \code{hvmf}, \code{logistic_gaussian},
\code{logistic_gaussian_ar1}, \code{cardioid}, \code{small_circle},
\code{uniform_beta_mixture}, and the joint \code{sunspots_joint_time_space}
model with shared hemisphere regression. Use \code{\link{gof_test}} for a common
interface or a \code{multiplier_bootstrap_*} function for model controls.
The joint sunspots model supports only the fast multiplier bootstrap.
}
\seealso{\code{\link{gof_test}}}
