#!/usr/bin/env Rscript

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("bootstrap/multiplier_bootstrap.R")

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("The information validation requires the installed `numDeriv` package.")
}

hvmf_validation_exact_sample <- function(n, mu, kappa) {
  q <- length(mu) - 1L
  m <- (q - 2) / 2
  proposal_gap <- if (m == 0) 0 else min(m / 2, kappa / 2)
  proposal_rate <- kappa - proposal_gap
  if (m == 0) {
    z <- stats::rexp(n, rate = kappa)
  } else {
    maximizer <- max(0, m / proposal_gap - 2)
    log_envelope <- m * log(maximizer + 2) -
      proposal_gap * maximizer
    z <- numeric(n)
    accepted <- 0L
    while (accepted < n) {
      remaining <- n - accepted
      proposal <- stats::rgamma(
        max(1024L, as.integer(ceiling(1.5 * remaining))),
        shape = m + 1,
        rate = proposal_rate
      )
      keep <- log(stats::runif(length(proposal))) <=
        pmin(
          m * log(proposal + 2) -
            proposal_gap * proposal - log_envelope,
          0
        )
      accepted_values <- proposal[keep]
      take <- min(length(accepted_values), remaining)
      if (take > 0L) {
        z[(accepted + 1L):(accepted + take)] <-
          accepted_values[seq_len(take)]
        accepted <- accepted + take
      }
    }
  }

  radial <- sqrt(z * (z + 2))
  directions <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
  directions <- directions / sqrt(rowSums(directions^2))
  base <- cbind(1 + z, radial * directions)
  spatial <- mu[-1L]
  spatial_norm <- sqrt(sum(spatial^2))
  if (spatial_norm > sqrt(.Machine$double.eps)) {
    direction <- spatial / spatial_norm
    boost <- matrix(0, q + 1L, q + 1L)
    boost[1L, 1L] <- mu[[1L]]
    boost[1L, -1L] <- spatial
    boost[-1L, 1L] <- spatial
    boost[-1L, -1L] <- diag(q) +
      (mu[[1L]] - 1) * tcrossprod(direction)
    base <- base %*% t(boost)
  }
  attr(base, "hvmf_sampler") <- "validation_exact_rejection"
  base
}

hvmf_validation_case <- function(label, q, kappa, rho, n_mc, seed) {
  spatial <- rep.int(0, q)
  spatial[[1L]] <- sinh(rho)
  mu <- c(cosh(rho), spatial)
  xi <- kappa * mu
  set.seed(seed)
  started <- proc.time()[["elapsed"]]
  sample <- hvmf_validation_exact_sample(n_mc, mu = mu, kappa = kappa)
  sampling_seconds <- proc.time()[["elapsed"]] - started

  score <- hvmf_canonical_score_matrix(sample, xi)
  information <- hvmf_canonical_information(xi)
  score_mean <- colMeans(score)
  score_outer <- crossprod(score) / nrow(score)
  score_mean_se <- sqrt(diag(information) / nrow(score))

  score_jacobian <- numDeriv::jacobian(
    func = function(value) {
      colMeans(hvmf_canonical_score_matrix(sample, value))
    },
    x = xi,
    method = "Richardson",
    method.args = list(eps = 1e-4, d = 1e-4)
  )
  negative_score_jacobian <- -score_jacobian

  J <- diag(c(-1, rep.int(1, q)))
  mean_log_likelihood <- function(value) {
    value <- as.numeric(value)
    value_kappa_sq <- -hvmf_minkowski_inner_product(value, value)
    if (!is.finite(value_kappa_sq) || value_kappa_sq <= 0 || value[[1L]] <= 0) {
      return(-Inf)
    }
    value_kappa <- sqrt(value_kappa_sq)
    hvmf_log_normalizing_constant(q, value_kappa) +
      mean(drop(sample %*% J %*% value))
  }
  hessian <- numDeriv::hessian(
    func = mean_log_likelihood,
    x = xi,
    method = "Richardson",
    method.args = list(eps = 1e-3, d = 1e-3)
  )
  negative_hessian <- -hessian

  matrix_metrics <- function(comparison) {
    difference <- comparison - information
    c(
      max_abs = max(abs(difference)),
      rmse = sqrt(mean(difference^2)),
      relative_frobenius =
        sqrt(sum(difference^2)) / sqrt(sum(information^2))
    )
  }
  outer_metrics <- matrix_metrics(score_outer)
  jacobian_metrics <- matrix_metrics(negative_score_jacobian)
  hessian_metrics <- matrix_metrics(negative_hessian)

  data.frame(
    case = label,
    q = q,
    kappa = kappa,
    rho = rho,
    n_mc = n_mc,
    sampler = attr(sample, "hvmf_sampler") %||% NA_character_,
    sampling_seconds = sampling_seconds,
    max_abs_score_mean = max(abs(score_mean)),
    max_abs_standardized_score_mean = max(abs(score_mean / score_mean_se)),
    score_outer_max_abs_error = outer_metrics[["max_abs"]],
    score_outer_rmse = outer_metrics[["rmse"]],
    score_outer_relative_frobenius = outer_metrics[["relative_frobenius"]],
    negative_score_jacobian_max_abs_error = jacobian_metrics[["max_abs"]],
    negative_score_jacobian_rmse = jacobian_metrics[["rmse"]],
    negative_score_jacobian_relative_frobenius =
      jacobian_metrics[["relative_frobenius"]],
    negative_loglik_hessian_max_abs_error = hessian_metrics[["max_abs"]],
    negative_loglik_hessian_rmse = hessian_metrics[["rmse"]],
    negative_loglik_hessian_relative_frobenius =
      hessian_metrics[["relative_frobenius"]],
    information_min_eigenvalue = min(eigen(information, symmetric = TRUE)$values),
    information_condition_number = base::kappa(information),
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) default else substring(hit[[1L]], nchar(prefix) + 1L)
}
output_path <- get_arg(
  "output",
  file.path("benchmarks", "hvmf_canonical_information_validation.csv")
)
n_h2 <- as.integer(get_arg("n_h2", 200000L))
n_h10 <- as.integer(get_arg("n_h10", 100000L))
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

results <- rbind(
  hvmf_validation_case(
    label = "H2_kappa200_paper",
    q = 2L,
    kappa = 200,
    rho = 0.5,
    n_mc = n_h2,
    seed = 2026073011L
  ),
  hvmf_validation_case(
    label = "H10_kappa10_paper",
    q = 10L,
    kappa = 10,
    rho = 0.35,
    n_mc = n_h10,
    seed = 2026073012L
  )
)
utils::write.csv(results, output_path, row.names = FALSE)
print(results)
cat(sprintf("Wrote %s\n", output_path))
