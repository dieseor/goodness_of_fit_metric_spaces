#!/usr/bin/env Rscript

# Audit the fitted objects in the paper's composite-null correction
#
#   y_{omega,t} + dot F_{omega}^{theta_hat}(t)^T
#                    Vhat^{-1} psi_{theta_hat}.
#
# This deliberately excludes the distance-profile derivative dot F, whose
# numerical audit is separate.  It checks theta_hat, psi_{theta_hat}, Vhat,
# Vhat^{-1}, and the equivalent fitted-MLE influence representation used by
# the fast Normal and logistic-Gaussian branches.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[pieces[[1L]]]] <- if (length(pieces) == 1L) "TRUE" else {
      paste(pieces[-1L], collapse = "=")
    }
  }
  out
}

parse_integer_csv <- function(value, default) {
  if (is.null(value) || !nzchar(value)) return(as.integer(default))
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
}

parse_character_csv <- function(value, default) {
  if (is.null(value) || !nzchar(value)) return(as.character(default))
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
M <- as.integer(args$M %||% 500L)
cores <- as.integer(args$cores %||% 2L)
dimensions <- parse_integer_csv(args$dimensions, c(2L, 5L, 10L))
n_values <- parse_integer_csv(args$n_values, c(50L, 100L, 200L, 400L))
families <- unique(tolower(parse_character_csv(args$families, c("normal", "lg", "vmf"))))
seed <- as.integer(args$seed %||% 20260804L)
output_dir <- args$output_dir %||% file.path(
  "simulation_results", "section6_new_scenarios",
  sprintf("paper_quantities_audit_M%d", M)
)

if (!is.finite(M) || M < 2L || !is.finite(cores) || cores < 1L) {
  stop("`M` must be at least two and `cores` must be positive.")
}
if (!all(families %in% c("normal", "lg", "vmf"))) {
  stop("`families` must be a subset of normal, lg, vmf.")
}

source(file.path("bootstrap", "multiplier_bootstrap.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

e1 <- function(d) c(1, rep.int(0, d - 1L))
section6_sigma_plus <- function(d) {
  Sigma <- diag(d)
  Sigma[1L, 2L] <- Sigma[2L, 1L] <- 0.75
  Sigma
}

gaussian_score_matrix <- function(x, mu, Sigma) {
  d <- ncol(x)
  Sigma_inv <- solve(Sigma)
  centered <- sweep(x, 2L, mu, FUN = "-")
  score_mu <- centered %*% t(Sigma_inv)
  score_sigma <- t(vapply(seq_len(nrow(x)), function(i) {
    rr <- centered[i, , drop = FALSE]
    fast_multiplier_sym_score_to_vech(
      0.5 * (Sigma_inv %*% crossprod(rr) %*% Sigma_inv - Sigma_inv)
    )
  }, numeric(d * (d + 1L) / 2L)))
  cbind(score_mu, score_sigma)
}

relative_frobenius_error <- function(x, y) {
  sqrt(sum((x - y)^2)) / sqrt(sum(y^2))
}

quantile95 <- function(x) stats::quantile(x, 0.95, names = FALSE, type = 8L)

run_replicates <- function(M, cores, worker, seed_offset) {
  seeds <- as.integer((seed + seed_offset + 1009L * seq_len(M)) %% 2147483646L + 1L)
  if (.Platform$OS.type == "unix" && cores > 1L) {
    return(parallel::mclapply(
      seq_len(M), function(i) worker(seeds[[i]]),
      mc.cores = min(cores, M), mc.set.seed = FALSE
    ))
  }
  lapply(seq_len(M), function(i) worker(seeds[[i]]))
}

one_gaussian_fit <- function(family, n, mu, Sigma, seed_value) {
  set.seed(seed_value)
  if (identical(family, "normal")) {
    x <- mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
    theta <- fit_mvnormal_theta(
      x, null = list(type = "composite"), unknown_param = "both"
    )
    z <- x
    mu_hat <- theta$mu
    Sigma_hat <- theta$Sigma
  } else {
    x <- rlogistic_gaussian_simplex(n, mu_ilr = mu, Sigma_ilr = Sigma)
    theta <- fit_logistic_gaussian_theta(
      x, null = list(type = "composite"), unknown_param = "both"
    )
    z <- logistic_gaussian_ilr_matrix(x)
    mu_hat <- theta$mu_ilr
    Sigma_hat <- theta$Sigma_ilr
  }

  score <- gaussian_score_matrix(z, mu_hat, Sigma_hat)
  paper_Vhat <- fast_multiplier_gaussian_paper_vhat(Sigma_hat, "both")
  paper_Vtrue <- fast_multiplier_gaussian_paper_vhat(Sigma, "both")
  paper_influence <- -score %*% t(solve(paper_Vhat))
  influence_closed_form <- cbind(
    sweep(z, 2L, mu_hat, FUN = "-"),
    t(vapply(seq_len(nrow(z)), function(i) {
      rr <- sweep(z[i, , drop = FALSE], 2L, mu_hat, FUN = "-")
      fast_multiplier_vech(crossprod(rr) - Sigma_hat)
    }, numeric(ncol(z) * (ncol(z) + 1L) / 2L)))
  )
  diagnostics <- fast_multiplier_matrix_condition_diagnostics(paper_Vhat)

  list(
    score_mean_max_abs = max(abs(colMeans(score))),
    influence_identity_max_abs = max(abs(paper_influence - influence_closed_form)),
    Vhat_relative_error = relative_frobenius_error(paper_Vhat, paper_Vtrue),
    Vhat_inverse_relative_error = relative_frobenius_error(
      solve(paper_Vhat), solve(paper_Vtrue)
    ),
    Vhat_condition_number = diagnostics$condition_number,
    Vhat_rcond = diagnostics$rcond,
    min_abs_Vhat_eigenvalue = min(abs(diagnostics$eigenvalues)),
    mu_error_l2 = sqrt(sum((mu_hat - mu)^2)),
    Sigma_relative_error = relative_frobenius_error(Sigma_hat, Sigma),
    xi_error_l2 = NA_real_
  )
}

one_vmf_fit <- function(n, mu, kappa, seed_value) {
  set.seed(seed_value)
  x <- normalize_vmf_data(rotasym::r_vMF(n, mu = mu, kappa = kappa))
  theta <- fit_vmf_theta(
    x, null = list(type = "composite"), unknown_param = "xi"
  )
  q <- length(mu) - 1L
  score <- t(vapply(seq_len(nrow(x)), function(i) {
    psi_xi(x[i, ], theta$xi, q)
  }, numeric(length(mu))))
  paper_Vhat <- dot_psi_xi(theta$xi, q)
  paper_Vtrue <- dot_psi_xi(kappa * mu, q)
  information_hat <- -paper_Vhat
  influence_from_paper <- -score %*% t(solve(paper_Vhat))
  influence_from_fast_metric <- score %*% t(solve(information_hat))
  diagnostics <- fast_multiplier_matrix_condition_diagnostics(paper_Vhat)

  list(
    score_mean_max_abs = max(abs(colMeans(score))),
    influence_identity_max_abs = max(abs(influence_from_paper - influence_from_fast_metric)),
    Vhat_relative_error = relative_frobenius_error(paper_Vhat, paper_Vtrue),
    Vhat_inverse_relative_error = relative_frobenius_error(
      solve(paper_Vhat), solve(paper_Vtrue)
    ),
    Vhat_condition_number = diagnostics$condition_number,
    Vhat_rcond = diagnostics$rcond,
    min_abs_Vhat_eigenvalue = min(abs(diagnostics$eigenvalues)),
    mu_error_l2 = sqrt(sum((theta$mu - mu)^2)),
    Sigma_relative_error = NA_real_,
    xi_error_l2 = sqrt(sum((theta$xi - kappa * mu)^2))
  )
}

rows <- list()
row_index <- 1L
cell_index <- 0L
total_cells <- length(families) * 2L * length(dimensions) * length(n_values)
for (family in c("normal", "lg")) {
  if (!family %in% families) next
  for (scenario in c("mixture_null", "t3_null")) {
    for (d in dimensions) {
      mu <- if (identical(scenario, "mixture_null")) 0.5 * e1(d) else rep.int(0, d)
      Sigma <- if (identical(scenario, "mixture_null")) section6_sigma_plus(d) else diag(d)
      for (n in n_values) {
        cell_index <- cell_index + 1L
        offset <- 100000L * match(family, c("normal", "lg")) +
          10000L * match(scenario, c("mixture_null", "t3_null")) + 100L * d + n
        fits <- run_replicates(M, cores, function(s) {
          one_gaussian_fit(family, n, mu, Sigma, s)
        }, offset)
        get <- function(name) vapply(fits, `[[`, numeric(1), name)
        rows[[row_index]] <- data.frame(
          family = family,
          scenario = scenario,
          d = d,
          n = n,
          M = M,
          paper_Vhat_method = "analytic_expected_score_jacobian",
          median_Vhat_relative_error = stats::median(get("Vhat_relative_error")),
          q95_Vhat_relative_error = quantile95(get("Vhat_relative_error")),
          median_Vhat_inverse_relative_error = stats::median(get("Vhat_inverse_relative_error")),
          q95_Vhat_inverse_relative_error = quantile95(get("Vhat_inverse_relative_error")),
          median_Vhat_condition_number = stats::median(get("Vhat_condition_number")),
          max_Vhat_condition_number = max(get("Vhat_condition_number")),
          min_Vhat_rcond = min(get("Vhat_rcond")),
          min_abs_Vhat_eigenvalue = min(get("min_abs_Vhat_eigenvalue")),
          max_abs_mean_score = max(get("score_mean_max_abs")),
          max_influence_identity_error = max(get("influence_identity_max_abs")),
          median_mu_error_l2 = stats::median(get("mu_error_l2")),
          median_Sigma_relative_error = stats::median(get("Sigma_relative_error")),
          median_xi_error_l2 = NA_real_,
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
        message(sprintf("completed %d/%d cells", cell_index, total_cells))
      }
    }
  }
}

if ("vmf" %in% families) for (scenario in c("antipodal_null", "projected_normal_null")) {
  for (d in dimensions) {
    mu <- e1(d + 1L)
    kappa <- if (identical(scenario, "antipodal_null")) d else 1.5 * d
    for (n in n_values) {
      cell_index <- cell_index + 1L
      offset <- 300000L + 10000L * match(
        scenario, c("antipodal_null", "projected_normal_null")
      ) + 100L * d + n
      fits <- run_replicates(M, cores, function(s) {
        one_vmf_fit(n, mu, kappa, s)
      }, offset)
      get <- function(name) vapply(fits, `[[`, numeric(1), name)
      rows[[row_index]] <- data.frame(
        family = "vmf",
        scenario = scenario,
        d = d,
        n = n,
        M = M,
        paper_Vhat_method = "analytic_expected_score_jacobian",
        median_Vhat_relative_error = stats::median(get("Vhat_relative_error")),
        q95_Vhat_relative_error = quantile95(get("Vhat_relative_error")),
        median_Vhat_inverse_relative_error = stats::median(get("Vhat_inverse_relative_error")),
        q95_Vhat_inverse_relative_error = quantile95(get("Vhat_inverse_relative_error")),
        median_Vhat_condition_number = stats::median(get("Vhat_condition_number")),
        max_Vhat_condition_number = max(get("Vhat_condition_number")),
        min_Vhat_rcond = min(get("Vhat_rcond")),
        min_abs_Vhat_eigenvalue = min(get("min_abs_Vhat_eigenvalue")),
        max_abs_mean_score = max(get("score_mean_max_abs")),
        max_influence_identity_error = max(get("influence_identity_max_abs")),
        median_mu_error_l2 = stats::median(get("mu_error_l2")),
        median_Sigma_relative_error = NA_real_,
        median_xi_error_l2 = stats::median(get("xi_error_l2")),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
      message(sprintf("completed %d/%d cells", cell_index, total_cells))
    }
  }
}

summary <- do.call(rbind, rows)
utils::write.csv(summary, file.path(output_dir, "paper_quantities_summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(
  M = M,
  cores = cores,
  families = paste(families, collapse = ","),
  dimensions = paste(dimensions, collapse = ","),
  n_values = paste(n_values, collapse = ","),
  seed = seed,
  dot_F_audited = FALSE,
  dot_F_reason = "distance-profile derivative deliberately excluded",
  stringsAsFactors = FALSE
), file.path(output_dir, "config.csv"), row.names = FALSE)
message("Wrote paper-quantity audit to: ", normalizePath(output_dir))
