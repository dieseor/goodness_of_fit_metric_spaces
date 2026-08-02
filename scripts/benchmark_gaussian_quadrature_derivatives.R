#!/usr/bin/env Rscript

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source(file.path("bootstrap", "multiplier_bootstrap.R"))

benchmark_arg <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  value <- args[startsWith(args, prefix)]
  if (!length(value)) default else substring(value[[1L]], nchar(prefix) + 1L)
}

benchmark_fixed_multipliers <- function(weight_matrix) {
  weight_matrix <- as.matrix(weight_matrix)
  index <- 0L
  list(
    name = "fixed paired Exp(1) multipliers",
    mean = 1,
    sd = 1,
    generator = function(n) {
      index <<- index + 1L
      if (index > nrow(weight_matrix) || n != ncol(weight_matrix)) {
        stop("Fixed multiplier matrix was consumed incompatibly.")
      }
      weight_matrix[index, ]
    }
  )
}

benchmark_extract_dot_f <- function(prep, sorted_thresholds) {
  derivative <- prep$D_cvm
  if (is.list(derivative) && !is.null(derivative$derivative_sorted)) {
    return(as.matrix(derivative$derivative_sorted))
  }
  if (!is.list(derivative) || is.null(derivative$aux_order_matrix)) {
    return(as.matrix(derivative))
  }
  score <- as.matrix(prep$Psi_aux)
  do.call(rbind, lapply(seq_len(nrow(derivative$aux_order_matrix)), function(i) {
    cumulative <- apply(
      score[derivative$aux_order_matrix[i, ], , drop = FALSE],
      2L, cumsum
    )
    if (is.null(dim(cumulative))) cumulative <- matrix(cumulative, ncol = ncol(score))
    index <- findInterval(
      sorted_thresholds[i, ], derivative$aux_sorted_distance_matrix[i, ]
    )
    output <- matrix(0, nrow = ncol(sorted_thresholds), ncol = ncol(score))
    positive <- index > 0L
    output[positive, ] <- cumulative[index[positive], , drop = FALSE] /
      nrow(score)
    output
  }))
}

benchmark_error <- function(reference, comparison, prefix) {
  difference <- as.numeric(comparison) - as.numeric(reference)
  setNames(c(max(abs(difference)), sqrt(mean(difference^2))),
           paste0(prefix, c("_max_abs", "_rmse")))
}

benchmark_case <- function(family, q, n, B, derivative_mc_size, seed) {
  Sigma <- toeplitz(0.3^(0:(q - 1L)))
  mu <- seq(-0.25, 0.25, length.out = q)
  set.seed(seed)
  z <- mvtnorm::rmvnorm(n, mu, Sigma)
  data <- if (identical(family, "normal")) {
    z
  } else {
    logistic_gaussian_ilr_to_simplex(z, ambient_dim = q + 1L)
  }
  spec <- if (identical(family, "normal")) {
    make_mvnormal_spec("both")
  } else {
    make_logistic_gaussian_spec("both")
  }

  fit_time <- system.time({
    theta <- spec$fit_theta(data, NULL, list(type = "composite"), list())
  })[["elapsed"]]
  quadrature_control <- list(
    derivative_method = "quadrature",
    gaussian_quadrature_abs_tol = 1e-5,
    gaussian_quadrature_max_terms = 1000000L,
    fast_multiplier_store_paper_quantities = TRUE,
    fast_multiplier_backend = "r"
  )
  score_control <- list(
    derivative_method = "score_mc",
    derivative_mc_size = derivative_mc_size,
    derivative_mc_seed = seed + 1L,
    fast_multiplier_store_paper_quantities = TRUE,
    fast_multiplier_backend = "r"
  )
  profile_time <- system.time({
    ks_prep <- prepare_ks_observed_data(
      data, spec, theta, make_sample_unique_distance_ks_grid(),
      control = quadrature_control, light = TRUE
    )
    cvm_prep <- prepare_cvm_observed_data_from_sample_ks(
      data, spec, theta, ks_prep, quadrature_control
    )
  })[["elapsed"]]

  prepare_one <- function(control) {
    if (identical(family, "normal")) {
      prepare_mvnormal_fast_multiplier(
        spec, data, theta, ks_prep, cvm_prep, control, "both"
      )
    } else {
      prepare_logistic_gaussian_fast_multiplier(
        spec, data, theta, ks_prep, cvm_prep, control, "both"
      )
    }
  }
  quadrature_time <- system.time({
    quadrature <- prepare_one(quadrature_control)
  })[["elapsed"]]
  score_time <- system.time({
    score_mc <- prepare_one(score_control)
  })[["elapsed"]]
  dot_quadrature <- benchmark_extract_dot_f(
    quadrature, ks_prep$sorted_distance_matrix
  )
  dot_score <- benchmark_extract_dot_f(score_mc, ks_prep$sorted_distance_matrix)

  correction_time <- system.time({
    operator_quadrature <- dot_quadrature %*%
      t(solve(quadrature$paper_Vhat))
    operator_score <- dot_score %*% t(solve(score_mc$paper_Vhat))
    correction_quadrature <- -quadrature$paper_score_obs %*%
      t(solve(quadrature$paper_Vhat)) %*% t(dot_quadrature)
    correction_score <- -score_mc$paper_score_obs %*%
      t(solve(score_mc$paper_Vhat)) %*% t(dot_score)
  })[["elapsed"]]

  set.seed(seed + 2L)
  multiplier_matrix <- matrix(rexp(B * n), nrow = B, ncol = n)
  run_one <- function(control) {
    arguments <- list(
      data = data,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = B,
      multipliers = benchmark_fixed_multipliers(multiplier_matrix),
      bootstrap_method = "fast_multiplier",
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = FALSE
      ),
      control = control
    )
    if (identical(family, "normal")) {
      do.call(multiplier_bootstrap_mvnormal, c(
        arguments,
        list(unknown_param = "both", fast_multiplier_backend = "r")
      ))
    } else {
      do.call(multiplier_bootstrap_logistic_gaussian, c(
        arguments, list(unknown_param = "both")
      ))
    }
  }
  quadrature_bootstrap_time <- system.time({
    quadrature_result <- run_one(quadrature_control)
  })[["elapsed"]]
  score_bootstrap_time <- system.time({
    score_result <- run_one(score_control)
  })[["elapsed"]]

  metrics <- c(
    benchmark_error(dot_quadrature, dot_score, "dot_f"),
    benchmark_error(quadrature$paper_Vhat, score_mc$paper_Vhat, "Vhat"),
    benchmark_error(operator_quadrature, operator_score, "dot_f_Vhat_inverse"),
    benchmark_error(correction_quadrature, correction_score, "fast_correction"),
    benchmark_error(
      quadrature_result$bootstrap$statistics$ks,
      score_result$bootstrap$statistics$ks,
      "ks_replicates"
    ),
    benchmark_error(
      quadrature_result$bootstrap$statistics$cvm,
      score_result$bootstrap$statistics$cvm,
      "cvm_replicates"
    )
  )
  row <- data.frame(
    family = family, q = q, simplex_parts = if (family == "lg") q + 1L else NA,
    n = n, B = B, derivative_mc_size = derivative_mc_size,
    fit_seconds = fit_time, profile_preparation_seconds = profile_time,
    quadrature_F_dotF_seconds = quadrature_time,
    score_mc_dotF_seconds = score_time,
    correction_seconds = correction_time,
    quadrature_bootstrap_total_seconds = quadrature_bootstrap_time,
    score_mc_bootstrap_total_seconds = score_bootstrap_time,
    quadrature_terms_max = quadrature$quadrature_diagnostics$max_terms_used,
    quadrature_residual_error_estimate =
      quadrature$quadrature_diagnostics$max_residual_error_estimate,
    t(metrics),
    ks_critical_quadrature = quadrature_result$inference$ks$critical_value,
    ks_critical_score_mc = score_result$inference$ks$critical_value,
    ks_p_quadrature = quadrature_result$inference$ks$p_value,
    ks_p_score_mc = score_result$inference$ks$p_value,
    cvm_critical_quadrature = quadrature_result$inference$cvm$critical_value,
    cvm_critical_score_mc = score_result$inference$cvm$critical_value,
    cvm_p_quadrature = quadrature_result$inference$cvm$p_value,
    cvm_p_score_mc = score_result$inference$cvm$p_value,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(row = row, quadrature = quadrature_result, score_mc = score_result)
}

n_values <- as.integer(strsplit(benchmark_arg("n_values", "50,200"), ",",
                                fixed = TRUE)[[1L]])
dimensions <- as.integer(strsplit(benchmark_arg("dimensions", "2,10"), ",",
                                  fixed = TRUE)[[1L]])
families <- strsplit(benchmark_arg("families", "normal,lg"), ",",
                     fixed = TRUE)[[1L]]
B <- as.integer(benchmark_arg("B", 199L))
derivative_mc_size <- as.integer(benchmark_arg("derivative_mc_size", 10000L))
n_cores <- as.integer(benchmark_arg("n_cores", 1L))
output_dir <- benchmark_arg(
  "output_dir",
  file.path("simulation_results", "gaussian_quadrature_validation")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

case_grid <- expand.grid(
  family = families, q = dimensions, n = n_values,
  stringsAsFactors = FALSE
)
run_index <- function(index) {
  item <- case_grid[index, ]
  message(sprintf(
    "[Gaussian quadrature benchmark] family=%s q=%d n=%d",
    item$family, item$q, item$n
  ))
  benchmark_case(
    item$family, item$q, item$n, B, derivative_mc_size,
    seed = 2026080100L + 1000L * item$n + 100L * item$q + index
  )
}
if (.Platform$OS.type != "windows" && n_cores > 1L) {
  cases <- parallel::mclapply(
    seq_len(nrow(case_grid)), run_index,
    mc.cores = min(n_cores, nrow(case_grid)), mc.set.seed = FALSE
  )
} else {
  cases <- lapply(seq_len(nrow(case_grid)), run_index)
}
summary <- do.call(rbind, lapply(cases, `[[`, "row"))
utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
saveRDS(cases, file.path(output_dir, "paired_results.rds"))
print(summary)
