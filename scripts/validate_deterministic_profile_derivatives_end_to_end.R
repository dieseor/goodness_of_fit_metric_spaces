#!/usr/bin/env Rscript

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("bootstrap/multiplier_bootstrap.R")

fixed_multiplier_spec <- function(weight_matrix) {
  weight_matrix <- as.matrix(weight_matrix)
  draw_index <- 0L
  list(
    name = "fixed paired Exp(1) draws",
    mean = 1,
    sd = 1,
    generator = function(n) {
      draw_index <<- draw_index + 1L
      if (draw_index > nrow(weight_matrix) || n != ncol(weight_matrix)) {
        stop("The fixed multiplier matrix was consumed with incompatible dimensions.")
      }
      weight_matrix[draw_index, ]
    }
  )
}

extract_complete_dot_f <- function(fast_prep, sorted_thresholds) {
  derivative <- fast_prep$D_cvm
  if (is.list(derivative) && !is.null(derivative$derivative_sorted)) {
    derivative <- derivative$derivative_sorted
  }
  if (is.list(derivative) && !is.null(derivative$aux_order_matrix)) {
    order_matrix <- derivative$aux_order_matrix
    sorted_aux_distances <- derivative$aux_sorted_distance_matrix
    score <- as.matrix(fast_prep$Psi_aux)
    derivative <- do.call(rbind, lapply(seq_len(nrow(order_matrix)), function(i) {
      cumulative_score <- apply(
        score[order_matrix[i, ], , drop = FALSE],
        2L,
        cumsum
      )
      if (is.null(dim(cumulative_score))) {
        cumulative_score <- matrix(cumulative_score, ncol = ncol(score))
      }
      threshold_index <- findInterval(
        sorted_thresholds[i, ],
        sorted_aux_distances[i, ]
      )
      output <- matrix(
        0,
        nrow = ncol(sorted_thresholds),
        ncol = ncol(score)
      )
      positive <- threshold_index > 0L
      output[positive, ] <- cumulative_score[
        threshold_index[positive],
        ,
        drop = FALSE
      ] / nrow(score)
      output
    }))
  }
  as.matrix(derivative)
}

paired_error_metrics <- function(reference, comparison) {
  difference <- as.numeric(comparison) - as.numeric(reference)
  c(
    max_abs = max(abs(difference)),
    rmse = sqrt(mean(difference^2)),
    mean_difference = mean(difference)
  )
}

hvmf_xi_jacobian_for_score_mc <- function(theta_hat) {
  q <- theta_hat$q
  kappa <- theta_hat$kappa
  mu <- theta_hat$mu
  if (q == 2L) {
    chi <- theta_hat$chi
    theta_angle <- theta_hat$theta
    dmu_dchi <- c(
      sinh(chi),
      cosh(chi) * cos(theta_angle),
      cosh(chi) * sin(theta_angle)
    )
    dmu_dtheta <- c(
      0,
      -sinh(chi) * sin(theta_angle),
      sinh(chi) * cos(theta_angle)
    )
    return(cbind(kappa * dmu_dchi, kappa * dmu_dtheta, mu))
  }

  eta <- mu[-1L]
  mu0 <- mu[[1L]]
  jacobian <- matrix(0, nrow = q + 1L, ncol = q + 1L)
  jacobian[1L, seq_len(q)] <- kappa * eta / mu0
  jacobian[-1L, seq_len(q)] <- kappa * diag(q)
  jacobian[, q + 1L] <- kappa * mu
  jacobian
}

prepare_complete_derivatives <- function(model, data, theta_hat, control) {
  if (identical(model, "vmf")) {
    spec <- make_vmf_spec(distance_type = "geodesic", unknown_param = "xi")
  } else {
    spec <- make_hvmf_spec(unknown_param = "both")
  }
  ks_prep <- prepare_ks_observed_data(
    data = data,
    spec = spec,
    theta_hat = theta_hat,
    ks_grid = make_sample_unique_distance_ks_grid(),
    control = control,
    light = TRUE
  )
  cvm_prep <- prepare_cvm_observed_data_from_sample_ks(
    data = data,
    spec = spec,
    theta_hat = theta_hat,
    ks_prep = ks_prep,
    control = control
  )
  fast_prep <- if (identical(model, "vmf")) {
    prepare_vmf_fast_multiplier(
      data = data,
      theta_hat = theta_hat,
      spec = spec,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      control = control,
      distance_type = "geodesic"
    )
  } else {
    prepare_hvmf_fast_multiplier(
      spec = spec,
      data = data,
      theta_hat = theta_hat,
      ks_prep = ks_prep,
      cvm_prep = cvm_prep,
      control = control
    )
  }
  list(
    dot_f = extract_complete_dot_f(
      fast_prep,
      sorted_thresholds = ks_prep$sorted_distance_matrix
    ),
    fast_prep = fast_prep
  )
}

run_paired_case <- function(model,
                            data,
                            B,
                            multiplier_matrix,
                            derivative_mc_size,
                            derivative_mc_seed,
                            grid_size) {
  n <- nrow(data)
  spec <- if (identical(model, "vmf")) {
    make_vmf_spec(distance_type = "geodesic", unknown_param = "xi")
  } else {
    make_hvmf_spec(unknown_param = "both")
  }
  theta_hat <- spec$fit_theta(
    data = data,
    weights = NULL,
    null = list(type = "composite"),
    control = list()
  )
  profile_grid_control <- if (identical(model, "vmf")) {
    list(vmf_profile_n_u = as.integer(grid_size))
  } else {
    list(hvmf_profile_n_y = as.integer(grid_size))
  }
  quadrature_control <- c(
    list(
      derivative_method = "quadrature",
      fast_multiplier_backend = "r",
      fast_multiplier_cvm_block_size = 50L
    ),
    profile_grid_control
  )
  score_mc_control <- c(
    list(
      derivative_method = "score_mc",
      derivative_mc_size = as.integer(derivative_mc_size),
      derivative_mc_seed = as.integer(derivative_mc_seed),
      fast_multiplier_backend = "r",
      fast_multiplier_cvm_block_size = 50L
    ),
    profile_grid_control
  )

  quadrature_derivative <- prepare_complete_derivatives(
    model, data, theta_hat, quadrature_control
  )
  score_mc_derivative <- prepare_complete_derivatives(
    model, data, theta_hat, score_mc_control
  )
  if (!identical(dim(quadrature_derivative$dot_f), dim(score_mc_derivative$dot_f))) {
    stop("Paired complete derivative matrices have incompatible dimensions.")
  }

  common <- list(
    data = data,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = as.integer(B),
    alpha = 0.05,
    n_cores = 1L,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    )
  )
  run_method <- function(method_control) {
    method_args <- c(
      common,
      list(
        multipliers = fixed_multiplier_spec(multiplier_matrix),
        control = method_control
      )
    )
    if (identical(model, "vmf")) {
      do.call(multiplier_bootstrap_vmf, c(
        method_args,
        list(distance_type = "geodesic", unknown_param = "xi")
      ))
    } else {
      do.call(multiplier_bootstrap_hvmf, c(
        method_args,
        list(
          unknown_param = "both",
          fast_multiplier_backend = "r",
          fuse_ks_cvm = TRUE
        )
      ))
    }
  }
  quadrature_result <- run_method(quadrature_control)
  score_mc_result <- run_method(score_mc_control)

  theta_quadrature <- as.numeric(quadrature_result$observed$theta_hat$xi)
  theta_score_mc <- as.numeric(score_mc_result$observed$theta_hat$xi)
  if (max(abs(theta_quadrature - theta_score_mc)) > 1e-12 ||
      max(abs(theta_quadrature - as.numeric(theta_hat$xi))) > 1e-12) {
    stop("The paired methods did not use the same canonical fitted parameter.")
  }

  quadrature_dot_f_score_coordinates <-
    if (identical(model, "hvmf")) {
      quadrature_derivative$dot_f %*%
        hvmf_xi_jacobian_for_score_mc(theta_hat)
    } else {
      quadrature_derivative$dot_f
    }
  dot_metrics <- paired_error_metrics(
    quadrature_dot_f_score_coordinates,
    score_mc_derivative$dot_f
  )
  dot_row <- data.frame(
    model = model,
    n = n,
    B = B,
    derivative_mc_size = derivative_mc_size,
    grid_size = grid_size,
    comparison = "complete_dot_F_in_score_mc_coordinates",
    statistic = "all",
    n_values_compared = length(quadrature_derivative$dot_f),
    max_abs_difference = dot_metrics[["max_abs"]],
    rmse_difference = dot_metrics[["rmse"]],
    mean_difference = dot_metrics[["mean_difference"]],
    observed_quadrature = NA_real_,
    observed_score_mc = NA_real_,
    critical_quadrature = NA_real_,
    critical_score_mc = NA_real_,
    p_value_quadrature = NA_real_,
    p_value_score_mc = NA_real_,
    stringsAsFactors = FALSE
  )

  bootstrap_rows <- lapply(c("ks", "cvm"), function(statistic) {
    quadrature_replicates <- quadrature_result$bootstrap$statistics[[statistic]]
    score_mc_replicates <- score_mc_result$bootstrap$statistics[[statistic]]
    metrics <- paired_error_metrics(quadrature_replicates, score_mc_replicates)
    data.frame(
      model = model,
      n = n,
      B = B,
      derivative_mc_size = derivative_mc_size,
      grid_size = grid_size,
      comparison = "bootstrap_replicates",
      statistic = statistic,
      n_values_compared = length(quadrature_replicates),
      max_abs_difference = metrics[["max_abs"]],
      rmse_difference = metrics[["rmse"]],
      mean_difference = metrics[["mean_difference"]],
      observed_quadrature = quadrature_result$inference[[statistic]]$observed,
      observed_score_mc = score_mc_result$inference[[statistic]]$observed,
      critical_quadrature = quadrature_result$inference[[statistic]]$critical_value,
      critical_score_mc = score_mc_result$inference[[statistic]]$critical_value,
      p_value_quadrature = quadrature_result$inference[[statistic]]$p_value,
      p_value_score_mc = score_mc_result$inference[[statistic]]$p_value,
      stringsAsFactors = FALSE
    )
  })

  list(
    summary = do.call(rbind, c(list(dot_row), bootstrap_rows)),
    dot_f = list(
      quadrature_xi = quadrature_derivative$dot_f,
      score_mc_historical_coordinates = score_mc_derivative$dot_f,
      quadrature_in_score_mc_coordinates =
        quadrature_dot_f_score_coordinates
    ),
    results = list(
      quadrature = quadrature_result,
      score_mc = score_mc_result
    ),
    theta_hat = theta_hat
  )
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) default else substring(hit[[1L]], nchar(prefix) + 1L)
}

n <- as.integer(get_arg("n", 200L))
B <- as.integer(get_arg("B", 999L))
derivative_mc_size <- as.integer(get_arg("derivative_mc_size", 1000L))
grid_size <- as.integer(get_arg("grid_size", 4097L))
output_path <- get_arg(
  "output",
  file.path("benchmarks", "deterministic_profile_derivative_end_to_end.csv")
)
details_path <- sub("\\.csv$", ".rds", output_path)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

set.seed(2026073001L)
vmf_data <- rotasym::r_vMF(n, mu = c(1, 0, 0), kappa = 2)
hvmf_mu <- c(cosh(0.5), sinh(0.5), 0)
hvmf_data <- rhvmf_h2_polar(n, mu = hvmf_mu, kappa = 200)
set.seed(2026073002L)
multiplier_matrix <- matrix(stats::rexp(B * n), nrow = B, ncol = n)

cases <- list(
  vmf = run_paired_case(
    "vmf", vmf_data, B, multiplier_matrix,
    derivative_mc_size, 2026073003L, grid_size
  ),
  hvmf = run_paired_case(
    "hvmf", hvmf_data, B, multiplier_matrix,
    derivative_mc_size, 2026073004L, grid_size
  )
)
summary <- do.call(rbind, lapply(cases, `[[`, "summary"))
utils::write.csv(summary, output_path, row.names = FALSE)
saveRDS(
  list(
    configuration = list(
      n = n,
      B = B,
      derivative_mc_size = derivative_mc_size,
      grid_size = grid_size,
      sample_seed = 2026073001L,
      multiplier_seed = 2026073002L,
      vmf_derivative_seed = 2026073003L,
      hvmf_derivative_seed = 2026073004L,
      multiplier_matrix = multiplier_matrix
    ),
    cases = cases,
    summary = summary
  ),
  details_path
)
print(summary)
cat(sprintf("Wrote %s and %s\n", output_path, details_path))
