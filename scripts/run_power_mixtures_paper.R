#!/usr/bin/env Rscript

resolve_paper_power_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_paper_power_path("scripts", "run_power_mixtures_pilot.R"))

paper_repo_root <- normalizePath(
  dirname(dirname(resolve_paper_power_path("scripts", "run_power_mixtures_pilot.R"))),
  winslash = "/",
  mustWork = TRUE
)

hyperboloid_point_h2 <- function(chi, theta_angle) {
  c(
    cosh(chi),
    sinh(chi) * cos(theta_angle),
    sinh(chi) * sin(theta_angle)
  )
}

paper_scenario_catalog <- function() {
  mu_hvmf_null <- c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2))
  sigma_case2_alt <- matrix(0.85, nrow = 3L, ncol = 3L)
  diag(sigma_case2_alt) <- 1

  list(
    vmf_s2_antipodal = list(
      model = "vmf",
      alternative = "antipodal_mixture",
      mu0 = c(1, 0, 0),
      mu1 = c(-1, 0, 0),
      kappa = 2
    ),
    hvmf_h2_case1 = list(
      model = "hvmf",
      alternative = "halfway_mixture_case1",
      mu0 = mu_hvmf_null,
      mu1 = c(sqrt(2), -1 / sqrt(2), 1 / sqrt(2)),
      kappa = 200
    ),
    hvmf_h2_case2 = list(
      model = "hvmf",
      alternative = "halfway_mixture_case2",
      mu0 = mu_hvmf_null,
      mu1 = c(sqrt(2), -1 / sqrt(2), -1 / sqrt(2)),
      kappa = 200
    ),
    mvnormal_d3_equal_sigma_opposite_mean = list(
      model = "mvnormal",
      alternative = "halfway_mixture_equal_sigma_opposite_mean",
      mu0 = c(1.25, -1.0, 0.75),
      Sigma0 = diag(3L),
      mu1 = c(-1.25, 1.0, -0.75),
      Sigma1 = diag(3L)
    ),
    mvnormal_d3_equal_mean_distinct_sigma = list(
      model = "mvnormal",
      alternative = "halfway_mixture_equal_mean_distinct_sigma",
      mu0 = c(0, 0, 0),
      Sigma0 = diag(3L),
      mu1 = c(0, 0, 0),
      Sigma1 = sigma_case2_alt
    ),
    mvnormal_d2_moderate_location_correlation = list(
      model = "mvnormal",
      alternative = "halfway_mixture_opposite_mean_correlation",
      mu0 = c(1, -1),
      Sigma0 = matrix(c(1, 0.25, 0.25, 1), nrow = 2L),
      mu1 = c(-1, 1),
      Sigma1 = matrix(c(1, -0.25, -0.25, 1), nrow = 2L)
    ),
    mvnormal_d2_small_location_rho050 = list(
      model = "mvnormal",
      alternative = "halfway_mixture_small_opposite_mean_rho_0_50",
      mu0 = c(0, 0.5),
      Sigma0 = matrix(c(1, 0.5, 0.5, 1), nrow = 2L),
      mu1 = c(0, -0.5),
      Sigma1 = matrix(c(1, -0.5, -0.5, 1), nrow = 2L)
    ),
    mvnormal_d2_small_location_rho075 = list(
      model = "mvnormal",
      alternative = "halfway_mixture_small_opposite_mean_rho_0_75",
      mu0 = c(0, 0.5),
      Sigma0 = matrix(c(1, 0.75, 0.75, 1), nrow = 2L),
      mu1 = c(0, -0.5),
      Sigma1 = matrix(c(1, -0.75, -0.75, 1), nrow = 2L)
    ),
    mvnormal_d2_small_location_rho095 = list(
      model = "mvnormal",
      alternative = "halfway_mixture_small_opposite_mean_rho_0_95",
      mu0 = c(0, 0.5),
      Sigma0 = matrix(c(1, 0.95, 0.95, 1), nrow = 2L),
      mu1 = c(0, -0.5),
      Sigma1 = matrix(c(1, -0.95, -0.95, 1), nrow = 2L)
    )
  )
}

paper_table_catalog <- function() {
  list(
    vmf_s2_antipodal_paper = list(
      caption = paste(
        "Empirical rejection rates at level 0.05 for the composite GOF test under",
        "$(1-\\beta/2)\\,\\mathrm{vMF}(\\mu,2)+(\\beta/2)\\,\\mathrm{vMF}(-\\mu,2)$",
        "on $\\mathbb{S}^2$, with $\\mu=(1,0,0)^\\top$, $M=1000$, and $B=5000$."
      ),
      label = "tab:power-vmf-s2-antipodal-paper"
    ),
    logistic_gaussian_simplex_d3_dirichlet_1_1_8_paper = list(
      caption = paste(
        "Empirical rejection rates at level 0.05 for the composite GOF test under",
        "$(1-\\beta)\\,\\mathrm{LG}(0,I_3)+\\beta\\,\\mathrm{Dir}(1,1,8)$",
        "on $\\Delta^2$, with $M=1000$, and $B=5000$."
      ),
      label = "tab:power-lg-dir-118-paper"
    ),
    hvmf_h2_case1 = list(
      caption = paste(
        "Empirical rejection rates at level 0.05 for the composite GOF test under",
        "$(1-\\beta/2)\\,\\mathrm{HvMF}(\\mu_0,200)+(\\beta/2)\\,\\mathrm{HvMF}(\\mu_1,200)$",
        "on $\\mathbb{H}^2$, with $\\mu_0=(\\sqrt{2},1/\\sqrt{2},1/\\sqrt{2})^\\top$ and",
        "$\\mu_1=(\\sqrt{2},-1/\\sqrt{2},1/\\sqrt{2})^\\top$, $M=1000$, and $B=5000$."
      ),
      label = "tab:power-hvmf-h2-case1"
    ),
    hvmf_h2_case2 = list(
      caption = paste(
        "Empirical rejection rates at level 0.05 for the composite GOF test under",
        "$(1-\\beta/2)\\,\\mathrm{HvMF}(\\mu_0,200)+(\\beta/2)\\,\\mathrm{HvMF}(\\mu_1,200)$",
        "on $\\mathbb{H}^2$, with $\\mu_0=(\\sqrt{2},1/\\sqrt{2},1/\\sqrt{2})^\\top$ and",
        "$\\mu_1=(\\sqrt{2},-1/\\sqrt{2},-1/\\sqrt{2})^\\top$, $M=1000$, and $B=5000$."
      ),
      label = "tab:power-hvmf-h2-case2"
    ),
    mvnormal_d3_equal_sigma_opposite_mean = list(
      caption = paste(
        "Empirical rejection rates at level 0.05 for the composite GOF test under",
        "$(1-\\beta/2)\\,\\mathcal{N}_3(\\mu,I_3)+(\\beta/2)\\,\\mathcal{N}_3(-\\mu,I_3)$",
        "in $\\mathbb{R}^3$, with $\\mu=(1.25,-1,0.75)^\\top$, $M=1000$, and $B=5000$."
      ),
      label = "tab:power-mvnormal-equal-sigma-opposite-mean"
    ),
    mvnormal_d3_equal_mean_distinct_sigma = list(
      caption = paste(
        "Empirical rejection rates at level 0.05 for the composite GOF test under",
        "$(1-\\beta/2)\\,\\mathcal{N}_3(0,I_3)+(\\beta/2)\\,\\mathcal{N}_3(0,\\Sigma_1)$",
        "in $\\mathbb{R}^3$, where $\\Sigma_1$ has unit marginals and pairwise correlation 0.85,",
        "with $M=1000$, and $B=5000$."
      ),
      label = "tab:power-mvnormal-equal-mean-distinct-sigma"
    )
  )
}

default_paper_output_dir <- function(block) {
  file.path(
    paper_repo_root,
    "simulation_results",
    switch(
      block,
      vmf_missing = "power_mixtures_vmf_beta0125_B5000",
      hvmf = "power_mixtures_hvmf_B5000",
      mvnormal = "power_mixtures_mvnormal_B5000",
      tables = "power_mixtures_paper_tables",
      stop(sprintf("Unsupported block '%s'.", block))
    )
  )
}

parse_paper_task <- function(value) {
  value <- tolower(trimws(as.character(value %||% "simulate")))
  if (!value %in% c("simulate", "tables", "all")) {
    stop("`task` must be one of 'simulate', 'tables', or 'all'.")
  }
  value
}

parse_paper_block <- function(value) {
  value <- tolower(trimws(as.character(value %||% "all_missing")))
  if (!value %in% c("all_missing", "vmf_missing", "hvmf", "mvnormal")) {
    stop("`block` must be one of 'all_missing', 'vmf_missing', 'hvmf', or 'mvnormal'.")
  }
  value
}

paper_block_scenarios <- function(block) {
  switch(
    block,
    vmf_missing = c("vmf_s2_antipodal"),
    hvmf = c("hvmf_h2_case1", "hvmf_h2_case2"),
    mvnormal = c(
      "mvnormal_d3_equal_sigma_opposite_mean",
      "mvnormal_d3_equal_mean_distinct_sigma",
      "mvnormal_d2_moderate_location_correlation",
      "mvnormal_d2_small_location_rho050",
      "mvnormal_d2_small_location_rho075",
      "mvnormal_d2_small_location_rho095"
    ),
    stop(sprintf("Unsupported block '%s'.", block))
  )
}

paper_block_beta_values <- function(block, beta_values_default) {
  if (identical(block, "vmf_missing")) {
    return(0.125)
  }
  as.numeric(beta_values_default)
}

paper_halfway_weight <- function(beta) {
  beta <- as.numeric(beta)
  if (length(beta) != 1L || !is.finite(beta) || beta < 0 || beta > 1) {
    stop("`beta` must be a finite scalar in [0, 1].")
  }
  beta / 2
}

generate_hvmf_halfway_sample <- function(n, beta, mu0, mu1, kappa) {
  alt_weight <- paper_halfway_weight(beta)
  choose_alt <- stats::runif(n) < alt_weight
  sample <- matrix(0, nrow = n, ncol = length(mu0))
  n_alt <- sum(choose_alt)
  n_null <- n - n_alt

  if (n_null > 0L) {
    sample[!choose_alt, ] <- rhvmf_h2_polar(n = n_null, mu = mu0, kappa = kappa)
  }
  if (n_alt > 0L) {
    sample[choose_alt, ] <- rhvmf_h2_polar(n = n_alt, mu = mu1, kappa = kappa)
  }

  sample
}

generate_mvnormal_halfway_sample <- function(n, beta, mu0, Sigma0, mu1, Sigma1) {
  alt_weight <- paper_halfway_weight(beta)
  choose_alt <- stats::runif(n) < alt_weight
  sample <- matrix(0, nrow = n, ncol = length(mu0))
  n_alt <- sum(choose_alt)
  n_null <- n - n_alt

  if (n_null > 0L) {
    sample[!choose_alt, ] <- mvtnorm::rmvnorm(n = n_null, mean = mu0, sigma = Sigma0)
  }
  if (n_alt > 0L) {
    sample[choose_alt, ] <- mvtnorm::rmvnorm(n = n_alt, mean = mu1, sigma = Sigma1)
  }

  sample
}

make_hvmf_power_ks_grid <- function(mu0,
                                    mu1,
                                    kappa,
                                    n_ref = 300L,
                                    n_omega = 10L,
                                    n_t = 10L,
                                    seed = 20260705L) {
  make_sample_unique_distance_ks_grid()
}

make_mvnormal_power_ks_grid <- function(theta_list,
                                        omega_scale = 1.5,
                                        n_t = 10L) {
  theta_list <- lapply(theta_list, normalize_mvnormal_theta)
  omega_rows <- list()
  idx <- 1L

  for (theta in theta_list) {
    eig <- eigen(theta$Sigma, symmetric = TRUE)
    omega_rows[[idx]] <- theta$mu
    idx <- idx + 1L

    for (j in seq_len(theta$ambient_dim)) {
      step <- sqrt(max(eig$values[[j]], 0))
      if (step <= 0) {
        next
      }
      direction <- eig$vectors[, j]
      omega_rows[[idx]] <- theta$mu + omega_scale * step * direction
      idx <- idx + 1L
      omega_rows[[idx]] <- theta$mu - omega_scale * step * direction
      idx <- idx + 1L
    }
  }

  omega_grid <- do.call(rbind, omega_rows)
  omega_grid <- as.matrix(unique.data.frame(as.data.frame(round(omega_grid, 12L))))

  t_max <- max(vapply(theta_list, function(theta) {
    shift_matrix <- matrix(
      rep(theta$mu, each = nrow(omega_grid)),
      nrow = nrow(omega_grid),
      ncol = theta$ambient_dim
    ) - omega_grid
    shift_norm_max <- max(sqrt(rowSums(shift_matrix^2)))
    lambda_max <- max(theta$eigenvalues_full)
    shift_norm_max + 5 * sqrt(max(lambda_max, 0) * theta$ambient_dim)
  }, numeric(1)))

  if (!is.finite(t_max) || t_max <= 0) {
    t_max <- 1
  }

  list(
    omega_grid = omega_grid,
    t_grid = seq(0, t_max, length.out = as.integer(n_t))
  )
}

make_paper_design_grid <- function(scenarios, n_values, beta_values) {
  catalog <- paper_scenario_catalog()
  rows <- list()
  idx <- 1L

  for (scenario in scenarios) {
    if (is.null(catalog[[scenario]])) {
      stop(sprintf("Unknown scenario '%s'.", scenario))
    }
    for (n in n_values) {
      for (beta in beta_values) {
        rows[[idx]] <- data.frame(
          scenario = scenario,
          alternative = as.character(catalog[[scenario]]$alternative),
          n = as.integer(n),
          beta = as.numeric(beta),
          gamma_deg = NA_real_,
          dirichlet_a = NA_real_,
          dirichlet_alpha1 = NA_real_,
          dirichlet_alpha2 = NA_real_,
          dirichlet_alpha3 = NA_real_,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  design <- do.call(rbind, rows)
  design$design_id <- seq_len(nrow(design))
  design
}

prepare_paper_ks_grids <- function(scenarios, base_seed = 20260705L) {
  grids <- vector("list", length(scenarios))
  names(grids) <- scenarios
  for (i in seq_along(scenarios)) {
    grids[[i]] <- paper_sample_ks_grid()
  }
  grids
}

paper_sample_ks_grid <- function() {
  make_sample_unique_distance_ks_grid()
}

run_single_paper_job <- function(job_row,
                                 design_row,
                                 B,
                                 alpha_nominal,
                                 statistics,
                                 bootstrap_method,
                                 bootstrap_n_cores,
                                 base_seed,
                                 ks_grid,
                                 derivative_mc_size,
                                 fast_multiplier_cvm_block_size,
                                 logistic_gaussian_quadform_method) {
  catalog <- paper_scenario_catalog()
  scenario_spec <- catalog[[as.character(design_row$scenario)]]
  if (is.null(scenario_spec)) {
    stop(sprintf("Unknown scenario '%s'.", design_row$scenario))
  }

  rep_id <- as.integer(job_row$rep)
  design_id <- as.integer(design_row$design_id)
  data_seed <- seed_from_job(base_seed, design_id = design_id, rep = rep_id, stream = 0L)
  bootstrap_seed <- seed_from_job(base_seed, design_id = design_id, rep = rep_id, stream = 1L)
  derivative_seed <- seed_from_job(base_seed, design_id = design_id, rep = rep_id, stream = 2L)

  start_time <- proc.time()[["elapsed"]]

  out <- data.frame(
    scenario = design_row$scenario,
    alternative = design_row$alternative,
    n = as.integer(design_row$n),
    beta = as.numeric(design_row$beta),
    gamma_deg = as.numeric(design_row$gamma_deg %||% NA_real_),
    dirichlet_a = as.numeric(design_row$dirichlet_a %||% NA_real_),
    dirichlet_alpha1 = as.numeric(design_row$dirichlet_alpha1 %||% NA_real_),
    dirichlet_alpha2 = as.numeric(design_row$dirichlet_alpha2 %||% NA_real_),
    dirichlet_alpha3 = as.numeric(design_row$dirichlet_alpha3 %||% NA_real_),
    rep = rep_id,
    ks_stat = NA_real_,
    cvm_stat = NA_real_,
    ks_pvalue = NA_real_,
    cvm_pvalue = NA_real_,
    status = "ok",
    error_message = NA_character_,
    bootstrap_method_requested = bootstrap_method,
    bootstrap_method_effective = NA_character_,
    fallback_to_reestimated = NA,
    seed_data = data_seed,
    seed_bootstrap = bootstrap_seed,
    seed_derivative = derivative_seed,
    elapsed_seconds = NA_real_,
    stringsAsFactors = FALSE
  )

  result <- tryCatch(
    {
      set.seed(data_seed)

      if (identical(scenario_spec$model, "vmf")) {
        x <- generate_vmf_antipodal_sample(
          n = as.integer(design_row$n),
          beta = as.numeric(design_row$beta),
          mu = scenario_spec$mu0,
          kappa = scenario_spec$kappa
        )

        multiplier_bootstrap_vmf(
          data = x,
          null = list(type = "composite"),
          statistics = statistics,
          ks_grid = ks_grid,
          B = as.integer(B),
          alpha = alpha_nominal,
          n_cores = as.integer(bootstrap_n_cores),
          seed = bootstrap_seed,
          bootstrap_method = bootstrap_method,
          keep = list(
            observed_process = FALSE,
            bootstrap_statistics = FALSE,
            bootstrap_thetas = FALSE
          ),
          control = list(
            derivative_method = "score_mc",
            derivative_mc_size = as.integer(derivative_mc_size),
            derivative_mc_seed = as.integer(derivative_seed),
            fast_multiplier_cvm_block_size = as.integer(fast_multiplier_cvm_block_size),
            vmf_profile_method = "tabulated",
            vmf_profile_n_u = 4097L
          ),
          distance_type = "geodesic",
          unknown_param = "xi"
        )
      } else if (identical(scenario_spec$model, "hvmf")) {
        x <- generate_hvmf_halfway_sample(
          n = as.integer(design_row$n),
          beta = as.numeric(design_row$beta),
          mu0 = scenario_spec$mu0,
          mu1 = scenario_spec$mu1,
          kappa = scenario_spec$kappa
        )

        multiplier_bootstrap_hvmf(
          data = x,
          null = list(type = "composite"),
          statistics = statistics,
          ks_grid = ks_grid,
          B = as.integer(B),
          alpha = alpha_nominal,
          n_cores = as.integer(bootstrap_n_cores),
          seed = bootstrap_seed,
          bootstrap_method = bootstrap_method,
          keep = list(
            observed_process = FALSE,
            bootstrap_statistics = FALSE,
            bootstrap_thetas = FALSE
          ),
          control = list(
            derivative_method = "score_mc",
            derivative_mc_size = as.integer(derivative_mc_size),
            derivative_mc_seed = as.integer(derivative_seed),
            fast_multiplier_cvm_block_size = as.integer(fast_multiplier_cvm_block_size),
            hvmf_profile_method = "tabulated",
            hvmf_profile_n_y = 4097L
          ),
          unknown_param = "both"
        )
      } else if (identical(scenario_spec$model, "mvnormal")) {
        x <- generate_mvnormal_halfway_sample(
          n = as.integer(design_row$n),
          beta = as.numeric(design_row$beta),
          mu0 = scenario_spec$mu0,
          Sigma0 = scenario_spec$Sigma0,
          mu1 = scenario_spec$mu1,
          Sigma1 = scenario_spec$Sigma1
        )

        multiplier_bootstrap_mvnormal(
          data = x,
          null = list(type = "composite"),
          statistics = statistics,
          ks_grid = ks_grid,
          B = as.integer(B),
          alpha = alpha_nominal,
          n_cores = as.integer(bootstrap_n_cores),
          seed = bootstrap_seed,
          bootstrap_method = bootstrap_method,
          keep = list(
            observed_process = FALSE,
            bootstrap_statistics = FALSE,
            bootstrap_thetas = FALSE
          ),
          control = list(
            derivative_method = "score_mc",
            derivative_mc_size = as.integer(derivative_mc_size),
            derivative_mc_seed = as.integer(derivative_seed),
            fast_multiplier_cvm_block_size = as.integer(fast_multiplier_cvm_block_size),
            mvnormal_quadform_method = logistic_gaussian_quadform_method
          ),
          unknown_param = "both"
        )
      } else {
        stop(sprintf("Unsupported model '%s'.", scenario_spec$model))
      }
    },
    error = function(e) e
  )

  out$elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  if (inherits(result, "error")) {
    out$status <- "error"
    out$error_message <- conditionMessage(result)
    out$bootstrap_method_effective <- NA_character_
    out$fallback_to_reestimated <- NA
    return(out)
  }

  out$ks_stat <- as.numeric(result$observed$ks$statistic %||% NA_real_)
  out$cvm_stat <- as.numeric(result$observed$cvm$statistic %||% NA_real_)
  out$ks_pvalue <- as.numeric(result$inference$ks$p_value %||% NA_real_)
  out$cvm_pvalue <- as.numeric(result$inference$cvm$p_value %||% NA_real_)
  out$bootstrap_method_effective <- as.character(result$diagnostics$effective_bootstrap_method %||% NA_character_)
  out$fallback_to_reestimated <- isTRUE(result$diagnostics$fallback_to_reestimated)
  out
}

run_power_mixtures_paper_block <- function(block,
                                           output_dir = default_paper_output_dir(block),
                                           scenarios = NULL,
                                           M = 1000L,
                                           B = 5000L,
                                           n_values = c(50L, 100L, 200L),
                                           beta_values = c(0, 0.25, 0.5, 1),
                                           statistics = c("ks", "cvm"),
                                           alpha_nominal = 0.05,
                                           alpha_grid = seq(0, 1, by = 0.01),
                                           bootstrap_method = "fast_multiplier",
                                           bootstrap_n_cores = 1L,
                                           n_cores_outer = 1L,
                                           base_seed = 20260705L,
                                           derivative_mc_size = 1000L,
                                           fast_multiplier_cvm_block_size = 50L,
                                           logistic_gaussian_quadform_method = "auto",
                                           show_progress = TRUE) {
  scenarios_default <- paper_block_scenarios(block)
  scenarios <- if (is.null(scenarios)) {
    scenarios_default
  } else {
    as.character(scenarios)
  }
  if (!all(scenarios %in% scenarios_default)) {
    stop(sprintf(
      "For block '%s', `scenarios` must be chosen from: %s.",
      block,
      paste(scenarios_default, collapse = ", ")
    ))
  }
  beta_values_block <- paper_block_beta_values(block, beta_values)
  n_values <- sort(unique(as.integer(n_values)))
  ks_grid <- paper_sample_ks_grid()

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(show_progress)) {
    cat(sprintf(
      "Running %d Monte Carlo jobs with up to %d outer core(s) for paper block '%s'.\n",
      length(scenarios) * length(n_values) * length(beta_values_block) * as.integer(M),
      max(1L, as.integer(n_cores_outer)),
      block
    ))
  }

  raw_rows <- list()
  raw_row_offset <- 0L
  n_block_seconds <- numeric(length(n_values))
  names(n_block_seconds) <- as.character(n_values)

  write_partial_outputs <- function(raw_rows_current,
                                    n_values_completed,
                                    final = FALSE) {
    if (length(raw_rows_current) == 0L) {
      return(invisible(NULL))
    }

    raw_results_current <- do.call(rbind, raw_rows_current)
    raw_results_current <- raw_results_current[order(
      raw_results_current$scenario,
      raw_results_current$alternative,
      raw_results_current$n,
      raw_results_current$beta,
      raw_results_current$rep
    ), , drop = FALSE]

    summary_df_current <- summarize_power_results(raw_results_current, alpha = alpha_nominal)
    alpha_curve_df_current <- summarize_alpha_curves(raw_results_current, alpha_grid = alpha_grid)

    utils::write.csv(raw_results_current, file.path(output_dir, "power_replications_long.csv"), row.names = FALSE)
    utils::write.csv(summary_df_current, file.path(output_dir, "power_summary_005.csv"), row.names = FALSE)
    utils::write.csv(alpha_curve_df_current, file.path(output_dir, "alpha_curve_summary.csv"), row.names = FALSE)

    if (nrow(summary_df_current) > 0L && length(unique(summary_df_current$beta)) > 1L) {
      save_all_plots(summary_df_current, alpha_curve_df_current, output_dir = output_dir)
    }
    write_metadata_files(
      output_dir = output_dir,
      M = M,
      B = B,
      n_values = n_values_completed,
      beta_values = beta_values_block,
      statistics = statistics,
      bootstrap_method = bootstrap_method,
      bootstrap_n_cores = bootstrap_n_cores,
      n_cores_outer = n_cores_outer,
      base_seed = base_seed,
      derivative_mc_size = derivative_mc_size,
      fast_multiplier_cvm_block_size = fast_multiplier_cvm_block_size,
      logistic_gaussian_quadform_method = logistic_gaussian_quadform_method,
      alpha_grid = alpha_grid,
      block = block
    )

    saveRDS(
      list(
        raw_results = raw_results_current,
        summary = summary_df_current,
        alpha_curves = alpha_curve_df_current,
        n_values_completed = as.integer(n_values_completed),
        final = isTRUE(final)
      ),
      file = file.path(output_dir, "checkpoint_latest.rds")
    )

    invisible(list(
      raw_results = raw_results_current,
      summary = summary_df_current,
      alpha_curves = alpha_curve_df_current
    ))
  }

  for (n_index in seq_along(n_values)) {
    n_value <- n_values[[n_index]]
    if (isTRUE(show_progress)) {
      if (n_index == 1L) {
        cat(sprintf("Starting n = %d.\n", n_value))
      } else {
        previous_n <- n_values[[n_index - 1L]]
        cat(sprintf(
          "Completed n = %d in %.2f seconds. Starting n = %d.\n",
          previous_n,
          n_block_seconds[[as.character(previous_n)]],
          n_value
        ))
      }
    }

    design_n <- make_paper_design_grid(
      scenarios = scenarios,
      n_values = n_value,
      beta_values = as.numeric(beta_values_block)
    )
    jobs_n <- expand.grid(
      design_id = design_n$design_id,
      rep = seq_len(as.integer(M)),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    block_started_at <- proc.time()[["elapsed"]]

    run_job_by_index <- function(i) {
      job_row <- jobs_n[i, , drop = FALSE]
      design_row <- design_n[design_n$design_id == job_row$design_id, , drop = FALSE]
      run_single_paper_job(
        job_row = job_row,
        design_row = design_row,
        B = B,
        alpha_nominal = alpha_nominal,
        statistics = statistics,
        bootstrap_method = bootstrap_method,
        bootstrap_n_cores = bootstrap_n_cores,
        base_seed = base_seed,
        ks_grid = ks_grid,
        derivative_mc_size = derivative_mc_size,
        fast_multiplier_cvm_block_size = fast_multiplier_cvm_block_size,
        logistic_gaussian_quadform_method = logistic_gaussian_quadform_method
      )
    }

    raw_rows_n <- run_jobs_with_progress(
      total_jobs = nrow(jobs_n),
      run_job_by_index = run_job_by_index,
      n_cores_outer = n_cores_outer,
      show_progress = show_progress
    )
    n_block_seconds[[as.character(n_value)]] <- proc.time()[["elapsed"]] - block_started_at
    raw_rows[(raw_row_offset + 1L):(raw_row_offset + length(raw_rows_n))] <- raw_rows_n
    raw_row_offset <- raw_row_offset + length(raw_rows_n)

    write_partial_outputs(
      raw_rows_current = raw_rows[seq_len(raw_row_offset)],
      n_values_completed = n_values[seq_len(n_index)],
      final = (n_index == length(n_values))
    )
  }

  if (isTRUE(show_progress) && length(n_values) > 0L) {
    final_n <- n_values[[length(n_values)]]
    cat(sprintf(
      "Completed n = %d in %.2f seconds.\n",
      final_n,
      n_block_seconds[[as.character(final_n)]]
    ))
  }

  final_outputs <- write_partial_outputs(
    raw_rows_current = raw_rows,
    n_values_completed = n_values,
    final = TRUE
  )

  invisible(list(
    raw_results = final_outputs$raw_results,
    summary = final_outputs$summary,
    alpha_curves = final_outputs$alpha_curves
  ))
}

read_metadata_value <- function(metadata_path, key) {
  if (!file.exists(metadata_path)) {
    stop(sprintf("Metadata file not found: %s", metadata_path))
  }

  lines <- readLines(metadata_path, warn = FALSE)
  prefix <- sprintf("%s:", key)
  matches <- lines[startsWith(lines, prefix)]
  if (length(matches) == 0L) {
    stop(sprintf("Could not find key '%s' in metadata file: %s", key, metadata_path))
  }

  trimws(sub(prefix, "", matches[[1L]], fixed = TRUE))
}

read_checked_summary <- function(output_dir, expected_M = 1000L, expected_B = 5000L) {
  summary_path <- file.path(output_dir, "power_summary_005.csv")
  metadata_path <- file.path(output_dir, "metadata.txt")

  if (!file.exists(summary_path)) {
    stop(sprintf("Summary CSV not found: %s", summary_path))
  }

  observed_M <- as.integer(read_metadata_value(metadata_path, "M"))
  observed_B <- as.integer(read_metadata_value(metadata_path, "B"))

  if (observed_M != as.integer(expected_M)) {
    stop(sprintf("Expected M = %d in %s, found %d.", expected_M, output_dir, observed_M))
  }
  if (observed_B != as.integer(expected_B)) {
    stop(sprintf("Expected B = %d in %s, found %d.", expected_B, output_dir, observed_B))
  }

  list(
    summary = utils::read.csv(summary_path, stringsAsFactors = FALSE),
    metadata = list(M = observed_M, B = observed_B)
  )
}

read_checked_alpha_curves <- function(output_dir, expected_M = 1000L, expected_B = 5000L) {
  alpha_path <- file.path(output_dir, "alpha_curve_summary.csv")
  metadata_path <- file.path(output_dir, "metadata.txt")

  if (!file.exists(alpha_path)) {
    stop(sprintf("Alpha-curve CSV not found: %s", alpha_path))
  }

  observed_M <- as.integer(read_metadata_value(metadata_path, "M"))
  observed_B <- as.integer(read_metadata_value(metadata_path, "B"))

  if (observed_M != as.integer(expected_M)) {
    stop(sprintf("Expected M = %d in %s, found %d.", expected_M, output_dir, observed_M))
  }
  if (observed_B != as.integer(expected_B)) {
    stop(sprintf("Expected B = %d in %s, found %d.", expected_B, output_dir, observed_B))
  }

  list(
    alpha_curves = utils::read.csv(alpha_path, stringsAsFactors = FALSE),
    metadata = list(M = observed_M, B = observed_B)
  )
}

map_vmf_old_beta_to_paper <- function(beta) {
  out <- rep(NA_real_, length(beta))
  out[abs(beta - 0) < 1e-12] <- 0
  out[abs(beta - 0.125) < 1e-12] <- 0.25
  out[abs(beta - 0.25) < 1e-12] <- 0.5
  out[abs(beta - 0.5) < 1e-12] <- 1

  if (any(is.na(out))) {
    stop("Encountered a vMF beta value that does not match the paper relabeling.")
  }

  out
}

assert_complete_beta_grid <- function(df,
                                      scenario,
                                      n_values = c(50L, 100L, 200L),
                                      beta_values = c(0, 0.25, 0.5, 1)) {
  subset_df <- df[df$scenario == scenario, , drop = FALSE]
  expected <- expand.grid(
    n = as.integer(n_values),
    beta = as.numeric(beta_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  observed <- unique(subset_df[, c("n", "beta")])
  key_expected <- paste(expected$n, format(expected$beta, trim = TRUE), sep = "|")
  key_observed <- paste(observed$n, format(observed$beta, trim = TRUE), sep = "|")
  missing_keys <- setdiff(key_expected, key_observed)

  if (length(missing_keys) > 0L) {
    stop(sprintf(
      "Scenario '%s' is missing the following (n, beta) combinations: %s",
      scenario,
      paste(missing_keys, collapse = ", ")
    ))
  }

  invisible(TRUE)
}

build_vmf_paper_summary <- function() {
  beta0_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_vmf_beta0_B5000")
  betarest_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_vmf_beta_rest_B5000")
  beta0125_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_vmf_beta0125_B5000")

  beta0 <- read_checked_summary(beta0_dir)$summary
  betarest <- read_checked_summary(betarest_dir)$summary
  beta0125 <- read_checked_summary(beta0125_dir)$summary

  keep0 <- beta0[beta0$scenario == "vmf_s2_antipodal" & abs(beta0$beta - 0) < 1e-12, , drop = FALSE]
  keeprest <- betarest[
    betarest$scenario == "vmf_s2_antipodal" &
      (abs(betarest$beta - 0.25) < 1e-12 | abs(betarest$beta - 0.5) < 1e-12),
    ,
    drop = FALSE
  ]
  keep0125 <- beta0125[
    beta0125$scenario == "vmf_s2_antipodal" & abs(beta0125$beta - 0.125) < 1e-12,
    ,
    drop = FALSE
  ]

  combined <- rbind(keep0, keep0125, keeprest)
  combined$beta <- map_vmf_old_beta_to_paper(combined$beta)
  combined$scenario <- "vmf_s2_antipodal_paper"
  combined$alternative <- "antipodal_halfway_mixture_paper"
  rownames(combined) <- NULL

  assert_complete_beta_grid(combined, "vmf_s2_antipodal_paper")
  combined
}

build_lg_paper_summary <- function() {
  lg_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_logistic_gaussian_B5000")
  summary_df <- read_checked_summary(lg_dir)$summary
  keep <- summary_df[
    summary_df$scenario == "logistic_gaussian_simplex_d3" &
      abs(summary_df$dirichlet_alpha1 - 1) < 1e-12 &
      abs(summary_df$dirichlet_alpha2 - 1) < 1e-12 &
      abs(summary_df$dirichlet_alpha3 - 8) < 1e-12,
    ,
    drop = FALSE
  ]
  keep$scenario <- "logistic_gaussian_simplex_d3_dirichlet_1_1_8_paper"
  keep$alternative <- "dirichlet_mixture_1_1_8_paper"
  rownames(keep) <- NULL

  assert_complete_beta_grid(keep, "logistic_gaussian_simplex_d3_dirichlet_1_1_8_paper")
  keep
}

build_vmf_paper_alpha_curves <- function() {
  beta0_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_vmf_beta0_B5000")
  betarest_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_vmf_beta_rest_B5000")
  beta0125_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_vmf_beta0125_B5000")

  beta0 <- read_checked_alpha_curves(beta0_dir)$alpha_curves
  betarest <- read_checked_alpha_curves(betarest_dir)$alpha_curves
  beta0125 <- read_checked_alpha_curves(beta0125_dir)$alpha_curves

  keep0 <- beta0[beta0$scenario == "vmf_s2_antipodal" & abs(beta0$beta - 0) < 1e-12, , drop = FALSE]
  keeprest <- betarest[
    betarest$scenario == "vmf_s2_antipodal" &
      (abs(betarest$beta - 0.25) < 1e-12 | abs(betarest$beta - 0.5) < 1e-12),
    ,
    drop = FALSE
  ]
  keep0125 <- beta0125[
    beta0125$scenario == "vmf_s2_antipodal" & abs(beta0125$beta - 0.125) < 1e-12,
    ,
    drop = FALSE
  ]

  combined <- rbind(keep0, keep0125, keeprest)
  combined$beta <- map_vmf_old_beta_to_paper(combined$beta)
  combined$scenario <- "vmf_s2_antipodal_paper"
  combined$alternative <- "antipodal_halfway_mixture_paper"
  rownames(combined) <- NULL
  combined
}

build_lg_paper_alpha_curves <- function() {
  lg_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_logistic_gaussian_B5000")
  alpha_df <- read_checked_alpha_curves(lg_dir)$alpha_curves
  keep <- alpha_df[
    alpha_df$scenario == "logistic_gaussian_simplex_d3" &
      abs(alpha_df$dirichlet_alpha1 - 1) < 1e-12 &
      abs(alpha_df$dirichlet_alpha2 - 1) < 1e-12 &
      abs(alpha_df$dirichlet_alpha3 - 8) < 1e-12,
    ,
    drop = FALSE
  ]
  keep$scenario <- "logistic_gaussian_simplex_d3_dirichlet_1_1_8_paper"
  keep$alternative <- "dirichlet_mixture_1_1_8_paper"
  rownames(keep) <- NULL
  keep
}

build_hvmf_paper_summary <- function() {
  hvmf_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_hvmf_B5000")
  summary_df <- read_checked_summary(hvmf_dir)$summary
  for (scenario in c("hvmf_h2_case1", "hvmf_h2_case2")) {
    assert_complete_beta_grid(summary_df, scenario)
  }
  summary_df
}

build_hvmf_paper_alpha_curves <- function() {
  hvmf_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_hvmf_B5000")
  read_checked_alpha_curves(hvmf_dir)$alpha_curves
}

build_mvnormal_paper_summary <- function() {
  mvnormal_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_mvnormal_B5000")
  summary_df <- read_checked_summary(mvnormal_dir)$summary
  for (scenario in c(
    "mvnormal_d3_equal_sigma_opposite_mean",
    "mvnormal_d3_equal_mean_distinct_sigma"
  )) {
    assert_complete_beta_grid(summary_df, scenario)
  }
  summary_df
}

build_mvnormal_paper_alpha_curves <- function() {
  mvnormal_dir <- file.path(paper_repo_root, "simulation_results", "power_mixtures_mvnormal_B5000")
  read_checked_alpha_curves(mvnormal_dir)$alpha_curves
}

paper_calibration_band <- function(alpha = 0.05, M = 1000L, z = 1.96) {
  se <- sqrt(alpha * (1 - alpha) / as.integer(M))
  c(alpha - z * se, alpha + z * se)
}

value_in_band <- function(x, band, tol = 1e-12) {
  is.finite(x) && x >= band[[1L]] - tol && x <= band[[2L]] + tol
}

extract_single_power_value <- function(df, n_value, beta_value, column) {
  hit <- df[
    df$n == as.integer(n_value) &
      abs(df$beta - as.numeric(beta_value)) < 1e-12,
    ,
    drop = FALSE
  ]
  if (nrow(hit) != 1L) {
    stop(sprintf(
      "Expected exactly one row for n = %s and beta = %s, found %d.",
      n_value,
      format(beta_value, trim = TRUE),
      nrow(hit)
    ))
  }
  as.numeric(hit[[column]][[1L]])
}

build_power_table_data <- function(summary_df,
                                   scenario,
                                   beta_values = c(0, 0.25, 0.5, 1),
                                   calibration_band = paper_calibration_band()) {
  subset_df <- summary_df[summary_df$scenario == scenario, , drop = FALSE]
  n_values <- sort(unique(subset_df$n))
  rows <- list()
  idx <- 1L

  for (n_value in n_values) {
    ks_beta0 <- extract_single_power_value(subset_df, n_value, 0, "power_ks_005")
    cvm_beta0 <- extract_single_power_value(subset_df, n_value, 0, "power_cvm_005")

    rows[[idx]] <- data.frame(
      n = as.integer(n_value),
      statistic = "KS",
      beta_0 = ks_beta0,
      beta_0_25 = extract_single_power_value(subset_df, n_value, 0.25, "power_ks_005"),
      beta_0_5 = extract_single_power_value(subset_df, n_value, 0.5, "power_ks_005"),
      beta_1 = extract_single_power_value(subset_df, n_value, 1, "power_ks_005"),
      h0_in_band = value_in_band(ks_beta0, calibration_band),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L

    rows[[idx]] <- data.frame(
      n = as.integer(n_value),
      statistic = "CvM",
      beta_0 = cvm_beta0,
      beta_0_25 = extract_single_power_value(subset_df, n_value, 0.25, "power_cvm_005"),
      beta_0_5 = extract_single_power_value(subset_df, n_value, 0.5, "power_cvm_005"),
      beta_1 = extract_single_power_value(subset_df, n_value, 1, "power_cvm_005"),
      h0_in_band = value_in_band(cvm_beta0, calibration_band),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  do.call(rbind, rows)
}

format_prob_entry <- function(value, bold = FALSE, digits = 3L) {
  if (!is.finite(value)) {
    return("NA")
  }

  entry <- formatC(value, digits = digits, format = "f")
  if (isTRUE(bold)) {
    sprintf("\\textbf{%s}", entry)
  } else {
    entry
  }
}

write_latex_power_table <- function(table_df,
                                    file_path,
                                    caption,
                                    label,
                                    calibration_band,
                                    M = 1000L) {
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\begin{tabular}{llcccc}",
    "\\toprule",
    "$n$ & Stat. & $\\beta=0$ & $\\beta=0.25$ & $\\beta=0.5$ & $\\beta=1$ \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(table_df))) {
    row <- table_df[i, , drop = FALSE]
    lines <- c(
      lines,
      sprintf(
        "%d & %s & %s & %s & %s & %s \\\\",
        as.integer(row$n),
        as.character(row$statistic),
        format_prob_entry(row$beta_0, bold = isTRUE(row$h0_in_band)),
        format_prob_entry(row$beta_0_25),
        format_prob_entry(row$beta_0_5),
        format_prob_entry(row$beta_1)
      )
    )
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    sprintf(
      "\\caption{%s Boldface at $\\beta=0$ indicates that the empirical size lies in the Monte Carlo 95\\%% interval $[%.4f, %.4f]$ under $H_0$, with $M=%d$.}",
      caption,
      calibration_band[[1L]],
      calibration_band[[2L]],
      as.integer(M)
    ),
    sprintf("\\label{%s}", label),
    "\\end{table}"
  )

  writeLines(lines, con = file_path)
}

compile_paper_power_tables <- function(output_dir = default_paper_output_dir("tables")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  calibration_band <- paper_calibration_band()

  combined_summary <- rbind(
    build_vmf_paper_summary(),
    build_lg_paper_summary(),
    build_hvmf_paper_summary(),
    build_mvnormal_paper_summary()
  )
  combined_alpha_curves <- rbind(
    build_vmf_paper_alpha_curves(),
    build_lg_paper_alpha_curves(),
    build_hvmf_paper_alpha_curves(),
    build_mvnormal_paper_alpha_curves()
  )
  utils::write.csv(
    combined_summary,
    file.path(output_dir, "combined_power_summary_005_paper.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    combined_alpha_curves,
    file.path(output_dir, "combined_alpha_curve_summary_paper.csv"),
    row.names = FALSE
  )
  save_all_plots(
    summary_df = combined_summary,
    alpha_curve_df = combined_alpha_curves,
    output_dir = output_dir
  )

  table_specs <- paper_table_catalog()
  table_scenarios <- names(table_specs)
  manifest_lines <- c(
    sprintf("created_at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("paper_repo_root: %s", paper_repo_root),
    sprintf("calibration_band_lower: %.8f", calibration_band[[1L]]),
    sprintf("calibration_band_upper: %.8f", calibration_band[[2L]])
  )

  for (scenario in table_scenarios) {
    table_df <- build_power_table_data(
      summary_df = combined_summary,
      scenario = scenario,
      calibration_band = calibration_band
    )
    csv_path <- file.path(output_dir, sprintf("%s_power_table.csv", safe_slug(scenario)))
    tex_path <- file.path(output_dir, sprintf("%s_power_table.tex", safe_slug(scenario)))

    utils::write.csv(table_df, csv_path, row.names = FALSE)
    write_latex_power_table(
      table_df = table_df,
      file_path = tex_path,
      caption = table_specs[[scenario]]$caption,
      label = table_specs[[scenario]]$label,
      calibration_band = calibration_band,
      M = 1000L
    )

    manifest_lines <- c(
      manifest_lines,
      sprintf("scenario: %s", scenario),
      sprintf("  csv: %s", csv_path),
      sprintf("  tex: %s", tex_path)
    )
  }

  writeLines(manifest_lines, con = file.path(output_dir, "manifest.txt"))

  invisible(list(
    combined_summary = combined_summary,
    output_dir = output_dir
  ))
}

run_selected_paper_blocks <- function(block,
                                      scenarios = NULL,
                                      M,
                                      B,
                                      n_values,
                                      beta_values,
                                      statistics,
                                      alpha_nominal,
                                      alpha_grid,
                                      bootstrap_method,
                                      bootstrap_n_cores,
                                      n_cores_outer,
                                      base_seed,
                                      derivative_mc_size,
                                      fast_multiplier_cvm_block_size,
                                      logistic_gaussian_quadform_method,
                                      show_progress) {
  blocks <- if (identical(block, "all_missing")) {
    c("vmf_missing", "hvmf", "mvnormal")
  } else {
    c(block)
  }

  results <- vector("list", length(blocks))
  names(results) <- blocks

  for (i in seq_along(blocks)) {
    block_i <- blocks[[i]]
    results[[block_i]] <- run_power_mixtures_paper_block(
      block = block_i,
      output_dir = default_paper_output_dir(block_i),
      scenarios = if (identical(block_i, block)) scenarios else NULL,
      M = M,
      B = B,
      n_values = n_values,
      beta_values = beta_values,
      statistics = statistics,
      alpha_nominal = alpha_nominal,
      alpha_grid = alpha_grid,
      bootstrap_method = bootstrap_method,
      bootstrap_n_cores = bootstrap_n_cores,
      n_cores_outer = n_cores_outer,
      base_seed = base_seed + 1000L * (i - 1L),
      derivative_mc_size = derivative_mc_size,
      fast_multiplier_cvm_block_size = fast_multiplier_cvm_block_size,
      logistic_gaussian_quadform_method = logistic_gaussian_quadform_method,
      show_progress = show_progress
    )
  }

  invisible(results)
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_power_mixtures(commandArgs(trailingOnly = TRUE))
  task <- parse_paper_task(args$task %||% "simulate")
  block <- parse_paper_block(args$block %||% "all_missing")

  if (task %in% c("simulate", "all")) {
    run_selected_paper_blocks(
      block = block,
      scenarios = if (!is.null(args$scenarios)) {
        parse_character_csv(args$scenarios, character(0))
      } else {
        NULL
      },
      M = as.integer(args$M %||% 1000L),
      B = as.integer(args$B %||% 5000L),
      n_values = parse_integer_csv(args$n_values, c(50L, 100L, 200L)),
      beta_values = parse_numeric_csv(args$beta_values, c(0, 0.25, 0.5, 1)),
      statistics = parse_character_csv(args$statistics, c("ks", "cvm")),
      alpha_nominal = as.numeric(args$alpha_nominal %||% 0.05),
      alpha_grid = parse_numeric_csv(args$alpha_grid, seq(0, 1, by = 0.01)),
      bootstrap_method = as.character(args$bootstrap_method %||% "fast_multiplier"),
      bootstrap_n_cores = as.integer(args$bootstrap_n_cores %||% 1L),
      n_cores_outer = as.integer(args$n_cores_outer %||% 1L),
      base_seed = as.integer(args$seed %||% 20260705L),
      derivative_mc_size = as.integer(args$derivative_mc_size %||% 1000L),
      fast_multiplier_cvm_block_size = as.integer(args$fast_multiplier_cvm_block_size %||% 50L),
      logistic_gaussian_quadform_method = as.character(args$logistic_gaussian_quadform_method %||% "auto"),
      show_progress = parse_logical_flag(args$show_progress, default = TRUE)
    )
  }

  if (task %in% c("tables", "all")) {
    compile_paper_power_tables(
      output_dir = as.character(args$tables_output_dir %||% default_paper_output_dir("tables"))
    )
  }
}
