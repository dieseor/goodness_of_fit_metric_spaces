resolve_fast_paper_path <- function(...) {
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

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

source(resolve_fast_paper_path("bootstrap", "calibration_study.R"))
source(resolve_fast_paper_path("scripts", "run_fast_multiplier_validation_all_models.R"))
source(resolve_fast_paper_path("real_data", "logistic_gaussian", "utils_logistic_gaussian_screening.R"))
source(resolve_fast_paper_path("real_data", "wind", "run_risoe_125m_screening_ks_cvm.R"))
source(resolve_fast_paper_path("real_data", "sunspots", "run_sunspots_weighted_mixture_rolling_windows_gof.R"))
source(resolve_fast_paper_path("scripts", "run_comets_distance_profile_cardioid.R"))
source(resolve_fast_paper_path("scripts", "run_comets_distance_profile_small_circle_benchmark.R"))
source(resolve_fast_paper_path("scripts", "run_comets_distance_profile_jp_short_benchmark.R"))
source(resolve_fast_paper_path("scripts", "run_comets_rotational_mixtures_short_long.R"))
source(resolve_fast_paper_path("real_data", "comets", "utils_comets_data.R"))

fast_multiplier_common_control <- function(derivative_mc_size = 1000L,
                                           derivative_mc_seed = 20260613L,
                                           cvm_block_size = 50L) {
  list(
    derivative_method = "score_mc",
    derivative_mc_size = as.integer(derivative_mc_size),
    derivative_mc_seed = as.integer(derivative_mc_seed),
    fast_multiplier_cvm_block_size = as.integer(cvm_block_size)
  )
}

prepare_fast_scenario <- function(scenario,
                                  derivative_mc_size = 1000L,
                                  derivative_mc_seed = 20260613L,
                                  cvm_block_size = 50L) {
  if (identical(scenario$null$type, "composite")) {
    scenario$bootstrap_method <- "fast_multiplier"
    scenario$control <- utils::modifyList(
      scenario$control %||% list(),
      fast_multiplier_common_control(
        derivative_mc_size = derivative_mc_size,
        derivative_mc_seed = derivative_mc_seed,
        cvm_block_size = cvm_block_size
      )
    )
  }
  scenario
}

paper_calibration_jobs <- function(derivative_mc_size = 1000L,
                                   derivative_mc_seed = 20260613L) {
  list(
    list(
      id = "normal_composite",
      scenarios = list(prepare_fast_scenario(make_normal_composite_calibration_scenario(), derivative_mc_size, derivative_mc_seed)),
      statistics = c("ks", "cvm"),
      n_values = c(50L, 100L, 200L),
      M_outer = 1000L,
      B = 1000L
    ),
    list(
      id = "vmf_composite",
      scenarios = lapply(list(
        make_vmf_composite_calibration_scenario(0.5),
        make_vmf_composite_calibration_scenario(2.0)
      ), prepare_fast_scenario, derivative_mc_size = derivative_mc_size, derivative_mc_seed = derivative_mc_seed),
      statistics = c("ks", "cvm"),
      n_values = c(50L, 100L, 200L),
      M_outer = 1000L,
      B = 1000L
    ),
    list(
      id = "hvmf_composite",
      scenarios = lapply(default_hvmf_composite_calibration_scenarios(), prepare_fast_scenario, derivative_mc_size = derivative_mc_size, derivative_mc_seed = derivative_mc_seed),
      statistics = c("cvm"),
      n_values = c(50L, 100L, 200L),
      M_outer = 1000L,
      B = 1000L
    ),
    list(
      id = "logistic_gaussian_composite",
      scenarios = lapply(default_logistic_gaussian_composite_calibration_scenarios(), prepare_fast_scenario, derivative_mc_size = derivative_mc_size, derivative_mc_seed = derivative_mc_seed),
      statistics = c("ks", "cvm"),
      n_values = c(50L, 100L, 200L),
      M_outer = 1000L,
      B = 1000L
    ),
    list(
      id = "small_circle_composite",
      scenarios = lapply(list(
        make_small_circle_composite_calibration_scenario(5, 0)
      ), prepare_fast_scenario, derivative_mc_size = derivative_mc_size, derivative_mc_seed = derivative_mc_seed),
      statistics = c("ks", "cvm"),
      n_values = c(50L, 100L, 200L),
      M_outer = 1000L,
      B = 1000L
    ),
    list(
      id = "spherical_cauchy_composite",
      scenarios = lapply(list(
        make_spherical_cauchy_composite_calibration_scenario(0.3)
      ), prepare_fast_scenario, derivative_mc_size = derivative_mc_size, derivative_mc_seed = derivative_mc_seed),
      statistics = c("ks", "cvm"),
      n_values = c(50L, 100L, 200L),
      M_outer = 1000L,
      B = 1000L
    ),
    list(
      id = "beta_mixture2_composite",
      scenarios = lapply(list(
        make_beta_mixture2_composite_calibration_scenario(weight1 = 0.4, alpha1 = 2, beta1 = 8, alpha2 = 8, beta2 = 2)
      ), prepare_fast_scenario, derivative_mc_size = derivative_mc_size, derivative_mc_seed = derivative_mc_seed),
      statistics = c("ks", "cvm"),
      n_values = c(50L, 100L, 200L),
      M_outer = 1000L,
      B = 1000L
    )
  )
}

run_fast_paper_calibration <- function(output_root = file.path("output", "fast_multiplier", "paper_results", "calibration"),
                                       n_cores = 12L,
                                       derivative_mc_size = 1000L,
                                       derivative_mc_seed = 20260613L) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  jobs <- paper_calibration_jobs(
    derivative_mc_size = derivative_mc_size,
    derivative_mc_seed = derivative_mc_seed
  )
  job_paths <- character(length(jobs))
  names(job_paths) <- vapply(jobs, `[[`, character(1), "id")

  for (i in seq_along(jobs)) {
    job <- jobs[[i]]
    job_dir <- file.path(output_root, job$id)
    message(sprintf("[fast paper calibration] %d/%d %s", i, length(jobs), job$id))
    result <- run_bootstrap_calibration_study(
      scenarios = job$scenarios,
      n_values = job$n_values,
      M_outer = as.integer(job$M_outer),
      B = as.integer(job$B),
      statistics = job$statistics,
      output_dir = job_dir,
      n_cores_outer = as.integer(n_cores),
      seed = 20260613L + i
    )
    saveRDS(result, file = file.path(job_dir, "result.rds"))
    job_paths[[i]] <- job_dir
  }

  writeLines(job_paths, con = file.path(output_root, "job_paths.txt"))
  invisible(job_paths)
}

run_fast_comets_jp_short_long <- function(output_root,
                                          B = 1000L,
                                          n_cores = 12L,
                                          bootstrap_method = "fast_multiplier",
                                          derivative_mc_size = 1000L,
                                          derivative_mc_seed = 20260613L) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  comets_data <- load_comets_real_data(finite_normals = "both")
  datasets <- list(
    short_period = as.matrix(comets_data$short$normal),
    long_period = as.matrix(comets_data$long$normal)
  )
  rows <- list()
  idx <- 1L
  for (dataset_name in names(datasets)) {
    for (stat_name in c("ks", "cvm")) {
      result <- run_single_jp_short_comet_model(
        data_matrix = datasets[[dataset_name]],
        statistic = stat_name,
        B = as.integer(B),
        M_value = as.integer(B),
        n_cores = as.integer(n_cores),
        seed = 20270000L + idx,
        bootstrap_method = bootstrap_method,
        distance_type = "geodesic",
        control = utils::modifyList(
          list(
            jp_profile_method = "tabulated",
            jp_profile_n_u = 1025L,
            jp_profile_n_delta = 257L,
            jp_vmf_switch_abs_kappa_psi = 1e-3,
            jp_mle_max_abs_kappa_psi = 6
          ),
          fast_multiplier_common_control(
            derivative_mc_size = derivative_mc_size,
            derivative_mc_seed = derivative_mc_seed + idx
          )
        )
      )
      inf <- result$inference[[stat_name]]
      theta_hat <- result$observed$theta_hat
      rows[[idx]] <- data.frame(
        dataset = dataset_name,
        statistic = stat_name,
        p_value = inf$p_value,
        observed_statistic = inf$observed,
        kappa_hat = theta_hat$kappa,
        psi_hat = theta_hat$psi,
        bootstrap_method = bootstrap_method,
        stringsAsFactors = FALSE
      )
      saveRDS(result, file = file.path(output_root, sprintf("%s_%s_result.rds", dataset_name, stat_name)))
      idx <- idx + 1L
    }
  }
  summary_df <- do.call(rbind, rows)
  utils::write.csv(summary_df, file.path(output_root, "jp_short_long_summary.csv"), row.names = FALSE)
  invisible(summary_df)
}

run_fast_paper_real_data <- function(output_root = file.path("output", "fast_multiplier", "paper_results", "real_data"),
                                     n_cores = 12L,
                                     derivative_mc_size = 1000L,
                                     derivative_mc_seed = 20260613L) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  message("[fast paper real_data] comets cardioid")
  run_comets_distance_profile_cardioid(
    output_root = file.path(output_root, "comets_cardioid"),
    stages = c("oort_cvm", "oort_ks", "short_cvm", "short_ks"),
    cvm_B = 1000L,
    ks_B = 1000L,
    n_cores = as.integer(n_cores),
    bootstrap_method = "fast_multiplier",
    control = utils::modifyList(
      list(cardioid_optim_control = list(maxit = 1000)),
      fast_multiplier_common_control(derivative_mc_size, derivative_mc_seed)
    )
  )

  message("[fast paper real_data] comets small circle")
  for (dataset_name in c("short", "long")) {
    for (stat_name in c("ks", "cvm")) {
      run_comets_distance_profile_small_circle_benchmark(
        output_root = file.path(output_root, sprintf("comets_small_circle_%s_%s", dataset_name, stat_name)),
        dataset = dataset_name,
        B_values = 1000L,
        statistic = stat_name,
        n_cores = as.integer(n_cores),
        bootstrap_method = "fast_multiplier",
        control = utils::modifyList(
          list(
            small_circle_profile_method = "legendre",
            small_circle_L_max = 200L,
            small_circle_quad_n = 400L,
            small_circle_tol = 1e-10,
            small_circle_optim_control = list(maxit = 300L, reltol = 1e-9)
          ),
          fast_multiplier_common_control(derivative_mc_size, derivative_mc_seed)
        )
      )
    }
  }

  message("[fast paper real_data] comets jp")
  run_fast_comets_jp_short_long(
    output_root = file.path(output_root, "comets_jp"),
    B = 1000L,
    n_cores = as.integer(n_cores),
    bootstrap_method = "fast_multiplier",
    derivative_mc_size = derivative_mc_size,
    derivative_mc_seed = derivative_mc_seed
  )

  message("[fast paper real_data] comets beta mixture")
  run_comets_mixtures_short_long(
    output_root = file.path(output_root, "comets_beta_mixture2"),
    datasets = c("short", "long"),
    models = "beta_mixture2",
    B = 1000L,
    statistics = c("ks", "cvm"),
    n_cores = as.integer(n_cores),
    bootstrap_method = "fast_multiplier"
  )

  message("[fast paper real_data] logistic Gaussian datasets")
  run_logistic_gaussian_screening_batch(
    B = 1000L,
    n_cores = as.integer(n_cores),
    bootstrap_method = "fast_multiplier",
    output_dir = file.path(output_root, "logistic_gaussian"),
    control = fast_multiplier_common_control(derivative_mc_size, derivative_mc_seed)
  )

  message("[fast paper real_data] Risoe 125m screening")
  run_risoe_125m_screening_ks_cvm(
    output_dir = file.path(output_root, "risoe_125m"),
    B = 1000L,
    n_cores = as.integer(n_cores),
    bootstrap_method = "fast_multiplier"
  )

  message("[fast paper real_data] sunspots rolling windows")
  run_sunspots_weighted_mixture_rolling_windows_gof(
    output_dir = file.path(output_root, "sunspots_weighted_mixture"),
    statistics = "ks",
    B = 1000L,
    n_cores = as.integer(n_cores),
    bootstrap_method = "fast_multiplier"
  )

  invisible(output_root)
}

parse_fast_paper_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

run_fast_multiplier_paper_results <- function(mode = c("all", "paper_only", "validation", "calibration", "real_data"),
                                              output_root = file.path("output", "fast_multiplier", "paper_results"),
                                              validation_output_root = file.path("output", "fast_multiplier", "validation_all_models_12cores_B1000_M5"),
                                              n_cores = 12L,
                                              derivative_mc_size = 1000L,
                                              derivative_mc_seed = 20260613L,
                                              validation_M_outer = 5L) {
  mode <- match.arg(mode)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  if (mode %in% c("all", "validation")) {
    run_fast_multiplier_validation_all_models(
      output_root = validation_output_root,
      B = 1000L,
      M_outer = as.integer(validation_M_outer),
      derivative_mc_size = as.integer(derivative_mc_size),
      base_seed = 20260613L,
      n_cores = as.integer(n_cores)
    )
  }

  if (mode %in% c("all", "paper_only", "calibration")) {
    run_fast_paper_calibration(
      output_root = file.path(output_root, "calibration"),
      n_cores = as.integer(n_cores),
      derivative_mc_size = as.integer(derivative_mc_size),
      derivative_mc_seed = as.integer(derivative_mc_seed)
    )
  }

  if (mode %in% c("all", "paper_only", "real_data")) {
    run_fast_paper_real_data(
      output_root = file.path(output_root, "real_data"),
      n_cores = as.integer(n_cores),
      derivative_mc_size = as.integer(derivative_mc_size),
      derivative_mc_seed = as.integer(derivative_mc_seed)
    )
  }

  invisible(output_root)
}

if (sys.nframe() == 0L) {
  args <- parse_fast_paper_args(commandArgs(trailingOnly = TRUE))
  run_fast_multiplier_paper_results(
    mode = args$mode %||% "all",
    output_root = args$output_root %||% file.path("output", "fast_multiplier", "paper_results"),
    validation_output_root = args$validation_output_root %||% file.path("output", "fast_multiplier", "validation_all_models_12cores_B1000_M5"),
    n_cores = as.integer(args$n_cores %||% 12L),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 1000L),
    derivative_mc_seed = as.integer(args$derivative_mc_seed %||% 20260613L),
    validation_M_outer = as.integer(args$validation_M_outer %||% 5L)
  )
}
