resolve_fast_problem_diag_path <- function(...) {
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

source(resolve_fast_problem_diag_path("scripts", "run_fast_multiplier_validation_all_models.R"))
source(resolve_fast_problem_diag_path("bootstrap", "calibration_study.R"))

problem_case_order <- function(model_name) {
  cases <- make_validation_cases()
  idx <- which(vapply(cases, `[[`, character(1), "model") == model_name)
  if (length(idx) != 1L) {
    stop(sprintf("Could not identify a unique validation case for model '%s'.", model_name))
  }
  idx
}

problem_case_seed_info <- function(model_name,
                                   replicate_id,
                                   M_outer = 5L,
                                   base_seed = 20260613L) {
  case_pos <- problem_case_order(model_name)
  global_idx <- (case_pos - 1L) * as.integer(M_outer) + as.integer(replicate_id)
  list(
    case_position = case_pos,
    global_idx = global_idx,
    sample_seed = as.integer(base_seed + 1000L * global_idx),
    bootstrap_seed = as.integer(base_seed + 2000L * global_idx),
    derivative_mc_seed = as.integer(base_seed + 3000L * global_idx)
  )
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

first_or_na <- function(x, default = NA) {
  if (length(x) == 0L) {
    return(default)
  }
  x[[1L]]
}

rbind_fill_frames <- function(frames) {
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(frames, names)))
  frames_aligned <- lapply(frames, function(df) {
    missing_names <- setdiff(all_names, names(df))
    for (nm in missing_names) {
      df[[nm]] <- NA
    }
    df[, all_names, drop = FALSE]
  })
  do.call(rbind, frames_aligned)
}

theta_diag_row <- function(case, theta_hat) {
  as.data.frame(
    validation_extract_theta_diagnostics(
      case = case,
      theta_hat = theta_hat,
      control = case$control
    ),
    stringsAsFactors = FALSE
  )
}

score_coord_summary <- function(Psi_aux) {
  Psi_aux <- as.matrix(Psi_aux)
  data.frame(
    score_coord = seq_len(ncol(Psi_aux)),
    score_mean = colMeans(Psi_aux),
    score_sd = apply(Psi_aux, 2L, stats::sd),
    stringsAsFactors = FALSE
  )
}

beta_score_block_summary <- function(Psi_aux) {
  Psi_aux <- as.matrix(Psi_aux)
  means <- colMeans(Psi_aux)
  sds <- apply(Psi_aux, 2L, stats::sd)
  data.frame(
    block = c("mu_1", "mu_2", "tau", "log_alpha1", "log_beta1", "log_alpha2", "log_beta2"),
    score_mean = means,
    score_sd = sds,
    stringsAsFactors = FALSE
  )
}

problem_prep_vhat_diagnostics <- function(prep,
                                          theta_hat) {
  if (!is.null(prep$vhat_diagnostics)) {
    return(prep$vhat_diagnostics)
  }

  par0 <- if (!is.null(theta_hat$xi)) {
    as.numeric(theta_hat$xi)
  } else if (!is.null(theta_hat$phi)) {
    as.numeric(theta_hat$phi)
  } else {
    rep(NA_real_, ncol(prep$Vhat))
  }

  fast_multiplier_vhat_diagnostics(
    S_obs = prep$S_obs,
    Psi_aux = prep$Psi_aux,
    Vhat = prep$Vhat,
    par0 = par0
  )
}

prepare_problem_fast_diag <- function(case,
                                      data,
                                      theta_hat,
                                      derivative_mc_size = 1000L,
                                      derivative_mc_seed = 20260613L) {
  spec <- case$spec_fn()
  control_fast <- utils::modifyList(
    case$control,
    list(
      derivative_method = "score_mc",
      derivative_mc_size = as.integer(derivative_mc_size),
      derivative_mc_seed = as.integer(derivative_mc_seed),
      fast_multiplier_allow_singular_vhat_diagnostics = TRUE
    )
  )
  if (identical(case$model, "cardioid")) {
    control_fast$cardioid_fast_boundary_eps <- -1
  }
  spec_fast_multiplier_prepare(
    spec = spec,
    data = data,
    theta_hat = theta_hat,
    ks_prep = NULL,
    cvm_prep = NULL,
    control = control_fast
  )
}

diagnose_validation_problem_case <- function(case,
                                             replicate_id,
                                             validation_raw,
                                             M_outer = 5L,
                                             base_seed = 20260613L,
                                             derivative_mc_size = 1000L) {
  seed_info <- problem_case_seed_info(
    model_name = case$model,
    replicate_id = replicate_id,
    M_outer = M_outer,
    base_seed = base_seed
  )
  set.seed(seed_info$sample_seed)
  data <- case$sample_fn(case$n)
  spec <- case$spec_fn()
  theta_hat <- spec$fit_theta(
    data = data,
    weights = NULL,
    null = case$wrapper_args$null,
    control = case$control
  )

  prep <- prepare_problem_fast_diag(
    case = case,
    data = data,
    theta_hat = theta_hat,
    derivative_mc_size = derivative_mc_size,
    derivative_mc_seed = seed_info$derivative_mc_seed
  )
  theta_diag <- theta_diag_row(case, theta_hat)
  raw_sub <- validation_raw[
    validation_raw$model == case$model &
      validation_raw$scenario == case$scenario &
      validation_raw$replicate_id == replicate_id,
    ,
    drop = FALSE
  ]
  raw_ks <- raw_sub[raw_sub$statistic == "ks", , drop = FALSE]
  raw_cvm <- raw_sub[raw_sub$statistic == "cvm", , drop = FALSE]
  prep <- tryCatch(
    prepare_problem_fast_diag(
      case = case,
      data = data,
      theta_hat = theta_hat,
      derivative_mc_size = derivative_mc_size,
      derivative_mc_seed = seed_info$derivative_mc_seed
    ),
    error = identity
  )
  prep_diag <- if (inherits(prep, "error")) NULL else problem_prep_vhat_diagnostics(prep, theta_hat)

  out <- data.frame(
    model = case$model,
    scenario = case$scenario,
    replicate_id = replicate_id,
    sample_seed = seed_info$sample_seed,
    bootstrap_seed = seed_info$bootstrap_seed,
    derivative_mc_seed = seed_info$derivative_mc_seed,
    old_p_value_ks = first_or_na(raw_ks$old_p_value, NA_real_),
    fast_p_value_ks = first_or_na(raw_ks$fast_p_value, NA_real_),
    old_p_value_cvm = first_or_na(raw_cvm$old_p_value, NA_real_),
    fast_p_value_cvm = first_or_na(raw_cvm$fast_p_value, NA_real_),
    abs_diff_ks = first_or_na(raw_ks$abs_p_value_diff, NA_real_),
    abs_diff_cvm = first_or_na(raw_cvm$abs_p_value_diff, NA_real_),
    fast_error_ks = first_or_na(raw_ks$fast_error, NA_character_),
    fast_error_cvm = first_or_na(raw_cvm$fast_error, NA_character_),
    prep_error = if (inherits(prep, "error")) conditionMessage(prep) else NA_character_,
    score_mean_norm = if (inherits(prep, "error")) NA_real_ else prep_diag$score_mean_aux_norm,
    Vhat_min_eigenvalue = if (inherits(prep, "error")) NA_real_ else safe_min(prep_diag$Vhat_eigenvalues),
    Vhat_rcond = if (inherits(prep, "error")) NA_real_ else prep_diag$Vhat_rcond,
    Vhat_condition_number = if (inherits(prep, "error")) NA_real_ else prep_diag$Vhat_condition_number,
    stringsAsFactors = FALSE
  )

  cbind(out, theta_diag)
}

diagnose_beta_paper_failure <- function(output_root,
                                        derivative_mc_size = 1000L,
                                        max_replicates = 1000L) {
  scenario <- make_beta_mixture2_composite_calibration_scenario(
    weight1 = 0.4,
    alpha1 = 2,
    beta1 = 8,
    alpha2 = 8,
    beta2 = 2
  )

  scenario_seed <- 20260620L + 1000L
  n_seed <- scenario_seed + 100000L
  task_grid <- build_calibration_tasks_single_n(
    n_value = 100L,
    M_outer = as.integer(max_replicates),
    seed = n_seed
  )

  spec <- make_beta_mixture2_spec(distance_type = "geodesic")
  rows <- list()
  first_failure <- NULL

  for (i in seq_len(nrow(task_grid))) {
    task_row <- task_grid[i, , drop = FALSE]
    set.seed(task_row$sample_seed)
    data <- simulate_h0_sample(
      scenario = scenario,
      n = task_row$n,
      replicate_id = task_row$replicate_id
    )
    theta_hat <- spec$fit_theta(
      data = data,
      weights = NULL,
      null = scenario$null,
      control = scenario$control
    )
    prep <- prepare_problem_fast_diag(
      case = list(
        model = "beta_mixture2",
        scenario = scenario$id,
        n = task_row$n,
        spec_fn = function() spec,
        wrapper_args = list(null = scenario$null),
        control = scenario$control
      ),
      data = data,
      theta_hat = theta_hat,
      derivative_mc_size = derivative_mc_size,
      derivative_mc_seed = 20260613L
    )
    prep_diag <- problem_prep_vhat_diagnostics(prep, theta_hat)
    theta_diag <- as.list(theta_diag_row(
      list(model = "beta_mixture2", control = scenario$control),
      theta_hat
    )[1, , drop = TRUE])
    row <- data.frame(
      paper_scenario = scenario$id,
      n = task_row$n,
      replicate_id = task_row$replicate_id,
      sample_seed = task_row$sample_seed,
      bootstrap_seed = task_row$bootstrap_seed,
      score_mean_norm = prep_diag$score_mean_aux_norm,
      Vhat_min_eigenvalue = safe_min(prep_diag$Vhat_eigenvalues),
      Vhat_rcond = prep_diag$Vhat_rcond,
      Vhat_condition_number = prep_diag$Vhat_condition_number,
      stringsAsFactors = FALSE
    )
    rows[[i]] <- cbind(row, as.data.frame(theta_diag, stringsAsFactors = FALSE))

    invalid <- !is.finite(prep_diag$Vhat_rcond) ||
      prep_diag$Vhat_rcond <= 1e-12 ||
      safe_min(prep_diag$Vhat_eigenvalues) <= 0
    if (isTRUE(invalid)) {
      first_failure <- rows[[i]]
      break
    }
  }

  out_df <- do.call(rbind, rows)
  utils::write.csv(out_df, file.path(output_root, "beta_mixture2_paper_scan.csv"), row.names = FALSE)
  if (!is.null(first_failure)) {
    utils::write.csv(first_failure, file.path(output_root, "beta_mixture2_paper_first_failure.csv"), row.names = FALSE)
  }

  list(scan = out_df, first_failure = first_failure)
}

run_fast_multiplier_problem_diagnostics <- function(output_root = file.path("output", "fast_multiplier", "problem_diagnostics"),
                                                    validation_raw_csv = file.path("output", "fast_multiplier", "validation_all_models_12cores_B1000_M5", "validation_raw.csv"),
                                                    M_outer = 5L,
                                                    base_seed = 20260613L,
                                                    derivative_mc_size = 1000L) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  validation_raw <- utils::read.csv(validation_raw_csv, stringsAsFactors = FALSE)
  cases <- make_validation_cases()
  cases <- cases[vapply(cases, function(case) case$model %in% c("beta_mixture2", "jp", "cardioid"), logical(1))]

  validation_rows <- list()
  beta_block_rows <- list()
  beta_coord_rows <- list()
  idx <- 1L

  for (case in cases) {
    for (replicate_id in seq_len(M_outer)) {
      row <- diagnose_validation_problem_case(
        case = case,
        replicate_id = replicate_id,
        validation_raw = validation_raw,
        M_outer = M_outer,
        base_seed = base_seed,
        derivative_mc_size = derivative_mc_size
      )
      validation_rows[[idx]] <- row

      if (identical(case$model, "beta_mixture2")) {
        seed_info <- problem_case_seed_info(case$model, replicate_id, M_outer = M_outer, base_seed = base_seed)
        set.seed(seed_info$sample_seed)
        data <- case$sample_fn(case$n)
        theta_hat <- case$spec_fn()$fit_theta(data = data, weights = NULL, null = case$wrapper_args$null, control = case$control)
        prep <- prepare_problem_fast_diag(
          case = case,
          data = data,
          theta_hat = theta_hat,
          derivative_mc_size = derivative_mc_size,
          derivative_mc_seed = seed_info$derivative_mc_seed
        )
        block_df <- beta_score_block_summary(prep$Psi_aux)
        block_df$model <- case$model
        block_df$scenario <- case$scenario
        block_df$replicate_id <- replicate_id
        beta_block_rows[[length(beta_block_rows) + 1L]] <- block_df

        coord_df <- score_coord_summary(prep$Psi_aux)
        coord_df$model <- case$model
        coord_df$scenario <- case$scenario
        coord_df$replicate_id <- replicate_id
        beta_coord_rows[[length(beta_coord_rows) + 1L]] <- coord_df
      }

      idx <- idx + 1L
    }
  }

  validation_df <- rbind_fill_frames(validation_rows)
  beta_block_df <- rbind_fill_frames(beta_block_rows)
  beta_coord_df <- rbind_fill_frames(beta_coord_rows)

  utils::write.csv(validation_df, file.path(output_root, "validation_problem_diagnostics.csv"), row.names = FALSE)
  utils::write.csv(beta_block_df, file.path(output_root, "beta_mixture2_score_block_diagnostics.csv"), row.names = FALSE)
  utils::write.csv(beta_coord_df, file.path(output_root, "beta_mixture2_score_coordinate_diagnostics.csv"), row.names = FALSE)

  paper_beta <- diagnose_beta_paper_failure(
    output_root = output_root,
    derivative_mc_size = derivative_mc_size,
    max_replicates = 1000L
  )

  list(
    validation = validation_df,
    beta_block = beta_block_df,
    beta_coord = beta_coord_df,
    paper_beta = paper_beta
  )
}

parse_fast_problem_diag_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

if (sys.nframe() == 0L) {
  args <- parse_fast_problem_diag_args(commandArgs(trailingOnly = TRUE))
  result <- run_fast_multiplier_problem_diagnostics(
    output_root = args$output_root %||% file.path("output", "fast_multiplier", "problem_diagnostics"),
    validation_raw_csv = args$validation_raw_csv %||% file.path("output", "fast_multiplier", "validation_all_models_12cores_B1000_M5", "validation_raw.csv"),
    M_outer = as.integer(args$M_outer %||% 5L),
    base_seed = as.integer(args$base_seed %||% 20260613L),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 1000L)
  )
  cat("Validation diagnostics:", file.path(args$output_root %||% file.path("output", "fast_multiplier", "problem_diagnostics"), "validation_problem_diagnostics.csv"), "\n")
  cat("Beta block diagnostics:", file.path(args$output_root %||% file.path("output", "fast_multiplier", "problem_diagnostics"), "beta_mixture2_score_block_diagnostics.csv"), "\n")
  cat("Beta paper scan:", file.path(args$output_root %||% file.path("output", "fast_multiplier", "problem_diagnostics"), "beta_mixture2_paper_scan.csv"), "\n")
}
