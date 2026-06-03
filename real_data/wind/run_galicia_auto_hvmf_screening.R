source(file.path("wind", "run_kolkata_monthly_hvmf_screening.R"))

month_label_slug <- function(months_vec) {
  paste(tolower(month.abb[as.integer(months_vec)]), collapse = "_")
}

choose_thinning_for_block <- function(clean_df,
                                      steps = 4:7,
                                      min_obs = 100L,
                                      max_obs = 300L,
                                      max_day = 31L) {
  candidates <- list()
  idx <- 1L

  for (step in as.integer(steps)) {
    for (start_day in seq_len(step)) {
      thinned_df <- thin_kolkata_by_day_pattern(
        df = clean_df,
        start_day = start_day,
        step = step,
        max_day = max_day
      )
      n_thinned <- nrow(thinned_df)

      if (n_thinned >= as.integer(min_obs) && n_thinned <= as.integer(max_obs)) {
        candidates[[idx]] <- list(
          step = as.integer(step),
          start_day = as.integer(start_day),
          n_thinned = as.integer(n_thinned),
          thinned_df = thinned_df
        )
        idx <- idx + 1L
      }
    }
  }

  if (length(candidates) == 0L) {
    return(NULL)
  }

  candidate_df <- do.call(rbind, lapply(candidates, function(x) {
    data.frame(
      step = x$step,
      start_day = x$start_day,
      n_thinned = x$n_thinned,
      distance_to_center = abs(x$n_thinned - 200L),
      stringsAsFactors = FALSE
    )
  }))

  best_idx <- order(
    candidate_df$step,
    candidate_df$distance_to_center,
    candidate_df$start_day
  )[[1L]]

  candidates[[best_idx]]
}

build_case_definition <- function(all_daily_df,
                                  months_vec,
                                  steps = 4:7,
                                  min_obs = 100L,
                                  max_obs = 300L) {
  months_vec <- sort(unique(as.integer(months_vec)))
  block_df <- all_daily_df[all_daily_df$month %in% months_vec, , drop = FALSE]
  block_df <- block_df[order(block_df$date), , drop = FALSE]
  rownames(block_df) <- NULL

  clean_df <- block_df[
    is.finite(block_df$WD10M) &
      is.finite(block_df$WS10M) &
      block_df$WS10M > 0,
    ,
    drop = FALSE
  ]
  rownames(clean_df) <- NULL

  thinning <- choose_thinning_for_block(
    clean_df = clean_df,
    steps = steps,
    min_obs = min_obs,
    max_obs = max_obs,
    max_day = 31L
  )

  if (is.null(thinning)) {
    return(NULL)
  }

  list(
    months_vec = months_vec,
    case_label = month_label_slug(months_vec),
    step = thinning$step,
    start_day = thinning$start_day,
    n_raw = nrow(block_df),
    n_clean = nrow(clean_df),
    n_thinned = thinning$n_thinned,
    block_df = block_df,
    clean_df = clean_df,
    thinned_df = thinning$thinned_df
  )
}

build_automatic_case_set <- function(all_daily_df,
                                     steps = 4:7,
                                     min_obs = 100L,
                                     max_obs = 300L,
                                     granularity_mode = "automatic") {
  granularity_mode <- match.arg(
    granularity_mode,
    choices = c("automatic", "monthly", "bimonthly")
  )

  monthly_cases <- lapply(1:12, function(m) {
    build_case_definition(
      all_daily_df = all_daily_df,
      months_vec = m,
      steps = steps,
      min_obs = min_obs,
      max_obs = max_obs
    )
  })

  monthly_ok <- vapply(monthly_cases, function(x) !is.null(x), logical(1))

  if (granularity_mode == "monthly") {
    selected <- monthly_cases[monthly_ok]
    design_mode <- "monthly"
  } else if (granularity_mode == "bimonthly") {
    pair_months <- split(1:12, ceiling((1:12) / 2))
    bimonthly_cases <- lapply(pair_months, function(block) {
      build_case_definition(
        all_daily_df = all_daily_df,
        months_vec = block,
        steps = steps,
        min_obs = min_obs,
        max_obs = max_obs
      )
    })
    selected <- bimonthly_cases[vapply(bimonthly_cases, function(x) !is.null(x), logical(1))]
    design_mode <- "bimonthly"
  } else {
    if (all(monthly_ok)) {
      selected <- monthly_cases
      design_mode <- "monthly"
    } else {
      pair_months <- split(1:12, ceiling((1:12) / 2))
      bimonthly_cases <- lapply(pair_months, function(block) {
        build_case_definition(
          all_daily_df = all_daily_df,
          months_vec = block,
          steps = steps,
          min_obs = min_obs,
          max_obs = max_obs
        )
      })
      bimonthly_ok <- vapply(bimonthly_cases, function(x) !is.null(x), logical(1))
      if (!all(bimonthly_ok)) {
        stop("Could not build a complete automatic case set with requested constraints.")
      }
      selected <- bimonthly_cases
      design_mode <- "bimonthly"
    }
  }

  if (length(selected) == 0L) {
    stop("No valid cases were found with the requested thinning constraints.")
  }

  selected <- selected[order(vapply(selected, function(x) x$months_vec[[1L]], integer(1)))]

  design_table <- do.call(rbind, lapply(seq_along(selected), function(i) {
    x <- selected[[i]]
    data.frame(
      case_id = as.integer(i),
      case_label = x$case_label,
      months = paste(month.abb[x$months_vec], collapse = "+"),
      month_start = min(x$months_vec),
      month_end = max(x$months_vec),
      step = x$step,
      start_day = x$start_day,
      n_raw = x$n_raw,
      n_clean = x$n_clean,
      n_thinned = x$n_thinned,
      design_mode = design_mode,
      stringsAsFactors = FALSE
    )
  }))

  list(
    cases = selected,
    design_mode = design_mode,
    design_table = design_table
  )
}

embed_case_h2 <- function(case_def) {
  wind_df_thinned <- data.frame(
    datetime = as.POSIXct(case_def$thinned_df$date, tz = "UTC"),
    ws10m = case_def$thinned_df$WS10M,
    wd10m = case_def$thinned_df$WD10M,
    stringsAsFactors = FALSE
  )

  embedded_df <- build_hvmf_wind_set(
    df = wind_df_thinned,
    speed_col = "ws10m",
    direction_col = "wd10m",
    height_m = 10L,
    fixed_tz = "UTC"
  )

  embedded_df$months_id <- case_def$case_label
  embedded_df$speed_mean <- embedded_df$speed_mean_height
  embedded_df$dataset_source <- "NASA POWER CERES/MERRA2 daily point"
  embedded_df$start_day <- as.integer(case_def$start_day)
  embedded_df$step <- as.integer(case_def$step)

  embedded_df
}

run_galicia_cases <- function(all_daily_df,
                              case_set,
                              location_name,
                              location_slug,
                              start_year,
                              end_year,
                              B,
                              n_cores,
                              profile_method,
                              output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  processed_dir <- file.path(output_dir, "processed")
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

  summary_rows <- list()

  for (i in seq_along(case_set$cases)) {
    case_def <- case_set$cases[[i]]
    dataset_id <- sprintf(
      "galicia_%s_%s_%d_%d_step%d_start%d",
      location_slug,
      case_def$case_label,
      as.integer(start_year),
      as.integer(end_year),
      as.integer(case_def$step),
      as.integer(case_def$start_day)
    )

    cat(sprintf("Running Galicia case %s\n", dataset_id))

    case_row <- tryCatch(
      {
        embedded_df <- embed_case_h2(case_def)

        raw_csv <- file.path(processed_dir, sprintf("%s_raw_daily.csv", dataset_id))
        thinned_csv <- file.path(processed_dir, sprintf("%s_thinned.csv", dataset_id))
        utils::write.csv(case_def$clean_df, raw_csv, row.names = FALSE)
        utils::write.csv(case_def$thinned_df, thinned_csv, row.names = FALSE)

        benchmark_row <- run_composite_benchmark_case(
          df = embedded_df,
          dataset_id = dataset_id,
          output_dir = output_dir,
          B = B,
          n_cores = n_cores,
          profile_method = profile_method
        )

        diag_case <- diagnose_power_daily_data(case_def$block_df)

        benchmark_row$location_name <- location_name
        benchmark_row$case_label <- case_def$case_label
        benchmark_row$months <- paste(month.abb[case_def$months_vec], collapse = "+")
        benchmark_row$month_start <- min(case_def$months_vec)
        benchmark_row$month_end <- max(case_def$months_vec)
        benchmark_row$step <- case_def$step
        benchmark_row$start_day <- case_def$start_day
        benchmark_row$n_raw_case <- case_def$n_raw
        benchmark_row$n_clean_case <- case_def$n_clean
        benchmark_row$n_thinned_case <- case_def$n_thinned
        benchmark_row$n_duplicate_dates <- as.integer(diag_case$n_duplicate_dates)
        benchmark_row$n_na_wd10m <- as.integer(diag_case$n_na_wd10m)
        benchmark_row$n_na_ws10m <- as.integer(diag_case$n_na_ws10m)
        benchmark_row$n_nonpositive_ws10m <- as.integer(diag_case$n_nonpositive_ws10m)
        benchmark_row$n_wd10m_outside_0_360 <- as.integer(diag_case$n_wd10m_outside_0_360)
        benchmark_row$raw_case_csv <- raw_csv
        benchmark_row$thinned_case_csv <- thinned_csv
        benchmark_row$status <- "ok"
        benchmark_row$error_message <- NA_character_
        benchmark_row
      },
      error = function(e) {
        failed <- data.frame(
          dataset_id = dataset_id,
          location_name = location_name,
          case_label = case_def$case_label,
          months = paste(month.abb[case_def$months_vec], collapse = "+"),
          month_start = min(case_def$months_vec),
          month_end = max(case_def$months_vec),
          step = case_def$step,
          start_day = case_def$start_day,
          n_raw_case = case_def$n_raw,
          n_clean_case = case_def$n_clean,
          n_thinned_case = case_def$n_thinned,
          kappa_hat = NA_real_,
          theta_deg_hat = NA_real_,
          statistic_CvM_composite = NA_real_,
          p_value_composite = NA_real_,
          B = B,
          n_cores = n_cores,
          processed_csv = NA_character_,
          result_rds = NA_character_,
          log_file = NA_character_,
          elapsed_seconds = NA_real_,
          status = "failed",
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
        failed
      }
    )

    summary_rows[[i]] <- case_row
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  summary_df
}

pilot_location_score <- function(summary_df) {
  if (!"p_value_composite" %in% names(summary_df)) {
    return(list(significant_rate = NA_real_, mean_kappa = NA_real_))
  }

  pvals <- suppressWarnings(as.numeric(summary_df$p_value_composite))
  kappas <- suppressWarnings(as.numeric(summary_df$kappa_hat))

  list(
    significant_rate = mean(pvals < 0.05, na.rm = TRUE),
    mean_kappa = mean(kappas, na.rm = TRUE)
  )
}

select_best_galicia_location <- function(start_year,
                                         end_year,
                                         B_select = 99L,
                                         n_cores_select = 2L,
                                         profile_method = "tabulated",
                                         granularity_mode = "automatic",
                                         output_dir = file.path("wind", "galicia_auto_pipeline")) {
  locations <- data.frame(
    slug = c("santiago", "acoruna", "vigo", "lugo", "ourense", "fisterra"),
    name = c("Santiago de Compostela", "A Coruna", "Vigo", "Lugo", "Ourense", "Fisterra"),
    latitude = c(42.8782, 43.3623, 42.2406, 43.0121, 42.3358, 42.9049),
    longitude = c(-8.5448, -8.4115, -8.7207, -7.5559, -7.8639, -9.2629),
    stringsAsFactors = FALSE
  )

  pilot_root <- file.path(output_dir, "pilot_selection")
  dir.create(pilot_root, recursive = TRUE, showWarnings = FALSE)

  score_rows <- list()

  for (i in seq_len(nrow(locations))) {
    loc <- locations[i, , drop = FALSE]
    cat(sprintf("Pilot selection for %s\n", loc$name[[1L]]))

    loc_result <- tryCatch({
      raw <- NULL
      for (attempt in 1:3) {
        raw <- tryCatch(
          load_kolkata_power_daily(
            start_year = start_year,
            end_year = end_year,
            latitude = loc$latitude[[1L]],
            longitude = loc$longitude[[1L]]
          ),
          error = function(e) {
            cat(sprintf("  Download attempt %d failed for %s: %s\n", attempt, loc$name[[1L]], conditionMessage(e)))
            if (attempt < 3L) Sys.sleep(5)
            NULL
          }
        )
        if (!is.null(raw)) break
      }
      if (is.null(raw)) stop(sprintf("All download attempts failed for %s.", loc$name[[1L]]))

      case_set <- build_automatic_case_set(
        all_daily_df = raw$data,
        steps = 4:7,
        min_obs = 100L,
        max_obs = 300L,
        granularity_mode = granularity_mode
      )

      pilot_dir <- file.path(pilot_root, loc$slug[[1L]])
      summary_df <- run_galicia_cases(
        all_daily_df = raw$data,
        case_set = case_set,
        location_name = loc$name[[1L]],
        location_slug = loc$slug[[1L]],
        start_year = start_year,
        end_year = end_year,
        B = B_select,
        n_cores = n_cores_select,
        profile_method = profile_method,
        output_dir = pilot_dir
      )

      score <- pilot_location_score(summary_df)

      design_csv <- file.path(pilot_dir, "case_design.csv")
      summary_csv <- file.path(pilot_dir, "pilot_summary.csv")
      utils::write.csv(case_set$design_table, design_csv, row.names = FALSE)
      utils::write.csv(summary_df, summary_csv, row.names = FALSE)

      data.frame(
        slug = loc$slug[[1L]],
        name = loc$name[[1L]],
        latitude = loc$latitude[[1L]],
        longitude = loc$longitude[[1L]],
        design_mode = case_set$design_mode,
        n_cases = nrow(case_set$design_table),
        significant_rate = score$significant_rate,
        mean_kappa = score$mean_kappa,
        download_ok = TRUE,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat(sprintf("  Skipping %s: %s\n", loc$name[[1L]], conditionMessage(e)))
      data.frame(
        slug = loc$slug[[1L]],
        name = loc$name[[1L]],
        latitude = loc$latitude[[1L]],
        longitude = loc$longitude[[1L]],
        design_mode = NA_character_,
        n_cases = NA_integer_,
        significant_rate = NA_real_,
        mean_kappa = NA_real_,
        download_ok = FALSE,
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })

    score_rows[[i]] <- loc_result
  }

  score_df <- do.call(rbind, score_rows)
  ok_rows <- !is.na(score_df$download_ok) & score_df$download_ok
  if (!any(ok_rows)) {
    stop("All location downloads failed. Check network connectivity.")
  }
  score_df_ok <- score_df[ok_rows, , drop = FALSE]
  score_df_ok <- score_df_ok[order(-score_df_ok$significant_rate, -score_df_ok$mean_kappa), , drop = FALSE]
  score_df <- rbind(score_df_ok, score_df[!ok_rows, , drop = FALSE])
  rownames(score_df) <- NULL

  score_csv <- file.path(pilot_root, "location_scores.csv")
  utils::write.csv(score_df, score_csv, row.names = FALSE)

  best <- score_df_ok[1L, , drop = FALSE]

  list(
    best = best,
    scores = score_df,
    score_csv = score_csv
  )
}

run_galicia_auto_pipeline <- function(start_year = 1982L,
                                      end_year = 2022L,
                                      granularity_mode = "automatic",
                                      B = 500L,
                                      n_cores = 6L,
                                      B_select = 99L,
                                      n_cores_select = 2L,
                                      profile_method = "tabulated",
                                      output_dir = file.path("wind", "galicia_auto_pipeline")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  selection <- select_best_galicia_location(
    start_year = start_year,
    end_year = end_year,
    B_select = B_select,
    n_cores_select = n_cores_select,
    profile_method = profile_method,
    granularity_mode = granularity_mode,
    output_dir = output_dir
  )

  best <- selection$best
  cat(sprintf("Selected location: %s (%s, %s)\n", best$name[[1L]], best$latitude[[1L]], best$longitude[[1L]]))

  raw <- load_kolkata_power_daily(
    start_year = start_year,
    end_year = end_year,
    latitude = best$latitude[[1L]],
    longitude = best$longitude[[1L]]
  )

  case_set <- build_automatic_case_set(
    all_daily_df = raw$data,
    steps = c(7L),
    min_obs = 100L,
    max_obs = 300L,
    granularity_mode = granularity_mode
  )

  final_dir <- file.path(output_dir, sprintf("final_%s", best$slug[[1L]]))
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

  all_daily_csv <- file.path(final_dir, sprintf("galicia_%s_daily_%d_%d_all_months.csv", best$slug[[1L]], as.integer(start_year), as.integer(end_year)))
  design_csv <- file.path(final_dir, "case_design.csv")
  utils::write.csv(raw$data, all_daily_csv, row.names = FALSE)
  utils::write.csv(case_set$design_table, design_csv, row.names = FALSE)

  summary_df <- run_galicia_cases(
    all_daily_df = raw$data,
    case_set = case_set,
    location_name = best$name[[1L]],
    location_slug = best$slug[[1L]],
    start_year = start_year,
    end_year = end_year,
    B = B,
    n_cores = n_cores,
    profile_method = profile_method,
    output_dir = final_dir
  )

  summary_csv <- file.path(
    final_dir,
    sprintf("galicia_hvmf_screening_%s_%d_%d_B%d.csv", best$slug[[1L]], as.integer(start_year), as.integer(end_year), as.integer(B))
  )
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  qa_global <- diagnose_power_daily_data(raw$data)
  qa_global_df <- data.frame(
    location_name = best$name[[1L]],
    location_slug = best$slug[[1L]],
    latitude = best$latitude[[1L]],
    longitude = best$longitude[[1L]],
    n_rows = qa_global$n_rows,
    n_duplicate_dates = qa_global$n_duplicate_dates,
    n_na_wd10m = qa_global$n_na_wd10m,
    n_na_ws10m = qa_global$n_na_ws10m,
    n_nonpositive_ws10m = qa_global$n_nonpositive_ws10m,
    n_wd10m_outside_0_360 = qa_global$n_wd10m_outside_0_360,
    max_identical_wd_run = qa_global$max_identical_wd_run,
    max_identical_ws_run = qa_global$max_identical_ws_run,
    max_identical_pair_run = qa_global$max_identical_pair_run,
    stringsAsFactors = FALSE
  )
  qa_global_csv <- file.path(final_dir, "qa_global.csv")
  utils::write.csv(qa_global_df, qa_global_csv, row.names = FALSE)

  list(
    selected_location = best,
    selection_scores_csv = selection$score_csv,
    design_csv = design_csv,
    all_daily_csv = all_daily_csv,
    summary_csv = summary_csv,
    qa_global_csv = qa_global_csv,
    source_url = raw$url,
    output_dir = final_dir
  )
}

if (sys.nframe() == 0L) {
  result <- run_galicia_auto_pipeline(
    start_year = 1982L,
    end_year = 2022L,
    granularity_mode = "automatic",
    B = 500L,
    n_cores = 6L,
    B_select = 99L,
    n_cores_select = 2L,
    profile_method = "tabulated",
    output_dir = file.path("wind", "galicia_auto_pipeline")
  )

  cat("selected_location=", result$selected_location$name[[1L]], "\n", sep = "")
  cat("selection_scores_csv=", result$selection_scores_csv, "\n", sep = "")
  cat("design_csv=", result$design_csv, "\n", sep = "")
  cat("all_daily_csv=", result$all_daily_csv, "\n", sep = "")
  cat("summary_csv=", result$summary_csv, "\n", sep = "")
  cat("qa_global_csv=", result$qa_global_csv, "\n", sep = "")
  cat("source_url=", result$source_url, "\n", sep = "")
}
