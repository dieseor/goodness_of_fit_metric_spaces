source(file.path("wind", "preprocess_risoe_modern_hvmf.R"))
source(file.path("wind", "run_hvmf_real_data_cvm.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))

jensen_like_day_patterns <- function() {
  list(
    set12 = c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L),
    set3 = c(1L, 5L, 9L, 13L, 17L, 21L, 25L, 29L)
  )
}

filter_month_scope <- function(df, month_scope) {
  month_scope <- match.arg(month_scope, c("nov_dec", "nov", "dec"))
  keep_months <- switch(
    month_scope,
    nov_dec = c(11L, 12L),
    nov = 11L,
    dec = 12L
  )

  output <- df[df$month %in% keep_months, , drop = FALSE]
  rownames(output) <- NULL
  output
}

filter_day_pattern <- function(df, day_pattern) {
  output <- df[df$day %in% as.integer(day_pattern), , drop = FALSE]
  rownames(output) <- NULL
  output
}

build_jensen_like_hvmf_set <- function(selected_df,
                                       speed_col,
                                       direction_col,
                                       height_m,
                                       day_pattern,
                                       month_scope,
                                       fixed_tz = "UTC") {
  scoped_df <- filter_month_scope(selected_df, month_scope = month_scope)
  pattern_df <- filter_day_pattern(scoped_df, day_pattern = day_pattern)

  if (nrow(pattern_df) == 0L) {
    stop("No rows remain after Jensen-like month/day filtering.")
  }

  build_hvmf_wind_set(
    df = pattern_df,
    speed_col = speed_col,
    direction_col = direction_col,
    height_m = height_m,
    fixed_tz = fixed_tz
  )
}

write_jensen_like_hvmf_csv <- function(df, output_csv, fixed_tz = "UTC") {
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write_hvmf_csv(df, output_csv, fixed_tz = fixed_tz)
  output_csv
}

run_simple_hvmf_screening_case <- function(dataset_name,
                                           output_csv,
                                           B = 5000L,
                                           n_cores = 10L,
                                           alpha = 0.05,
                                           seed = 52000L) {
  prepared <- prepare_existing_h2_csv_dataset(output_csv, tol = 1e-8)
  theta_hat <- fit_hvmf_theta(
    data = prepared$data_matrix,
    weights = NULL,
    null = list(type = "composite"),
    unknown_param = "both",
    control = list()
  )

  start_time <- Sys.time()
  result <- multiplier_bootstrap_hvmf(
    data = prepared$data_matrix,
    null = list(type = "simple", theta = theta_hat),
    statistics = "cvm",
    B = B,
    alpha = alpha,
    n_cores = n_cores,
    seed = seed,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE
    ),
    unknown_param = "both"
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  data.frame(
    dataset = dataset_name,
    n = prepared$n,
    kappa_hat = theta_hat$kappa,
    mu_hat_x0 = theta_hat$mu[[1L]],
    mu_hat_x1 = theta_hat$mu[[2L]],
    mu_hat_x2 = theta_hat$mu[[3L]],
    sinh_chi_hat = theta_hat$sinh_chi,
    theta_deg_hat = theta_hat$theta_deg,
    statistic_CvM = result$observed$cvm$statistic,
    p_value = result$inference$cvm$p_value,
    B = B,
    n_cores = n_cores,
    elapsed_seconds = elapsed_seconds,
    output_csv = output_csv,
    stringsAsFactors = FALSE
  )
}

make_jensen_like_screening_configs <- function(output_dir = file.path("wind", "jensen_like_screening")) {
  patterns <- jensen_like_day_patterns()
  month_scopes <- c("nov_dec", "nov", "dec")
  heights <- list(
    list(id = "77m", speed_col = "ws77", direction_col = "wd77", height_m = 77L),
    list(id = "125m", speed_col = "ws125", direction_col = "wd125", height_m = 125L)
  )

  configs <- list()
  idx <- 1L

  for (pattern_name in names(patterns)) {
    for (month_scope in month_scopes) {
      for (height in heights) {
        dataset_name <- sprintf("risoe_%s_%s_%s", height$id, month_scope, pattern_name)
        output_csv <- file.path(output_dir, paste0(dataset_name, "_hvmf.csv"))
        configs[[idx]] <- list(
          dataset_name = dataset_name,
          pattern_name = pattern_name,
          day_pattern = patterns[[pattern_name]],
          month_scope = month_scope,
          speed_col = height$speed_col,
          direction_col = height$direction_col,
          height_m = height$height_m,
          output_csv = output_csv
        )
        idx <- idx + 1L
      }
    }
  }

  configs
}

run_risoe_jensen_like_screening <- function(input_nc = "wind/risoe_m_concurent.nc",
                                            output_dir = file.path("wind", "jensen_like_screening"),
                                            B = 5000L,
                                            n_cores = 10L,
                                            fixed_tz = "UTC",
                                            tie_break = "earliest") {
  concurrent_df <- load_risoe_concurrent(input_nc, fixed_tz = fixed_tz)
  selected_df <- select_noon_nov_dec(concurrent_df, tie_break = tie_break, fixed_tz = fixed_tz)
  configs <- make_jensen_like_screening_configs(output_dir = output_dir)

  summary_rows <- vector("list", length(configs))

  for (i in seq_along(configs)) {
    config <- configs[[i]]
    cat("\n=== Jensen-like screening ", i, "/", length(configs), ": ", config$dataset_name, " ===\n", sep = "")
    df_case <- build_jensen_like_hvmf_set(
      selected_df = selected_df,
      speed_col = config$speed_col,
      direction_col = config$direction_col,
      height_m = config$height_m,
      day_pattern = config$day_pattern,
      month_scope = config$month_scope,
      fixed_tz = fixed_tz
    )
    write_jensen_like_hvmf_csv(df_case, config$output_csv, fixed_tz = fixed_tz)
    summary_rows[[i]] <- cbind(
      data.frame(
        pattern_name = config$pattern_name,
        month_scope = config$month_scope,
        height_m = config$height_m,
        stringsAsFactors = FALSE
      ),
      run_simple_hvmf_screening_case(
        dataset_name = config$dataset_name,
        output_csv = config$output_csv,
        B = B,
        n_cores = n_cores,
        seed = 52000L + i
      )
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_csv <- file.path(output_dir, "risoe_jensen_like_screening_summary.csv")
  write.csv(summary_df, summary_csv, row.names = FALSE)
  summary_df <- summary_df[order(summary_df$p_value, decreasing = TRUE), , drop = FALSE]
  rownames(summary_df) <- NULL

  cat("\nTop Jensen-like simple-fit candidates by p-value:\n")
  print(summary_df[, c("dataset", "n", "pattern_name", "month_scope", "height_m", "statistic_CvM", "p_value")], row.names = FALSE)

  invisible(list(
    summary = summary_df,
    summary_csv = summary_csv
  ))
}

if (sys.nframe() == 0L) {
  run_risoe_jensen_like_screening()
}
