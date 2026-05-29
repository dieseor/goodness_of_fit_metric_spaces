source(file.path("wind", "preprocess_risoe_modern_hvmf.R"))

download_power_daily_point <- function(latitude,
                                       longitude,
                                       parameters,
                                       start,
                                       end,
                                       community = "RE") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Install it with install.packages('jsonlite').", call. = FALSE)
  }

  parameter_string <- paste(parameters, collapse = ",")
  url <- sprintf(
    paste0(
      "https://power.larc.nasa.gov/api/temporal/daily/point",
      "?parameters=%s&community=%s&longitude=%s&latitude=%s&start=%s&end=%s&format=JSON"
    ),
    utils::URLencode(parameter_string, reserved = TRUE),
    utils::URLencode(community, reserved = TRUE),
    utils::URLencode(as.character(longitude), reserved = TRUE),
    utils::URLencode(as.character(latitude), reserved = TRUE),
    utils::URLencode(as.character(start), reserved = TRUE),
    utils::URLencode(as.character(end), reserved = TRUE)
  )

  payload <- jsonlite::fromJSON(url)
  list(url = url, payload = payload)
}

power_parameter_to_numeric <- function(parameter_block, name) {
  if (is.null(parameter_block[[name]])) {
    stop(sprintf("NASA POWER response does not contain `%s`.", name))
  }

  values <- parameter_block[[name]]
  data.frame(
    date = as.Date(names(values), format = "%Y%m%d"),
    value = as.numeric(unname(values)),
    stringsAsFactors = FALSE
  )
}

load_kolkata_power_daily <- function(start_year = 1982L,
                                     end_year = 2022L,
                                     latitude = 22.57,
                                     longitude = 88.36) {
  start_year <- as.integer(start_year)
  end_year <- as.integer(end_year)

  if (!is.finite(start_year) || !is.finite(end_year) || start_year > end_year) {
    stop("`start_year` and `end_year` must be finite integers with start_year <= end_year.")
  }

  download <- download_power_daily_point(
    latitude = latitude,
    longitude = longitude,
    parameters = c("WD10M", "WS10M"),
    start = sprintf("%04d0101", start_year),
    end = sprintf("%04d1231", end_year)
  )

  parameter_block <- download$payload$properties$parameter
  wd_df <- power_parameter_to_numeric(parameter_block, "WD10M")
  ws_df <- power_parameter_to_numeric(parameter_block, "WS10M")
  merged <- merge(wd_df, ws_df, by = "date", suffixes = c("_WD10M", "_WS10M"), all = TRUE)
  names(merged)[names(merged) == "value_WD10M"] <- "WD10M"
  names(merged)[names(merged) == "value_WS10M"] <- "WS10M"

  calendar_fields <- extract_calendar_fields(as.POSIXct(merged$date, tz = "UTC"), fixed_tz = "UTC")
  merged$year <- calendar_fields$year
  merged$month <- calendar_fields$month
  merged$day <- calendar_fields$day
  merged <- merged[
    merged$year >= start_year & merged$year <= end_year,
    c("date", "year", "month", "day", "WD10M", "WS10M"),
    drop = FALSE
  ]
  rownames(merged) <- NULL

  list(
    data = merged,
    url = download$url,
    geometry = download$payload$geometry,
    header = download$payload$header
  )
}

filter_power_month <- function(df, month = 8L) {
  month <- as.integer(month)
  if (!month %in% 1:12) {
    stop("`month` must be an integer in 1:12.")
  }
  if (!is.data.frame(df) || !"month" %in% names(df)) {
    stop("`df` must contain a `month` column.")
  }

  output <- df[df$month == month, , drop = FALSE]
  rownames(output) <- NULL
  output
}

run_length_table <- function(x) {
  r <- rle(as.character(x))
  if (length(r$lengths) == 0L) {
    return(data.frame(value = character(0), run_length = integer(0), stringsAsFactors = FALSE))
  }

  data.frame(
    value = r$values,
    run_length = as.integer(r$lengths),
    stringsAsFactors = FALSE
  )
}

diagnose_power_daily_data <- function(df) {
  if (!is.data.frame(df) || !all(c("date", "WD10M", "WS10M") %in% names(df))) {
    stop("`df` must contain at least `date`, `WD10M`, and `WS10M`.")
  }

  wd_runs <- run_length_table(round(df$WD10M, 6))
  ws_runs <- run_length_table(round(df$WS10M, 6))
  repeated_pairs <- paste(round(df$WD10M, 6), round(df$WS10M, 6), sep = "|")
  pair_counts <- sort(table(repeated_pairs), decreasing = TRUE)

  list(
    n_rows = nrow(df),
    date_min = if (nrow(df) == 0L) as.Date(NA) else min(df$date),
    date_max = if (nrow(df) == 0L) as.Date(NA) else max(df$date),
    n_missing_dates = sum(is.na(df$date)),
    n_duplicate_dates = sum(duplicated(df$date)),
    n_na_wd10m = sum(!is.finite(df$WD10M)),
    n_na_ws10m = sum(!is.finite(df$WS10M)),
    n_nonpositive_ws10m = sum(is.finite(df$WS10M) & df$WS10M <= 0),
    wd10m_min = suppressWarnings(min(df$WD10M, na.rm = TRUE)),
    wd10m_max = suppressWarnings(max(df$WD10M, na.rm = TRUE)),
    ws10m_min = suppressWarnings(min(df$WS10M, na.rm = TRUE)),
    ws10m_max = suppressWarnings(max(df$WS10M, na.rm = TRUE)),
    n_wd10m_outside_0_360 = sum(is.finite(df$WD10M) & (df$WD10M < 0 | df$WD10M > 360)),
    max_identical_wd_run = if (nrow(wd_runs) == 0L) 0L else max(wd_runs$run_length),
    max_identical_ws_run = if (nrow(ws_runs) == 0L) 0L else max(ws_runs$run_length),
    max_identical_pair_run = if (length(pair_counts) == 0L) 0L else as.integer(pair_counts[[1L]]),
    top_repeated_pairs = utils::head(data.frame(
      pair = names(pair_counts),
      count = as.integer(pair_counts),
      stringsAsFactors = FALSE
    ), 10L)
  )
}

thin_kolkata_every_four_days <- function(df) {
  thin_kolkata_by_day_pattern(df = df, start_day = 1L, step = 4L, max_day = 31L)
}

make_day_pattern <- function(start_day = 1L, step = 4L, max_day = 31L) {
  start_day <- as.integer(start_day)
  step <- as.integer(step)
  max_day <- as.integer(max_day)

  if (!is.finite(start_day) || start_day < 1L || start_day > max_day) {
    stop("`start_day` must be an integer in [1, max_day].")
  }
  if (!is.finite(step) || step < 1L) {
    stop("`step` must be a strictly positive integer.")
  }
  if (!is.finite(max_day) || max_day < 1L) {
    stop("`max_day` must be a strictly positive integer.")
  }

  seq.int(start_day, max_day, by = step)
}

thin_kolkata_by_day_pattern <- function(df,
                                        start_day = 1L,
                                        step = 4L,
                                        max_day = 31L) {
  if (!is.data.frame(df) || !all(c("date", "day", "WD10M", "WS10M") %in% names(df))) {
    stop("`df` must contain at least `date`, `day`, `WD10M`, and `WS10M`.")
  }

  pattern <- make_day_pattern(start_day = start_day, step = step, max_day = max_day)
  output <- df[df$day %in% pattern, , drop = FALSE]
  output <- output[order(output$date), , drop = FALSE]
  rownames(output) <- NULL
  output
}

monthly_power_quality_summary <- function(df,
                                          start_day = 1L,
                                          step = 4L,
                                          max_day = 31L) {
  if (!is.data.frame(df) || !all(c("date", "year", "month", "day", "WD10M", "WS10M") %in% names(df))) {
    stop("`df` must contain date/year/month/day/WD10M/WS10M columns.")
  }

  months <- sort(unique(as.integer(df$month)))
  rows <- lapply(months, function(month_value) {
    month_df <- df[df$month == month_value, , drop = FALSE]
    clean_df <- month_df[
      is.finite(month_df$WD10M) &
        is.finite(month_df$WS10M) &
        month_df$WS10M > 0,
      ,
      drop = FALSE
    ]
    thinned_df <- thin_kolkata_by_day_pattern(
      df = clean_df,
      start_day = start_day,
      step = step,
      max_day = max_day
    )
    diag_month <- diagnose_power_daily_data(month_df)

    data.frame(
      month = as.integer(month_value),
      month_label = month.abb[as.integer(month_value)],
      n_raw_month = nrow(month_df),
      n_clean_month = nrow(clean_df),
      n_thinned_month = nrow(thinned_df),
      n_duplicate_dates = as.integer(diag_month$n_duplicate_dates),
      n_na_wd10m = as.integer(diag_month$n_na_wd10m),
      n_na_ws10m = as.integer(diag_month$n_na_ws10m),
      n_nonpositive_ws10m = as.integer(diag_month$n_nonpositive_ws10m),
      n_wd10m_outside_0_360 = as.integer(diag_month$n_wd10m_outside_0_360),
      max_identical_wd_run = as.integer(diag_month$max_identical_wd_run),
      max_identical_ws_run = as.integer(diag_month$max_identical_ws_run),
      max_identical_pair_count = as.integer(diag_month$max_identical_pair_run),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

prepare_kolkata_power_hvmf <- function(start_year = 1982L,
                                       end_year = 2022L,
                                       month = 8L,
                                       start_day = 1L,
                                       step = 4L,
                                       latitude = 22.57,
                                       longitude = 88.36,
                                       output_dir = file.path("wind", "processed")) {
  raw_download <- load_kolkata_power_daily(
    start_year = start_year,
    end_year = end_year,
    latitude = latitude,
    longitude = longitude
  )
  raw_df_all <- raw_download$data
  diagnostics_all <- diagnose_power_daily_data(raw_df_all)
  raw_df <- filter_power_month(raw_df_all, month = month)

  valid <- is.finite(raw_df$WD10M) & is.finite(raw_df$WS10M) & raw_df$WS10M > 0
  clean_df <- raw_df[valid, , drop = FALSE]
  rownames(clean_df) <- NULL
  thinned_df <- thin_kolkata_by_day_pattern(
    df = clean_df,
    start_day = start_day,
    step = step,
    max_day = 31L
  )

  wind_df <- data.frame(
    datetime = as.POSIXct(clean_df$date, tz = "UTC"),
    ws10m = clean_df$WS10M,
    wd10m = clean_df$WD10M,
    stringsAsFactors = FALSE
  )
  wind_df_thinned <- data.frame(
    datetime = as.POSIXct(thinned_df$date, tz = "UTC"),
    ws10m = thinned_df$WS10M,
    wd10m = thinned_df$WD10M,
    stringsAsFactors = FALSE
  )

  embedded_df <- build_hvmf_wind_set(
    df = wind_df_thinned,
    speed_col = "ws10m",
    direction_col = "wd10m",
    height_m = 10,
    fixed_tz = "UTC"
  )

  embedded_df$dataset_id <- sprintf(
    "kolkata_%s_%d_%d_start%d_step%d",
    tolower(month.abb[as.integer(month)]),
    as.integer(start_year),
    as.integer(end_year),
    as.integer(start_day),
    as.integer(step)
  )
  embedded_df$source <- "NASA POWER CERES/MERRA2 daily point"
  embedded_df$station <- "Kolkata"
  embedded_df$latitude <- latitude
  embedded_df$longitude <- longitude

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  all_daily_csv <- file.path(
    output_dir,
    sprintf("kolkata_daily_%d_%d_all_months.csv", as.integer(start_year), as.integer(end_year))
  )
  raw_csv <- file.path(
    output_dir,
    sprintf("kolkata_%s_%d_%d_raw_daily.csv", tolower(month.abb[as.integer(month)]), as.integer(start_year), as.integer(end_year))
  )
  thinned_csv <- file.path(
    output_dir,
    sprintf(
      "kolkata_%s_%d_%d_start%d_step%d.csv",
      tolower(month.abb[as.integer(month)]),
      as.integer(start_year),
      as.integer(end_year),
      as.integer(start_day),
      as.integer(step)
    )
  )
  hvmf_csv <- file.path(
    output_dir,
    sprintf(
      "kolkata_%s_%d_%d_start%d_step%d_hvmf.csv",
      tolower(month.abb[as.integer(month)]),
      as.integer(start_year),
      as.integer(end_year),
      as.integer(start_day),
      as.integer(step)
    )
  )

  utils::write.csv(raw_df_all, all_daily_csv, row.names = FALSE)
  utils::write.csv(clean_df, raw_csv, row.names = FALSE)
  utils::write.csv(thinned_df, thinned_csv, row.names = FALSE)
  write_hvmf_csv(embedded_df, hvmf_csv, fixed_tz = "UTC")

  list(
    all_daily_data = raw_df_all,
    raw_data = clean_df,
    thinned_data = thinned_df,
    embedded_data = embedded_df,
    diagnostics_all = diagnostics_all,
    all_daily_csv = all_daily_csv,
    raw_csv = raw_csv,
    thinned_csv = thinned_csv,
    hvmf_csv = hvmf_csv,
    source_url = raw_download$url,
    has_intensity = all(is.finite(raw_df_all$WS10M) & raw_df_all$WS10M > 0),
    has_direction = all(is.finite(raw_df_all$WD10M)),
    month = as.integer(month),
    month_label = month.abb[as.integer(month)],
    start_day = as.integer(start_day),
    step = as.integer(step),
    n_raw_all = nrow(raw_df_all),
    n_raw_month = nrow(raw_df),
    n_clean_month = nrow(clean_df),
    n_every4days = nrow(thinned_df)
  )
}

if (sys.nframe() == 0L) {
  result <- prepare_kolkata_power_hvmf()
  cat("source_url=", result$source_url, "\n", sep = "")
  cat("has_intensity=", result$has_intensity, "\n", sep = "")
  cat("has_direction=", result$has_direction, "\n", sep = "")
  cat("n_raw_all=", result$n_raw_all, "\n", sep = "")
  cat("n_raw_month=", result$n_raw_month, "\n", sep = "")
  cat("n_clean_month=", result$n_clean_month, "\n", sep = "")
  cat("n_every4days=", result$n_every4days, "\n", sep = "")
  cat("n_duplicate_dates_all=", result$diagnostics_all$n_duplicate_dates, "\n", sep = "")
  cat("n_na_wd10m_all=", result$diagnostics_all$n_na_wd10m, "\n", sep = "")
  cat("n_na_ws10m_all=", result$diagnostics_all$n_na_ws10m, "\n", sep = "")
  cat("n_nonpositive_ws10m_all=", result$diagnostics_all$n_nonpositive_ws10m, "\n", sep = "")
  cat("n_wd10m_outside_0_360_all=", result$diagnostics_all$n_wd10m_outside_0_360, "\n", sep = "")
  cat("max_identical_wd_run_all=", result$diagnostics_all$max_identical_wd_run, "\n", sep = "")
  cat("max_identical_ws_run_all=", result$diagnostics_all$max_identical_ws_run, "\n", sep = "")
  cat("max_identical_pair_count_all=", result$diagnostics_all$max_identical_pair_run, "\n", sep = "")
  cat("all_daily_csv=", result$all_daily_csv, "\n", sep = "")
  cat("raw_csv=", result$raw_csv, "\n", sep = "")
  cat("thinned_csv=", result$thinned_csv, "\n", sep = "")
  cat("hvmf_csv=", result$hvmf_csv, "\n", sep = "")
}