check_ncdf4_available <- function() {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required. Install it with install.packages('ncdf4').", call. = FALSE)
  }
}

parse_netcdf_time_units <- function(units_string) {
  if (!is.character(units_string) || length(units_string) != 1L || !nzchar(units_string)) {
    stop("NetCDF time units must be a non-empty character scalar.")
  }

  match <- regexec("^\\s*([A-Za-z]+)\\s+since\\s+(.+)\\s*$", units_string)
  parts <- regmatches(units_string, match)[[1L]]

  if (length(parts) != 3L) {
    stop(sprintf("Unsupported NetCDF time units string: %s", units_string))
  }

  unit_name <- tolower(parts[[2L]])
  origin_text <- parts[[3L]]
  seconds_per_unit <- switch(
    unit_name,
    second = 1,
    seconds = 1,
    minute = 60,
    minutes = 60,
    hour = 3600,
    hours = 3600,
    day = 86400,
    days = 86400,
    stop(sprintf("Unsupported NetCDF time unit: %s", unit_name))
  )

  list(
    unit_name = unit_name,
    origin_text = origin_text,
    seconds_per_unit = seconds_per_unit
  )
}

extract_calendar_fields <- function(datetime, fixed_tz = "UTC") {
  datetime_lt <- as.POSIXlt(datetime, tz = fixed_tz)

  data.frame(
    year = datetime_lt$year + 1900L,
    month = datetime_lt$mon + 1L,
    day = datetime_lt$mday,
    hour = datetime_lt$hour,
    minute = datetime_lt$min,
    stringsAsFactors = FALSE
  )
}

format_datetime_for_csv <- function(datetime, fixed_tz = "UTC") {
  format(datetime, "%Y-%m-%d %H:%M:%S", tz = fixed_tz, usetz = FALSE)
}

format_datetime_range <- function(datetime, fixed_tz = "UTC") {
  sprintf(
    "%s to %s",
    format(min(datetime), "%Y-%m-%d %H:%M:%S", tz = fixed_tz, usetz = FALSE),
    format(max(datetime), "%Y-%m-%d %H:%M:%S", tz = fixed_tz, usetz = FALSE)
  )
}

read_risoe_nc_metadata <- function(path, fixed_tz = "UTC") {
  check_ncdf4_available()

  if (!file.exists(path)) {
    stop(sprintf("Input NetCDF not found: %s", path), call. = FALSE)
  }

  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  required_vars <- c("ws77", "wd77", "ws125", "wd125")
  missing_vars <- setdiff(required_vars, names(nc$var))
  if (length(missing_vars) > 0L) {
    stop(sprintf("Missing required NetCDF variables: %s", paste(missing_vars, collapse = ", ")))
  }

  if (is.null(nc$dim[["time"]])) {
    stop("NetCDF file does not contain a `time` dimension.")
  }

  time_units <- nc$dim[["time"]]$units
  if (is.null(time_units) || !nzchar(time_units)) {
    stop("NetCDF `time` dimension does not provide units.")
  }

  time_calendar <- ncdf4::ncatt_get(nc, "time", "calendar")$value
  if (is.null(time_calendar) || identical(time_calendar, 0) || !nzchar(time_calendar)) {
    time_calendar <- NA_character_
  }

  time_info <- parse_netcdf_time_units(time_units)
  heights_m <- c(
    ws77 = as.numeric(ncdf4::ncatt_get(nc, "ws77", "height")$value),
    wd77 = as.numeric(ncdf4::ncatt_get(nc, "wd77", "height")$value),
    ws125 = as.numeric(ncdf4::ncatt_get(nc, "ws125", "height")$value),
    wd125 = as.numeric(ncdf4::ncatt_get(nc, "wd125", "height")$value)
  )

  list(
    path = path,
    fixed_tz = fixed_tz,
    required_vars = required_vars,
    time_units = time_units,
    time_calendar = time_calendar,
    time_origin_text = time_info$origin_text,
    time_unit_name = time_info$unit_name,
    time_seconds_per_unit = time_info$seconds_per_unit,
    time_length = nc$dim[["time"]]$len,
    heights_m = heights_m
  )
}

load_risoe_concurrent <- function(path, fixed_tz = "UTC") {
  metadata <- read_risoe_nc_metadata(path, fixed_tz = fixed_tz)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  time_values <- as.numeric(nc$dim[["time"]]$vals)
  origin_time <- as.POSIXct(metadata$time_origin_text, tz = fixed_tz)

  if (is.na(origin_time)) {
    stop(sprintf("Could not parse NetCDF time origin: %s", metadata$time_origin_text))
  }

  datetime <- origin_time + time_values * metadata$time_seconds_per_unit
  attr(datetime, "tzone") <- fixed_tz

  ws77 <- as.numeric(ncdf4::ncvar_get(nc, "ws77"))
  wd77 <- as.numeric(ncdf4::ncvar_get(nc, "wd77"))
  ws125 <- as.numeric(ncdf4::ncvar_get(nc, "ws125"))
  wd125 <- as.numeric(ncdf4::ncvar_get(nc, "wd125"))

  n_time <- length(datetime)
  series_lengths <- c(
    ws77 = length(ws77),
    wd77 = length(wd77),
    ws125 = length(ws125),
    wd125 = length(wd125)
  )

  if (any(series_lengths != n_time)) {
    bad <- names(series_lengths)[series_lengths != n_time]
    stop(sprintf("Length mismatch between `time` and variables: %s", paste(bad, collapse = ", ")))
  }

  calendar_fields <- extract_calendar_fields(datetime, fixed_tz = fixed_tz)

  data.frame(
    datetime = datetime,
    year = calendar_fields$year,
    month = calendar_fields$month,
    day = calendar_fields$day,
    hour = calendar_fields$hour,
    minute = calendar_fields$minute,
    ws77 = ws77,
    wd77 = wd77,
    ws125 = ws125,
    wd125 = wd125,
    stringsAsFactors = FALSE
  )
}

select_noon_nov_dec <- function(df, tie_break = "earliest", fixed_tz = "UTC") {
  tie_break <- match.arg(tie_break, c("earliest", "latest"))

  if (!is.data.frame(df) || !"datetime" %in% names(df)) {
    stop("`df` must be a data.frame containing a `datetime` column.")
  }
  if (!inherits(df$datetime, "POSIXct")) {
    stop("`df$datetime` must be a POSIXct vector.")
  }

  calendar_fields <- extract_calendar_fields(df$datetime, fixed_tz = fixed_tz)
  keep <- calendar_fields$month %in% c(11L, 12L)
  df_nov_dec <- df[keep, , drop = FALSE]

  if (nrow(df_nov_dec) == 0L) {
    stop("No November-December records are available after filtering.")
  }

  date_key <- as.Date(df_nov_dec$datetime, tz = fixed_tz)
  noon_datetime <- as.POSIXct(paste(date_key, "12:00:00"), tz = fixed_tz)
  noon_distance_sec <- abs(as.numeric(difftime(df_nov_dec$datetime, noon_datetime, units = "secs")))
  datetime_numeric <- as.numeric(df_nov_dec$datetime)

  order_index <- switch(
    tie_break,
    earliest = order(date_key, noon_distance_sec, datetime_numeric),
    latest = order(date_key, noon_distance_sec, -datetime_numeric)
  )

  ordered_df <- df_nov_dec[order_index, , drop = FALSE]
  ordered_date_key <- as.Date(ordered_df$datetime, tz = fixed_tz)
  selected_df <- ordered_df[!duplicated(ordered_date_key), , drop = FALSE]
  rownames(selected_df) <- NULL

  selected_df
}

build_hvmf_wind_set <- function(df, speed_col, direction_col, height_m, fixed_tz = "UTC") {
  required_cols <- c("datetime", speed_col, direction_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!inherits(df$datetime, "POSIXct")) {
    stop("`df$datetime` must be POSIXct.")
  }

  speed <- as.numeric(df[[speed_col]])
  direction <- as.numeric(df[[direction_col]])
  valid <- is.finite(speed) & is.finite(direction) & speed > 0

  filtered_df <- df[valid, "datetime", drop = FALSE]
  speed <- speed[valid]
  direction <- direction[valid]

  if (nrow(filtered_df) == 0L) {
    stop(sprintf("No valid rows remain for `%s` and `%s`.", speed_col, direction_col))
  }

  direction_deg <- direction %% 360
  speed_mean_height <- mean(speed, na.rm = TRUE)
  if (!is.finite(speed_mean_height) || speed_mean_height <= 0) {
    stop(sprintf("Invalid mean speed for `%s`.", speed_col))
  }

  speed_scaled <- speed / speed_mean_height
  angle_rad <- direction_deg * pi / 180
  sinh_scaled <- sinh(speed_scaled)
  x0 <- cosh(speed_scaled)
  x1 <- sinh_scaled * cos(angle_rad)
  x2 <- sinh_scaled * sin(angle_rad)
  minkowski_norm <- -x0^2 + x1^2 + x2^2

  calendar_fields <- extract_calendar_fields(filtered_df$datetime, fixed_tz = fixed_tz)

  result <- data.frame(
    datetime = filtered_df$datetime,
    year = calendar_fields$year,
    month = calendar_fields$month,
    day = calendar_fields$day,
    hour = calendar_fields$hour,
    minute = calendar_fields$minute,
    height_m = rep(as.numeric(height_m), length(speed)),
    speed = speed,
    direction_deg = direction_deg,
    speed_mean_height = rep(speed_mean_height, length(speed)),
    speed_scaled = speed_scaled,
    angle_rad = angle_rad,
    x0 = x0,
    x1 = x1,
    x2 = x2,
    minkowski_norm = minkowski_norm,
    stringsAsFactors = FALSE
  )

  attr(result, "dropped_days") <- nrow(df) - nrow(result)
  attr(result, "fixed_tz") <- fixed_tz
  attr(result, "speed_col") <- speed_col
  attr(result, "direction_col") <- direction_col

  result
}

plot_risoe_hvmf_scatter <- function(df, output_png) {
  png(filename = output_png, width = 900, height = 600, res = 120)
  on.exit(dev.off(), add = TRUE)

  plot(
    x = df$direction_deg,
    y = df$speed_scaled,
    xlim = c(0, 360),
    xaxs = "i",
    pch = 16,
    cex = 0.7,
    col = grDevices::rgb(0.1, 0.3, 0.6, 0.65),
    xlab = "direction_deg",
    ylab = "speed_scaled",
    main = sprintf("Risoe HvMF preprocessing scatter (%sm)", unique(df$height_m))
  )
  grid(col = "grey85")
}

write_hvmf_csv <- function(df, output_csv, fixed_tz = "UTC") {
  output_df <- df
  output_df$datetime <- format_datetime_for_csv(output_df$datetime, fixed_tz = fixed_tz)
  write.csv(output_df, file = output_csv, row.names = FALSE)
}

print_hvmf_set_summary <- function(label, df, fixed_tz = "UTC") {
  cat(label, ":\n", sep = "")
  cat("  n = ", nrow(df), "\n", sep = "")
  cat("  date range = ", format_datetime_range(df$datetime, fixed_tz = fixed_tz), "\n", sep = "")
  cat("  mean speed = ", sprintf("%.10f", df$speed_mean_height[[1L]]), "\n", sep = "")
  cat(
    "  min/max speed_scaled = ",
    sprintf("%.10f / %.10f", min(df$speed_scaled), max(df$speed_scaled)),
    "\n",
    sep = ""
  )
  cat(
    "  max |minkowski_norm + 1| = ",
    sprintf("%.10e", max(abs(df$minkowski_norm + 1))),
    "\n",
    sep = ""
  )
  cat("  days dropped after height-specific filtering = ", attr(df, "dropped_days"), "\n\n", sep = "")
}

write_risoe_hvmf_sets <- function(input_nc = "real_data/wind/risoe_m_concurent.nc",
                                  output_dir = "real_data/wind/processed",
                                  tie_break = "earliest",
                                  make_plots = TRUE,
                                  fixed_tz = "UTC") {
  metadata <- read_risoe_nc_metadata(input_nc, fixed_tz = fixed_tz)
  concurrent_df <- load_risoe_concurrent(input_nc, fixed_tz = fixed_tz)
  selected_df <- select_noon_nov_dec(concurrent_df, tie_break = tie_break, fixed_tz = fixed_tz)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  df_A <- build_hvmf_wind_set(
    df = selected_df,
    speed_col = "ws77",
    direction_col = "wd77",
    height_m = 77,
    fixed_tz = fixed_tz
  )
  df_B <- build_hvmf_wind_set(
    df = selected_df,
    speed_col = "ws125",
    direction_col = "wd125",
    height_m = 125,
    fixed_tz = fixed_tz
  )

  output_csv_A <- file.path(output_dir, "risoe_modern_set_A_77m_hvmf.csv")
  output_csv_B <- file.path(output_dir, "risoe_modern_set_B_125m_hvmf.csv")
  output_png_A <- file.path(output_dir, "risoe_modern_set_A_77m_scatter.png")
  output_png_B <- file.path(output_dir, "risoe_modern_set_B_125m_scatter.png")

  write_hvmf_csv(df_A, output_csv_A, fixed_tz = fixed_tz)
  write_hvmf_csv(df_B, output_csv_B, fixed_tz = fixed_tz)

  if (isTRUE(make_plots)) {
    plot_risoe_hvmf_scatter(df_A, output_png_A)
    plot_risoe_hvmf_scatter(df_B, output_png_B)
  }

  cat("Input NetCDF metadata:\n")
  cat("  file = ", metadata$path, "\n", sep = "")
  cat("  time units = ", metadata$time_units, "\n", sep = "")
  cat("  calendar = ", ifelse(is.na(metadata$time_calendar), "NA", metadata$time_calendar), "\n", sep = "")
  cat("  fixed_tz = ", metadata$fixed_tz, " (timestamps used exactly as stored; no local-time conversion)\n\n", sep = "")

  cat("Selected daily timestamps near noon (HH:MM):\n")
  print(table(format(selected_df$datetime, "%H:%M", tz = fixed_tz)))
  cat("\n")

  print_hvmf_set_summary("Modern set A, 77 m", df_A, fixed_tz = fixed_tz)
  print_hvmf_set_summary("Modern set B, 125 m", df_B, fixed_tz = fixed_tz)

  stopifnot(file.exists(output_csv_A))
  stopifnot(file.exists(output_csv_B))
  stopifnot(all(abs(df_A$minkowski_norm + 1) < 1e-8))
  stopifnot(all(abs(df_B$minkowski_norm + 1) < 1e-8))
  stopifnot(all(df_A$month %in% c(11L, 12L)))
  stopifnot(all(df_B$month %in% c(11L, 12L)))
  stopifnot(!anyDuplicated(as.Date(df_A$datetime, tz = fixed_tz)))
  stopifnot(!anyDuplicated(as.Date(df_B$datetime, tz = fixed_tz)))

  invisible(list(
    metadata = metadata,
    selected_df = selected_df,
    df_A = df_A,
    df_B = df_B,
    output_csv_A = output_csv_A,
    output_csv_B = output_csv_B,
    output_png_A = output_png_A,
    output_png_B = output_png_B
  ))
}

if (sys.nframe() == 0L) {
  write_risoe_hvmf_sets()
}
