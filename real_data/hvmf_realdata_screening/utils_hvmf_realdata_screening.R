resolve_hvmf_realdata_screening_path <- function(...) {
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

  file.path(...)
}

bootstrap_script_path_hvmf_realdata <- resolve_hvmf_realdata_screening_path("bootstrap", "multiplier_bootstrap.R")
if (!exists("multiplier_bootstrap_gof", mode = "function")) {
  source(bootstrap_script_path_hvmf_realdata)
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) {
    rhs
  } else {
    lhs
  }
}

capture_warnings_screening <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = unique(warnings))
}

format_timestamp_screening <- function(x) {
  format(as.POSIXct(x, tz = "UTC"), "%Y-%m-%d %H:%M:%S", tz = "UTC")
}

normalize_column_names <- function(x) {
  output <- trimws(as.character(x))
  output <- tolower(output)
  output <- gsub("[^a-z0-9]+", "_", output)
  output <- gsub("^_+|_+$", "", output)
  output <- gsub("_+", "_", output)
  output
}

trim_character_columns <- function(df) {
  for (name in names(df)) {
    if (is.character(df[[name]])) {
      df[[name]] <- trimws(df[[name]])
    }
  }
  df
}

screening_directories <- function(base_dir = "hvmf_realdata_screening") {
  list(
    base = base_dir,
    raw = file.path(base_dir, "raw"),
    processed = file.path(base_dir, "processed"),
    results = file.path(base_dir, "results"),
    logs = file.path(base_dir, "logs")
  )
}

ensure_screening_directories <- function(base_dir = "hvmf_realdata_screening") {
  dirs <- screening_directories(base_dir = base_dir)
  for (path in unname(unlist(dirs))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  dirs
}

download_file_if_needed <- function(url, destfile, overwrite = FALSE) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  if (!overwrite && file.exists(destfile)) {
    return(invisible(destfile))
  }

  utils::download.file(url = url, destfile = destfile, mode = "wb", quiet = TRUE)
  invisible(destfile)
}

download_coops_currents <- function(station,
                                    begin_date,
                                    end_date,
                                    bin = NULL,
                                    units = "metric",
                                    dest_dir = file.path("hvmf_realdata_screening", "raw"),
                                    overwrite = FALSE) {
  query <- c(
    product = "currents",
    station = station,
    begin_date = begin_date,
    end_date = end_date,
    units = units,
    time_zone = "gmt",
    format = "csv"
  )
  if (!is.null(bin)) {
    query <- c(query, bin = as.character(bin))
  }

  query_string <- paste(
    paste(names(query), utils::URLencode(unname(query), reserved = TRUE), sep = "="),
    collapse = "&"
  )
  url <- paste0("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?", query_string)
  filename <- sprintf(
    "coops_currents_%s_%s_%s%s.csv",
    station,
    begin_date,
    end_date,
    if (is.null(bin)) "" else paste0("_bin", bin)
  )
  destfile <- file.path(dest_dir, filename)

  download_result <- tryCatch(
    {
      download_file_if_needed(url = url, destfile = destfile, overwrite = overwrite)
      list(ok = TRUE, error_message = NA_character_)
    },
    error = function(e) {
      list(ok = FALSE, error_message = conditionMessage(e))
    }
  )

  if (!download_result$ok) {
    return(list(
      ok = FALSE,
      url = url,
      file = destfile,
      error_message = download_result$error_message
    ))
  }

  lines <- readLines(destfile, warn = FALSE, n = 10L)
  error_line <- lines[grepl("(^|\\s)(Error:|Wrong Bin Number)", lines)]
  if (length(error_line) > 0L) {
    return(list(
      ok = FALSE,
      url = url,
      file = destfile,
      error_message = paste(trimws(error_line), collapse = " ")
    ))
  }

  raw_df <- utils::read.csv(destfile, stringsAsFactors = FALSE, check.names = FALSE)
  names(raw_df) <- normalize_column_names(names(raw_df))
  raw_df <- trim_character_columns(raw_df)

  list(
    ok = TRUE,
    url = url,
    file = destfile,
    data = raw_df,
    error_message = NA_character_
  )
}

download_ndbc_stdmet <- function(station,
                                 year,
                                 dest_dir = file.path("hvmf_realdata_screening", "raw"),
                                 overwrite = FALSE) {
  filename <- sprintf("%sh%s.txt.gz", station, year)
  url <- sprintf("https://www.ndbc.noaa.gov/data/historical/stdmet/%s", filename)
  destfile <- file.path(dest_dir, filename)

  download_result <- tryCatch(
    {
      download_file_if_needed(url = url, destfile = destfile, overwrite = overwrite)
      list(ok = TRUE, error_message = NA_character_)
    },
    error = function(e) {
      list(ok = FALSE, error_message = conditionMessage(e))
    }
  )

  list(
    ok = download_result$ok,
    url = url,
    file = destfile,
    error_message = download_result$error_message
  )
}

parse_coops_datetime <- function(x) {
  as.POSIXct(x, format = "%Y-%m-%d %H:%M", tz = "UTC")
}

parse_ndbc_stdmet <- function(file) {
  con <- gzfile(file, open = "rt")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1L, warn = FALSE)
  units_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(header_line) == 0L) {
    stop(sprintf("NDBC file is empty: %s", file))
  }
  if (length(units_line) == 0L) {
    stop(sprintf("NDBC file is missing the units line: %s", file))
  }

  tokens <- strsplit(trimws(sub("^#", "", header_line)), "[[:space:]]+")[[1]]
  if (length(tokens) < 5L) {
    stop(sprintf("Could not parse NDBC header from %s", file))
  }

  tokens[1:5] <- c("year", "month", "day", "hour", "minute")
  col_names <- normalize_column_names(tokens)

  raw_df <- utils::read.table(
    con,
    header = FALSE,
    col.names = col_names,
    stringsAsFactors = FALSE,
    fill = TRUE
  )

  raw_df <- trim_character_columns(raw_df)
  raw_df$datetime <- as.POSIXct(
    sprintf(
      "%04d-%02d-%02d %02d:%02d",
      as.integer(raw_df$year),
      as.integer(raw_df$month),
      as.integer(raw_df$day),
      as.integer(raw_df$hour),
      as.integer(raw_df$minute)
    ),
    tz = "UTC"
  )

  raw_df
}

choose_first_existing <- function(names_vec, candidates = character(0), patterns = character(0)) {
  for (candidate in candidates) {
    if (candidate %in% names_vec) {
      return(candidate)
    }
  }

  for (pattern in patterns) {
    hits <- names_vec[grepl(pattern, names_vec)]
    if (length(hits) > 0L) {
      return(hits[[1L]])
    }
  }

  NULL
}

clean_speed_direction <- function(df,
                                  speed_col,
                                  direction_col,
                                  sentinel_values = numeric(0),
                                  direction_upper_exclusive = 360) {
  speed <- suppressWarnings(as.numeric(df[[speed_col]]))
  direction_deg <- suppressWarnings(as.numeric(df[[direction_col]]))

  keep <- is.finite(speed) & is.finite(direction_deg)
  keep <- keep & (speed > 0)
  keep <- keep & (direction_deg >= 0) & (direction_deg < direction_upper_exclusive)

  if (length(sentinel_values) > 0L) {
    sentinel_speed <- Reduce(`|`, lapply(sentinel_values, function(value) abs(speed - value) < 1e-12))
    sentinel_direction <- Reduce(`|`, lapply(sentinel_values, function(value) abs(direction_deg - value) < 1e-12))
    keep <- keep & !sentinel_speed & !sentinel_direction
  }

  output <- df[keep, , drop = FALSE]
  output$raw_speed <- speed[keep]
  output$raw_direction_deg <- direction_deg[keep]
  output
}

components_from_speed_direction <- function(speed, direction_deg) {
  theta <- as.numeric(direction_deg) * pi / 180
  speed <- as.numeric(speed)

  data.frame(
    v1 = speed * sin(theta),
    v2 = speed * cos(theta)
  )
}

split_phases_by_principal_axis <- function(v1, v2) {
  velocity_matrix <- cbind(as.numeric(v1), as.numeric(v2))
  if (nrow(velocity_matrix) < 2L) {
    stop("At least two observations are required to split phases.")
  }
  if (!all(is.finite(velocity_matrix))) {
    stop("Velocity components must be finite for the phase split.")
  }

  pca <- stats::prcomp(velocity_matrix, center = TRUE, scale. = FALSE)
  axis <- as.numeric(pca$rotation[, 1L])
  axis_norm <- sqrt(sum(axis^2))
  if (!is.finite(axis_norm) || axis_norm <= 0) {
    stop("Could not determine a valid principal axis.")
  }

  axis <- axis / axis_norm
  projection <- drop(velocity_matrix %*% axis)

  list(
    axis = axis,
    projection = projection,
    positive = projection > 0,
    negative = projection < 0,
    method = "pca"
  )
}

embed_components_to_h2 <- function(v1, v2, scale) {
  scale <- as.numeric(scale)
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("`scale` must be a strictly positive finite scalar.")
  }

  z1 <- as.numeric(v1) / scale
  z2 <- as.numeric(v2) / scale
  x0 <- sqrt(1 + z1^2 + z2^2)
  x1 <- z1
  x2 <- z2
  minkowski_norm <- -x0^2 + x1^2 + x2^2

  data.frame(
    scale = rep_len(scale, length(z1)),
    z1 = z1,
    z2 = z2,
    x0 = x0,
    x1 = x1,
    x2 = x2,
    minkowski_norm = minkowski_norm
  )
}

thin_to_max_n <- function(df, max_n = 500L, time_col = "datetime") {
  n <- nrow(df)
  if (n <= max_n) {
    return(df)
  }

  if (!time_col %in% names(df)) {
    stop(sprintf("Time column `%s` is missing.", time_col))
  }

  ord <- order(df[[time_col]])
  df_ordered <- df[ord, , drop = FALSE]
  idx <- unique(round(seq(1, nrow(df_ordered), length.out = min(max_n, nrow(df_ordered)))))
  df_ordered[idx, , drop = FALSE]
}

has_daily_coverage <- function(datetime_vec) {
  datetime_vec <- as.POSIXct(datetime_vec, tz = "UTC")
  if (length(datetime_vec) == 0L || all(is.na(datetime_vec))) {
    return(FALSE)
  }

  dates <- sort(unique(as.Date(datetime_vec)))
  if (length(dates) < 2L) {
    return(FALSE)
  }

  gaps <- diff(as.integer(dates))
  length(gaps) > 0L && all(gaps <= 1L)
}

thin_by_day_stride_if_daily_coverage <- function(df,
                                                 day_stride = 5L,
                                                 time_col = "datetime") {
  if (!time_col %in% names(df)) {
    stop(sprintf("Time column `%s` is missing.", time_col))
  }

  day_stride <- as.integer(day_stride)
  if (!is.finite(day_stride) || day_stride <= 1L) {
    return(list(
      data = df,
      applied = FALSE,
      unique_days_before = length(unique(as.Date(df[[time_col]]))),
      unique_days_after = length(unique(as.Date(df[[time_col]]))),
      day_stride = day_stride
    ))
  }

  ordered_df <- df[order(df[[time_col]]), , drop = FALSE]
  datetime_vec <- as.POSIXct(ordered_df[[time_col]], tz = "UTC")

  if (!has_daily_coverage(datetime_vec)) {
    return(list(
      data = ordered_df,
      applied = FALSE,
      unique_days_before = length(unique(as.Date(datetime_vec))),
      unique_days_after = length(unique(as.Date(datetime_vec))),
      day_stride = day_stride
    ))
  }

  dates <- as.Date(datetime_vec)
  unique_dates <- sort(unique(dates))
  keep_dates <- unique_dates[seq(1L, length(unique_dates), by = day_stride)]
  keep <- dates %in% keep_dates

  list(
    data = ordered_df[keep, , drop = FALSE],
    applied = TRUE,
    unique_days_before = length(unique_dates),
    unique_days_after = length(keep_dates),
    day_stride = day_stride
  )
}

verify_h2_matrix_screening <- function(X, tol = 1e-8) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  if (!all(dim(X) == c(nrow(X), 3L))) {
    stop("`X` must be an n x 3 matrix.")
  }
  if (!all(is.finite(X))) {
    stop("`X` contains non-finite values.")
  }
  if (any(X[, 1L] <= 0)) {
    stop("`X` contains rows with x0 <= 0.")
  }

  minkowski_norm <- -X[, 1L]^2 + X[, 2L]^2 + X[, 3L]^2
  max_error <- max(abs(minkowski_norm + 1))
  if (!is.finite(max_error) || max_error >= tol) {
    stop(sprintf("The embedding is not on H^2 within tolerance %.1e.", tol))
  }

  list(
    X = X,
    minkowski_norm = minkowski_norm,
    max_error = max_error
  )
}

compute_hvmf_auxiliary_parameters <- function(mu) {
  mu <- as.numeric(mu)
  if (length(mu) != 3L) {
    stop("`mu` must have length 3.")
  }

  sinh_chi_hat <- sqrt(mu[2]^2 + mu[3]^2)
  chi_hat <- asinh(sinh_chi_hat)
  theta_hat <- atan2(mu[3], mu[2])
  theta_deg_hat <- (theta_hat * 180 / pi) %% 360

  list(
    sinh_chi_hat = sinh_chi_hat,
    chi_hat = chi_hat,
    theta_hat = theta_hat,
    theta_deg_hat = theta_deg_hat
  )
}

run_hvmf_simple_plugin_cvm <- function(X,
                                       B = 5000L,
                                       n_cores = 3L,
                                       seed = NULL,
                                       control = list(
                                         hvmf_profile_method = "tabulated",
                                         hvmf_profile_n_y = 4097L
                                       )) {
  X_checked <- verify_h2_matrix_screening(X)$X

  fit_result <- capture_warnings_screening(
    hvmf_mle_h2(X_checked)
  )
  fit <- fit_result$value

  if (!is.finite(fit$kappa) || fit$kappa <= 0) {
    stop("The HvMF MLE returned a nonpositive or non-finite kappa.")
  }

  mu_error <- max(abs(-fit$mu[1]^2 + fit$mu[2]^2 + fit$mu[3]^2 + 1))
  if (!is.finite(mu_error) || mu_error >= 1e-8) {
    stop("The fitted mu is not on H^2 within tolerance 1e-8.")
  }

  theta0 <- list(mu = fit$mu, kappa = fit$kappa)
  spec <- make_hvmf_spec(unknown_param = "both")
  bootstrap_result <- capture_warnings_screening(
    multiplier_bootstrap_gof(
      data = X_checked,
      spec = spec,
      null = list(type = "simple", theta = theta0),
      statistics = "cvm",
      B = as.integer(B),
      n_cores = as.integer(n_cores),
      seed = seed,
      keep = list(
        observed_process = TRUE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = FALSE
      ),
      control = control
    )
  )

  list(
    fit = fit,
    theta0 = theta0,
    result = bootstrap_result$value,
    warnings = unique(c(fit_result$warnings, bootstrap_result$warnings))
  )
}

summarize_hvmf_result <- function(dataset_id,
                                  source,
                                  station = NA_character_,
                                  window_start = NA_character_,
                                  window_end = NA_character_,
                                  phase = NA_character_,
                                  n_raw = NA_integer_,
                                  n_clean = NA_integer_,
                                  n_final = NA_integer_,
                                  scale = NA_real_,
                                  fit = NULL,
                                  bootstrap_result = NULL,
                                  B = NA_integer_,
                                  n_cores = NA_integer_,
                                  elapsed_seconds = NA_real_,
                                  processed_csv = NA_character_,
                                  result_rds = NA_character_,
                                  log_file = NA_character_,
                                  status = "ok",
                                  error_message = NA_character_) {
  if (!identical(status, "ok") || is.null(fit) || is.null(bootstrap_result)) {
    return(data.frame(
      dataset_id = dataset_id,
      source = source,
      station = station,
      window_start = window_start,
      window_end = window_end,
      phase = phase,
      n_raw = n_raw,
      n_clean = n_clean,
      n_final = n_final,
      scale = scale,
      kappa_hat = NA_real_,
      mu_hat_x0 = NA_real_,
      mu_hat_x1 = NA_real_,
      mu_hat_x2 = NA_real_,
      sinh_chi_hat = NA_real_,
      theta_deg_hat = NA_real_,
      statistic_CvM = NA_real_,
      p_value = NA_real_,
      B = as.integer(B),
      n_cores = as.integer(n_cores),
      elapsed_seconds = elapsed_seconds,
      processed_csv = processed_csv,
      result_rds = result_rds,
      log_file = log_file,
      status = status,
      error_message = error_message,
      stringsAsFactors = FALSE
    ))
  }

  aux <- compute_hvmf_auxiliary_parameters(fit$mu)

  data.frame(
    dataset_id = dataset_id,
    source = source,
    station = station,
    window_start = window_start,
    window_end = window_end,
    phase = phase,
    n_raw = n_raw,
    n_clean = n_clean,
    n_final = n_final,
    scale = scale,
    kappa_hat = fit$kappa,
    mu_hat_x0 = fit$mu[1],
    mu_hat_x1 = fit$mu[2],
    mu_hat_x2 = fit$mu[3],
    sinh_chi_hat = aux$sinh_chi_hat,
    theta_deg_hat = aux$theta_deg_hat,
    statistic_CvM = bootstrap_result$observed$cvm$statistic,
    p_value = bootstrap_result$inference$cvm$p_value,
    B = as.integer(B),
    n_cores = as.integer(n_cores),
    elapsed_seconds = elapsed_seconds,
    processed_csv = processed_csv,
    result_rds = result_rds,
    log_file = log_file,
    status = status,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

extract_currents_components <- function(df) {
  names_df <- names(df)

  along_col <- choose_first_existing(
    names_df,
    candidates = c("along_channel_velocity", "along_channel", "along_channel_vel"),
    patterns = c("along.*channel.*vel", "along.*channel")
  )
  cross_col <- choose_first_existing(
    names_df,
    candidates = c("cross_channel_velocity", "cross_channel", "cross_channel_vel"),
    patterns = c("cross.*channel.*vel", "cross.*channel")
  )

  if (!is.null(along_col) && !is.null(cross_col)) {
    v1 <- suppressWarnings(as.numeric(df[[along_col]]))
    v2 <- suppressWarnings(as.numeric(df[[cross_col]]))

    return(list(
      method = "along_cross",
      v1 = v1,
      v2 = v2,
      raw_speed = sqrt(v1^2 + v2^2),
      raw_direction_deg = rep_len(NA_real_, length(v1)),
      along_col = along_col,
      cross_col = cross_col
    ))
  }

  speed_col <- choose_first_existing(
    names_df,
    candidates = c("speed", "velocity", "spd"),
    patterns = c("^speed$", "speed", "velocity")
  )
  direction_col <- choose_first_existing(
    names_df,
    candidates = c("direction", "dir", "direction_deg"),
    patterns = c("^direction$", "^dir$", "direction")
  )

  if (is.null(speed_col) || is.null(direction_col)) {
    stop("Could not identify either along/cross-channel columns or speed/direction columns in the currents dataset.")
  }

  cleaned <- clean_speed_direction(
    df = df,
    speed_col = speed_col,
    direction_col = direction_col,
    sentinel_values = c(999, 9999),
    direction_upper_exclusive = 360
  )
  components <- components_from_speed_direction(cleaned$raw_speed, cleaned$raw_direction_deg)

  list(
    method = "speed_direction",
    df = cleaned,
    v1 = components$v1,
    v2 = components$v2,
    raw_speed = cleaned$raw_speed,
    raw_direction_deg = cleaned$raw_direction_deg,
    speed_col = speed_col,
    direction_col = direction_col
  )
}

order_screening_summary <- function(summary_df) {
  if (nrow(summary_df) == 0L) {
    return(summary_df)
  }

  status_rank <- ifelse(summary_df$status == "ok", 0L, 1L)
  p_value_order <- ifelse(is.na(summary_df$p_value), -Inf, summary_df$p_value)
  n_final_order <- ifelse(is.na(summary_df$n_final), -Inf, summary_df$n_final)
  summary_df[order(status_rank, -p_value_order, -n_final_order, summary_df$dataset_id), , drop = FALSE]
}

write_processed_dataset_csv <- function(df, path) {
  output <- df
  if ("datetime" %in% names(output)) {
    output$datetime <- format(as.POSIXct(output$datetime, tz = "UTC"), "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  utils::write.csv(output, path, row.names = FALSE)
  invisible(path)
}

coerce_numeric_columns <- function(df, columns) {
  for (column in intersect(columns, names(df))) {
    df[[column]] <- suppressWarnings(as.numeric(df[[column]]))
  }
  df
}
