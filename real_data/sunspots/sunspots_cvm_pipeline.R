#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stats)
  library(utils)
  library(grDevices)
  library(graphics)
})

repo_root_default <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

parse_bool <- function(x, key) {
  x <- tolower(trimws(x))
  if (x %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (x %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop(sprintf("Option '--%s' must be TRUE or FALSE.", key))
}

parse_args <- function(args) {
  out <- list(
    external_repo = repo_root_default,
    B = 1000L,
    cores = 12L,
    rerun_ks = TRUE,
    rerun_cvm = TRUE,
    out_dir = file.path(
      repo_root_default,
      "real_data", "sunspots", "output", "fmgp_small_circle_window_comparison_with_cvm"
    ),
    work_dir = file.path(
      repo_root_default,
      "real_data", "sunspots", "output", "sunspots_cvm_pipeline_cache"
    )
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) > 1L) parts[[2L]] else NULL
    if (is.null(value) || identical(value, "")) {
      stop(sprintf("Option '--%s' requires the form --%s=value.", key, key))
    }
    if (key %in% c("external_repo", "out_dir", "work_dir")) {
      out[[key]] <- value
    } else if (key %in% c("B", "cores")) {
      out[[key]] <- as.integer(value)
    } else if (key %in% c("rerun_ks", "rerun_cvm")) {
      out[[key]] <- parse_bool(value, key)
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  out$external_repo <- normalizePath(out$external_repo, winslash = "/", mustWork = FALSE)
  out$out_dir <- normalizePath(out$out_dir, winslash = "/", mustWork = FALSE)
  out$work_dir <- normalizePath(out$work_dir, winslash = "/", mustWork = FALSE)
  out
}

safe_log10 <- function(x) {
  log10(pmax(as.numeric(x), 1e-12))
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) {
    return(NA_real_)
  }
  x_ok <- x[ok]
  y_ok <- y[ok]
  if (stats::sd(x_ok) < 1e-12 || stats::sd(y_ok) < 1e-12) {
    return(NA_real_)
  }
  stats::cor(x_ok, y_ok)
}

as_utc <- function(x) {
  as.POSIXct(x, tz = "UTC")
}

run_checked <- function(cmd, args, wd = NULL) {
  old_wd <- NULL
  if (!is.null(wd)) {
    old_wd <- getwd()
    setwd(wd)
    on.exit(setwd(old_wd), add = TRUE)
  }
  status <- system2(cmd, args = args, stdout = "", stderr = "", wait = TRUE)
  if (!identical(status, 0L)) {
    stop("Command failed: ", paste(c(cmd, args), collapse = " "))
  }
}

read_csv_base <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv_base <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
}

infer_group_from_filename <- function(path) {
  b <- basename(path)
  if (grepl("_N_", b, ignore.case = FALSE)) {
    return("N")
  }
  if (grepl("_S_", b, ignore.case = FALSE)) {
    return("S")
  }
  if (grepl("_all_", b, ignore.case = TRUE) || grepl("_both_", b, ignore.case = TRUE)) {
    return("all")
  }
  stop("Could not infer FMGP group from filename: ", b)
}

infer_cycle_from_dates <- function(dates) {
  mid <- as.POSIXct(stats::median(as.numeric(dates)), origin = "1970-01-01", tz = "UTC")
  y <- as.integer(format(mid, "%Y", tz = "UTC"))
  if (y >= 1976 && y <= 1986) {
    return(21L)
  }
  if (y >= 1987 && y <= 1996) {
    return(22L)
  }
  if (y >= 1997 && y <= 2009) {
    return(23L)
  }
  stop("Could not infer cycle from FMGP date range. Median year: ", y)
}

read_fmgp_rdata <- function(path) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!"test.results" %in% loaded && !exists("test.results", envir = env, inherits = FALSE)) {
    stop("FMGP RData does not contain `test.results`: ", path)
  }
  tr <- as.data.frame(get("test.results", envir = env))
  required <- c("date", "PAD_2")
  missing <- setdiff(required, names(tr))
  if (length(missing) > 0L) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }

  dates <- as_utc(tr$date)
  raw <- as.numeric(tr$PAD_2)

  data.frame(
    cycle = infer_cycle_from_dates(dates),
    group = infer_group_from_filename(path),
    fmgp_window_id = seq_len(nrow(tr)),
    fmgp_date = dates,
    fmgp_PAD_2_raw = raw,
    fmgp_PAD_2_BY = stats::p.adjust(raw, method = "BY"),
    fmgp_PAD = if ("PAD" %in% names(tr)) as.numeric(tr$PAD) else NA_real_,
    source_file = basename(path),
    stringsAsFactors = FALSE
  )
}

read_all_fmgp <- function(fmgp_dir) {
  if (!dir.exists(fmgp_dir)) {
    stop("FMGP directory not found: ", fmgp_dir)
  }
  files <- list.files(fmgp_dir, pattern = "\\.RData$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No .RData files found in FMGP directory: ", fmgp_dir)
  }

  out <- do.call(rbind, lapply(files, read_fmgp_rdata))
  out <- out[out$group == "all", , drop = FALSE]
  if (nrow(out) == 0L) {
    stop("No FMGP all/both .RData rows were found in: ", fmgp_dir)
  }

  key <- paste(out$cycle, out$group, out$fmgp_window_id, out$fmgp_date, out$fmgp_PAD_2_raw)
  out <- out[!duplicated(key), , drop = FALSE]
  out[order(out$cycle, out$fmgp_window_id), , drop = FALSE]
}

build_gof_paths <- function(repo, B, statistic) {
  stat_key <- tolower(statistic)
  if (identical(stat_key, "ks") && identical(as.integer(B), 1000L)) {
    output_dir <- file.path(
      repo,
      "real_data", "sunspots", "output",
      "cycles21_23_weighted_mixture_rolling_windows_10cr_B1000", "fast"
    )
  } else {
    stat_label <- if (identical(stat_key, "ks")) "KS" else "CvM"
    folder <- sprintf(
      "cycles21_23_weighted_mixture_rolling_windows_10cr_B%d_all_%s",
      as.integer(B),
      stat_label
    )
    output_dir <- file.path(repo, "real_data", "sunspots", "output", folder)
  }
  results_csv <- file.path(
    output_dir,
    "sunspots_cycles21_23_weighted_mixture_rolling_10cr_gof_results.csv"
  )
  list(output_dir = output_dir, results_csv = results_csv)
}

gof_results_complete <- function(path, statistic) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  dat <- tryCatch(read_csv_base(path), error = function(e) NULL)
  if (is.null(dat) || !("cycle" %in% names(dat))) {
    return(FALSE)
  }
  if ("statistic_type" %in% names(dat)) {
    dat <- dat[tolower(dat$statistic_type) == tolower(statistic), , drop = FALSE]
  }
  identical(sort(unique(as.integer(dat$cycle))), c(21L, 22L, 23L))
}

run_gof_if_needed <- function(repo, B, cores, statistic, rerun, work_dir) {
  paths <- build_gof_paths(repo, B, statistic)
  if (!rerun && gof_results_complete(paths$results_csv, statistic)) {
    message(sprintf("Reusing existing %s results: %s", statistic, paths$results_csv))
    return(paths)
  }

  script <- file.path(
    repo, "real_data", "sunspots", "run_sunspots_weighted_mixture_rolling_windows_gof.R"
  )
  if (!file.exists(script)) {
    stop("Missing GOF script: ", script)
  }

  stat_key <- tolower(statistic)
  paths$output_dir <- file.path(work_dir, sprintf("sunspots_%s_fast_B%d", stat_key, as.integer(B)))
  paths$results_csv <- file.path(
    paths$output_dir,
    "sunspots_cycles21_23_weighted_mixture_rolling_10cr_gof_results.csv"
  )
  dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(paths$output_dir, "sunspots_cycles21_23_weighted_mixture_rolling_10cr.log")
  existing_rows <- 0L
  if (file.exists(paths$results_csv)) {
    existing_rows <- tryCatch(nrow(read_csv_base(paths$results_csv)), error = function(e) 0L)
  }
  message(sprintf(
    "Launching %s fast-multiplier run with B=%d and n_cores=%d.",
    statistic, as.integer(B), as.integer(cores)
  ))
  message(sprintf("Resuming directory: %s", paths$output_dir))
  message(sprintf("Progress log: %s", log_path))
  message(sprintf("Existing completed rows in that directory: %d", existing_rows))

  expr <- sprintf(
    paste0(
      "source(%s); ",
      "run_sunspots_weighted_mixture_rolling_windows_gof(",
      "cycles = c(21L, 22L, 23L), ",
      "grid_groups = 'all', ",
      "statistics = %s, ",
      "B = %dL, ",
      "n_cores = %dL, ",
      "bootstrap_method = 'fast_multiplier', ",
      "output_dir = %s)"
    ),
    shQuote(script),
    shQuote(stat_key),
    as.integer(B),
    as.integer(cores),
    shQuote(paths$output_dir)
  )

  launch_script <- tempfile(
    pattern = sprintf("run_sunspots_%s_fast_", stat_key),
    fileext = ".R"
  )
  writeLines(expr, con = launch_script, useBytes = TRUE)
  on.exit(unlink(launch_script), add = TRUE)

  run_checked("Rscript", c("--vanilla", launch_script), wd = repo)

  if (!file.exists(paths$results_csv)) {
    stop("Expected GOF CSV was not created: ", paths$results_csv)
  }

  paths
}

read_sc_results <- function(path, statistic_name) {
  if (!file.exists(path)) {
    stop("Small-circle rolling results CSV not found: ", path)
  }
  sc <- read_csv_base(path)
  required <- c(
    "cycle", "window_id", "global_window_id",
    "start_date", "end_date", "center_date",
    "n", "n_north", "n_south",
    "p_value_raw", "test_statistic"
  )
  missing <- setdiff(required, names(sc))
  if (length(missing) > 0L) {
    stop("Missing columns in small-circle CSV: ", paste(missing, collapse = ", "))
  }

  if ("grid_group" %in% names(sc)) {
    sc <- sc[sc$grid_group == "all", , drop = FALSE]
  }
  if ("statistic_type" %in% names(sc)) {
    sc <- sc[tolower(sc$statistic_type) == tolower(statistic_name), , drop = FALSE]
  }

  sc$cycle <- as.integer(sc$cycle)
  sc$window_id <- as.integer(sc$window_id)
  sc$global_window_id <- as.integer(sc$global_window_id)
  sc$start_date <- as_utc(sc$start_date)
  sc$end_date <- as_utc(sc$end_date)
  sc$center_date <- as_utc(sc$center_date)
  sc$p_SC <- as.numeric(sc$p_value_raw)
  sc$p_SC_BY <- NA_real_

  by_groups <- interaction(sc$cycle, drop = TRUE)
  for (grp in levels(by_groups)) {
    idx <- which(by_groups == grp)
    ok <- is.finite(sc$p_SC[idx])
    if (any(ok)) {
      sc$p_SC_BY[idx[ok]] <- stats::p.adjust(sc$p_SC[idx[ok]], method = "BY")
    }
  }

  sc[order(sc$cycle, sc$window_id), , drop = FALSE]
}

nearest_fmgp_for_row <- function(sc_row, fmgp_group, tolerance_days = 14) {
  candidates <- fmgp_group[fmgp_group$cycle == sc_row$cycle, , drop = FALSE]
  if (nrow(candidates) == 0L) {
    return(NULL)
  }
  diffs <- as.numeric(difftime(candidates$fmgp_date, sc_row$end_date, units = "days"))
  idx <- which.min(abs(diffs))
  if (length(idx) == 0L || !is.finite(diffs[idx]) || abs(diffs[idx]) > tolerance_days) {
    return(NULL)
  }
  candidates[idx, , drop = FALSE]
}

make_matched_table <- function(sc, fmgp, tolerance_days = 14) {
  out <- sc
  out$FMGP_all_raw <- NA_real_
  out$FMGP_all_BY <- NA_real_
  out$FMGP_all_PAD <- NA_real_
  out$FMGP_all_date <- as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")
  out$FMGP_all_window_id <- NA_integer_
  out$FMGP_all_diff_end_days <- NA_real_

  fg <- fmgp[fmgp$group == "all", , drop = FALSE]
  for (i in seq_len(nrow(out))) {
    m <- nearest_fmgp_for_row(out[i, , drop = FALSE], fg, tolerance_days = tolerance_days)
    if (!is.null(m)) {
      out$FMGP_all_raw[i] <- m$fmgp_PAD_2_raw[1L]
      out$FMGP_all_BY[i] <- m$fmgp_PAD_2_BY[1L]
      out$FMGP_all_PAD[i] <- m$fmgp_PAD[1L]
      out$FMGP_all_date[i] <- m$fmgp_date[1L]
      out$FMGP_all_window_id[i] <- m$fmgp_window_id[1L]
      out$FMGP_all_diff_end_days[i] <- as.numeric(difftime(m$fmgp_date[1L], out$end_date[i], units = "days"))
    }
  }

  out
}

build_combined_table <- function(ks_matched, cvm_matched) {
  key_cols <- c(
    "cycle", "window_id", "global_window_id",
    "start_date", "end_date", "center_date",
    "n", "n_north", "n_south"
  )

  ks_keep <- c(
    key_cols, "p_SC", "p_SC_BY", "test_statistic",
    "FMGP_all_raw", "FMGP_all_BY", "FMGP_all_PAD",
    "FMGP_all_date", "FMGP_all_window_id", "FMGP_all_diff_end_days"
  )
  cvm_keep <- c(key_cols, "p_SC", "p_SC_BY", "test_statistic")

  ks_df <- ks_matched[, ks_keep, drop = FALSE]
  cvm_df <- cvm_matched[, cvm_keep, drop = FALSE]

  names(ks_df)[names(ks_df) == "p_SC"] <- "p_value_ks"
  names(ks_df)[names(ks_df) == "p_SC_BY"] <- "p_value_ks_by"
  names(ks_df)[names(ks_df) == "test_statistic"] <- "stat_ks"

  names(cvm_df)[names(cvm_df) == "p_SC"] <- "p_value_cvm"
  names(cvm_df)[names(cvm_df) == "p_SC_BY"] <- "p_value_cvm_by"
  names(cvm_df)[names(cvm_df) == "test_statistic"] <- "stat_cvm"

  out <- merge(ks_df, cvm_df, by = key_cols, sort = FALSE)
  names(out)[names(out) == "FMGP_all_raw"] <- "p_value_uniformity"
  names(out)[names(out) == "FMGP_all_BY"] <- "p_value_uniformity_by"

  out <- out[order(out$cycle, out$window_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

compute_correlations <- function(df) {
  cycles <- sort(unique(as.integer(df$cycle)))
  rows <- vector("list", length(cycles) + 1L)

  for (i in seq_along(cycles)) {
    cy <- cycles[i]
    d <- df[df$cycle == cy, , drop = FALSE]
    rows[[i]] <- data.frame(
      cycle = as.character(cy),
      n_windows = nrow(d),
      cor_log10_ks_uniformity = safe_cor(safe_log10(d$p_value_ks), safe_log10(d$p_value_uniformity)),
      cor_log10_cvm_uniformity = safe_cor(safe_log10(d$p_value_cvm), safe_log10(d$p_value_uniformity)),
      stringsAsFactors = FALSE
    )
  }

  rows[[length(rows)]] <- data.frame(
    cycle = "21-23",
    n_windows = nrow(df),
    cor_log10_ks_uniformity = safe_cor(safe_log10(df$p_value_ks), safe_log10(df$p_value_uniformity)),
    cor_log10_cvm_uniformity = safe_cor(safe_log10(df$p_value_cvm), safe_log10(df$p_value_uniformity)),
    stringsAsFactors = FALSE
  )

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

draw_series <- function(x, y, col, pch, lty, lwd = 1.4) {
  ord <- order(x)
  graphics::lines(x[ord], y[ord], col = col, lty = lty, lwd = lwd)
  graphics::points(x[ord], y[ord], col = col, pch = pch, cex = 0.8)
}

plot_pvalues <- function(df, file, corrected = FALSE) {
  cycles <- sort(unique(as.integer(df$cycle)))
  y_ks <- if (corrected) "p_value_ks_by" else "p_value_ks"
  y_cvm <- if (corrected) "p_value_cvm_by" else "p_value_cvm"
  y_unif <- if (corrected) "p_value_uniformity_by" else "p_value_uniformity"

  grDevices::png(file, width = 3600, height = 2400, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(length(cycles), 1L), mar = c(4, 5, 3, 1), oma = c(0, 0, 1, 0))

  for (i in seq_along(cycles)) {
    cy <- cycles[i]
    d <- df[df$cycle == cy, , drop = FALSE]
    x <- as.Date(d$end_date, tz = "UTC")
    y1 <- as.numeric(d[[y_ks]])
    y2 <- as.numeric(d[[y_cvm]])
    y3 <- as.numeric(d[[y_unif]])

    graphics::plot(
      x, y1,
      type = "n",
      ylim = c(0, 1),
      xlab = "Window end date",
      ylab = "p-value",
      main = sprintf("Cycle %d", cy)
    )
    graphics::abline(h = 0.05, col = "grey40", lty = "dashed", lwd = 1)
    graphics::abline(h = 0.10, col = "grey40", lty = "dotted", lwd = 1)
    draw_series(x, y1, col = "#1f77b4", pch = 16, lty = "dotdash")
    draw_series(x, y2, col = "#228B22", pch = 15, lty = "solid")
    draw_series(x, y3, col = "#d62728", pch = 17, lty = "dashed")

    if (i == 1L) {
      graphics::legend(
        "top",
        legend = c("Small circle (KS)", "Small circle (CvM)", "Longitudinal uniformity"),
        col = c("#1f77b4", "#228B22", "#d62728"),
        pch = c(16, 15, 17),
        lty = c("dotdash", "solid", "dashed"),
        horiz = TRUE,
        bty = "n",
        inset = 0.01
      )
    }
  }
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(args$work_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(args$external_repo)) {
    stop("External repo not found: ", args$external_repo)
  }

  message("Using external repo: ", args$external_repo)

  ks_paths <- run_gof_if_needed(
    repo = args$external_repo,
    B = args$B,
    cores = args$cores,
    statistic = "KS",
    rerun = args$rerun_ks,
    work_dir = args$work_dir
  )
  cvm_paths <- run_gof_if_needed(
    repo = args$external_repo,
    B = args$B,
    cores = args$cores,
    statistic = "CvM",
    rerun = args$rerun_cvm,
    work_dir = args$work_dir
  )

  fmgp_dir <- file.path(args$external_repo, "real_data", "sunspots", "fmgp_rdata")
  fmgp <- read_all_fmgp(fmgp_dir)
  ks_matched <- make_matched_table(read_sc_results(ks_paths$results_csv, "ks"), fmgp)
  cvm_matched <- make_matched_table(read_sc_results(cvm_paths$results_csv, "cvm"), fmgp)
  combined <- build_combined_table(ks_matched, cvm_matched)
  correlations <- compute_correlations(combined)

  combined_csv <- file.path(args$out_dir, "sc_vs_fmgp_pvalues.csv")
  corr_csv <- file.path(args$out_dir, "sunspots_rolling_correlations.csv")
  raw_png <- file.path(args$out_dir, "sc_vs_fmgp_pvalue_timeplot_raw.png")
  by_png <- file.path(args$out_dir, "sc_vs_fmgp_pvalue_timeplot_BY.png")

  write_csv_base(combined, combined_csv)
  write_csv_base(correlations, corr_csv)
  plot_pvalues(combined, raw_png, corrected = FALSE)
  plot_pvalues(combined, by_png, corrected = TRUE)

  message("Saved:")
  message("  ", combined_csv)
  message("  ", corr_csv)
  message("  ", raw_png)
  message("  ", by_png)
}

main()
