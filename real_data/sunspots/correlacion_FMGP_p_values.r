

#!/usr/bin/env Rscript

# Compare rolling-window small-circle GOF p-values with the longitudinal
# uniformity p-values of Fernández-de-Marcos and García-Portugués.
#
# Expected inputs:
# 1. Our rolling-window GOF results CSV, by default:
#    real_data/sunspots/output/cycles21_23_weighted_mixture_rolling_windows_10cr_B1000/
#      sunspots_cycles21_23_weighted_mixture_rolling_10cr_gof_results.csv
#
# 2. FMGP .RData files in, by default:
#    real_data/sunspots/fmgp_rdata/
#
#    The directory should contain the RData files for all/both hemispheres,
#    10-Carrington-rotation windows, cycles 21--23. If the directory also
#    contains N/S files, they are ignored by this script.
#
# Outputs:
#   real_data/sunspots/output/fmgp_small_circle_window_comparison/
#
# Main outputs:
#   - fmgp_pvalues_long.csv
#   - sc_vs_fmgp_window_matched_cycles21_23.csv
#   - sc_vs_fmgp_correlations.csv
#   - sc_vs_fmgp_contingency_tables.csv
#   - sc_vs_fmgp_block_summary.csv
#   - sc_vs_fmgp_pvalue_timeplot_raw.png
#   - sc_vs_fmgp_pvalue_timeplot_BY.png
#   - sc_vs_fmgp_scatter_all.png
#   - sc_vs_fmgp_sessionInfo.txt
#   - sc_vs_fmgp_alignment_summary.csv
#   - sc_vs_fmgp_observed_date_alignment.csv, if a sunspots CSV is available

suppressPackageStartupMessages({
  library(stats)
  library(utils)
  library(grDevices)
  library(graphics)
})

plot_col_sc <- "#1f77b4"
plot_col_fmgp <- "#d62728"
plot_col_ref <- "gray40"

repo_dir <- "/Users/Diego/Desktop/Codigo/goodness_of_fit_metric_spaces"
if (dir.exists(repo_dir)) {
  setwd(repo_dir)
}

#
# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------
# The script can be run with default paths, but the paths can also be overridden
# from the terminal. This is useful because the rolling-window GOF results and
# the FMGP .RData files may live in different folders while we are experimenting.
parse_args <- function() {
  out <- list(
    sc_csv = file.path(
      "real_data", "sunspots", "output",
      "cycles21_23_weighted_mixture_rolling_windows_10cr_B1000",
      "sunspots_cycles21_23_weighted_mixture_rolling_10cr_gof_results.csv"
    ),
    fmgp_dir = file.path("real_data", "sunspots", "fmgp_rdata"),
    sunspots_csv = file.path("real_data", "sunspots", "output", "sunspots_cycles21_23_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "fmgp_small_circle_window_comparison"),
    tolerance_days = 14,
    alpha = 0.05
  )

  args <- commandArgs(trailingOnly = TRUE)
  if (any(args %in% c("--help", "-h"))) {
    cat(
      "Usage:\n",
      "  Rscript real_data/sunspots/correlacion_FMGP_p_values.r [options]\n\n",
      "Options:\n",
      "  --sc_csv=PATH          Rolling-window small-circle GOF results CSV.\n",
      "  --fmgp_dir=PATH        Directory containing FMGP .RData files.\n",
      "  --sunspots_csv=PATH    Optional sunspots CSV used to diagnose FMGP date alignment.\n",
      "  --output_dir=PATH      Output directory.\n",
      "  --tolerance_days=NUM   Maximum allowed date mismatch for nearest-date matching. Default: 14.\n",
      "  --alpha=NUM            Significance level for binary summaries. Default: 0.05.\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }

  for (arg in args) {
    if (grepl("^--sc_csv=", arg)) {
      out$sc_csv <- sub("^--sc_csv=", "", arg)
    } else if (grepl("^--fmgp_dir=", arg)) {
      out$fmgp_dir <- sub("^--fmgp_dir=", "", arg)
    } else if (grepl("^--sunspots_csv=", arg)) {
      out$sunspots_csv <- sub("^--sunspots_csv=", "", arg)
    } else if (grepl("^--output_dir=", arg)) {
      out$output_dir <- sub("^--output_dir=", "", arg)
    } else if (grepl("^--tolerance_days=", arg)) {
      out$tolerance_days <- as.numeric(sub("^--tolerance_days=", "", arg))
    } else if (grepl("^--alpha=", arg)) {
      out$alpha <- as.numeric(sub("^--alpha=", "", arg))
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  out
}

#
# -----------------------------------------------------------------------------
# Small date/name helpers
# -----------------------------------------------------------------------------
# All comparisons are done in UTC to avoid accidental shifts from local time
# zones. The FMGP files are identified in two ways:
#   - the group, N/S/all, is inferred from the file name;
#   - the solar cycle is inferred from the range of dates inside the file.
# The main analysis keeps only the all/both group.
as_utc <- function(x) {
  as.POSIXct(x, tz = "UTC")
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

# -----------------------------------------------------------------------------
# Diagnose what FMGP's stored `date` represents
# -----------------------------------------------------------------------------
# Fernández-de-Marcos suggested that their stored window date might have been
# computed automatically from the first/last sunspot actually present in the
# window, rather than from the theoretical Carrington grid endpoint. This helper
# checks that hypothesis.
#
# For each matched window and group all, we take the sunspots in our window
# [start_date, end_date), and compare FMGP's
# stored date with:
#   - the theoretical end_date;
#   - the first observed date in that group/window;
#   - the last observed date in that group/window.
# If FMGP dates are much closer to the last observed date than to end_date, then
# the group-specific shifts are explained by gaps with no
# sunspots near the theoretical endpoint.
read_sunspots_for_alignment <- function(path) {
  candidate_paths <- trimws(strsplit(path, ",", fixed = TRUE)[[1L]])

  if (length(candidate_paths) == 1L && !file.exists(candidate_paths)) {
    default_dir <- dirname(candidate_paths)
    separate_candidates <- file.path(
      default_dir,
      sprintf("sunspots_cycle%d_s2_all.csv", 21:23)
    )
    if (all(file.exists(separate_candidates))) {
      candidate_paths <- separate_candidates
    }
  }

  existing_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing_paths) == 0L) {
    return(NULL)
  }

  data_list <- lapply(existing_paths, function(one_path) {
    dat <- utils::read.csv(one_path, stringsAsFactors = FALSE)
    required <- c("cycle", "date", "hemisphere")
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
      warning("Cannot use sunspots CSV for observed-date alignment diagnostic; file `", one_path,
              "` is missing columns: ", paste(missing, collapse = ", "))
      return(NULL)
    }
    dat$cycle <- as.integer(dat$cycle)
    dat$date <- as_utc(dat$date)
    dat
  })

  data_list <- Filter(Negate(is.null), data_list)
  if (length(data_list) == 0L) {
    return(NULL)
  }

  out <- do.call(rbind, data_list)
  out[order(out$cycle, out$date), , drop = FALSE]
}

make_observed_date_alignment <- function(matched, sunspots) {
  if (is.null(sunspots)) {
    return(NULL)
  }

  rows <- list()
  k <- 0L
  groups <- c("all")

  for (i in seq_len(nrow(matched))) {
    row_i <- matched[i, , drop = FALSE]
    for (g in groups) {
      fmgp_date_col <- paste0("FMGP_", g, "_date")
      if (!fmgp_date_col %in% names(row_i) || is.na(row_i[[fmgp_date_col]][1L])) {
        next
      }

      in_window <- sunspots$cycle == row_i$cycle[1L] &
        sunspots$date >= row_i$start_date[1L] &
        sunspots$date < row_i$end_date[1L]

      dates_g <- sunspots$date[in_window]
      if (length(dates_g) == 0L) {
        next
      }

      first_obs <- min(dates_g)
      last_obs <- max(dates_g)
      fmgp_date <- row_i[[fmgp_date_col]][1L]

      diff_to_end <- as.numeric(difftime(fmgp_date, row_i$end_date[1L], units = "days"))
      diff_to_first_obs <- as.numeric(difftime(fmgp_date, first_obs, units = "days"))
      diff_to_last_obs <- as.numeric(difftime(fmgp_date, last_obs, units = "days"))
      gap_last_obs_to_end <- as.numeric(difftime(row_i$end_date[1L], last_obs, units = "days"))
      gap_start_to_first_obs <- as.numeric(difftime(first_obs, row_i$start_date[1L], units = "days"))

      candidate <- c(
        end_date = abs(diff_to_end),
        first_observed = abs(diff_to_first_obs),
        last_observed = abs(diff_to_last_obs)
      )
      closest <- names(candidate)[which.min(candidate)]

      k <- k + 1L
      rows[[k]] <- data.frame(
        cycle = row_i$cycle[1L],
        window_id = row_i$window_id[1L],
        group = g,
        n_group_window = length(dates_g),
        start_date = row_i$start_date[1L],
        end_date = row_i$end_date[1L],
        fmgp_date = fmgp_date,
        first_observed_date = first_obs,
        last_observed_date = last_obs,
        diff_fmgp_to_end_days = diff_to_end,
        diff_fmgp_to_first_observed_days = diff_to_first_obs,
        diff_fmgp_to_last_observed_days = diff_to_last_obs,
        abs_diff_fmgp_to_end_days = abs(diff_to_end),
        abs_diff_fmgp_to_first_observed_days = abs(diff_to_first_obs),
        abs_diff_fmgp_to_last_observed_days = abs(diff_to_last_obs),
        gap_start_to_first_observed_days = gap_start_to_first_obs,
        gap_last_observed_to_end_days = gap_last_obs_to_end,
        closest_reference = closest,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(NULL)
  }

  do.call(rbind, rows)
}

summarize_observed_date_alignment <- function(observed_alignment) {
  if (is.null(observed_alignment) || nrow(observed_alignment) == 0L) {
    return(NULL)
  }

  combos <- unique(observed_alignment[, c("cycle", "group")])
  combos <- combos[order(combos$cycle, combos$group), , drop = FALSE]
  rows <- vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    cy <- combos$cycle[i]
    g <- combos$group[i]
    df <- observed_alignment[observed_alignment$cycle == cy & observed_alignment$group == g, , drop = FALSE]
    rows[[i]] <- data.frame(
      cycle = cy,
      group = g,
      n_windows = nrow(df),
      median_abs_diff_to_end_days = stats::median(df$abs_diff_fmgp_to_end_days, na.rm = TRUE),
      median_abs_diff_to_first_observed_days = stats::median(df$abs_diff_fmgp_to_first_observed_days, na.rm = TRUE),
      median_abs_diff_to_last_observed_days = stats::median(df$abs_diff_fmgp_to_last_observed_days, na.rm = TRUE),
      proportion_closest_to_end = mean(df$closest_reference == "end_date", na.rm = TRUE),
      proportion_closest_to_first_observed = mean(df$closest_reference == "first_observed", na.rm = TRUE),
      proportion_closest_to_last_observed = mean(df$closest_reference == "last_observed", na.rm = TRUE),
      median_gap_last_observed_to_end_days = stats::median(df$gap_last_observed_to_end_days, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

#
# -----------------------------------------------------------------------------
# Read one FMGP .RData file
# -----------------------------------------------------------------------------
# Each FMGP file contains a data.frame called `test.results`. The important
# column for us is `PAD_2`, the p-value for their projected Anderson--Darling
# longitudinal uniformity test in 10-Carrington-rotation rolling windows.
#
# The values stored in `PAD_2` are treated as raw p-values. We also compute
# Benjamini--Yekutieli adjusted p-values within each file, matching the multiple
# testing correction used for the sequential plots in FMGP.
read_fmgp_rdata <- function(path) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!"test.results" %in% loaded && !exists("test.results", envir = env, inherits = FALSE)) {
    stop("FMGP RData does not contain `test.results`: ", path)
  }
  tr <- get("test.results", envir = env)
  tr <- as.data.frame(tr)

  required <- c("date", "PAD_2")
  missing <- setdiff(required, names(tr))
  if (length(missing) > 0L) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }

  dates <- as_utc(tr$date)
  group <- infer_group_from_filename(path)
  cycle <- infer_cycle_from_dates(dates)

  raw <- as.numeric(tr$PAD_2)
  by <- stats::p.adjust(raw, method = "BY")

  data.frame(
    cycle = cycle,
    group = group,
    fmgp_window_id = seq_len(nrow(tr)),
    fmgp_date = dates,
    fmgp_PAD_2_raw = raw,
    fmgp_PAD_2_BY = by,
    fmgp_PAD = if ("PAD" %in% names(tr)) as.numeric(tr$PAD) else NA_real_,
    source_file = basename(path),
    stringsAsFactors = FALSE
  )
}

#
# -----------------------------------------------------------------------------
# Read and stack all FMGP files
# -----------------------------------------------------------------------------
#
# The expected folder contains the all/both .RData files for cycles 21--23.
# If the directory also contains N/S files, they are read and then ignored.
# During manual downloads it is easy to duplicate a file, so after reading
# everything we remove exact duplicate rows.
read_all_fmgp <- function(fmgp_dir) {
  if (!dir.exists(fmgp_dir)) {
    stop(
      "FMGP directory not found: ", fmgp_dir, "\n",
      "Create it and place the FMGP all/both .RData files there, e.g.\n",
      "  sunspots_cir_test_all_10_results_cycle21.RData\n",
      "  sunspots_cir_test_all_10_results_cycle22.RData\n",
      "  sunspots_cir_test_all_10_results_cycle23.RData\n"
    )
  }
  files <- list.files(fmgp_dir, pattern = "\\.RData$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No .RData files found in FMGP directory: ", fmgp_dir)
  }

  out <- do.call(rbind, lapply(files, read_fmgp_rdata))

  # The main comparison is only against FMGP all/both windows. FMGP N/S were
  # useful to diagnose their grid construction, but they are not part of the
  # joint small-circle mixture comparison.
  out <- out[out$group == "all", , drop = FALSE]
  if (nrow(out) == 0L) {
    stop("No FMGP all/both .RData rows were found in: ", fmgp_dir)
  }

  # Remove exact duplicated files if the same cycle/group/date/raw appears twice.
  key <- paste(out$cycle, out$group, out$fmgp_window_id, out$fmgp_date, out$fmgp_PAD_2_raw)
  out <- out[!duplicated(key), , drop = FALSE]

  out[order(out$cycle, out$group, out$fmgp_window_id), , drop = FALSE]
}

#
# -----------------------------------------------------------------------------
# Read our small-circle rolling-window GOF results
# -----------------------------------------------------------------------------
# These are the p-values from our composite KS test for the weighted mixture of
# two small circles, computed on rolling windows of 10 Carrington rotations.
# We store p_SC and -log10(p_SC), the latter being convenient for correlation
# plots because small p-values become large positive values.
read_sc_results <- function(path) {
  if (!file.exists(path)) {
    stop("Small-circle rolling results CSV not found: ", path)
  }
  sc <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("cycle", "window_id", "start_date", "end_date", "center_date", "n", "p_value_raw", "test_statistic")
  missing <- setdiff(required, names(sc))
  if (length(missing) > 0L) {
    stop("Missing columns in small-circle CSV: ", paste(missing, collapse = ", "))
  }
  sc$cycle <- as.integer(sc$cycle)
  sc$window_id <- as.integer(sc$window_id)
  sc$start_date <- as_utc(sc$start_date)
  sc$end_date <- as_utc(sc$end_date)
  sc$center_date <- as_utc(sc$center_date)
  sc$p_SC <- as.numeric(sc$p_value_raw)
  sc$neglog10_p_SC <- -log10(pmax(sc$p_SC, .Machine$double.xmin))
  if ("grid_group" %in% names(sc)) {
    sc <- sc[sc$grid_group == "all", , drop = FALSE]
  }

  # Match FMGP's multiple-testing correction: BY/FDR adjustment over the
  # rolling windows within each displayed cycle. If several statistics are
  # present in the same CSV, adjust them separately.
  sc$p_SC_BY <- NA_real_
  by_groups <- interaction(
    sc$cycle,
    if ("statistic_type" %in% names(sc)) sc$statistic_type else rep("stat", nrow(sc)),
    drop = TRUE
  )
  for (grp in levels(by_groups)) {
    idx <- which(by_groups == grp)
    ok <- is.finite(sc$p_SC[idx])
    if (any(ok)) {
      sc$p_SC_BY[idx[ok]] <- stats::p.adjust(sc$p_SC[idx[ok]], method = "BY")
    }
  }
  sc$neglog10_p_SC_BY <- -log10(pmax(sc$p_SC_BY, .Machine$double.xmin))

  sc[order(sc$cycle, sc$window_id), , drop = FALSE]
}

#
# -----------------------------------------------------------------------------
# Match one of our windows to one FMGP window
# -----------------------------------------------------------------------------
# The natural comparison is window-by-window. In principle, our rolling windows
# and the FMGP windows should be the same 10-Carrington-rotation windows. To be
# robust to tiny differences in how dates were stored, we match by the closest
# FMGP `date` to our window `end_date`, and reject matches that are farther away
# than `tolerance_days`.
nearest_fmgp_for_row <- function(sc_row, fmgp_group, tolerance_days) {
  candidates <- fmgp_group[fmgp_group$cycle == sc_row$cycle, , drop = FALSE]
  if (nrow(candidates) == 0L) {
    return(NULL)
  }
  diffs <- as.numeric(difftime(candidates$fmgp_date, sc_row$end_date, units = "days"))
  idx <- which.min(abs(diffs))
  if (length(idx) == 0L || !is.finite(diffs[idx])) {
    return(NULL)
  }
  if (abs(diffs[idx]) > tolerance_days) {
    return(NULL)
  }
  candidates[idx, , drop = FALSE]
}

#
# -----------------------------------------------------------------------------
# Build the main comparison table
# -----------------------------------------------------------------------------
# Starting from our small-circle table, this adds the corresponding FMGP all/both
# p-values. The main analysis intentionally ignores FMGP N/S: our GOF statistic
# is for the joint small-circle mixture fitted to all sunspots in each window.
make_matched_table <- function(sc, fmgp, tolerance_days) {
  out <- sc

  out$FMGP_all_raw <- NA_real_
  out$FMGP_all_BY <- NA_real_
  out$FMGP_all_PAD <- NA_real_
  out$FMGP_all_date <- as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")
  out$FMGP_all_window_id <- NA_integer_
  out$FMGP_all_diff_end_days <- NA_real_

  fg <- fmgp[fmgp$group == "all", , drop = FALSE]
  for (i in seq_len(nrow(out))) {
    m <- nearest_fmgp_for_row(out[i, , drop = FALSE], fg, tolerance_days)
    if (!is.null(m)) {
      out$FMGP_all_raw[i] <- m$fmgp_PAD_2_raw[1L]
      out$FMGP_all_BY[i] <- m$fmgp_PAD_2_BY[1L]
      out$FMGP_all_PAD[i] <- m$fmgp_PAD[1L]
      out$FMGP_all_date[i] <- m$fmgp_date[1L]
      out$FMGP_all_window_id[i] <- m$fmgp_window_id[1L]
      out$FMGP_all_diff_end_days[i] <- as.numeric(difftime(m$fmgp_date[1L], out$end_date[i], units = "days"))
    }
  }

  out$neglog10_FMGP_all_raw <- -log10(pmax(out$FMGP_all_raw, .Machine$double.xmin))
  out$neglog10_FMGP_all_BY <- -log10(pmax(out$FMGP_all_BY, .Machine$double.xmin))

  out
}

#
# -----------------------------------------------------------------------------
# Numerical summaries: correlations and contingency tables
# -----------------------------------------------------------------------------
# The main continuous diagnostic is the correlation between -log10(p_SC) and
# -log10(p_FMGP) window by window. Large values of both mean strong evidence
# against the corresponding nulls.
#
# The contingency tables are a cruder binary version: they compare which windows
# reject at a chosen alpha level for our small-circle GOF and for FMGP's
# longitudinal uniformity test.
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor(x[ok], y[ok])
}

make_correlations <- function(matched) {
  targets <- c("FMGP_all_raw", "FMGP_all_BY")
  cycles <- c(sort(unique(matched$cycle)), NA_integer_)
  rows <- list()
  k <- 0L
  for (cy in cycles) {
    df <- if (is.na(cy)) matched else matched[matched$cycle == cy, , drop = FALSE]
    label <- if (is.na(cy)) "all_cycles" else as.character(cy)
    for (target in targets) {
      p_col <- target
      nl_col <- paste0("neglog10_", target)
      ok <- is.finite(df$p_SC) & is.finite(df[[p_col]])
      k <- k + 1L
      rows[[k]] <- data.frame(
        cycle = label,
        SC_version = "raw",
        comparator = target,
        n_matched = sum(ok),
        cor_neglog10 = safe_cor(df$neglog10_p_SC, df[[nl_col]]),
        median_p_SC_when_FMGP_lt_005 = if (any(ok & df[[p_col]] < 0.05)) stats::median(df$p_SC[ok & df[[p_col]] < 0.05]) else NA_real_,
        median_p_SC_when_FMGP_ge_005 = if (any(ok & df[[p_col]] >= 0.05)) stats::median(df$p_SC[ok & df[[p_col]] >= 0.05]) else NA_real_,
        stringsAsFactors = FALSE
      )

      k <- k + 1L
      rows[[k]] <- data.frame(
        cycle = label,
        SC_version = "BY",
        comparator = target,
        n_matched = sum(ok),
        cor_neglog10 = safe_cor(df$neglog10_p_SC_BY, df[[nl_col]]),
        median_p_SC_when_FMGP_lt_005 = if (any(ok & df[[p_col]] < 0.05)) stats::median(df$p_SC_BY[ok & df[[p_col]] < 0.05]) else NA_real_,
        median_p_SC_when_FMGP_ge_005 = if (any(ok & df[[p_col]] >= 0.05)) stats::median(df$p_SC_BY[ok & df[[p_col]] >= 0.05]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

make_contingencies <- function(matched, alpha) {
  comparators <- c("FMGP_all_raw", "FMGP_all_BY")
  cycles <- c(sort(unique(matched$cycle)), NA_integer_)
  rows <- list()
  k <- 0L
  for (cy in cycles) {
    df <- if (is.na(cy)) matched else matched[matched$cycle == cy, , drop = FALSE]
    cycle_label <- if (is.na(cy)) "all_cycles" else as.character(cy)
    for (comp in comparators) {
      ok <- is.finite(df$p_SC) & is.finite(df[[comp]])
      if (!any(ok)) next
      fm_rej <- df[[comp]][ok] < alpha
      for (sc_version in c("raw", "BY")) {
        sc_p <- if (sc_version == "raw") df$p_SC[ok] else df$p_SC_BY[ok]
        sc_rej <- sc_p < alpha
        tab <- table(
          SC_reject = factor(sc_rej, levels = c(FALSE, TRUE)),
          FMGP_reject = factor(fm_rej, levels = c(FALSE, TRUE))
        )
        k <- k + 1L
        rows[[k]] <- data.frame(
          cycle = cycle_label,
          SC_version = sc_version,
          comparator = comp,
          alpha = alpha,
          n = sum(ok),
          SC_no_FMGP_no = as.integer(tab["FALSE", "FALSE"]),
          SC_no_FMGP_yes = as.integer(tab["FALSE", "TRUE"]),
          SC_yes_FMGP_no = as.integer(tab["TRUE", "FALSE"]),
          SC_yes_FMGP_yes = as.integer(tab["TRUE", "TRUE"]),
          FMGP_reject_implies_SC_reject_rate = if (sum(fm_rej) > 0) mean(sc_rej[fm_rej]) else NA_real_,
          SC_reject_implies_FMGP_reject_rate = if (sum(sc_rej) > 0) mean(fm_rej[sc_rej]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

#
# -----------------------------------------------------------------------------
# Consecutive rejection/non-rejection blocks
# -----------------------------------------------------------------------------
# Since the windows are overlapping and ordered in time, isolated windows are
# less interpretable than runs of consecutive windows. This summarizes temporal
# blocks where a chosen FMGP p-value is below alpha or above alpha, and reports
# the corresponding small-circle p-values inside each block.
make_blocks <- function(matched, p_col, alpha) {
  rows <- list()
  k <- 0L
  for (cy in sort(unique(matched$cycle))) {
    df <- matched[matched$cycle == cy, , drop = FALSE]
    df <- df[order(df$window_id), , drop = FALSE]
    flag <- is.finite(df[[p_col]]) & df[[p_col]] < alpha
    if (length(flag) == 0L) next
    r <- rle(flag)
    ends <- cumsum(r$lengths)
    starts <- c(1L, head(ends, -1L) + 1L)
    for (j in seq_along(r$lengths)) {
      idx <- starts[j]:ends[j]
      k <- k + 1L
      rows[[k]] <- data.frame(
        cycle = cy,
        comparator = p_col,
        alpha = alpha,
        block_id = j,
        rejected = r$values[j],
        first_window_id = df$window_id[min(idx)],
        last_window_id = df$window_id[max(idx)],
        start_date = df$start_date[min(idx)],
        end_date = df$end_date[max(idx)],
        n_windows = length(idx),
        median_p_SC = stats::median(df$p_SC[idx], na.rm = TRUE),
        min_p_SC = min(df$p_SC[idx], na.rm = TRUE),
        median_p_FMGP = stats::median(df[[p_col]][idx], na.rm = TRUE),
        min_p_FMGP = min(df[[p_col]][idx], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#
# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
# The time plot overlays our p-values and the FMGP p-values along the rolling
# windows. The scatter plot compares the strength of evidence directly via
# -log10 p-values.
plot_time_pvalues <- function(matched, output_file, adjusted = FALSE) {
  png(output_file, width = 1800, height = 1200, res = 180)
  on.exit(dev.off(), add = TRUE)
  cycles <- sort(unique(matched$cycle))
  op <- par(mfrow = c(3, 1), mar = c(3.8, 4.2, 2.5, 1.2), oma = c(0, 0, 2, 0))
  on.exit(par(op), add = TRUE)

  sc_col <- if (adjusted) "p_SC_BY" else "p_SC"
  fmgp_col <- if (adjusted) "FMGP_all_BY" else "FMGP_all_raw"
  legend_labels <- if (adjusted) {
    c("Small circle (BY)", "Longitudinal uniformity (BY)", "0.05", "0.10")
  } else {
    c("Small circle", "Longitudinal uniformity", "0.05", "0.10")
  }

  for (cy in cycles) {
    df <- matched[matched$cycle == cy, , drop = FALSE]
    df <- df[order(df$end_date), , drop = FALSE]
    plot(
      df$end_date, df[[sc_col]],
      type = "b", pch = 16, cex = 0.55, col = plot_col_sc,
      lwd = 1.5,
      ylim = c(0, 1),
      xlab = if (identical(cy, cycles[length(cycles)])) "Window end date" else "",
      ylab = "p-value",
      main = paste("Cycle", cy)
    )
    lines(
      df$end_date, df[[fmgp_col]],
      type = "b", pch = 17, cex = 0.5,
      lty = 2, col = plot_col_fmgp, lwd = 1.5
    )
    abline(h = c(0.05, 0.10), lty = c(2, 3), col = plot_col_ref)
    if (identical(cy, cycles[1L])) {
      legend(
        "topright",
        legend = legend_labels,
        lty = c(1, 2, 2, 3),
        pch = c(16, 17, NA, NA),
        col = c(plot_col_sc, plot_col_fmgp, plot_col_ref, plot_col_ref),
        bty = "n",
        cex = 0.75
      )
    }
  }
}

plot_scatter_all <- function(matched, output_file) {
  png(output_file, width = 1400, height = 1200, res = 180)
  on.exit(dev.off(), add = TRUE)

  x <- matched$neglog10_FMGP_all_raw
  y <- matched$neglog10_p_SC
  ok <- is.finite(x) & is.finite(y)
  plot(
    x[ok], y[ok],
    pch = 16,
    col = plot_col_sc,
    xlab = expression(-log[10](p[FMGP~all~raw])),
    ylab = expression(-log[10](p[SC])),
    main = "Window-by-window comparison: joint small-circle GOF vs FMGP all"
  )
  grid()
  if (sum(ok) >= 3L) {
    fit <- stats::lm(y[ok] ~ x[ok])
    abline(fit, lwd = 2, col = plot_col_fmgp)
    legend(
      "topleft",
      legend = sprintf("cor = %.3f", stats::cor(x[ok], y[ok])),
      bty = "n"
    )
  }
}

#
# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------
# Workflow:
#   1. read our rolling-window small-circle p-values;
#   2. read the FMGP longitudinal p-values from the .RData files;
#   3. match windows by date;
#   4. compute correlations, binary rejection tables and temporal blocks;
#   5. save CSVs and plots for inspection.
args <- parse_args()
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading small-circle results: ", args$sc_csv)
sc <- read_sc_results(args$sc_csv)

message("Reading FMGP RData files from: ", args$fmgp_dir)
fmgp <- read_all_fmgp(args$fmgp_dir)

message("Reading sunspots data for observed-date alignment diagnostic, if available: ", args$sunspots_csv)
sunspots_for_alignment <- read_sunspots_for_alignment(args$sunspots_csv)

message("Matching windows by nearest FMGP date to small-circle end_date, tolerance = ", args$tolerance_days, " days")
matched <- make_matched_table(sc, fmgp, tolerance_days = args$tolerance_days)

correlations <- make_correlations(matched)
contingencies_005 <- make_contingencies(matched, alpha = 0.05)
contingencies_010 <- make_contingencies(matched, alpha = 0.10)
contingencies <- rbind(contingencies_005, contingencies_010)
blocks <- make_blocks(matched, p_col = "FMGP_all_raw", alpha = 0.05)

alignment_summary <- do.call(rbind, lapply(c("all"), function(g) {
  diff_col <- paste0("FMGP_", g, "_diff_end_days")
  rows <- lapply(sort(unique(matched$cycle)), function(cy) {
    d <- matched[matched$cycle == cy, diff_col]
    data.frame(
      cycle = cy,
      group = g,
      n_matched = sum(is.finite(d)),
      median_abs_diff_end_days = stats::median(abs(d), na.rm = TRUE),
      max_abs_diff_end_days = max(abs(d), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}))

observed_date_alignment <- make_observed_date_alignment(matched, sunspots_for_alignment)
observed_date_alignment_summary <- summarize_observed_date_alignment(observed_date_alignment)

utils::write.csv(
  fmgp,
  file.path(args$output_dir, "fmgp_pvalues_long.csv"),
  row.names = FALSE
)
utils::write.csv(
  matched,
  file.path(args$output_dir, "sc_vs_fmgp_window_matched_cycles21_23.csv"),
  row.names = FALSE
)
utils::write.csv(
  correlations,
  file.path(args$output_dir, "sc_vs_fmgp_correlations.csv"),
  row.names = FALSE
)
utils::write.csv(
  contingencies,
  file.path(args$output_dir, "sc_vs_fmgp_contingency_tables.csv"),
  row.names = FALSE
)
utils::write.csv(
  blocks,
  file.path(args$output_dir, "sc_vs_fmgp_block_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  alignment_summary,
  file.path(args$output_dir, "sc_vs_fmgp_alignment_summary.csv"),
  row.names = FALSE
)

if (!is.null(observed_date_alignment)) {
  utils::write.csv(
    observed_date_alignment,
    file.path(args$output_dir, "sc_vs_fmgp_observed_date_alignment.csv"),
    row.names = FALSE
  )
}
if (!is.null(observed_date_alignment_summary)) {
  utils::write.csv(
    observed_date_alignment_summary,
    file.path(args$output_dir, "sc_vs_fmgp_observed_date_alignment_summary.csv"),
    row.names = FALSE
  )
}

plot_time_pvalues(
  matched,
  file.path(args$output_dir, "sc_vs_fmgp_pvalue_timeplot_raw.png"),
  adjusted = FALSE
)
plot_time_pvalues(
  matched,
  file.path(args$output_dir, "sc_vs_fmgp_pvalue_timeplot_BY.png"),
  adjusted = TRUE
)
plot_scatter_all(
  matched,
  file.path(args$output_dir, "sc_vs_fmgp_scatter_all.png")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(args$output_dir, "sc_vs_fmgp_sessionInfo.txt")
)

message("Done. Outputs written to: ", args$output_dir)
message("Key files:")
message("- ", file.path(args$output_dir, "sc_vs_fmgp_window_matched_cycles21_23.csv"))
message("- ", file.path(args$output_dir, "sc_vs_fmgp_correlations.csv"))
message("- ", file.path(args$output_dir, "sc_vs_fmgp_contingency_tables.csv"))
message("- ", file.path(args$output_dir, "sc_vs_fmgp_alignment_summary.csv"))
if (!is.null(observed_date_alignment_summary)) {
  message("- ", file.path(args$output_dir, "sc_vs_fmgp_observed_date_alignment_summary.csv"))
}