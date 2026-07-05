#!/usr/bin/env Rscript

resolve_hvmf_generation_check_path <- function(...) {
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

utils_path_hvmf_generation_check <- resolve_hvmf_generation_check_path("utils.R")
calibration_path_hvmf_generation_check <- resolve_hvmf_generation_check_path(
  "bootstrap",
  "calibration_study.R"
)

source(utils_path_hvmf_generation_check)
source(calibration_path_hvmf_generation_check)

parse_named_args_hvmf_generation_check <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(pieces) < 2L) {
      next
    }
    key <- trimws(pieces[[1L]])
    value <- trimws(paste(pieces[-1L], collapse = "="))
    if (nzchar(key)) {
      out[[key]] <- value
    }
  }
  out
}

parse_flag_hvmf_generation_check <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_integer_hvmf_generation_check <- function(x, default) {
  if (is.null(x) || !nzchar(x)) {
    return(as.integer(default))
  }
  value <- suppressWarnings(as.integer(x))
  if (!is.finite(value) || value < 1L) {
    as.integer(default)
  } else {
    value
  }
}

extract_hvmf_metadata_from_path <- function(path) {
  filename <- basename(path)
  parent_dir <- basename(dirname(path))
  grandparent_dir <- basename(dirname(dirname(path)))

  kappa_match <- regmatches(filename, regexpr("kappa[0-9]+", filename))
  if (length(kappa_match) == 0L || !nzchar(kappa_match)) {
    kappa_match <- regmatches(grandparent_dir, regexpr("kappa[0-9]+", grandparent_dir))
  }

  n_match <- regmatches(filename, regexpr("N_[0-9]+", filename))
  if (length(n_match) == 0L || !nzchar(n_match)) {
    n_match <- regmatches(parent_dir, regexpr("n[0-9]+", parent_dir))
  }

  sample_match <- regmatches(filename, regexpr("samp_[0-9]+", filename))

  list(
    kappa = suppressWarnings(as.numeric(sub("^kappa", "", kappa_match))),
    n = suppressWarnings(as.integer(sub("^N_|^n", "", n_match))),
    sample_id = suppressWarnings(as.integer(sub("^samp_", "", sample_match)))
  )
}

find_hvmf_generation_file <- function(data_file = NULL,
                                      data_dir = default_hvmf_calibration_data_dir(),
                                      kappa = 50,
                                      n = 50,
                                      sample_id = 1L) {
  if (!is.null(data_file) && nzchar(data_file)) {
    resolved_file <- resolve_hvmf_generation_check_path(data_file)
    if (!file.exists(resolved_file)) {
      stop(sprintf("The requested data file does not exist: %s", data_file))
    }
    return(normalizePath(resolved_file, winslash = "/", mustWork = TRUE))
  }

  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  candidates <- c(
    file.path(
      data_dir,
      sprintf("kappa%s", format(kappa, scientific = FALSE, trim = TRUE)),
      sprintf("n%s", as.integer(n)),
      sprintf("samp_%s_N_%s_kappa%s.csv", as.integer(sample_id), as.integer(n), format(kappa, scientific = FALSE, trim = TRUE))
    ),
    file.path(
      data_dir,
      sprintf("kappa%s", format(kappa, scientific = FALSE, trim = TRUE)),
      sprintf("samp_%s_N_%s_kappa%s.csv", as.integer(sample_id), as.integer(n), format(kappa, scientific = FALSE, trim = TRUE))
    )
  )

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  pattern <- sprintf(
    "^samp_%s_N_%s_kappa%s\\.csv$",
    as.integer(sample_id),
    as.integer(n),
    format(kappa, scientific = FALSE, trim = TRUE)
  )
  matches <- list.files(data_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(matches) >= 1L) {
    return(normalizePath(matches[[1L]], winslash = "/", mustWork = TRUE))
  }

  stop(sprintf(
    paste(
      "No HvMF dataset found for kappa = %s, n = %s, sample_id = %s.",
      "Expected, for example, `%s`."
    ),
    format(kappa, scientific = FALSE, trim = TRUE),
    as.integer(n),
    as.integer(sample_id),
    candidates[[1L]]
  ))
}

load_or_generate_hvmf_sample <- function(data_file = NULL,
                                         data_dir = default_hvmf_calibration_data_dir(),
                                         kappa = 50,
                                         n = 50,
                                         sample_id = 1L,
                                         regenerate = FALSE,
                                         seed = 123,
                                         mu = hvmf_typeiv_fixed_mu()) {
  resolved_path <- NULL
  sample_source <- "file"
  file_error <- NULL

  resolved_path <- tryCatch(
    find_hvmf_generation_file(
      data_file = data_file,
      data_dir = data_dir,
      kappa = kappa,
      n = n,
      sample_id = sample_id
    ),
    error = function(e) {
      file_error <<- conditionMessage(e)
      NULL
    }
  )

  if (!is.null(resolved_path)) {
    message(sprintf("Loaded HvMF dataset: %s", resolved_path))
    dat <- utils::read.csv(resolved_path)
    xyz <- as.matrix(dat[, c("V1", "V2", "V3"), drop = FALSE])
    xyz <- normalize_hvmf_h2_data(xyz)
    metadata <- extract_hvmf_metadata_from_path(resolved_path)
    if (is.finite(metadata$kappa)) {
      kappa <- metadata$kappa
    }
    if (is.finite(metadata$n)) {
      n <- metadata$n
    }
  } else {
    if (!isTRUE(regenerate)) {
      stop(file_error %||% paste(
        "No existing HvMF dataset was found.",
        "Set `regenerate=TRUE` to allow a small fallback simulation."
      ))
    }
    set.seed(as.integer(seed))
    xyz <- rhvmf_h2_gig(n = as.integer(n), mu = mu, kappa = kappa, check = TRUE)
    sample_source <- "generated"
    message(sprintf(
      "Generated fallback HvMF sample with n = %s, kappa = %s, seed = %s.",
      as.integer(n),
      format(kappa, scientific = FALSE, trim = TRUE),
      as.integer(seed)
    ))
  }

  list(
    data = xyz,
    source = sample_source,
    data_file = resolved_path,
    mu = as.numeric(mu),
    kappa = as.numeric(kappa),
    n = as.integer(n),
    sample_id = as.integer(sample_id)
  )
}

compute_hvmf_mu_projection <- function(data, mu) {
  x <- normalize_hvmf_h2_data(data)
  mu <- as.numeric(normalize_hvmf_h2_data(mu)[1L, , drop = TRUE])
  y <- as.numeric(x %*% c(mu[[1L]], -mu[[2L]], -mu[[3L]]))
  pmax(y, 1)
}

build_default_hvmf_omega_grid <- function(data,
                                          mu,
                                          n_omega = 6L) {
  omega_grid_obj <- make_hvmf_ks_grid(
    data = data,
    mu = mu,
    n_omega = max(2L, as.integer(n_omega)),
    n_t = 2L
  )
  omega_grid <- omega_grid_obj$omega_grid
  mu <- as.numeric(normalize_hvmf_h2_data(mu)[1L, , drop = TRUE])

  has_mu <- apply(omega_grid, 1L, function(omega) {
    max(abs(as.numeric(omega) - mu)) < 1e-10
  })

  if (!any(has_mu)) {
    omega_grid <- rbind(mu, omega_grid)
  } else if (!isTRUE(has_mu[[1L]])) {
    omega_grid <- rbind(mu, omega_grid[!has_mu, , drop = FALSE])
  }

  omega_grid[seq_len(min(nrow(omega_grid), as.integer(n_omega))), , drop = FALSE]
}

compute_single_hvmf_projection_panel <- function(omega,
                                                 data,
                                                 mu,
                                                 kappa,
                                                 density_grid_size = 1000L) {
  omega <- as.numeric(normalize_hvmf_h2_data(omega)[1L, , drop = TRUE])
  mu <- as.numeric(normalize_hvmf_h2_data(mu)[1L, , drop = TRUE])
  x <- normalize_hvmf_h2_data(data)

  alpha <- as.numeric(-minkowski_inner_product(mu, omega))
  alpha <- max(alpha, 1)
  y <- as.numeric(x %*% c(omega[[1L]], -omega[[2L]], -omega[[3L]]))
  y <- pmax(y, 1)
  y <- y[is.finite(y)]
  if (length(y) == 0L) {
    stop("The projected sample contains no finite values in [1, infinity).")
  }

  y_max <- max(y)
  y_upper <- max(y_max * 1.02, 1.1)
  y_grid <- seq(1, y_upper, length.out = as.integer(density_grid_size))
  density_grid <- hvmf_projection_density_h2(y_grid, alpha = alpha, kappa = kappa)

  list(
    omega = omega,
    alpha = alpha,
    y = y,
    y_grid = y_grid,
    density_grid = density_grid
  )
}

compute_hvmf_projection_panels <- function(data,
                                           omega_grid,
                                           mu,
                                           kappa,
                                           n_cores = 12L,
                                           density_grid_size = 1000L) {
  omega_grid <- normalize_hvmf_h2_data(omega_grid)
  n_cores <- max(1L, min(as.integer(n_cores), 12L, nrow(omega_grid)))

  panels <- parallel::mclapply(
    X = seq_len(nrow(omega_grid)),
    FUN = function(i) {
      compute_single_hvmf_projection_panel(
        omega = omega_grid[i, , drop = FALSE],
        data = data,
        mu = mu,
        kappa = kappa,
        density_grid_size = density_grid_size
      )
    },
    mc.cores = n_cores
  )

  names(panels) <- sprintf("omega_%02d", seq_along(panels))
  panels
}

plot_hvmf_generation_check <- function(panel_data,
                                       kappa,
                                       output_file,
                                       dataset_label,
                                       source_label,
                                       hist_breaks = 10) {
  if (length(panel_data) == 0L) {
    stop("`panel_data` must contain at least one omega panel.")
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  n_panels <- length(panel_data)
  n_col <- if (n_panels <= 4L) 2L else 3L
  n_row <- ceiling(n_panels / n_col)
  hist_breaks <- max(20L, as.integer(hist_breaks))

  grDevices::pdf(output_file, width = 4.8 * n_col, height = 3.8 * n_row)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(n_row, n_col),
    mar = c(4.2, 4.2, 3.2, 1.2),
    oma = c(0, 0, 3, 0)
  )

  for (panel in panel_data) {
    breaks <- seq(1, max(panel$y), length.out = hist_breaks + 1L)
    hist(
      panel$y,
      breaks = breaks,
      freq = FALSE,
      col = "#d9e6f2",
      border = "white",
      xlab = "Y = -<omega, X>_M",
      ylab = "Density",
      main = sprintf("alpha = %.3f", panel$alpha)
    )
    lines(panel$y_grid, panel$density_grid, lwd = 2.2, col = "#c23b23")
    grid(col = "#d9d9d9")
  }

  mtext(
    sprintf(
      "%s | HvMF generation check (%s), kappa = %s",
      dataset_label,
      source_label,
      format(kappa, scientific = FALSE, trim = TRUE)
    ),
    outer = TRUE,
    cex = 1.1,
    line = 1
  )
  grDevices::dev.off()

  invisible(output_file)
}

run_hvmf_generation_check <- function(data_file = NULL,
                                      data_dir = default_hvmf_calibration_data_dir(),
                                      output_file = file.path(
                                        "output",
                                        "hvmf_generation_check",
                                        "hvmf_typeiv_generation_check.pdf"
                                      ),
                                      kappa = 50,
                                      n = 50,
                                      sample_id = 1L,
                                      regenerate = FALSE,
                                      seed = 123,
                                      n_omega = 6L,
                                      n_cores = 12L,
                                      hist_breaks = 45L) {
  mu <- hvmf_typeiv_fixed_mu()
  sample_obj <- load_or_generate_hvmf_sample(
    data_file = data_file,
    data_dir = data_dir,
    kappa = kappa,
    n = n,
    sample_id = sample_id,
    regenerate = regenerate,
    seed = seed,
    mu = mu
  )

  dataset_label <- if (!is.null(sample_obj$data_file)) {
    basename(sample_obj$data_file)
  } else {
    sprintf(
      "generated_sample_N_%s_kappa%s",
      sample_obj$n,
      format(sample_obj$kappa, scientific = FALSE, trim = TRUE)
    )
  }

  omega_grid <- build_default_hvmf_omega_grid(
    data = sample_obj$data,
    mu = sample_obj$mu,
    n_omega = n_omega
  )
  panel_data <- compute_hvmf_projection_panels(
    data = sample_obj$data,
    omega_grid = omega_grid,
    mu = sample_obj$mu,
    kappa = sample_obj$kappa,
    n_cores = n_cores
  )

  plot_hvmf_generation_check(
    panel_data = panel_data,
    kappa = sample_obj$kappa,
    output_file = output_file,
    dataset_label = dataset_label,
    source_label = sample_obj$source,
    hist_breaks = hist_breaks
  )

  message(sprintf("Saved HvMF generation check plot: %s", normalizePath(output_file, winslash = "/", mustWork = TRUE)))

  invisible(list(
    output_file = normalizePath(output_file, winslash = "/", mustWork = TRUE),
    data_file = sample_obj$data_file,
    source = sample_obj$source,
    n = sample_obj$n,
    kappa = sample_obj$kappa,
    omega_grid = omega_grid,
    panels = panel_data
  ))
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_hvmf_generation_check(commandArgs(trailingOnly = TRUE))

  data_file <- if (!is.null(args$data_file) && nzchar(args$data_file)) {
    args$data_file
  } else {
    NULL
  }

  data_dir <- if (!is.null(args$data_dir) && nzchar(args$data_dir)) {
    args$data_dir
  } else {
    default_hvmf_calibration_data_dir()
  }

  output_file <- if (!is.null(args$output_file) && nzchar(args$output_file)) {
    args$output_file
  } else {
    file.path("output", "hvmf_generation_check", "hvmf_typeiv_generation_check.pdf")
  }

  kappa <- if (!is.null(args$kappa) && nzchar(args$kappa)) {
    suppressWarnings(as.numeric(args$kappa))
  } else {
    50
  }
  if (!is.finite(kappa) || kappa <= 0) {
    stop("`kappa` must be a strictly positive finite scalar.")
  }

  n <- parse_integer_hvmf_generation_check(args$n, default = 50L)
  sample_id <- parse_integer_hvmf_generation_check(args$sample_id, default = 1L)
  regenerate <- parse_flag_hvmf_generation_check(args$regenerate, default = FALSE)
  seed <- parse_integer_hvmf_generation_check(args$seed, default = 123L)
  n_omega <- parse_integer_hvmf_generation_check(args$n_omega, default = 6L)
  n_cores <- parse_integer_hvmf_generation_check(args$n_cores, default = 12L)
  hist_breaks <- parse_integer_hvmf_generation_check(args$hist_breaks, default = 45L)

  invisible(run_hvmf_generation_check(
    data_file = data_file,
    data_dir = data_dir,
    output_file = output_file,
    kappa = kappa,
    n = n,
    sample_id = sample_id,
    regenerate = regenerate,
    seed = seed,
    n_omega = n_omega,
    n_cores = n_cores,
    hist_breaks = hist_breaks
  ))
}
