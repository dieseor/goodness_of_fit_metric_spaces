#!/usr/bin/env Rscript

resolve_simplex_contour_root <- function() {
  candidates <- c(".", "..")

  for (candidate in candidates) {
    utils_path <- file.path(
      candidate,
      "real_data",
      "logistic_gaussian",
      "utils_logistic_gaussian_screening.R"
    )
    if (file.exists(utils_path)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "Could not locate the project root containing real_data/logistic_gaussian/utils_logistic_gaussian_screening.R.",
    call. = FALSE
  )
}

parse_simplex_contour_args <- function(args) {
  parsed <- list()

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    parsed[[key]] <- value
  }

  parsed
}

project_logistic_gaussian_contour <- function(fit,
                                              probability,
                                              n_path = 401L) {
  angles <- seq(0, 2 * pi, length.out = as.integer(n_path))
  radius <- sqrt(stats::qchisq(probability, df = length(fit$mu_hat)))
  unit_circle <- cbind(cos(angles), sin(angles))

  eig <- eigen(fit$Sigma_hat, symmetric = TRUE)
  sqrt_sigma <- eig$vectors %*%
    diag(sqrt(pmax(eig$values, 0)), nrow = ncol(fit$Sigma_hat)) %*%
    t(eig$vectors)

  z_path <- t(
    matrix(fit$mu_hat, nrow = length(fit$mu_hat), ncol = length(angles)) +
      radius * sqrt_sigma %*% t(unit_circle)
  )
  x_path <- inverse_ilr_to_closed(z_path, V = fit$ilr_basis)
  projected <- project_simplex_to_ternary(x_path)

  data.frame(
    x = projected$x,
    y = projected$y,
    level = sprintf("%d%%", round(100 * probability)),
    stringsAsFactors = FALSE
  )
}

build_simplex_contour_plot <- function(dataset_name) {
  data_prep <- prepare_composition_dataset(dataset_name)
  if (!identical(data_prep$status %||% "ok", "ok")) {
    stop(sprintf("Dataset %s could not be prepared.", dataset_name), call. = FALSE)
  }
  if (data_prep$D != 3L) {
    stop(sprintf("Dataset %s does not have D = 3.", dataset_name), call. = FALSE)
  }
  if (isTRUE(data_prep$has_zeros)) {
    stop(sprintf("Dataset %s contains zeros; the standard logistic Gaussian fit was not attempted.", dataset_name), call. = FALSE)
  }

  fit <- fit_logistic_gaussian_plugin(data_prep$X_closed)

  observed <- project_simplex_to_ternary(data_prep$X_closed)
  contour_probs <- c(0.50, 0.80, 0.95)
  contour_data <- do.call(
    rbind,
    lapply(contour_probs, function(probability) {
      project_logistic_gaussian_contour(
        fit = fit,
        probability = probability
      )
    })
  )
  contour_data$level <- factor(
    contour_data$level,
    levels = c("50%", "80%", "95%")
  )

  simplex_boundary <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, sqrt(3) / 2, 0)
  )
  vertex_labels <- data.frame(
    x = c(-0.028, 1.028, 0.5),
    y = c(-0.026, -0.026, sqrt(3) / 2 + 0.028),
    label = data_prep$component_names,
    hjust = c(1, 0, 0.5),
    vjust = c(1, 1, 0)
  )

  # These are fitted logistic-Gaussian/Aitchison probability contours:
  # Gaussian ellipses in ilr coordinates transported back to the simplex,
  # not Euclidean density contours drawn directly on the triangle.
  plt <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = simplex_boundary,
      ggplot2::aes(x = x, y = y),
      colour = "grey20",
      linewidth = 0.7,
      lineend = "round"
    ) +
    ggplot2::geom_path(
      data = contour_data,
      ggplot2::aes(x = x, y = y, colour = level, group = level),
      linewidth = 1.1,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = observed,
      ggplot2::aes(x = x, y = y),
      colour = "#2b2b2b",
      size = 2.6,
      alpha = 0.88
    ) +
    ggplot2::geom_text(
      data = vertex_labels,
      ggplot2::aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
      colour = "grey20",
      size = 4.8
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "50%" = "#d95f02",
        "80%" = "#1b9e77",
        "95%" = "#7570b3"
      )
    ) +
    ggplot2::coord_equal(
      xlim = c(-0.05, 1.05),
      ylim = c(-0.045, sqrt(3) / 2 + 0.05),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        title = NULL,
        override.aes = list(linewidth = 1.3)
      )
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = c(0.84, 0.84),
      legend.justification = c(0, 1),
      legend.text = ggplot2::element_text(size = 15),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.86),
        colour = NA
      ),
      legend.margin = ggplot2::margin(2, 2, 2, 2),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(10, 28, 10, 10)
    )

  list(
    plot = plt,
    data_prep = data_prep,
    fit = fit
  )
}

run_simplex_contours <- function(args = commandArgs(trailingOnly = TRUE)) {
  root_dir <- resolve_simplex_contour_root()
  source(file.path(
    root_dir,
    "real_data",
    "logistic_gaussian",
    "utils_logistic_gaussian_screening.R"
  ))

  parsed_args <- parse_simplex_contour_args(args)
  dataset_name <- parsed_args$dataset %||% "Sediments"
  dataset_slug <- slugify_dataset_name(dataset_name)
  default_output_dir <- file.path(
    root_dir,
    "real_data",
    "logistic_gaussian",
    "screening",
    "fast",
    "isnps_20260621_simplex_contours",
    "plots"
  )
  output_dir <- parsed_args$output_dir %||% default_output_dir
  output_dir <- normalizePath(output_dir, mustWork = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  plot_result <- build_simplex_contour_plot(dataset_name = dataset_name)
  pdf_file <- file.path(output_dir, sprintf("%s_simplex_logistic_gaussian_contours.pdf", dataset_slug))
  png_file <- file.path(output_dir, sprintf("%s_simplex_logistic_gaussian_contours.png", dataset_slug))

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot_result$plot,
    device = grDevices::pdf,
    width = 7.2,
    height = 5.1,
    bg = "white",
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_result$plot,
    width = 7.2,
    height = 5.1,
    dpi = 300,
    bg = "white"
  )

  message(sprintf("Saved %s simplex contour plot to:", dataset_name))
  message("  ", pdf_file)
  message("  ", png_file)

  invisible(list(
    pdf = pdf_file,
    png = png_file,
    output_dir = output_dir,
    dataset_name = dataset_name,
    data_prep = plot_result$data_prep,
    fit = plot_result$fit
  ))
}

if (sys.nframe() == 0L) {
  run_simplex_contours()
}
