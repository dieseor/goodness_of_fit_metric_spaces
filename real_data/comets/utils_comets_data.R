load_comets_real_data <- function(finite_normals = c("none", "short", "long", "both"),
                                  warn_incomplete = FALSE,
                                  warn_nonfinite = FALSE) {
  finite_normals <- match.arg(finite_normals)

  if (!requireNamespace("sphunif", quietly = TRUE)) {
    stop("Package `sphunif` is required for the comet analyses.")
  }

  data("comets", package = "sphunif")
  comets$normal <- cbind(
    sin(comets$i) * sin(comets$om),
    -sin(comets$i) * cos(comets$om),
    cos(comets$i)
  )

  valid_rows <-
    !is.na(comets$class) &
    is.finite(comets$per_y) &
    is.finite(comets$i) &
    is.finite(comets$om) &
    !is.na(comets$frag)

  dropped_incomplete <- sum(!valid_rows)
  if (warn_incomplete && dropped_incomplete > 0L) {
    warning(
      sprintf(
        "Dropping %d comet rows with incomplete orbital elements before period filtering.",
        dropped_incomplete
      ),
      call. = FALSE
    )
  }

  comets_valid <- comets[valid_rows, , drop = FALSE]

  short_selector <-
    !(comets_valid$class %in% c("HYP", "PAR")) &
    comets_valid$per_y < 200 &
    !comets_valid$frag
  long_selector <-
    !(comets_valid$class %in% c("HYP", "PAR")) &
    comets_valid$per_y >= 200 &
    !comets_valid$frag

  comets_short <- comets_valid[short_selector, , drop = FALSE]
  comets_long <- comets_valid[long_selector, , drop = FALSE]

  filter_finite_rows <- function(df, label) {
    normal_matrix <- as.matrix(df$normal)
    finite_rows <- apply(normal_matrix, 1L, function(r) all(is.finite(r)))
    dropped_nonfinite <- sum(!finite_rows)
    if (warn_nonfinite && dropped_nonfinite > 0L) {
      warning(
        sprintf(
          "Dropping %d %s comet rows with non-finite normal coordinates after filtering.",
          dropped_nonfinite,
          label
        ),
        call. = FALSE
      )
    }
    df[finite_rows, , drop = FALSE]
  }

  if (finite_normals %in% c("short", "both")) {
    comets_short <- filter_finite_rows(comets_short, "short-period")
  }
  if (finite_normals %in% c("long", "both")) {
    comets_long <- filter_finite_rows(comets_long, "long-period")
  }

  list(
    raw = comets,
    valid = comets_valid,
    short = comets_short,
    long = comets_long
  )
}

write_comets_dataset_summary <- function(path, dataset_label, data_matrix) {
  utils::write.csv(
    data.frame(
      dataset = dataset_label,
      n = nrow(data_matrix),
      ambient_dim = ncol(data_matrix),
      stringsAsFactors = FALSE
    ),
    file = path,
    row.names = FALSE
  )
}
