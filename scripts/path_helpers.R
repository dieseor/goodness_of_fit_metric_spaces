resolve_project_path_helper <- function(...) {
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

canonical_comets_root <- function(...) {
  file.path("real_data", "comets", ...)
}

canonical_comets_cardioid_dir <- function(run_name = NULL, speed = NULL, ...) {
  parts <- c(canonical_comets_root("cardioid"), if (!is.null(run_name)) run_name, if (!is.null(speed)) speed, list(...))
  do.call(file.path, as.list(parts))
}

canonical_comets_jp_dir <- function(run_name = NULL, speed = NULL, ...) {
  parts <- c(canonical_comets_root("jp"), if (!is.null(run_name)) run_name, if (!is.null(speed)) speed, list(...))
  do.call(file.path, as.list(parts))
}

canonical_comets_small_circle_dir <- function(run_name, speed = NULL, ...) {
  parts <- c(canonical_comets_root("small_circle", run_name), if (!is.null(speed)) speed, list(...))
  do.call(file.path, as.list(parts))
}

canonical_comets_spherical_cauchy_dir <- function(run_name = NULL, speed = NULL, ...) {
  parts <- c(canonical_comets_root("spherical_cauchy"), if (!is.null(run_name)) run_name, if (!is.null(speed)) speed, list(...))
  do.call(file.path, as.list(parts))
}

canonical_comets_mixture_dir <- function(run_name = NULL, speed = NULL, ...) {
  parts <- c(canonical_comets_root("mixture"), if (!is.null(run_name)) run_name, if (!is.null(speed)) speed, list(...))
  do.call(file.path, as.list(parts))
}

canonical_comets_logs_dir <- function(...) {
  file.path(canonical_comets_root("logs"), ...)
}

canonical_logistic_gaussian_root <- function(...) {
  file.path("real_data", "logistic_gaussian", ...)
}

canonical_logistic_gaussian_dataset_downloads_dir <- function(...) {
  file.path(canonical_logistic_gaussian_root("dataset_downloads"), ...)
}

canonical_logistic_gaussian_screening_dir <- function(speed = c("slow", "fast"), run_name = NULL, ...) {
  speed <- match.arg(speed)
  parts <- c(canonical_logistic_gaussian_root("screening", speed), if (!is.null(run_name)) run_name, list(...))
  do.call(file.path, as.list(parts))
}

canonical_calibration_bootstrap_dir <- function(model, speed = NULL, ...) {
  parts <- c("output", "calibration", "bootstrap", model, if (!is.null(speed)) speed, list(...))
  do.call(file.path, as.list(parts))
}

canonical_fast_multiplier_validation_dir <- function(...) {
  file.path("output", "validation", "fast_multiplier", ...)
}

canonical_fast_multiplier_diagnostics_dir <- function(...) {
  file.path("output", "diagnostics", "fast_multiplier", ...)
}

canonical_fast_multiplier_logs_dir <- function(...) {
  file.path("output", "logs", "fast_multiplier", ...)
}

canonical_wind_dir <- function(...) {
  file.path("real_data", "wind", ...)
}

canonical_wind_screening_dir <- function(speed = c("slow", "fast"), ...) {
  speed <- match.arg(speed)
  file.path(canonical_wind_dir("month_diagnostics", "screening_125m_ks_cvm_b1000", speed), ...)
}

canonical_sunspots_dir <- function(...) {
  file.path("real_data", "sunspots", ...)
}

canonical_sunspots_weighted_windows_dir <- function(speed = c("slow", "fast"), ...) {
  speed <- match.arg(speed)
  file.path(canonical_sunspots_dir("output", "cycles21_23_weighted_mixture_rolling_windows_10cr_B1000", speed), ...)
}
