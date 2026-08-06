# Optional compiled kernels for distance-profile evaluation.
#
# The project is script-based rather than an installed R package, so the
# compiled backend is loaded lazily with Rcpp::sourceCpp().  Nothing is
# compiled when the default R backend is used.

.distance_profile_cpp_state <- new.env(parent = emptyenv())
.distance_profile_cpp_state$loaded <- FALSE
.distance_profile_cpp_state$exports <- new.env(parent = globalenv())
.distance_profile_cpp_state$active_backend <- "r"

resolve_distance_profile_backend_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  matched <- candidates[file.exists(candidates)]
  if (length(matched) == 0L) {
    stop(sprintf("Could not resolve distance-profile backend path: %s", file.path(...)))
  }
  normalizePath(matched[[1L]], winslash = "/", mustWork = TRUE)
}

normalize_distance_profile_backend <- function(backend = c("r", "cpp")) {
  backend <- tolower(as.character(backend))
  if (length(backend) > 1L) {
    backend <- backend[[1L]]
  }
  if (length(backend) != 1L || is.na(backend) || !backend %in% c("r", "cpp")) {
    stop("`distance_profile_backend` must be either 'r' or 'cpp'.")
  }
  backend
}

distance_profile_backend_from_control <- function(control = list()) {
  value <- control$distance_profile_backend
  if (is.null(value)) value <- "r"
  normalize_distance_profile_backend(value)
}

distance_profile_backend_current <- function() {
  .distance_profile_cpp_state$active_backend
}

# Rcpp::sourceCpp() uses a cache directory shared by independent R sessions.
# Serialise compilation/loading for this project: without a lock, two sessions
# can observe the cache while its shared object is being replaced.
with_distance_profile_cpp_cache_lock <- function(cache_dir, code,
                                                 timeout_seconds = 1200,
                                                 stale_seconds = 900) {
  lock_dir <- file.path(cache_dir, ".distance_profile_sourcecpp.lock")
  deadline <- Sys.time() + timeout_seconds
  acquired <- FALSE
  repeat {
    if (isTRUE(dir.create(lock_dir, showWarnings = FALSE))) {
      acquired <- TRUE
      break
    }
    lock_info <- suppressWarnings(file.info(lock_dir))
    lock_age <- as.numeric(difftime(
      Sys.time(), lock_info$mtime[[1L]], units = "secs"
    ))
    if (is.finite(lock_age) && lock_age > stale_seconds) {
      unlink(lock_dir, recursive = TRUE, force = TRUE)
      next
    }
    if (Sys.time() >= deadline) {
      stop(sprintf(
        "Timed out after %d seconds waiting for the C++ distance-profile cache lock.",
        timeout_seconds
      ), call. = FALSE)
    }
    Sys.sleep(0.05)
  }
  on.exit({
    if (acquired) unlink(lock_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  force(code)
}

with_distance_profile_backend <- function(backend, expr) {
  backend <- normalize_distance_profile_backend(backend)
  previous <- .distance_profile_cpp_state$active_backend
  on.exit({
    .distance_profile_cpp_state$active_backend <- previous
  }, add = TRUE)
  if (identical(backend, "cpp")) {
    ensure_distance_profile_cpp_loaded()
  }
  .distance_profile_cpp_state$active_backend <- backend
  eval.parent(substitute(expr))
}

ensure_distance_profile_cpp_loaded <- function() {
  if (isTRUE(.distance_profile_cpp_state$loaded)) {
    return(invisible(.distance_profile_cpp_state$exports))
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("The C++ distance-profile backend requires the Rcpp package.")
  }

  source_file <- resolve_distance_profile_backend_path(
    "cpp",
    "distance_profile_backend.cpp"
  )
  cache_dir <- file.path(
    tools::R_user_dir("goodness_of_fit_metric_spaces", which = "cache"),
    "sourceCpp"
  )
  suppressWarnings(
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  )
  cache_probe <- if (dir.exists(cache_dir)) {
    tempfile("write-probe-", tmpdir = cache_dir)
  } else {
    NA_character_
  }
  cache_writable <- !is.na(cache_probe) &&
    suppressWarnings(file.create(cache_probe))
  if (isTRUE(cache_writable)) {
    unlink(cache_probe)
  } else {
    cache_dir <- file.path(
      tempdir(),
      "goodness_of_fit_metric_spaces",
      "sourceCpp"
    )
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  source_cpp <- function(rebuild) {
    Rcpp::sourceCpp(
      file = source_file,
      env = .distance_profile_cpp_state$exports,
      cacheDir = cache_dir,
      rebuild = rebuild,
      showOutput = FALSE,
      verbose = FALSE
    )
  }
  with_distance_profile_cpp_cache_lock(cache_dir, {
    first_attempt <- tryCatch(
      source_cpp(rebuild = FALSE),
      error = function(error) error
    )
    if (inherits(first_attempt, "error")) {
      tryCatch(
        source_cpp(rebuild = TRUE),
        error = function(error) {
          stop(sprintf(
            "Failed to compile or load the C++ distance-profile backend: %s",
            conditionMessage(error)
          ), call. = FALSE)
        }
      )
    }
  })

  .distance_profile_cpp_state$loaded <- TRUE
  invisible(.distance_profile_cpp_state$exports)
}

distance_profile_cpp_is_loaded <- function() {
  isTRUE(.distance_profile_cpp_state$loaded)
}

distance_profile_cpp_supports_spec <- function(spec_name) {
  spec_name <- as.character(spec_name)
  length(spec_name) == 1L && !is.na(spec_name) && (
    identical(spec_name, "normal") ||
      grepl("^small_circle_weighted_mixture2_(chordal|geodesic)$", spec_name) ||
      grepl("^sunspots_joint_time_space_(asymmetric|shared)$", spec_name)
  )
}

assert_distance_profile_cpp_spec_available <- function(spec_name) {
  spec_name <- as.character(spec_name)
  if (length(spec_name) == 1L && !is.na(spec_name) && grepl("^jp_", spec_name)) {
    stop(
      "The Jones-Pewsey model is intentionally excluded from the C++ distance-profile backend.",
      call. = FALSE
    )
  }
  if (!distance_profile_cpp_supports_spec(spec_name)) {
    stop(sprintf(
      "The C++ distance-profile backend was not retained for model '%s' because it did not pass the exactness and end-to-end performance gates.",
      paste(spec_name, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

distance_profile_cpp_call <- function(name, ...) {
  exports <- ensure_distance_profile_cpp_loaded()
  if (!exists(name, envir = exports, mode = "function", inherits = FALSE)) {
    stop(sprintf("C++ distance-profile kernel '%s' is not available.", name))
  }
  get(name, envir = exports, inherits = FALSE)(...)
}

make_distance_profile_backend_wrapper <- function(implementation,
                                                  cpp_supported = TRUE,
                                                  public_name = NULL) {
  stopifnot(is.function(implementation))
  wrapper_env <- new.env(parent = environment(implementation))
  wrapper_env$.distance_profile_implementation <- implementation
  wrapper_env$.distance_profile_cpp_supported <- isTRUE(cpp_supported)
  wrapper_env$.distance_profile_public_name <- public_name

  wrapper <- function() {
    backend_was_supplied <- !missing(backend)
    call <- match.call(expand.dots = TRUE)
    args <- as.list(call)[-1L]
    args$backend <- NULL
    args <- lapply(args, eval, envir = parent.frame())
    if (!backend_was_supplied &&
        !exists("distance_profile_backend_current", mode = "function", inherits = TRUE)) {
      return(do.call(.distance_profile_implementation, args))
    }
    backend_value <- if (backend_was_supplied) backend else distance_profile_backend_current()
    backend_value <- normalize_distance_profile_backend(backend_value)
    if (backend_was_supplied && identical(backend_value, "cpp") &&
        !.distance_profile_cpp_supported) {
      stop(sprintf(
        "The C++ distance-profile backend was not retained for '%s' because it did not pass the exactness and end-to-end performance gates.",
        .distance_profile_public_name
      ), call. = FALSE)
    }
    with_distance_profile_backend(
      backend_value,
      do.call(.distance_profile_implementation, args)
    )
  }
  formals(wrapper) <- c(
    formals(implementation),
    alist(backend = c("r", "cpp"))
  )
  environment(wrapper) <- wrapper_env
  attr(wrapper, "distance_profile_backend_wrapper") <- TRUE
  wrapper
}

install_distance_profile_backend_wrappers <- function(function_names,
                                                       envir = globalenv(),
                                                       cpp_supported = TRUE) {
  for (name in function_names) {
    if (!exists(name, envir = envir, mode = "function", inherits = FALSE)) next
    implementation <- get(name, envir = envir, inherits = FALSE)
    if (isTRUE(attr(implementation, "distance_profile_backend_wrapper"))) next
    assign(
      name,
      make_distance_profile_backend_wrapper(
        implementation = implementation,
        cpp_supported = cpp_supported,
        public_name = name
      ),
      envir = envir
    )
  }
  invisible(NULL)
}
