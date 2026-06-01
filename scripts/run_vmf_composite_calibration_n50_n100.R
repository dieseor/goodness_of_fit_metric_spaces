resolve_vmf_composite_calibration_path <- function(...) {
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

calibration_study_script_vmf <- resolve_vmf_composite_calibration_path(
  "bootstrap",
  "calibration_study.R"
)
source(calibration_study_script_vmf)

run_vmf_composite_calibration_n50_n100 <- function(
  output_root = file.path(
    "output",
    "bootstrap_calibration",
    "vmf_composite_M1000_B1000_n50_100"
  ),
  n_values = c(50L, 100L),
  M_outer = 1000L,
  B = 1000L,
  n_cores_outer = 12L,
  statistics = c("ks", "cvm"),
  alphas = c(0.01, 0.05, 0.10),
  seed = 20260602L,
  show_progress = TRUE,
  verbose = TRUE
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  scenarios_all <- default_bootstrap_composite_calibration_scenarios()
  scenarios <- scenarios_all[vapply(scenarios_all, function(s) identical(s$model, "vmf"), logical(1))]
  if (length(scenarios) == 0L) {
    stop("No vMF composite scenarios found in default_bootstrap_composite_calibration_scenarios().")
  }

  result <- run_bootstrap_calibration_study(
    scenarios = scenarios,
    n_values = as.integer(n_values),
    M_outer = as.integer(M_outer),
    B = as.integer(B),
    alpha_nominal = 0.05,
    alphas = as.numeric(alphas),
    statistics = statistics,
    n_cores_outer = as.integer(n_cores_outer),
    seed = as.integer(seed),
    output_dir = output_root,
    show_progress = isTRUE(show_progress),
    verbose = isTRUE(verbose)
  )

  saveRDS(result, file = file.path(output_root, "run_result.rds"))
  result
}

parse_named_args_vmf_composite <- function(args) {
  if (length(args) == 0L) {
    return(list())
  }

  args <- args[startsWith(args, "--")]
  if (length(args) == 0L) {
    return(list())
  }

  output <- vector("list", length(args))
  names(output) <- rep("", length(args))

  for (i in seq_along(args)) {
    arg <- substring(args[[i]], 3L)
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    output[[i]] <- value
    names(output)[[i]] <- key
  }

  output
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_vmf_composite(commandArgs(trailingOnly = TRUE))

  output_root <- if (!is.null(args$output_root)) {
    args$output_root
  } else {
    file.path("output", "bootstrap_calibration", "vmf_composite_M1000_B1000_n50_100")
  }

  n_values <- if (!is.null(args$n_values)) {
    as.integer(strsplit(args$n_values, ",", fixed = TRUE)[[1L]])
  } else {
    c(50L, 100L)
  }

  M_outer <- if (!is.null(args$M_outer)) as.integer(args$M_outer) else 1000L
  B <- if (!is.null(args$B)) as.integer(args$B) else 1000L
  n_cores_outer <- if (!is.null(args$n_cores_outer)) as.integer(args$n_cores_outer) else 12L
  seed <- if (!is.null(args$seed)) as.integer(args$seed) else 20260602L

  run_vmf_composite_calibration_n50_n100(
    output_root = output_root,
    n_values = n_values,
    M_outer = M_outer,
    B = B,
    n_cores_outer = n_cores_outer,
    seed = seed,
    show_progress = TRUE,
    verbose = TRUE
  )
}
