source(file.path("bootstrap", "calibration_study.R"))

args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  match <- args[startsWith(args, prefix)]
  if (length(match) == 0L) {
    return(default)
  }
  sub(prefix, "", match[[1L]], fixed = TRUE)
}

output_dir <- get_arg_value("output_dir", file.path("output", "calibration", "bootstrap", "hvmf_simple_ks_cvm_full"))
n_cores_outer <- as.integer(get_arg_value("n_cores_outer", "1"))
seed <- as.integer(get_arg_value("seed", "123"))
M_outer <- as.integer(get_arg_value("M_outer", "1000"))
B <- as.integer(get_arg_value("B", "1000"))
mode <- tolower(get_arg_value("mode", "full"))
show_progress <- tolower(get_arg_value("show_progress", "false")) %in% c("true", "1", "yes")
verbose <- !tolower(get_arg_value("verbose", "true")) %in% c("false", "0", "no")

if (identical(mode, "smoke")) {
  res <- run_smoke_hvmf_simple_ks_cvm_calibration_study(
    output_dir = output_dir,
    n_cores_outer = n_cores_outer,
    seed = seed,
    show_progress = show_progress,
    verbose = verbose
  )
} else if (identical(mode, "full")) {
  res <- run_full_hvmf_simple_ks_cvm_calibration_study(
    output_dir = output_dir,
    n_cores_outer = n_cores_outer,
    seed = seed,
    M_outer = M_outer,
    B = B,
    show_progress = show_progress,
    verbose = verbose
  )
} else {
  stop("`mode` must be either 'smoke' or 'full'.")
}

summary(res)
