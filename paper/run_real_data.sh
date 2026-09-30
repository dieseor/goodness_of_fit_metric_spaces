#!/usr/bin/env bash
# Real-data GOF tables of the paper (sunspots, comets, wind, compositions),
# with the settings and seeds of the runs used in the paper. Those runs were
# made in three steps: this script up to step 6 (2026-08-31), the sunspot run
# with Nderiv = 10000 (2026-09-08) and steps 7-8 (2026-09-08 and 2026-09-16).
# The wind step needs real_data/wind/risoe_m_all.nc from the DTU Risø archive.
set -euo pipefail

cd "$(dirname "$0")/.."
git rev-parse --is-inside-work-tree >/dev/null

# Avoid hidden BLAS/OpenMP parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

RUN_TAG="paper_main_realdata_B1000_3cores_$(date +%Y%m%d_%H%M%S)"
OUT="real_data/reruns/${RUN_TAG}"
mkdir -p "$OUT/logs"

run() {
  local label="$1"
  shift
  printf '\n===== %s =====\n' "$label"
  "$@" 2>&1 | tee "$OUT/logs/${label}.log"
}

# 1. Sunspots: cycle 23, joint time-space model.
run sunspots Rscript --vanilla real_data/sunspots/run_sunspots_cycle23_joint_time_space_gof.R \
  --input_csv=real_data/sunspots/output/sunspots_cycle23_s2_all.csv \
  --statistics=ks,cvm \
  --hemisphere_regression=shared \
  --start_date=1996-08-06 \
  --end_date=2008-11-28 \
  --n_sample_centers=all \
  --B=1000 \
  --derivative_mc_size=10000 \
  --n_cores=3 \
  --observed_profile_n_cores=3 \
  --center_seed=20260711 \
  --dequantization_seed=20260712 \
  --derivative_mc_seed=20260713 \
  --bootstrap_seed=20260714 \
  --time_quad_n=64 \
  --allow_boundary_fast=true \
  --output_dir="$OUT/sunspots"

# 2. Comets: cardioid C2, long and short periods, KS and CvM.
# 20260526 keeps the effective C2 seeds of the original run, where C2 was the
# third model in the list.
export OUT_CARDIOID="$OUT/comets/cardioid_C2"
run comets_cardioid_C2 Rscript -e '
source("scripts/run_comets_distance_profile_cardioid.R")
run_comets_distance_profile_cardioid(
  output_root = Sys.getenv("OUT_CARDIOID"),
  stages = c("oort_cvm", "oort_ks", "short_cvm", "short_ks"),
  model_ids = "C2",
  cvm_B = 1000L,
  ks_B = 1000L,
  ks_grid_mode = "sample_points_unique_distances",
  n_cores = 3L,
  base_seed = 20260526L,
  bootstrap_method = "fast_multiplier",
  distance_type = "geodesic",
  control = list(cardioid_optim_control = list(maxit = 1000L))
)
'

# 3. Comets: small circle. Only the automatic one-core retry is replaced by an
# error, so the run never continues with fewer than three workers.
export OUT_SMALL_CIRCLE="$OUT/comets/small_circle"
run comets_small_circle Rscript -e '
source("scripts/run_comets_distance_profile_small_circle_benchmark.R")

run_small_circle_benchmark_stage <- function(
    data_matrix, statistic, B_value, n_cores, seed,
    bootstrap_method = "reestimated", distance_type, ks_grid_mode,
    manual_ks_omega_points, manual_ks_t_points, control) {
  multiplier_bootstrap_gof(
    data = data_matrix,
    spec = make_small_circle_spec(distance_type = distance_type),
    null = list(type = "composite"),
    statistics = statistic,
    ks_grid = if (identical(statistic, "ks")) {
      make_small_circle_ks_grid_comets(
        data_matrix = data_matrix,
        distance_type = distance_type,
        ks_grid_mode = ks_grid_mode,
        manual_ks_omega_points = manual_ks_omega_points,
        manual_ks_t_points = manual_ks_t_points
      )
    } else {
      NULL
    },
    B = B_value,
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    bootstrap_method = bootstrap_method,
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = control
  )
}

tasks <- data.frame(
  dataset = c("short", "short", "long", "long"),
  statistic = c("ks", "cvm", "ks", "cvm"),
  run_id = c("short_ks", "short_cvm", "long_ks", "long_cvm"),
  stringsAsFactors = FALSE
)

for (j in seq_len(nrow(tasks))) {
  task <- tasks[j, , drop = FALSE]
  run_comets_distance_profile_small_circle_benchmark(
    output_root = file.path(Sys.getenv("OUT_SMALL_CIRCLE"), task$run_id),
    dataset = task$dataset[[1L]],
    B_values = 1000L,
    statistic = task$statistic[[1L]],
    n_cores = 3L,
    ks_grid_mode = "sample_points_unique_distances",
    base_seed = 20260531L,
    bootstrap_method = "fast_multiplier",
    distance_type = "geodesic"
  )
}
'

# 4. Comets: uniform-beta mixture, short and long periods.
export OUT_UNIFORM_BETA="$OUT/comets/uniform_beta"
run comets_uniform_beta Rscript -e '
source("scripts/run_comets_rotational_mixtures_short_long.R")
run_comets_mixtures_short_long(
  output_root = Sys.getenv("OUT_UNIFORM_BETA"),
  datasets = c("short", "long"),
  models = c("uniform_beta_mixture"),
  B = 1000L,
  statistics = c("ks", "cvm"),
  n_cores = 3L,
  seed = 20260710L,
  bootstrap_method = "fast_multiplier",
  distance_type = "geodesic",
  control_uniform = list(
    uniform_beta_mixture_profile_method = "legendre",
    uniform_beta_mixture_quad_n = 100L,
    uniform_beta_mixture_optim_control = list(maxit = 350L, reltol = 1e-9)
  )
)
'

# 5. Risø wind: the two windows in the table.
export OUT_WIND="$OUT/wind"
run wind Rscript -e '
source("real_data/wind/run_risoe_mjj77_plot_and_ndj_screening_mu_hat.R")

out <- Sys.getenv("OUT_WIND")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

selected_df <- select_noon_all_months(
  load_risoe_concurrent(
    file.path(repo_root, "real_data", "wind", "risoe_m_all.nc"),
    fixed_tz = "UTC"
  ),
  fixed_tz = "UTC"
)
configs <- make_configs()

results <- run_case(
  selected_df, configs[[1L]],
  B = 1000L, n_cores = 3L, seed = 2026080320L
)
utils::write.csv(results, file.path(out, "wind_B1000_checkpoint.csv"), row.names = FALSE)

results <- rbind(
  results,
  run_case(
    selected_df, configs[[2L]],
    B = 1000L, n_cores = 3L, seed = 2026080404L
  )
)
utils::write.csv(results, file.path(out, "wind_B1000_checkpoint.csv"), row.names = FALSE)
utils::write.csv(
  make_summary(results),
  file.path(out, "wind_B1000_summary.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
'

# 6. Compositions: the eight datasets in the table. Their positions in the
# original list of 29 are kept so each dataset keeps its seed.
export OUT_COMPOSITIONS="$OUT/compositions"
run compositions Rscript -e '
expressions <- parse(file = "scripts/run_logistic_gaussian_paper_defaults_mc.R")
runner_env <- new.env(parent = .GlobalEnv)
for (expr in expressions[-length(expressions)]) eval(expr, envir = runner_env)

fn <- function(name) get(name, envir = runner_env, inherits = TRUE)
all_items <- get("paper_datasets", envir = runner_env)
active <- c(
  "Activity31", "ArcticLake", "ClamEast", "HouseholdExp",
  "PogoJump", "Sediments", "WhiteCells_microscopic", "Yatquat_panel"
)
indices <- match(active, all_items$dataset)
if (anyNA(indices) || anyDuplicated(indices)) stop("Invalid active-dataset selection.")

out <- Sys.getenv("OUT_COMPOSITIONS")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
settings <- list(
  output_dir = out,
  B = 1000L,
  derivative_mc_size = 10000L,
  n_cores = 3L,
  seed = 20260803L
)

saveRDS(
  list(
    datasets = all_items[indices, , drop = FALSE],
    original_row_indices = indices,
    settings = settings,
    fixed_control = list(
      derivative_method = "score_mc",
      mvnormal_quadform_method = "auto",
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_multiplier_cache_corrections = "auto",
      compute_auxiliary_diagnostics = FALSE
    ),
    created_at = Sys.time(),
    R_version = R.version.string
  ),
  file.path(out, "run_manifest.rds")
)

rows <- vector("list", length(indices))
for (j in seq_along(indices)) {
  i <- indices[[j]]
  item <- all_items[i, , drop = FALSE]
  dataset_seed <- as.integer(settings$seed + 100L * (i - 1L))
  derivative_seed <- as.integer(settings$seed + 50000L + i)
  prepared <- fn("prepare_composition_dataset")(item$source_dataset)
  if (!identical(prepared$status, "ok")) stop(sprintf("Could not prepare %s.", item$source_dataset))

  result <- fn("run_logistic_gaussian_screening")(
    dataset_name = item$source_dataset,
    B = settings$B,
    max_centers = prepared$n,
    bootstrap_mode = "composite_multiplier",
    seed = dataset_seed,
    alpha = 0.05,
    ridge = 1e-8,
    n_cores = settings$n_cores,
    bootstrap_method = "fast_multiplier",
    bootstrap_keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    compute_auxiliary_diagnostics = FALSE,
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = settings$derivative_mc_size,
      derivative_mc_seed = derivative_seed,
      mvnormal_quadform_method = "auto",
      fast_multiplier_backend = "cpp",
      fast_multiplier_cpp_kernel = "contiguous_double",
      fast_multiplier_fuse_ks_cvm = TRUE,
      fast_multiplier_cache_corrections = "auto"
    ),
    omega_grid_type = "sample_points",
    t_grid_type = "sample_distances",
    make_plots = FALSE,
    save_outputs = TRUE,
    output_dir = settings$output_dir,
    run_seed_sensitivity = FALSE,
    verbose = TRUE
  )

  fn("validate_result")(
    result, item, settings, dataset_seed, derivative_seed,
    expected_n = prepared$n
  )
  rows[[j]] <- fn("result_row")(result, item, dataset_seed, derivative_seed)
  utils::write.csv(
    do.call(rbind, rows[seq_len(j)]),
    file.path(out, "ks_cvm_pvalues.csv"),
    row.names = FALSE
  )
}
'

# 7. Comets: joint C2 and SC tests with Nderiv = 10000, using the fits above.
run comets_joint Rscript --vanilla scripts/run_paper_comets_joint_score_mc.R \
  --B=1000 --n_cores=3 \
  --reference_dir="$OUT" \
  --output_dir="$OUT/comets/joint"

# 8. Comets: uniform-beta mixture with the joint kernel.
run comets_uniform_beta_joint Rscript --vanilla scripts/run_comets_uniform_beta_joint_kernel.R \
  --B=1000 --n_cores=2 \
  --output_root="$OUT/comets/uniform_beta_joint"

printf '\nDone. Results in: %s\n' "$OUT"
