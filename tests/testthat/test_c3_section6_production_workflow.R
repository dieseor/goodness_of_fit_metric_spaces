library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

run_command <- function(command, args = character(), env = character()) {
  if (length(env)) {
    env_names <- sub("=.*", "", env)
    env_values <- sub("^[^=]*=", "", env)
    old_values <- Sys.getenv(env_names, unset = NA_character_)
    do.call(Sys.setenv, stats::setNames(as.list(env_values), env_names))
    on.exit({
      for (index in seq_along(env_names)) {
        if (is.na(old_values[[index]])) {
          Sys.unsetenv(env_names[[index]])
        } else {
          do.call(
            Sys.setenv,
            stats::setNames(list(old_values[[index]]), env_names[[index]])
          )
        }
      }
    }, add = TRUE)
  }
  output <- system2(
    command,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  list(
    output = output,
    status = if (is.null(status)) 0L else as.integer(status)
  )
}

write_executable <- function(path, lines) {
  writeLines(lines, path)
  Sys.chmod(path, mode = "0755")
}

make_mock_cluster_commands <- function(mock_bin) {
  dir.create(mock_bin, recursive = TRUE)
  for (command in c("module", "R", "scontrol")) {
    write_executable(
      file.path(mock_bin, command),
      c("#!/bin/sh", "exit 0")
    )
  }
  write_executable(
    file.path(mock_bin, "date"),
    c(
      "#!/bin/sh",
      "if [ \"${1:-}\" = \"+%s\" ]; then",
      "  echo 1700000000",
      "else",
      "  echo 2026-07-29T12:00:00+02:00",
      "fi"
    )
  )
  write_executable(
    file.path(mock_bin, "Rscript"),
    c(
      "#!/bin/sh",
      "printf '%s\\n' 'CALL' \"$@\" >> \"$WORKFLOW_CAPTURE\"",
      "exit 0"
    )
  )
}

test_that("production sbatch keeps d=2,10 as its default and isolates d=5", {
  sbatch <- normalizePath("scripts/run_c3_section6_production.sbatch")
  submit_dir <- tempfile("section6-sbatch-submit-")
  mock_bin <- tempfile("section6-sbatch-bin-")
  capture <- tempfile("section6-sbatch-capture-")
  dir.create(submit_dir, recursive = TRUE)
  make_mock_cluster_commands(mock_bin)
  on.exit(unlink(
    c(submit_dir, mock_bin, capture),
    recursive = TRUE,
    force = TRUE
  ), add = TRUE)

  existing_d2 <- file.path(
    submit_dir,
    "simulation_results",
    "section6_new_scenarios",
    "final_normal_d2_d10_n50_100_200_400_M1000_B5000_quadrature"
  )
  dir.create(existing_d2, recursive = TRUE)
  sentinel <- file.path(existing_d2, "sentinel.txt")
  writeLines("do not modify", sentinel)
  sentinel_before <- tools::md5sum(sentinel)

  common_env <- c(
    sprintf("PATH=%s:%s", mock_bin, Sys.getenv("PATH")),
    sprintf("WORKFLOW_CAPTURE=%s", capture),
    sprintf("SLURM_SUBMIT_DIR=%s", submit_dir),
    "SLURM_CPUS_PER_TASK=32",
    "SLURM_JOB_ID=test-job",
    "SLURM_JOB_NAME=test-job",
    "SLURM_NTASKS=1",
    "SLURM_JOB_NUM_NODES=1",
    "SLURM_MEM_PER_NODE=6900"
  )

  default_run <- run_command(
    "/bin/bash",
    c(sbatch, "normal"),
    env = common_env
  )
  expect_identical(default_run$status, 0L)
  default_calls <- readLines(capture, warn = FALSE)
  expect_true("--dimensions=2,10" %in% default_calls)
  expect_true(paste0(
    "--output_dir=simulation_results/section6_new_scenarios/",
    "final_normal_d2_d10_n50_100_200_400_M1000_B5000_quadrature"
  ) %in% default_calls)
  expect_true("--derivative_method=quadrature" %in% default_calls)
  expect_identical(unname(tools::md5sum(sentinel)), unname(sentinel_before))

  writeLines(character(), capture)
  d5_run <- run_command(
    "/bin/bash",
    c(sbatch, "normal", "5"),
    env = common_env
  )
  expect_identical(d5_run$status, 0L)
  d5_calls <- readLines(capture, warn = FALSE)
  expect_true("--dimensions=5" %in% d5_calls)
  expect_true(paste0(
    "--output_dir=simulation_results/section6_new_scenarios/",
    "final_normal_d5_n50_100_200_400_M1000_B5000_quadrature"
  ) %in% d5_calls)
  expect_identical(unname(tools::md5sum(sentinel)), unname(sentinel_before))
})

test_that("production validator derives 48,000 and 24,000 expected rows", {
  validator <- "scripts/check_c3_section6_production.R"
  d2_output <- tempfile("section6-new-d2-")
  d5_output <- tempfile("section6-new-d5-")
  on.exit(unlink(c(d2_output, d5_output), recursive = TRUE, force = TRUE),
          add = TRUE)

  d2 <- run_command(
    "Rscript",
    c(
      "--vanilla",
      validator,
      "--mode=preflight",
      "--family=normal",
      paste0("--output_dir=", d2_output),
      "--seed=20260728",
      "--allow_new=true"
    )
  )
  expect_identical(d2$status, 0L)
  expect_true(any(grepl("completed: 0/48000", d2$output, fixed = TRUE)))
  expect_false(dir.exists(d2_output))

  d5 <- run_command(
    "Rscript",
    c(
      "--vanilla",
      validator,
      "--mode=preflight",
      "--family=normal",
      paste0("--output_dir=", d5_output),
      "--seed=20260728",
      "--dimensions=5",
      "--allow_new=true"
    )
  )
  expect_identical(d5$status, 0L)
  expect_true(any(grepl("completed: 0/24000", d5$output, fixed = TRUE)))
  expect_false(dir.exists(d5_output))
})

test_that("new-output preflight refuses a nonempty directory", {
  validator <- "scripts/check_c3_section6_production.R"
  output_dir <- tempfile("section6-nonempty-d5-")
  dir.create(output_dir)
  writeLines("existing", file.path(output_dir, "unrelated.txt"))
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  result <- suppressWarnings(run_command(
    "Rscript",
    c(
      "--vanilla",
      validator,
      "--mode=preflight",
      "--family=normal",
      paste0("--output_dir=", output_dir),
      "--seed=20260728",
      "--dimensions=5",
      "--allow_new=true"
    )
  ))
  expect_identical(result$status, 1L)
  expect_true(any(grepl(
    "Refusing to initialize it",
    result$output,
    fixed = TRUE
  )))
})
