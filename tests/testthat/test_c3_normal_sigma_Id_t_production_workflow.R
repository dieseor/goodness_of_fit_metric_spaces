library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

runner_env <- new.env(parent = globalenv())
sys.source("scripts/run_normal_sigma_Id_t_pilot.R", envir = runner_env)

test_that("each beta-specific production invocation contains 4000 jobs", {
  for (d in c(2L, 5L)) {
    for (beta in c(0, 0.25, 0.5, 1)) {
      design <- runner_env$normal_sigma_Id_t_design(
        dimensions = d,
        n_values = c(50L, 100L, 200L, 400L),
        beta_values = beta,
        M = 1000L
      )
      expect_equal(nrow(design), 4000L)
      expect_equal(unique(design$beta), beta)
      expect_equal(as.integer(table(design$design_id)), rep(1000L, 4L))
    }
  }
})

test_that("design ids and seeds are stable when splitting by beta", {
  for (d in c(2L, 5L)) {
    full <- runner_env$normal_sigma_Id_t_design(
      dimensions = d,
      n_values = c(50L, 100L, 200L, 400L),
      beta_values = c(0, 0.25, 0.5, 1),
      M = 1L
    )

    for (beta in c(0, 0.25, 0.5, 1)) {
      one <- runner_env$normal_sigma_Id_t_design(
        dimensions = d,
        n_values = c(50L, 100L, 200L, 400L),
        beta_values = beta,
        M = 1L
      )
      reference <- full[full$beta == beta, ]

      expect_equal(one$design_id, reference$design_id)

      expect_equal(
        runner_env$normal_sigma_Id_t_seed(
          20260728L, one$design_id, one$replication, 0L
        ),
        runner_env$normal_sigma_Id_t_seed(
          20260728L, reference$design_id, reference$replication, 0L
        )
      )
    }
  }
})

test_that("the C3 wrapper uses beta-specific production resources", {
  lines <- readLines(
    "scripts/run_c3_normal_sigma_Id_t_production.sbatch",
    warn = FALSE
  )

  expect_true(any(grepl("^#SBATCH --cpus-per-task=16$", lines)))
  expect_true(any(grepl("^#SBATCH --time=01:00:00$", lines)))
  expect_true(any(grepl("^#SBATCH --mem=4G$", lines)))
  expect_true(any(grepl("--M=1000 --B=5000", lines, fixed = TRUE)))
  expect_true(any(grepl("--n_values=50,100,200,400", lines, fixed = TRUE)))
  expect_true(any(grepl("--beta_values=${beta_value}", lines, fixed = TRUE)))
  expect_true(any(grepl("--checkpoint_results=2000", lines, fixed = TRUE)))
  expect_true(any(grepl("--derivative_mc_size=10000", lines, fixed = TRUE)))
  expect_true(any(grepl("--cvm_block_size=50", lines, fixed = TRUE)))

  first <- grep('run_invocation 2 3 "$output_d2"', lines, fixed = TRUE)
  second <- grep('run_invocation 5 6 "$output_d5"', lines, fixed = TRUE)

  expect_length(first, 1L)
  expect_length(second, 1L)
  expect_lt(first, second)
})

test_that("old spherical_normal identifiers are absent", {
  paths <- c(
    "bootstrap/normal_sigma_Id_model_spec.R",
    "bootstrap/normal_sigma_Id_bootstrap.R",
    "scripts/run_normal_sigma_Id_t_pilot.R",
    "scripts/check_c3_normal_sigma_Id_t_production.R",
    "scripts/preflight_c3_normal_sigma_Id_t_production.sh",
    "scripts/run_c3_normal_sigma_Id_t_production.sbatch",
    "tests/testthat/test_normal_sigma_Id_model_spec.R"
  )

  contents <- unlist(lapply(paths, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("spherical_normal", contents, fixed = TRUE)))
})
