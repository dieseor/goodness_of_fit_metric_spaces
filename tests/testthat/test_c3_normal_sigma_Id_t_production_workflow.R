library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

runner_env <- new.env(parent = globalenv())
sys.source("scripts/run_normal_sigma_Id_t_pilot.R", envir = runner_env)

test_that("each production invocation contains exactly 16,000 jobs", {
  for (d in c(2L, 5L)) {
    design <- runner_env$normal_sigma_Id_t_design(
      dimensions = d,
      n_values = c(50L, 100L, 200L, 400L),
      beta_values = c(0, 0.25, 0.5, 1),
      M = 1000L
    )
    expect_equal(nrow(design), 16000L)
    expect_setequal(unique(design$n), c(50L, 100L, 200L, 400L))
    expect_setequal(unique(design$beta), c(0, 0.25, 0.5, 1))
    expect_equal(as.integer(table(design$design_id)), rep(1000L, 16L))
  }
})

test_that("the C3 wrapper uses the requested resources and sequential pairings", {
  lines <- readLines(
    "scripts/run_c3_normal_sigma_Id_t_production.sbatch", warn = FALSE
  )
  expect_true(any(grepl("^#SBATCH --cpus-per-task=32$", lines)))
  expect_true(any(grepl("^#SBATCH --time=07:15:00$", lines)))
  expect_true(any(grepl("^#SBATCH --mem=6900M$", lines)))
  expect_true(any(grepl("--M=1000 --B=5000", lines, fixed = TRUE)))
  expect_true(any(grepl("--n_values=50,100,200,400", lines, fixed = TRUE)))
  expect_true(any(grepl("--beta_values=0,0.25,0.5,1", lines, fixed = TRUE)))
  expect_true(any(grepl("--derivative_mc_size=10000", lines, fixed = TRUE)))
  expect_true(any(grepl("--cvm_block_size=50", lines, fixed = TRUE)))
  expect_true(any(grepl("--show_progress=true", lines, fixed = TRUE)))
  first <- grep("run_invocation 2 3", lines, fixed = TRUE)
  second <- grep("run_invocation 5 6", lines, fixed = TRUE)
  expect_length(first, 1L)
  expect_length(second, 1L)
  expect_lt(first, second)
})

test_that("old spherical_normal identifiers are absent from the renamed workflow", {
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
