#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys

MODEL = Path(
    "real_data/sunspots/"
    "sunspots_cycle23_joint_time_models_parsimonious.R"
)
TEST = Path(
    "tests/testthat/"
    "test_sunspots_temporal_uniform_beta_gof.R"
)

SELECT_OLD = '''  converged <- Filter(is_converged, fits)
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1L), "value"))]]
  selected_converged <- is_converged(best)
'''

SELECT_NEW = '''  converged <- Filter(is_converged, fits)
  candidate_pool <- if (length(converged) > 0L) converged else fits
  best <- candidate_pool[[
    which.min(vapply(candidate_pool, `[[`, numeric(1L), "value"))
  ]]
  selected_converged <- is_converged(best)
'''

REFINE_OLD = '''    lapply(candidates, function(candidate) {
      try(
        stats::optim(
          candidate,
          objective,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          control = lbfgsb_control
        ),
        silent = TRUE
      )
    })
'''

REFINE_NEW = '''    refined <- lapply(candidates, function(candidate) {
      try(
        stats::optim(
          candidate,
          objective,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          control = lbfgsb_control
        ),
        silent = TRUE
      )
    })

    if (!inherits(exploratory, "try-error") &&
        is.list(exploratory) &&
        length(exploratory$value) == 1L &&
        is.finite(exploratory$value)) {
      c(list(exploratory), refined)
    } else {
      refined
    }
'''

TEST_APPEND = r'''

test_that("a converged exploratory fit survives failed refinements", {
  set.seed(20260816)
  eta <- temporal_uniform_beta_test_eta()
  s <- sample_sunspots_joint_parsimonious_time(
    350L,
    eta,
    time_model = "uniform_beta"
  )

  fit <- suppressWarnings(
    fit_sunspots_joint_parsimonious_time(
      s,
      time_model = "uniform_beta",
      control = list(
        parsimonious_time_n_starts = 4L,
        parsimonious_time_nelder_mead_control =
          list(maxit = 2500L, reltol = 1e-10),
        parsimonious_time_optim_control =
          list(maxit = 0L)
      )
    )
  )

  expect_true(is.finite(fit$loglik))
  expect_true(fit$selected_converged)
  expect_identical(fit$opt$convergence, 0L)
})
'''

def run(*args):
    return subprocess.run(
        list(args),
        check=True,
        text=True,
        capture_output=True,
    )

def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)

def main():
    if not Path(".git").exists():
        fail("Run this script from the repository root.")

    for path in (MODEL, TEST):
        if not path.exists():
            fail(f"Required file not found: {path}")

    model = MODEL.read_text(encoding="utf-8")
    test = TEST.read_text(encoding="utf-8")

    if "candidate_pool <- if (length(converged) > 0L)" in model:
        fail("The robust temporal-fit selection patch is already installed.")
    if model.count(SELECT_OLD) != 1:
        fail("Could not identify the fit-selection block uniquely.")
    if model.count(REFINE_OLD) != 1:
        fail("Could not identify the refinement block uniquely.")
    if "a converged exploratory fit survives failed refinements" in test:
        fail("The regression test is already present.")

    original_model = model
    original_test = test

    try:
        model = model.replace(SELECT_OLD, SELECT_NEW, 1)
        model = model.replace(REFINE_OLD, REFINE_NEW, 1)
        test = test.rstrip() + TEST_APPEND + "\n"

        MODEL.write_text(model, encoding="utf-8")
        TEST.write_text(test, encoding="utf-8")

        run("Rscript", "-e", f"parse(file={str(MODEL)!r})")
        run("Rscript", "-e", f"parse(file={str(TEST)!r})")
        run("git", "diff", "--check")
    except Exception:
        MODEL.write_text(original_model, encoding="utf-8")
        TEST.write_text(original_test, encoding="utf-8")
        print(
            "Patch failed; both files were restored.",
            file=sys.stderr,
        )
        raise

    print("Robust uniform+Beta temporal fitting patch applied.")
    print("No commit or push was performed.")
    print()
    print("Modified:")
    print(f"  {MODEL}")
    print(f"  {TEST}")
    print()
    print("Next command:")
    print(
        "  Rscript -e "
        "'testthat::test_file("
        '"tests/testthat/test_sunspots_temporal_uniform_beta_gof.R"'
        ")'"
    )

if __name__ == "__main__":
    main()
