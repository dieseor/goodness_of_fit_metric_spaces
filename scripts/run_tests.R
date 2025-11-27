# Run all test scripts in 'tests' that follow naming convention test_*.R
scripts <- list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE)
if (length(scripts) == 0) {
  cat("No test scripts found in tests/ to run.\n")
} else {
  for (s in scripts) {
    cat("==== Running", s, "====\n")
    try(source(s))
  }
  cat("Done running tests (scripts).\n")
}
