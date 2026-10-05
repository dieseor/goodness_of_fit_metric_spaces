# Compatibility bindings for historical research scripts; implementation stays
# in the installed package. The final release pin is pending CRAN acceptance.
if (!requireNamespace("gofmetric", quietly = TRUE)) {
  stop("Install gofmetric before loading research scripts; see README.md.")
}
# Retain the dependency attachment expected by historical research scripts.
suppressPackageStartupMessages(
  for (research_pkg in c("parallel", "movMF", "sphunif", "pracma", "rotasym",
                         "mvtnorm", "gbutils", "ggplot2")) {
    library(research_pkg, character.only = TRUE)
  }
)
rm(research_pkg)
research_namespace <- asNamespace("gofmetric")
for (research_name in ls(research_namespace, all.names = TRUE)) {
  if (!startsWith(research_name, ".")) {
    assign(research_name, get(research_name, envir = research_namespace),
           envir = environment())
  }
}
research_roots <- c(".", "..", "../..")
research_root <- research_roots[file.exists(file.path(
  research_roots, "scripts", "research_only_helpers.R"))][1L]
if (is.na(research_root)) stop("Could not locate research_only_helpers.R.")
source(file.path(research_root, "scripts", "research_only_helpers.R"),
       local = environment())
rm(research_namespace, research_name, research_roots, research_root)
