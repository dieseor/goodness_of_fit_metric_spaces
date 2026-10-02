ns <- asNamespace("gofmetric")
for (name in ls(ns, all.names = TRUE)) {
  value <- get(name, envir = ns)
  if (is.function(value)) assign(name, value)
}
rm(ns, name, value)
