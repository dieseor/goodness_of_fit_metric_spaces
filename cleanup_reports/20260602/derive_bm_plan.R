root <- 'output/bootstrap_calibration'
dirs <- list.dirs(root, full.names = FALSE, recursive = FALSE)

parse_mb <- function(name) {
  m <- regexec('_M([0-9]+)_B([0-9]+)', name)
  g <- regmatches(name, m)[[1]]
  if (length(g) == 3) return(c(M = as.integer(g[2]), B = as.integer(g[3])))
  m <- regexec('_MB([0-9]+)', name)
  g <- regmatches(name, m)[[1]]
  if (length(g) == 2) {
    v <- as.integer(g[2])
    return(c(M = v, B = v))
  }
  m <- regexec('_B([0-9]+)', name)
  g <- regmatches(name, m)[[1]]
  if (length(g) == 2) return(c(M = 0L, B = as.integer(g[2])))
  c(M = NA_integer_, B = NA_integer_)
}

norm_key <- function(name) {
  key <- name
  key <- gsub('_M[0-9]+_B[0-9]+', '', key)
  key <- gsub('_MB[0-9]+', '', key)
  key <- gsub('_B[0-9]+', '', key)
  key <- gsub('_(launch|run[0-9]+|clean|c[0-9]+|pilot|queue|onecore|debug|new)$', '', key)
  key <- gsub('__+', '_', key)
  key <- sub('_$', '', key)
  key
}

parsed <- lapply(dirs, function(d) {
  mb <- parse_mb(d)
  if (anyNA(mb)) return(NULL)
  data.frame(
    key = norm_key(d),
    name = d,
    M = as.integer(mb['M']),
    B = as.integer(mb['B']),
    qual = as.integer(mb['M'] >= 500 & mb['B'] >= 500),
    stringsAsFactors = FALSE
  )
})
parsed <- do.call(rbind, parsed)
if (is.null(parsed) || nrow(parsed) == 0) {
  write.table(data.frame(), 'cleanup_reports/20260602/bootstrap_bm_plan.tsv', sep='\t', row.names=FALSE, quote=FALSE)
  quit(save='no')
}

plan_rows <- list()
for (k in unique(parsed$key)) {
  g <- parsed[parsed$key == k, , drop = FALSE]
  ord <- order(-g$qual, -g$M, -g$B, g$name)
  g <- g[ord, , drop = FALSE]
  g$action <- 'DELETE'
  g$action[1] <- 'KEEP'
  plan_rows[[length(plan_rows) + 1L]] <- g[, c('action','key','name','M','B','qual')]
}
plan <- do.call(rbind, plan_rows)
plan <- plan[order(plan$action, plan$key, plan$name), ]
write.table(plan, 'cleanup_reports/20260602/bootstrap_bm_plan.tsv', sep='\t', row.names=FALSE, quote=FALSE)

cat('Resumen plan B/M:\n')
cat('KEEP', sum(plan$action == 'KEEP'), 'DELETE', sum(plan$action == 'DELETE'), '\n\n')
cat('DELETE propuestos (preview):\n')
print(head(plan[plan$action == 'DELETE', ], 60), row.names = FALSE)
