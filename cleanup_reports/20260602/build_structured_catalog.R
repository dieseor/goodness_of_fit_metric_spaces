base <- 'output/structured'
if (!dir.exists(base)) stop('Missing output/structured')

areas <- list.dirs(base, full.names = FALSE, recursive = FALSE)
rows <- list()
for (a in sort(areas)) {
  subdirs <- list.dirs(file.path(base, a), full.names = FALSE, recursive = FALSE)
  for (nm in sort(subdirs)) {
    p <- file.path(base, a, nm)
    size <- suppressWarnings(as.numeric(system(sprintf("du -sm %s 2>/dev/null | awk '{print $1}'", shQuote(p)), intern = TRUE)[1]))
    rows[[length(rows)+1L]] <- data.frame(area=a, name=nm, canonical_path=p, size_mb=size, stringsAsFactors = FALSE)
  }
}

df <- do.call(rbind, rows)
if (is.null(df)) df <- data.frame(area=character(), name=character(), canonical_path=character(), size_mb=numeric())
write.csv(df, 'output/_catalog/structured_output_index_20260602.csv', row.names = FALSE)

agg <- aggregate(size_mb ~ area, data=df, FUN=function(x) sum(x, na.rm = TRUE))
agg <- agg[order(-agg$size_mb), , drop=FALSE]
write.csv(agg, 'output/_catalog/structured_output_area_summary_20260602.csv', row.names = FALSE)

legacy <- list.dirs('output', full.names = FALSE, recursive = FALSE)
legacy <- setdiff(legacy, c('_catalog','structured'))
legacy_links <- legacy[file.info(file.path('output', legacy))$isdir & !is.na(Sys.readlink(file.path('output', legacy))) ]

md <- c(
  '# Output Structure (Canonical)',
  '',
  'Canonical root: output/structured',
  '',
  'Legacy paths at output/<name> are symlink aliases for compatibility.',
  '',
  '## Areas',
  '- calibration',
  '- real_data',
  '- distance_profiles',
  '- convergence',
  '- diagnostics',
  '- comets',
  '- other',
  '',
  '## Area Sizes (MB)',
  ''
)
for (i in seq_len(nrow(agg))) {
  md <- c(md, sprintf('- %s: %.0f MB', agg$area[i], agg$size_mb[i]))
}
md <- c(md, '', sprintf('## Legacy alias count: %d', length(legacy_links)), '')
if (length(legacy_links) > 0) {
  md <- c(md, 'Aliases (old -> canonical):')
  for (nm in sort(legacy_links)) {
    target <- Sys.readlink(file.path('output', nm))
    md <- c(md, sprintf('- output/%s -> %s', nm, target))
  }
}
writeLines(md, 'output/_catalog/STRUCTURE_20260602.md')

cat('Structured entries:', nrow(df), '\n')
cat('Legacy aliases:', length(legacy_links), '\n')
