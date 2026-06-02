base <- 'output'
entries <- list.dirs(base, full.names = FALSE, recursive = FALSE)
entries <- entries[!startsWith(entries, '_')]

classify_area <- function(x) {
  lx <- tolower(x)
  if (grepl('bootstrap|calibration', lx)) return('calibration')
  if (grepl('screening|realdata|real_data|dataset', lx)) return('real_data')
  if (grepl('distance|geodesic|chordal|poor_mans|plots', lx)) return('distance_profiles')
  if (grepl('gaussian_process|convergence', lx)) return('convergence')
  if (grepl('diagnostic|validation|mle', lx)) return('diagnostics')
  if (grepl('comets', lx)) return('comets')
  'other'
}

df <- data.frame(
  path = file.path(base, entries),
  name = entries,
  area = vapply(entries, classify_area, character(1)),
  stringsAsFactors = FALSE
)

dir_size_mb <- function(path) {
  cmd <- sprintf("du -sm %s 2>/dev/null | awk '{print $1}'", shQuote(path))
  out <- suppressWarnings(system(cmd, intern = TRUE))
  if (length(out) == 0) return(NA_real_)
  as.numeric(out[[1]])
}

df$size_mb <- vapply(df$path, dir_size_mb, numeric(1))
df <- df[order(df$area, -df$size_mb, df$name), ]
write.csv(df, file = 'output/_catalog/output_index_20260602.csv', row.names = FALSE)

summary_df <- aggregate(size_mb ~ area, data = df, FUN = function(x) sum(x, na.rm = TRUE))
summary_df <- summary_df[order(-summary_df$size_mb), ]
write.csv(summary_df, file = 'output/_catalog/output_area_summary_20260602.csv', row.names = FALSE)

md <- c(
  '# Output Catalog (2026-06-02)',
  '',
  'Auto-generated after cleanup pass.',
  '',
  '## Area Summary (MB)',
  ''
)
for (i in seq_len(nrow(summary_df))) {
  md <- c(md, sprintf('- %s: %.0f MB', summary_df$area[i], summary_df$size_mb[i]))
}
md <- c(md, '', '## Top-Level Paths', '')
for (i in seq_len(nrow(df))) {
  mb <- if (is.na(df$size_mb[i])) 'NA' else sprintf('%.0f', df$size_mb[i])
  md <- c(md, sprintf('- %s | area=%s | size_mb=%s', df$path[i], df$area[i], mb))
}
writeLines(md, con = 'output/_catalog/README.md')

cat('Catalog generated:', nrow(df), 'entries\n')
