files <- list.files(".", recursive = T, pattern = 'comple')

diseases <- c()
c <- 0
for (file in files) {
  c <- c+1
  disease <- strsplit(file, "/")[[1]][1]
  diseases <- append(diseases, disease)
  df <- read.csv(file, sep = '\t')
  df <- df[df$gene == 'T_cell', , drop = F]
  df <- df[, c('gene', 'effectSize', 'effectSizeStandardError')]
 
  if (c == 1) {
    merged <- df
  } else {
    merged <- rbind(merged, df)
  }
}

rownames(merged) <- diseases

fit <- rmeta::meta.summaries(merged$effectSize, merged$effectSizeStandardError, method = 'random')
fit$summary
fit$se.summary

pval     <- 2*stats::pnorm( abs(fit$summary/fit$se.summary), lower.tail=FALSE )
pval
