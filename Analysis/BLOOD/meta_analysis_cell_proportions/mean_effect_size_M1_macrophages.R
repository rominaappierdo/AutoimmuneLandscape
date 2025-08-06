setwd("Analysis/BLOOD/meta_analysis_cell_proportions")

files <- list.files(".", full.names = T, pattern = 'complete', recursive = T)

es <- c()
for (file in files) {
  df <- read.csv(file, sep = '\t')
  #df <- df[df$effectSizeFDR < 0.05, ]
 
  m <- df[df$gene == 'macrophage_m0', , drop = F]
  if (nrow(m) == 0) {
    next
  }
  
  es <- append(es, m$effectSize)
  
}
mean(es)
print(sd(es))

plot(density(es))
