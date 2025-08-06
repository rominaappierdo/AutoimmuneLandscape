setwd("Supplementary/supp_Files/cellDeconv_metaAnalysis_signatures/")

for (file in list.files('.', pattern = 'comple')) {
  print(file)
  df <- read.csv(file, sep = '\t')
  
  df2 <- df[, c('gene', 'numStudies', 'effectSize', 'effectSizeStandardError', 'effectSizePval', 'effectSizeFDR', 'heterogeneityPval', 'tauSquared', 'cochranesQ')]
  colnames(df2) <- c('cellType', 'numStudies', 'effectSize', 'effectSizeStandardError', 'effectSizePval', 'effectSizeFDR', 'heterogeneityPval', 'tauSquared', 'cochranesQ')
  df2 <- df2[order(df2$effectSize, decreasing = T), ]
  write.table(df2, file = file, sep = '\t', quote = F, col.names = T, row.names = F)
  
}