setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE_sex_differences/meta_analysis_pathways_activity")

library(ggplot2)

files <- list.files(".", recursive = T, pattern = ".*results.tsv")

merged_df <- data.frame()
c <- 0
for (file in files) {
  c <- c+1
  
  pathway <- strsplit(file, "/")[[1]][1]
  disease <- strsplit(file, "/")[[1]][2]
  sample_source <- 'tissue'
  
  df <- read.csv(file, sep = '\t')
  df <- df[df$effectSizeFDR < 0.05, ]
  df <- df[abs(df$effectSize) > 0.4, ]
  df$directionality <- ifelse(df$effectSize > 0, 'UPREG_MALES', "UPREG_FEMALES")
  
  df <- df[, c(1, 15, 16)]
  df$gene_set <- rep(pathway, nrow(df))
  df$disease <- rep(disease, nrow(df))
  df$sample_source <- rep(sample_source, nrow(df))
  
  if (c == 1) {
    merged_df <- df
  } else {
    
    merged_df <- rbind(merged_df, df)
  }
  
  
}



ggplot(merged_df, aes(x = effectSize, y = reorder(gene, effectSize), fill = disease)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ disease, ncol = 3, scales = "free_y")
