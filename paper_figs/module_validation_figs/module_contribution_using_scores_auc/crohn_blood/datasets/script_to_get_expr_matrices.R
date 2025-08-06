setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/logistic_regression/crohn_blood/datasets")

gses <- c("GSE100833",
          "GSE119600",
          "GSE126124",
          "GSE169568",
          "GSE86434",
          "GSE94648",
          "GSE112057",
          "GSE186507")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/crohn/downloaded_microarray/microarray_dataset_gene_symbol/", pattern = '_gs.tsv', full.names = T)

for (gse in gses){
  file <- files[grepl(gse, files)]
  if (length(file) == 0) {
    print(gse)
    next
  } else {
    d <- read.csv(file, sep = '\t', row.names = 1)
    write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  }
  
}

rnaseq_gses <- c("GSE112057","GSE186507")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/crohn/downloaded_RNAseq/", pattern = 'rds', full.names = T, recursive = T)

for (gse in rnaseq_gses){
  file <- files[grepl(gse, files)]
  d <- readRDS(file)$expr
  write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  
}
