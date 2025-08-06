setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/logistic_regression/sjogren_tissue/datasets")

gses <- c("GSE127952",
          "GSE23117",
          "GSE40568",
          "GSE80805",
          "GSE154926",
          "GSE159574")


files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/sjogren/downloaded_microarray/microarray_dataset_gene_symbol/", pattern = '_gs.tsv', full.names = T)

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



rnaseq_gses <- c("GSE154926","GSE159574")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/sjogren/downloaded_RNAseq/", pattern = 'rds', full.names = T, recursive = T)

for (gse in rnaseq_gses){
  file <- files[grepl(gse, files)]
  d <- readRDS(file)$expr
  write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  
}
