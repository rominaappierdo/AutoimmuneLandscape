setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/logistic_regression/lupus_blood/datasets")

gses <- c("GSE100163", "GSE110169", "GSE110396", "GSE112087","GSE154851", "GSE22098", "GSE39088",
          "GSE45291", "GSE49454", "GSE50635", "GSE61635", "GSE65391", "GSE72326", "GSE72754",
          "GSE80183", "GSE99967")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/lupus/downloaded_microarray/microarray_dataset_gene_symbol/", pattern = '_gs.tsv', full.names = T)

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

rnaseq_gses <- c( "GSE112087", "GSE80183")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/lupus/downloaded_RNAseq/", pattern = 'rds', full.names = T, recursive = T)

for (gse in rnaseq_gses){
  file <- files[grepl(gse, files)]
  d <- readRDS(file)$expr
  write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  
}
