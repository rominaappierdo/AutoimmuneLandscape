setwd("paper_figs/module_validation_figs/module_contribution_using_scores_auc/lupus_tissue_kidney/datasets")

gses <- c("GSE112943", "GSE113342", "GSE200306", "GSE32591")


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

rnaseq_gses <- c("GSE175759")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/lupus/downloaded_RNAseq/", pattern = 'rds', full.names = T, recursive = T)

for (gse in rnaseq_gses){
  file <- files[grepl(gse, files)]
  d <- readRDS(file)$expr
  write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  
}
