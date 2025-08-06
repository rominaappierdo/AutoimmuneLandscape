setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/logistic_regression/crohn_tissue/datasets")

gses <- c("GSE102133", "GSE112366", "GSE126124", "GSE16879", "GSE207022", "GSE26305", "GSE36807", "GSE48634", 
                      "GSE52746", "GSE75214", "GSE66407", "GSE6731", "GSE68570", "GSE69762", "GSE77387", "GSE83448", 
                      "GSE95095", "GSE9686", "GSE101794", "GSE111889", "GSE117993", "GSE137344", "GSE179128", "GSE193677", 
                      "GSE208303", "GSE66207", "GSE83687", "GSE85499")


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

rnaseq_gses <- c("GSE101794", "GSE111889", "GSE117993", "GSE137344", "GSE179128", "GSE193677", 
                                        "GSE208303", "GSE66207", "GSE83687", "GSE85499")


files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/crohn/downloaded_RNAseq/", pattern = 'rds', full.names = T, recursive = T)

for (gse in rnaseq_gses){
  file <- files[grepl(gse, files)]
  d <- readRDS(file)$expr
  write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  
}
