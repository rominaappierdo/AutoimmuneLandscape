setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/logistic_regression/psoriasis_tissue/datasets")

gses <- c("GSE109248", "GSE13355", "GSE145304", "GSE14905", "GSE166388", "GSE181318", "GSE52471", 
          "GSE53431", "GSE68923", "GSE68924", "GSE68937", "GSE68939", "GSE75343", "GSE75890", 
          "GSE78097", "GSE79704", "GSE80047", "GSE83582", "GSE114286", "GSE121212", "GSE178197", 
          "GSE66511", "GSE74697")


files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/psoriasis/downloaded_microarray/microarray_dataset_gene_symbol/", pattern = '_gs.tsv', full.names = T)

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

rnaseq_gses <- c("GSE114286",
                 "GSE121212",
                 "GSE178197",
                 "GSE66511",
                 "GSE74697")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/psoriasis//downloaded_RNAseq/", pattern = 'rds', full.names = T, recursive = T)

for (gse in rnaseq_gses){
  file <- files[grepl(gse, files)]
  d <- readRDS(file)$expr
  write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  
}
