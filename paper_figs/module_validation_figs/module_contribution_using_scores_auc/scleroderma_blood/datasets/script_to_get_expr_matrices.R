setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/logistic_regression/scleroderma_blood//datasets")

gses <- c("GSE130953",
          "GSE45536")

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/scleroderma//downloaded_microarray/microarray_dataset_gene_symbol/", pattern = '_gs.tsv', full.names = T)

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
