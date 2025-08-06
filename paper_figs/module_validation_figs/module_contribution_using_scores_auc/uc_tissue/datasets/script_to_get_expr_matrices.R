setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/logistic_regression/uc_tissue/datasets")

gses <- c("GSE10191",
          "GSE105074",
          "GSE126124",
          "GSE134025",
          "GSE178753",
          "GSE22619",
          "GSE36807",
          "GSE37283",
          "GSE38713",
          "GSE42911",
          "GSE48634",
          "GSE51785",
          "GSE53306",
          "GSE65114",
          "GSE66407",
          "GSE6731",
          "GSE73661",
          "GSE75214",
          "GSE92415",
          "GSE9686",
          "GSE48958",
          "GSE16879",
          "GSE206285",
          "GSE111889",
          "GSE128682",
          "GSE137344",
          "GSE159008",
          "GSE179128",
          "GSE193677",
          "GSE83687"
)


files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/uc/downloaded_microarray/microarray_dataset_gene_symbol/", pattern = '_gs.tsv', full.names = T)

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

rnaseq_gses <- c("GSE111889","GSE128682",
                 "GSE137344","GSE159008", "GSE179128", "GSE193677","GSE83687", "GSE105074") #non ho il file microarray di GSE105074

files <- list.files("/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/uc/downloaded_RNAseq/", pattern = 'rds', full.names = T, recursive = T)

for (gse in rnaseq_gses){
  print(gse)
  file <- files[grepl(gse, files)]
  print(file)
  d <- readRDS(file)$expr
  write.table(d, file = paste0(gse, ".tsv"), row.names = T, col.names = NA, sep = '\t', quote = F)
  
}
