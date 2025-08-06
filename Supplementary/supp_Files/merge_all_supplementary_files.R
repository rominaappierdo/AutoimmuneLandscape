setwd("Supplementary/supp_Files")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")

files <- list.files(".", pattern = 'tsv', recursive = T, full.names = T)

c <- 0
for (file in files) {
  c <- c+1
  f <- read.csv(file, sep = '\t')
  feature_type <- colnames(f)[1]
  f <- rename_single_column(f, colnames(f)[1], 'Feature')
  f$Feature_type <- rep(feature_type, nrow(f))
  
  disease <- strsplit(file, "_|/|\\.")[[1]][6]
  tissue_type <- strsplit(file, "_|/|\\.")[[1]][9]
  
  f$Disease <- rep(disease, nrow(f))
  f$Tissue_type <- rep(tissue_type, nrow(f))
  
  if (c == 1) {
    merged <- f
  } else {
    merged <- rbind(merged, f)
  }
 
}


write.table(merged, "SupplementaryFile1.tsv", quote = F, sep = '\t', row.names = F, col.names = T)
