setwd("paper_figs/module_validation_figs/module_contribution_using_scores_auc/lupus_tissue_kidney/features_matrices")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/function_compute_all_features_for_module.R")
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")

mod1 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod1/mod1.txt")
mod2 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod2/mod2.txt")
mod3 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod3/mod3.txt")
mod4 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod4/mod4.txt")
mod5 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod5/mod5.txt")
mod6 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod6/mod6.txt")
mod7 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod7/mod7.txt")
mod8 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod8/mod8.txt")
mod9 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod9/mod9.txt")
mod10 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod10/mod10.txt")
mod11 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod11/mod11.txt")
mod12 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod12/mod12.txt")
mod13 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod13/mod13.txt")
mod14 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod14/mod14.txt")
mod15 <- readLines("../../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod15/mod15.txt")

list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <- c("Th1/Th2 balance",  "Macrophages", "TCR activation", "Autoimmune suppression", "IFN-g", "Transendothelial migration", "Cell cycle (celle expansion)", "Endothelial dysfunction", "TLR activation", "T cell differentiation","Regulation of hematopoiesis", "Wnt beta-catenin", "NF-kB and Apoptosis", "Inflammation", "Tregs")
mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")

metadata <- read.csv("/mnt/data3/appierdo/PhD-project/Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv", sep = '\t')
metadata <- metadata[metadata$timepoint %in% c(" ", "", 'baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),]
metadata <- metadata[!grepl("non|unaff", metadata$sample_source),]




gses <- c("GSE112943", "GSE113342", "GSE200306", "GSE32591", "GSE175759")


files <- list.files("../datasets/", recursive = T, pattern = 'tsv', full.names = T)


for (file in files) {
  expr <- read.csv(file, sep = '\t', row.names = 1)
  gse <- strsplit(file, "/|_|\\.")[[1]][6]

  if (gse %in% c('GSE113342', 'GSE200306')) {
    next
  }
  print(gse)
  sub_meta <- metadata[metadata$series_id == gse,]
  sub_meta <- sub_meta[grepl('health|contr|lupus', sub_meta$disease),]
  sub_meta <- sort_by_col(sub_meta, 'geo_accession')
  sub_meta <- sub_meta[grepl('kidn', sub_meta$sample_source),]
  expr <- expr[, colnames(expr) %in% sub_meta$geo_accession]
  expr <- order_by_col_names(expr)
  
  if(!identical(colnames(expr), sub_meta$geo_accession)) {
    print('different samples')
    stop()
  }
  
  counter_mod <- 0
  for (module in list_modules) {
    counter_mod <- counter_mod+1
    mod_name <- mods[counter_mod]
    features <- compute_features(expr = as.data.frame(expr), module = module)
    write.table(features, file = paste0(gse, "_features_", mod_name, ".tsv"), sep = '\t', row.names = T, col.names = NA, quote = F)
  }
  
}

print('done')
