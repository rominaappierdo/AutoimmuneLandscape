setwd("paper_figs/module_validation_figs/module_contribution_using_scores_auc/jia_blood")
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/function_compute_all_features_for_module.R")
library(glmnet)
library(caret)
library(pROC)
library(ggpubr)
library(survcomp)
set.seed(123) #set seed for reproducibility ***

disease_to_select_regex <- 'health|contr|juvenile idiopathic arth'
alpha_value <- 0 #ridge regression
class_one_regex <- "juvenile idiopathic arth"
metadata <- read.csv("/mnt/data3/appierdo/PhD-project/Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/JIA_metadata.tsv", sep = '\t')
metadata <- metadata[metadata$timepoint %in% c(" ", "", 'baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"),]

gses <- readLines("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/actually_used_gses/jia/jia.gses.txt")
gses <- gses[!gses == 'GSE44314'] #removed error gsva

regex_gses <- paste(gses, collapse = "|")


mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")
module_names <- c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                  "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/transendothelial_migration", 
                  "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/innate_immune_response", 
                  "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired cell proliferation", "Lipid metabolism")



c_mod_name <- 0

list_for_computing_simpson_index <- list()

for (mod in mods) {
  c_mod_name <- c_mod_name+1
  print(module_names[[c_mod_name]])
  mod_interpretation <- module_names[c_mod_name]
  
  module_features <- readLines(paste0("../../../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files//", mod, "/", mod, ".txt"))
  module_features <- gsub('_tf|_cyto', '', module_features)
  
  files <- list.files("features_matrices/", recursive = T, pattern = paste0(".*", mod, '.tsv'), full.names = T) #all features datasets files
  files <- files[grepl(regex_gses, files)]
  
  module_pvals_on_all_datasets <- c()
  module_aucs_on_all_datasets <- c()
  for (file in files) {
    train_features <- read.csv(file, sep = '\t', row.names = 1)
    
    gse <- strsplit(file, "/|_|\\.")[[1]][4]
    
    sub_meta_training <- metadata[metadata$series_id == gse,] #metadata for current training dataset
    sub_meta_training <- sub_meta_training[grepl(disease_to_select_regex, sub_meta_training$disease),]  #select samples based on disease
    sub_meta_training <- sort_by_col(sub_meta_training, 'geo_accession')
    
    train_features <- train_features[, colnames(train_features) %in% sub_meta_training$geo_accession] #select samples from dataset based on disease_state of interest
    train_features <- order_by_col_names(train_features)
    
    train_features <- na.omit(train_features) #removing features that contains NA *******
    
    
    if(!identical(colnames(train_features), sub_meta_training$geo_accession)) {
      print("Training samples and metadata samples don't match")
      stop()
    }
    
    train_features_scores <- as.data.frame(colSums(train_features))
    colnames(train_features_scores) <- "module"
    anno_dis <- ifelse(grepl(class_one_regex, sub_meta_training$disease), 'case', 'control')
    train_features_scores$disease <- anno_dis
    
    labels <- ifelse(grepl("case", train_features_scores$disease), 1, 0)
    scores <- train_features_scores$module
    
    roc_object <- roc(labels, scores, direction  = '<')
    auc_value <- auc(roc_object)
    auc_value <- round(as.numeric(gsub("Area under the curve: ", "", auc_value)), digits = 3)
    
    
    p <- ggplot(train_features_scores, aes(x = disease, y = module, fill =disease)) +
      geom_boxplot() +
      theme_classic() +
      stat_compare_means()
    
    pval <- as.numeric(ggplot_build(p)$data[[2]]['p.adj'])
    
    module_pvals_on_all_datasets <- append(module_pvals_on_all_datasets, pval)
    module_aucs_on_all_datasets <- append(module_aucs_on_all_datasets, auc_value)
  }     
  list_for_computing_simpson_index[[mod_interpretation]] <- module_aucs_on_all_datasets
  
  int_pvals<- combine.test(module_pvals_on_all_datasets, method = 'fisher')
  mean_auc <- mean(module_aucs_on_all_datasets)
  
  module_row <- data.frame("module" =  module_names[[c_mod_name]], 
                           "integrated_pvals" = int_pvals, 
                           "mean_aucs" = mean_auc)
  
  
  if (c_mod_name == 1) {
    merged <- module_row
  } else {
    merged <- rbind(merged, module_row)
    
  }
  
}


saveRDS(merged, file = 'scores.rds')
saveRDS(list_for_computing_simpson_index, file = 'input_simpson_index_jia_blood.rds')

