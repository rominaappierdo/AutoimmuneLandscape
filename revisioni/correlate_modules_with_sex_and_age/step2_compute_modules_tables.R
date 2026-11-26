setwd("/mnt/data3/appierdo/AutoimmuneLandscape/revisioni/correlate_modules_with_sex_and_age")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/PhD-project/utils_functions/function_df_all_metadata.R")
source("/mnt/data3/appierdo/PhD-project/utils_functions/function_compute_all_features_for_module.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")

meta <- complete_metadata() 
sub_meta <- meta[meta$gender %in% c("M", "F"), ]
sub_meta <- sub_meta[!sub_meta$age == "", ]
sub_meta <- sub_meta[!grepl("month|M", sub_meta$age), ]
sub_meta$age <- as.numeric(sub_meta$age)
sub_meta <- sub_meta[sub_meta$cell_type == "", ]
sub_meta <- sub_meta[sub_meta$timepoint %in% c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"), ]
gses_to_get <- unique(sub_meta$series_id)


list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <-  c("Th1_response",  "Monocyte_macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                   "Anti-inflammatory_response", "IFN_response", "IL6_trans_signaling_Transendothelial_migration", 
                   "Lymphoid_cell_development", "Endothelial_dysfunction", "TLR_NFKB response_Innate_immune_response", 
                   "B cell","Integrins", "ECM remodelling", "Inflammation_Wound healing", "Impaired_cell_proliferation", "Lipid_metabolism")

mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")

gses_to_get <- c(
  "GSE80342", "GSE74761", "GSE141646", "GSE139179", "GSE158952",
  "GSE46239", "GSE60424", "GSE9006", "GSE132607",
  "GSE93606", "GSE17590", "GSE116006",
  "GSE11908", "GSE200306", "GSE49454", "GSE65391", "GSE72326",
  "GSE72747", "GSE72754", "GSE88885", "GSE88886", "GSE99967",
  "GSE39088", "GSE19285", "GSE24427", "GSE26104", "GSE46280",
  "GSE182038", "GSE136757", "GSE53431", "GSE57225", "GSE57376",
  "GSE11827", "GSE126476", "GSE13026", "GSE151161", "GSE172188",
  "GSE198520", "GSE20690", "GSE24742", "GSE45867", "GSE47726",
  "GSE68689", "GSE74143", "GSE97165", "GSE130955", "GSE181228",
  "GSE181549", "GSE38713", "GSE112594", "GSE8650", "GSE113342"
)

gses_to_get <- c(
  "GSE9006", "GSE17590", "GSE49454","GSE65391",  "GSE72326",  "GSE72747",  "GSE72754",
  "GSE99967",  "GSE39088",  "GSE19285",  "GSE24427", 
  "GSE26104", "GSE46280", "GSE136757", "GSE53431", "GSE57225", 
  "GSE57376", "GSE11827", "GSE126476", "GSE13026", "GSE172188",
  "GSE198520", "GSE20690", "GSE24742", "GSE45867", "GSE47726",
  "GSE68689", "GSE74143", "GSE97165", "GSE130955", "GSE181228",
  "GSE181549", "GSE38713", "GSE112594", "GSE8650", "GSE113342", "GSE200306", "GSE151161"
)

effective_gse <- c()
for (gse in gses_to_get) {
  
  print(gse)
  
  if (dir.exists(gse)) {
    next
  }
 
  #get gene expr matrix
  
  files <- list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Data_Collection/downloaded_datasets/", recursive = T, full.names = T, pattern = gse)
  
  if (!length(files) == 0) {
    
    effective_gse <- append(effective_gse, gse)
    print(files)
    
    if (any(grepl('rds$', files))) {
      expr_file <- files[grepl('rds$', files)]
      expr <- readRDS(expr_file)
      expr <- expr$expr
      
    } else if (any(grepl('tsv$', files))) {
      expr_file <- files[grepl('tsv$', files)][1]
      expr <- read.csv(expr_file, sep = '\t', row.names = 1)
      
    } else {
      next
    }
    
    counter_mod <- 0
    for (module in list_modules) {
      counter_mod <- counter_mod+1
      total_df <- compute_features(expr = expr, module = module, deconvolution = T)
      #total_df <- t(scale(t(total_df))) 
      dir.create(paste0("modules_df/", gse))
      name_file <- paste0("modules_df/", gse, "/", gse, "_", module_names[counter_mod], ".tsv")
      
      write.table(total_df, file = name_file, sep = '\t', row.names = T, col.names = NA, quote = F)
  
      
    }
    
    
  }
 
}


print('DONE')
