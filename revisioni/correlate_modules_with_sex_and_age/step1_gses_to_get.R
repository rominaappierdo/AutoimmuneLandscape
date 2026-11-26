setwd("/mnt/data3/appierdo/AutoimmuneLandscape/revisioni/correlate_modules_with_sex_and_age")
source("/mnt/data3/appierdo/PhD-project/utils_functions/function_df_all_metadata.R")


meta <- complete_metadata() 
sub_meta <- meta[meta$gender %in% c("M", "F"), ]
sub_meta <- sub_meta[!sub_meta$age == "", ]
sub_meta <- sub_meta[!grepl("month|M", sub_meta$age), ]
sub_meta$age <- as.numeric(sub_meta$age)
sub_meta <- sub_meta[sub_meta$cell_type == "", ]
sub_meta <- sub_meta[sub_meta$timepoint %in% c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early"), ]
gses_to_get <- unique(sub_meta$series_id)
