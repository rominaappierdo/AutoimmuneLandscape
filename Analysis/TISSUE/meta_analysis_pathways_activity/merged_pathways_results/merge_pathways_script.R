setwd("Analysis/TISSUE/meta_analysis_pathways_activity/merged_pathways_results")


diseases <- c('crohn', 'psoriasis', 'ra', 'scleroderma', 'sjogren', 'uc')
#diseases <- c('KIDNEY_gsva_lupus', 'SKIN_gsva_lupus')

all_files <- list.files("../", pattern = 'complete', full.names = T, recursive = T)

for (disease in diseases) {
  disease_files <- all_files[grepl(disease, all_files)]
  
  df_disease_merged <- data.frame()
  c <- 0
  for (file in disease_files) {
    c <- c+1
    df <- read.csv(file, sep = '\t')
    if (c ==1) {
      df_disease_merged <- df
    } else {
      df_disease_merged <- rbind(df_disease_merged, df)
    }
  }
  
  file_name <- paste0(disease, '_complete_results_meged.tsv')
  
  write.table(df_disease_merged, file = file_name, sep = '\t', quote = F, col.names = T, row.names = F)
  
}