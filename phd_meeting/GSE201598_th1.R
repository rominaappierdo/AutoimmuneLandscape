setwd("/mnt/data3/appierdo/AutoimmuneLandscape/phd_meeting")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")
library(factoextra)
library(ggpubr)
source("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/function_compute_all_features_for_module.R")


gse <- 'GSE201598'

df <- read.csv("GSE201598_raw_reads_quant_table_hTcells.txt", sep = '\t', row.names = 1)
df$names <- rownames(df)
df <- new_mapping_ids(df, var = 'refseq', colname = 'names')
df$names <- NULL
df <- na.omit(df)
df <- sum_expression_gene_rowise(df, 'mapped_symbol')


expr <- as.data.frame(conte(df, gse = gse)$expr)



pca.res <- prcomp(t(expr))
fviz_pca_ind(pca.res, habillage = score_df$condition, mean.point = F)+
  theme_classic()



list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <-  c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                   "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/Transendothelial_migration", 
                   "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/Innate_immune_response", 
                   "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired_cell_proliferation", "Lipid metabolism")

mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")

for (module in list_modules[11]) {
  
  total_df <- compute_features(expr = expr, module = module, deconvolution = T)
  
  total_df <- remove_all_zeros_rows(total_df)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  
  score_df$condition <- c("CTRL", "CTRL", "CTRL", "LDL", "LDL", "LDL", "VDL", "VDL", "VDL")
  
  my_comparisons <- list( c("CTRL", "LDL"), c("CTRL", "VDL"), c("LDL", "VDL") )
  
  box <- ggplot(score_df, aes(x = condition, y = module, fill = condition))+
    geom_boxplot() +
    theme_classic() +
    #scale_fill_tableau(palette = "Classic Green-Orange 12", name = 'Response') +
    #scale_fill_manual(values = c('#add98e', 'bisque', '#ff800e'), name = 'Response', labels = c('F','control', 'T')) +
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.format") +
    ggtitle('Th1 module') +
    ylab("module") +
    xlab("") +
    geom_jitter()
  
  box
  
}