setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Supplementary/supp_figures")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
source("../../utils_functions/function_compute_all_features_for_module.R")
source("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")
set.seed(1234)


gse<-'GSE78193'

########## set module
module_name <- "mod5"
module <- readLines("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod5/mod5.txt")
module_interpretation <- 'IFN-g module'

########## get expression data
data <- getGEOData(gse)
expr <- data$originalData$GSE78193_GPL6848$expr
map_probes <- (as.data.frame(data$originalData$GSE78193_GPL6848$keys))
colnames(map_probes) <- 'keys'
expr <- merge(expr, map_probes, by = 'row.names')
expr$Row.names <- NULL
expr <- na.omit(expr)
expr <- mean_expression_gene_rowise(expr, colname = 'keys')
expr <- order_by_col_names(expr)



########## get metadata
meta <- fread("../../Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv")
meta <- meta[meta$series_id==gse]
meta <- meta[meta$cell_type=='cultured sample', ]
meta <- sort_by_col(meta, colname = 'geo_accession')

if (!identical(meta$geo_accession, colnames(expr))) {
  stop("different samples from expression and metadata")
}

original_res <- readRDS("../original_ifn_mod5.rds")
#T0
final_dis_treat <- c()
final_dis_untreat <- c()
#module has 49 features
for (i in 1:2) {
  module <- c(mod1, mod2, mod3, mod4, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
  module <- sample(module)
  module <- sample(module, size = 49)
  
  
  ########## compute features 
  total_df <- compute_features(expr = expr, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  rownames(total_df) <- gsub("KEGG_|HALLMARK_|WP_", "", rownames(total_df))
  
  
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  score_df$Condition <- paste0(meta$timepoint, "_", meta$treatment)
  
  score_df$Condition <- factor(score_df$Condition, levels = c("t0_unstimulated", "t0_IFN-gamma",
                                                              "24h_unstimulated", "24h_IFN-gamma",
                                                              "48h_unstimulated", "48h_IFN-gamma"))
  
  score_df$timepoint <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][1]}))
  score_df$stimulation <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][2]}))
  score_df$facet = factor(score_df$timepoint, levels = c("t0", "24h", "48h"))
  score_df <- score_df[score_df$timepoint == 't0', ]
  
  original_res <- original_res[original_res$timepoint == 't0', ]
  mean_original_treated <- mean(original_res[original_res$stimulation == 'IFN-gamma', ]$module)
  mean_original_untreated <- mean(original_res[!original_res$stimulation == 'IFN-gamma', ]$module)
  
  score_df_treat <- score_df[score_df$stimulation == 'IFN-gamma', ]
  score_df_untreat <- score_df[!score_df$stimulation == 'IFN-gamma', ]
  
  final_dis_treat <- append(final_dis_treat, mean(score_df_treat$module))
  final_dis_untreat <- append(final_dis_untreat, mean(score_df_untreat$module))
  
  
}


final <- data.frame(matrix(nrow = 20, ncol = 0))
final$treatment_distribution <- c(final_dis_treat, final_dis_untreat)
final$treatment <- c(rep('treated', length(final_dis_treat)), rep('untreated', length(final_dis_untreat)))
final$original <- c(rep(mean_original_treated, length(final_dis_treat)), rep(mean_original_untreated, length(final_dis_untreat)))


ggplot(final, aes(x = treatment_distribution, fill= treatment)) +
  geom_density() +
  geom_vline(xintercept = mean_original_treated)+
  geom_vline(xintercept = mean_original_untreated) +
  ggtitle('t0')