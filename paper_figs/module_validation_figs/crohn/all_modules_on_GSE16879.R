setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/module_validation_figs/crohn")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/function_compute_all_features_for_module.R")


gse <- 'GSE16879'

#read all modules
source("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")

list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <-  c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                   "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/Transendothelial_migration", 
                   "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/Innate_immune_response", 
                   "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired_cell_proliferation", "Lipid metabolism")

mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")

########## get gene expression data
data <- getGEOData(gse)

########## get metadata
metadata <- GEOquery::getGEO(gse)
meta <- metadata$GSE16879_series_matrix.txt.gz@phenoData@data
meta <- sort_by_col(meta, colname = 'geo_accession')

#select after treatment
meta <- meta[grepl("Colon", meta$`tissue:ch1`), ]
meta <- meta[grepl("CO|after", meta$title), ]



counter_mod <- 0 #keep track of module number
list_plot <- list() # list of heatmaps plots


for (module in list_modules) { #for each module
  counter_mod <- counter_mod+1
  
  ########## get gene expression data
  expr <- data$originalData$GSE16879$expr
  map_probes <- (as.data.frame(data$originalData$GSE16879$keys))
  colnames(map_probes) <- 'keys'
  expr <- merge(expr, map_probes, by = 'row.names')
  expr$Row.names <- NULL
  expr <- na.omit(expr)
  expr <- mean_expression_gene_rowise(expr, colname = 'keys')
  expr <- order_by_col_names(expr)
  
  expr <- expr[, colnames(expr) %in% meta$geo_accession]
  
  if (!identical(meta$geo_accession, colnames(expr))) {
    stop("different samples from expr and metadata")
  }
  
  
  total_df <- compute_features(expr = expr, module = module, deconvolution = T)
  total_df <- remove_all_zeros_rows(total_df)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  rownames(total_df) <- gsub("KEGG_|HALLMARK_|WP_", "", rownames(total_df))
  
  
  
  #compute module scores
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  
  score_df$Condition <- meta$`disease:ch1`
  
  score_df$Condition <- factor(score_df$Condition, levels = c("Control", "CD", "UC"))
  # 
  #score_df$timepoint <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][1]}))
  #score_df$stimulation <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][2]}))
  #sscore_df$facet = factor(score_df$timepoint, levels = c("t0", "24h", "48h"))
  
  
  tableau_colors <- ggthemes::tableau_color_pal("Tableau 10")(10)  # Get the full palette of 10 colors
  # Subset and reverse the colors (e.g., use the first 5 colors and reverse them)
  custom_colors <- tableau_colors[6:10]
  
  #my_comparisons <- list(c("IFN-gamma", "unstimulated"))
  
  
  score_df$response <- meta$`response to infliximab:ch1`
  
  #score_df$response <- factor(score_df$response, levels= c('Control', 'Non responder', 'Responder'))
  
  box <- ggplot(score_df, aes(x = Condition, y = module, fill = response))+
    geom_boxplot() +
    theme_classic() +
    #scale_fill_tableau(palette = "Classic Green-Orange 12", name = 'Response') +
    scale_fill_manual(values = c('#add98e', 'bisque', '#ff800e'), name = 'Response', labels = c('F','control', 'T')) +
    stat_compare_means(aes(group = response), 
                       method = "t.test",
                       label = "p.format") +
    ggtitle(module_names[[counter_mod]]) +
    ylab("module") +
    xlab("")
  
  
  
  list_plot[[counter_mod]] <- box
  stop()
  
  
}
do.call(grid.arrange,c(list_plot, ncol = 5))
