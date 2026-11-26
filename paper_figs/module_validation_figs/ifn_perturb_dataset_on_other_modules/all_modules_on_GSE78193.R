setwd("paper_figs/module_validation_figs/ifn_perturb_dataset_on_other_modules")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(gridExtra)
library(ggpubr)
library(mclust)
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/PhD-project/utils_functions/function_compute_all_features_for_module.R")


gse <- 'GSE78193'

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
meta <- GEOquery::getGEO(gse)
meta <- meta$`GSE78193-GPL6848_series_matrix.txt.gz`@phenoData@data
meta <- sort_by_col(meta, colname = 'geo_accession')
condition_annotation <- ifelse(grepl("Unstimulated", meta$title), 'unstimulated', 'IFN-gamma') #****
timepoint_annotation <- ifelse(grepl("Time 0", meta$title), 't0', ifelse(grepl("48", meta$title), '48h', '24h'))


counter_mod <- 0 #keep track of module number
list_plot <- list() # list of heatmaps plots


for (module in list_modules) { #for each module
  counter_mod <- counter_mod+1
  
  ########## get gene expression data
  expr <- data$originalData$GSE78193_GPL6848$expr
  map_probes <- (as.data.frame(data$originalData$GSE78193_GPL6848$keys))
  colnames(map_probes) <- 'keys'
  expr <- merge(expr, map_probes, by = 'row.names')
  expr$Row.names <- NULL
  expr <- na.omit(expr)
  expr <- mean_expression_gene_rowise(expr, colname = 'keys')
  expr <- order_by_col_names(expr)
  
  if (!identical(meta$geo_accession, colnames(expr))) {
    stop("different samples from expr and metadata")
  }
  
  
  total_df <- compute_features(expr = expr, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  rownames(total_df) <- gsub("KEGG_|HALLMARK_|WP_", "", rownames(total_df))
  
  
  #compute module scores
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  
  score_df$Condition <- paste0(timepoint_annotation, "_", condition_annotation)
  
  score_df$Condition <- factor(score_df$Condition, levels = c("t0_unstimulated", "t0_IFN-gamma",
                                                              "24h_unstimulated", "24h_IFN-gamma",
                                                              "48h_unstimulated", "48h_IFN-gamma"))
  
  score_df$timepoint <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][1]}))
  score_df$stimulation <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][2]}))
  score_df$facet = factor(score_df$timepoint, levels = c("t0", "24h", "48h"))
  
  
  tableau_colors <- ggthemes::tableau_color_pal("Tableau 10")(10)  # Get the full palette of 10 colors
  # Subset and reverse the colors (e.g., use the first 5 colors and reverse them)
  custom_colors <- tableau_colors[6:10]
  
  my_comparisons <- list(c("IFN-gamma", "unstimulated"))
  
  
  score_df$stimulation <- factor(score_df$stimulation, levels = c("unstimulated", "IFN-gamma"))
  
  box <- ggplot(score_df, aes(x = stimulation, y = module)) +
    geom_boxplot(aes(fill = stimulation)) +
    stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test") +
    ggtitle(module_names[counter_mod])+
    ylab("module score")+
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          strip.text = element_text(
            size = 14, color = "white", face = 'bold')
    ) +
    #scale_fill_tableau(direction = -1) +
    scale_fill_manual(values = custom_colors) +
    facet_wrap(~ facet)
  
  box
  
  
  g <- ggplot_gtable(ggplot_build(box))
  strip_both <- which(grepl('strip-', g$layout$name))
  fills <- c("#377eb8","#e41a1b",  "#4eaf4a")
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- fills[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    
    k <- k+1
  }
  grid.draw(g)
  
  list_plot[[counter_mod]] <- g
  

}
do.call(grid.arrange,c(list_plot, ncol = 5))
