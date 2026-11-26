setwd("paper_figs/module_validation_figs/psoriasis_GSE249936")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(gridExtra)
library(ggpubr)
library(mclust)
library(ggplot2)
library(ggthemes)
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/module_validation_figs/psoriasis_GSE117239/Untitled.R")


gse <- 'GSE249936'

#read all modules
source("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")

list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <-  c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                   "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/Transendothelial_migration", 
                   "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/Innate_immune_response", 
                   "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired_cell_proliferation", "Lipid metabolism")

mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")

########## get gene expression data
data <- read.csv("GSE249936_counts_genes_filtered_RNASeq_Benezeder_2023.txt", sep = '\t')
data <- new_mapping_ids(data, "ensembl_gene_id")
data$ensembl_gene_id <- NULL
data <- na.omit(data)
data <- sum_expression_gene_rowise(data, 'mapped_symbol')
data <- as.data.frame(conte(data, gse = gse)$expr)
data <- order_by_col_names(data)

########## get metadata
meta <- GEOquery::getGEO(gse)
meta <- meta$GSE249936_series_matrix.txt.gz@phenoData@data
meta <- sort_by_col(meta, colname = 'title')



colnames(data) <- meta$geo_accession

data <- order_by_col_names(data)
meta <- sort_by_col(meta, colname = 'geo_accession')


condition_annotation <- meta$`disease state:ch1`
timepoint_annotation <- meta$`treatment:ch1`


counter_mod <- 0 #keep track of module number
list_plot <- list() # list of heatmaps plots
plot_aucs <- list() # list of auc plots
auc_values <- c()
rand_indices <- c()
list_boxplot <- list()

counter_mod <- 0
counter_boxplot <- 0

#select significant ones a priori
#significant ones are

for (module in list_modules) {#for each module
  counter_mod <- counter_mod+1
  
  ########## get gene expression data
  expr <- data
  
  ########## compute features for current module
  total_df <- compute_features(expr, module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  total_df <- na.omit(total_df)
  
  palette_for_anno <- list(
    #timepoint = c('none' = '#377eb8', 'dithranol' = '#e41a1b'),
    condition = c("psoriasis skin lesion" = "blue", 
                  "psoriasis skin lesion after treatment" = "orange",
                  "non-lesional skin" = 'pink')
  )
  
  anno <- HeatmapAnnotation(
    condition = condition_annotation,
    # timepoint = timepoint_annotation,
    col = palette_for_anno
  )
  
  
  p <-  Heatmap(
    as.matrix(total_df),
    bottom_annotation = anno,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 6),
    row_names_gp = grid::gpar(fontsize = 6),
    #heatmap_height = unit(0.2, "cm")*nrow(total_df),
    #heatmap_width = unit(0.3, "cm")*ncol(total_df),
    heatmap_legend_param = list(title = "Zscores"),
    na_col = 'black',
    show_column_names = F,
    show_heatmap_legend = F,
    show_row_names = F,
    column_title = module_names[counter_mod]
  )
  
  
  list_plot[[counter_mod]] <- grid.grabExpr(draw(p))
  
  
  ######## Adjusted Rand Index
  
  df_dendro_clustering <- as.data.frame(cutree(as.hclust(column_dend(draw(p))), k = 2))
  df_dendro_clustering$geo_accession <- rownames(df_dendro_clustering)
  colnames(df_dendro_clustering) <- c('dendro_clustering', 'geo_accession')
  df_dendro_clustering <- df_dendro_clustering[order(df_dendro_clustering$geo_accession),]
  
  if(!identical(meta$geo_accession, df_dendro_clustering$geo_accession)) {
    stop("different order")
  }
  
  labels <- ifelse(grepl("non-les|treat", condition_annotation), "control", "case")
  
  true_labels <- ifelse(grepl('control', labels), 1, 2)
  predicted_labels <- df_dendro_clustering$dendro_clustering
  
  rand_index <- adjustedRandIndex(true_labels, predicted_labels)
  rand_indices <- append(rand_indices, rand_index)
  
  
  #boxplot
  
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  score_df$Condition <- meta$`disease state:ch1`
  
  my_comparisons <- list(c("psoriasis skin lesion after treatment", "psoriasis skin lesion"),
                         c("psoriasis skin lesion after treatment", "non-lesional skin"),
                         c("non-lesional skin", "psoriasis skin lesion")
  )
  
  tableau_colors <- ggthemes::tableau_color_pal("Tableau 10")(10)  # Get the full palette of 10 colors
  # Subset and reverse the colors (e.g., use the first 5 colors and reverse them)
  custom_colors <- tableau_colors[2:10]
  
  
  box <- ggplot(score_df, aes(x = Condition, y = module)) +
    geom_boxplot(aes(fill = Condition)) +
    ggtitle(paste0(module_names[counter_mod], "\n"))+
    ylab("module score")+
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_blank(), 
          #legend.text=element_text(size=18),
          #legend.title=element_text(size=18)+#, legend.position="none") +
    )+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
    scale_fill_manual(values = custom_colors)
  
  #for legend
  # box <- box +
  #   theme(legend.text=element_text(size=16), legend.title = element_text(size=20), legend.position="bottom")
  # 
  # jpeg(file = 'legend.jpeg', width = 30 , units = 'cm', height = 20, res = 600)
  # box
  # dev.off()
  
  
  if (counter_mod %in% c(1, 2, 3, 4, 5, 10, 13, 15)) { #plot only significant ones
    counter_boxplot <- counter_boxplot+1
    list_boxplot[[counter_boxplot]] <- box
  }
  
  
  
}

#saveRDS(list_boxplot, file = 'list_boxplot.rds')
#do.call(grid.arrange,c(list_boxplot, ncol = 4))



list_boxplot <- readRDS("list_boxplot.rds")
list_boxplot_new <- list()
list_boxplot_new[[1]] <- list_boxplot[[7]]
list_boxplot_new[[2]] <- list_boxplot[[4]]
list_boxplot_new[[3]] <- list_boxplot[[5]]
list_boxplot_new[[4]] <- list_boxplot[[8]]
list_boxplot_new[[5]] <- list_boxplot[[1]]
list_boxplot_new[[6]] <- list_boxplot[[2]]
list_boxplot_new[[7]] <- list_boxplot[[3]]
list_boxplot_new[[8]] <- list_boxplot[[6]]


for (i in 1:length(list_boxplot_new)) {
  plot <- list_boxplot_new[[i]]
  
  new_plot <-  plot + 
    theme(
      axis.title.y = element_text(size = 18),
      plot.title = ggtext::element_markdown(size = 25, hjust = 0.5, face = 'bold', margin=margin(0,0,20,0))
      
    )
  
  new_plot$layers[[2]]$aes_params$textsize <- 7
  
  
  list_boxplot_new[[i]] <- new_plot
  
}


# jpeg(file = 'Fig7C_boxplot_psoriasis.jpeg', width = 68 , units = 'cm', height = 35, res = 600)
# do.call(grid.arrange,c(list_boxplot_new, ncol = 4))
# dev.off()

jpeg(file = 'Fig7C_boxplot_psoriasis_phd_thesis.jpeg', width = 45 , units = 'cm', height = 65, res = 600)
do.call(grid.arrange,c(list_boxplot_new, ncol = 2))
dev.off()

