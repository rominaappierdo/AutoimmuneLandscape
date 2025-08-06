setwd("Supplementary/supp_figures")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
library(ggplot2)
library(ggrepel)
library(gridExtra)

diseases <- c('pso', 'sjog', 'sclero', 'uc', 'crohn', "ra_", "lupus", "lupus")
plot_titles <- c('Ps', 'SS', 'SSc', 'UC', 'CD', "RA", "SLE Kidney", "SLE Skin")

c <- 0
list_plot <- list()
for (disease in diseases ){
  c <- c+1
  file_blood <- list.files("../../Analysis/BLOOD/meta_analysis_gene_expression/", pattern = paste0('.*', disease, '.*complet.*'), recursive = T, full.names = T)
  
  if (!disease == 'lupus') {
    file_tissue <- list.files("../../Analysis/TISSUE/meta_analysis_gene_expression/", pattern = paste0('.*', disease, '.*complet.*'), recursive = T, full.names = T)
    
    if(!length(file_blood) + length(file_tissue) == 2) {
      stop('error')
    }
  } else if (disease == 'lupus' & c == 7) {
    file_tissue <- list.files("../../Analysis/TISSUE/meta_analysis_gene_expression/", pattern = paste0('.*', disease, '.*complet.*'), recursive = T, full.names = T)
    file_tissue <- file_tissue[1]
    disease <- 'lupsu_kidney'
    
  } else {
    file_tissue <- list.files("../../Analysis/TISSUE/meta_analysis_gene_expression/", pattern = paste0('.*', disease, '.*complet.*'), recursive = T, full.names = T)
    file_tissue <- file_tissue[2]
    disease <- 'lupsu_skin'
  }
  
  blood <- read.csv(file_blood, sep = '\t')
  blood <- blood[blood$numStudies >=2, ]
  blood <- blood[blood$effectSizeFDR < 0.05, ]
  blood <- rename_single_column(blood, 'effectSize', 'effectSizeBlood')
  tissue <- read.csv(file_tissue, sep = '\t')
  tissue <- tissue[tissue$numStudies >=2, ]
  tissue <- tissue[tissue$effectSizeFDR < 0.05, ]
  tissue <- rename_single_column(tissue, 'effectSize', 'effectSizeTissue')
  
  merged <- merge(blood, tissue, by = 'gene')
  
  directionality_vector <- ifelse((merged$effectSizeBlood > 0 & merged$effectSizeTissue > 0), 'Opposite direction', ifelse((merged$effectSizeBlood < 0 & merged$effectSizeTissue < 0), 'Opposite direction', 'Same direction'))
  merged$Directionality <- directionality_vector
  
  p <- ggplot(merged, aes(x =effectSizeBlood, y = effectSizeTissue)) +
    geom_point(aes(color = Directionality)) +
    scale_color_manual(values = c('red', 'blue'))+
    theme_classic() +
    ggtitle(plot_titles[c]) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed')
  
  list_plot[[c]] <- p
  
}

jpeg(file = 'SuppFig5.jpeg', width = 21, height = 22, units = 'cm', res = 600)
do.call(grid.arrange,c(list_plot, ncol = 2))
dev.off()


