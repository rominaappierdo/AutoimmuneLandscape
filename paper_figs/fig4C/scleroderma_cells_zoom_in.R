setwd("paper_figs/fig4C")
library(ggplot2)
library(ggpubr)
library(ggrepel)

blood <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/meta_analysis_cell_proportions/scleroderma/cell_prop_scleroderma_vs_healthy_whole blood_complete_results.tsv"
tissue <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_cell_proportions/scleroderma/SKIN_scleroderma_vs_healthy_tissue_complete_results.tsv"

red <- readLines("/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/redundant_cell_types.txt")

b <- read.csv(blood, sep = '\t')
t <- read.csv(tissue, sep = '\t')

b <- b[, c('gene', 'effectSize')]
t <- t[, c('gene', 'effectSize')]

colnames(b) <- c('gene', 'Blood')
colnames(t) <- c('gene', 'Tissue')

merged <- merge(b, t, by = 'gene')
merged <- merged[!merged$gene %in% red, ]


merged[merged$gene == 'basophil', ]$gene <- "basophils"
merged[merged$gene == 'CD14_positive_monocyte', ]$gene <- "CD14+ monocytes"
merged[merged$gene == 'CD16_positive_monocyte', ]$gene <- "CD16+ monocytes"
merged[merged$gene == 'CD4_positive_alpha_beta_T_cell', ]$gene <- "CD4+ T cells"
merged[merged$gene == 'CD8_positive_alpha_beta_T_cell', ]$gene <- "CD8+ T cells"
merged[merged$gene == 'hematopoietic_progenitor', ]$gene <- "hematopoietic progenitors"
merged[merged$gene == 'macrophage_m0', ]$gene <- "M0 macrophages"
merged[merged$gene == 'macrophage_m1', ]$gene <- "M1 macrophages"
merged[merged$gene == 'macrophage_m2', ]$gene <- "M2 macrophages"
merged[merged$gene == 'MAST_cell', ]$gene <- "MAST cells"
merged[merged$gene == 'memory_B_cell', ]$gene <- "Memory B cells"
merged[merged$gene == 'naive_B_cell', ]$gene <- "Naive B cells"
merged[merged$gene == 'plasma_cell', ]$gene <- "Plasma cells"
merged[merged$gene == 'plasmacytoid_dendritic_cell', ]$gene <- "Plasmacytoid dendritic cells"



p <- ggplot(merged, aes(x = Blood, y = Tissue)) +
  theme_classic()+
  stat_cor(size = 8) +
  geom_smooth(method = 'lm', fill = 'white', col = '#528ab9', linewidth = 1)+
  geom_text_repel(aes(label = gene), 
                  size = 6
  ) +
  geom_point(size = 2) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(color = "#528ab9", face = "bold", size = 19)) +
  ggtitle("Scleroderma | Cell types")



jpeg("fig4C_scleroderma_cells.jpeg", width = 24,height = 18, units = "cm", res = 600)
p
dev.off()