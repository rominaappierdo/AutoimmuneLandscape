setwd("paper_figs/fig4C")
library(ggplot2)
library(ggpubr)
library(ggrepel)

blood <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/TFs/meta_analysis_tf_activity/sjogren/tf_gsva_sjogren_vs_healthy_whole blood_complete_results.tsv"
tissue <- "/mnt/data3/appierdo/AutoimmuneLandscape//Analysis/TISSUE/TFs/meta_analysis_tf_activity/sjogren/SALIVARY.GLAND_tf_gsva_sjogren_vs_healthy_tissue_complete_results.tsv"

b <- read.csv(blood, sep = '\t')
t <- read.csv(tissue, sep = '\t')

b <- b[, c('gene', 'effectSize')]
t <- t[, c('gene', 'effectSize')]

colnames(b) <- c('gene', 'Blood')
colnames(t) <- c('gene', 'Tissue')

merged <- merge(b, t, by = 'gene')

merged$selected_label <- ifelse(merged$gene %in% c("IRF1", "IRF2", "STAT2", "STAT2", "IRF9", "KLF5"), merged$gene, NA)


p <- ggplot(merged, aes(x = Blood, y = Tissue)) +
  theme_classic()+
  stat_cor(size = 7) +
  geom_smooth(method = 'lm', fill = 'white', col = '#9b18d6', linewidth = 1)+
  geom_text_repel(aes(label = selected_label), size = 5) +
  geom_point(size = 2) +
  ggtitle("Sjogren's syndrome | Transcription factors")+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(color = "#9b18d6", face = "bold", size = 19))


jpeg("fig4C_sjogren_tf.jpeg", width = 24,height = 18, units = "cm", res = 600)
p
dev.off()