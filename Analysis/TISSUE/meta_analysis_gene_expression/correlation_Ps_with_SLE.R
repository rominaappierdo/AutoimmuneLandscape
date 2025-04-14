setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_gene_expression")


library(ggplot2)

#pso and sle
pso <- read.csv("psoriasis/SKIN_pso_vs_healthy_tissue_complete_results.tsv", sep = '\t')
pso <- pso[pso$numStudies >=2, ]
pso <- pso[pso$effectSizeFDR < 0.05, ]

sle <- read.csv("lupus/SKIN_lupus_vs_healthy_tissue_complete_results.tsv", sep = '\t')
sle <- sle[sle$numStudies >=2, ]
sle <- sle[sle$effectSizeFDR < 0.05, ]

pso <- pso[, c('effectSize', 'gene')]
sle <- sle[, c('effectSize', 'gene')]

colnames(pso) <- c('ES_pso', 'gene')
colnames(sle) <- c('ES_sle', 'gene')

merged <- merge(pso, sle, by = 'gene')

cor(merged$ES_pso, merged$ES_sle, method = 'spearman')


ggplot(merged, aes(x = ES_sle, y = ES_pso)) +
  geom_point()+
  stat_cor(method = 'spearman') +
  geom_smooth(method = 'lm')
