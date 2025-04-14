setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_gene_expression")


library(ggplot2)

#SS and SSc
sjogren <- read.csv("sjogren/SALIVARY.GLAND_sjogren_vs_healthy__tissue_complete_results.tsv", sep = '\t')
sjogren <- sjogren[sjogren$numStudies >=2, ]
sjogren <- sjogren[sjogren$effectSizeFDR < 0.05, ]

scleroderma <- read.csv("scleroderma/SKIN_scleroderma_vs_healthy_tissue_complete_results.tsv", sep = '\t')
scleroderma <- scleroderma[scleroderma$numStudies >=2, ]
scleroderma <- scleroderma[scleroderma$effectSizeFDR < 0.05, ]

sjogren <- sjogren[, c('effectSize', 'gene')]
scleroderma <- scleroderma[, c('effectSize', 'gene')]

colnames(sjogren) <- c('ES_sjogren', 'gene')
colnames(scleroderma) <- c('ES_scleroderma', 'gene')

merged <- merge(sjogren, scleroderma, by = 'gene')

cor(merged$ES_sjogren, merged$ES_scleroderma, method = 'spearman')



