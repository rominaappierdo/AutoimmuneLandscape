setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_gene_expression")


library(ggplot2)

#pso and RA
pso <- read.csv("psoriasis/SKIN_pso_vs_healthy_tissue_complete_results.tsv", sep = '\t')
pso <- pso[pso$numStudies >=2, ]
pso <- pso[pso$effectSizeFDR < 0.05, ]

ra <- read.csv("ra/SYNOVIUM_ra_vs_healthy_tissue_complete_results.tsv", sep = '\t')
ra <- ra[ra$numStudies >=2, ]
ra <- ra[ra$effectSizeFDR < 0.05, ]

pso <- pso[, c('effectSize', 'gene')]
ra <- ra[, c('effectSize', 'gene')]

colnames(pso) <- c('ES_pso', 'gene')
colnames(ra) <- c('ES_ra', 'gene')

merged <- merge(pso, ra, by = 'gene')

cor(merged$ES_pso, merged$ES_ra, method = 'spearman')

