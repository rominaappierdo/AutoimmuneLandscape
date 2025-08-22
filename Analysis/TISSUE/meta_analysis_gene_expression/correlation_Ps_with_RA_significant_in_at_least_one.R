setwd("Analysis/TISSUE/meta_analysis_gene_expression")
library(ggplot2)

#pso and RA
pso <- read.csv("psoriasis/SKIN_pso_vs_healthy_tissue_complete_results.tsv", sep = '\t')
pso <- pso[pso$numStudies >=2, ]
pso <- pso[pso$effectSizeFDR < 0.05, ]
pso_genes <- pso$gene

ra <- read.csv("ra/SYNOVIUM_ra_vs_healthy_tissue_complete_results.tsv", sep = '\t')
ra <- ra[ra$numStudies >=2, ]
ra <- ra[ra$effectSizeFDR < 0.05, ]
ra_genes <- ra$gene

pso <- read.csv("psoriasis/SKIN_pso_vs_healthy_tissue_complete_results.tsv", sep = '\t')
ra <- read.csv("ra/SYNOVIUM_ra_vs_healthy_tissue_complete_results.tsv", sep = '\t')

all_genes <- c(pso_genes, ra_genes)

pso <- pso[pso$gene %in% all_genes, ]
ra <- ra[ra$gene %in% all_genes, ]


pso <- pso[, c('effectSize', 'gene')]
ra <- ra[, c('effectSize', 'gene')]

colnames(pso) <- c('ES_pso', 'gene')
colnames(ra) <- c('ES_ra', 'gene')

merged <- merge(pso, ra, by = 'gene')

cor(merged$ES_pso, merged$ES_ra, method = 'spearman')


ggplot(merged, aes(x = ES_pso, y = ES_ra)) +
  geom_point() +
  xlim(-6, 6) +
  ylim(-3, 3)
