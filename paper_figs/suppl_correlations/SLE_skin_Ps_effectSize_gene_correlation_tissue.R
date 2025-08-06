setwd("paper_figs/suppl_correlations")

pso <- read.csv("../../Analysis/TISSUE/meta_analysis_gene_expression/psoriasis/SKIN_pso_vs_healthy_tissue_complete_results.tsv", sep = '\t')
pso <- pso[pso$numStudies >=2, ]
pso <- pso[pso$effectSizeFDR < 0.05, ]

sle <- read.csv("../../Analysis/TISSUE/meta_analysis_gene_expression/lupus//SKIN_lupus_vs_healthy_tissue_complete_results.tsv", sep = '\t')
sle <- sle[sle$numStudies >=2, ]
sle <- sle[sle$effectSizeFDR < 0.05, ]

merged <- merge(pso, sle, by = 'gene')

correlation <- psych::corr.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')
