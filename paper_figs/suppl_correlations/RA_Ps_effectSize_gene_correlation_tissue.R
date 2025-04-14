setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/suppl_correlations")

pso <- read.csv("../../Analysis/TISSUE/meta_analysis_gene_expression/psoriasis/SKIN_pso_vs_healthy_tissue_complete_results.tsv", sep = '\t')
pso <- pso[pso$numStudies >=2, ]
pso <- pso[pso$effectSizeFDR < 0.05, ]

ra <- read.csv("../../Analysis/TISSUE/meta_analysis_gene_expression/ra/SYNOVIUM_ra_vs_healthy_tissue_complete_results.tsv", sep = '\t')
ra <- ra[ra$numStudies >=2, ]
ra <- ra[ra$effectSizeFDR < 0.05, ]

merged <- merge(pso, ra, by = 'gene')

correlation <- psych::corr.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')
