setwd("paper_figs/suppl_correlations")

scleroderma <- read.csv("../../Analysis/TISSUE/meta_analysis_gene_expression/scleroderma/SKIN_scleroderma_vs_healthy_tissue_complete_results.tsv", sep = '\t')
scleroderma <- scleroderma[scleroderma$numStudies >=2, ]
scleroderma <- scleroderma[scleroderma$effectSizeFDR < 0.05, ]

sjogren <- read.csv("../../Analysis/TISSUE/meta_analysis_gene_expression/sjogren/SALIVARY.GLAND_sjogren_vs_healthy__tissue_complete_results.tsv", sep = '\t')
sjogren <- sjogren[sjogren$numStudies >=2, ]
sjogren <- sjogren[sjogren$effectSizeFDR < 0.05, ]

merged <- merge(scleroderma, sjogren, by = 'gene')

correlation <- psych::corr.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')
