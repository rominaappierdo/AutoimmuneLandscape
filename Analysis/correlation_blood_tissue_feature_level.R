setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis")

#SS - TF
sjogren_blood <- read.csv("BLOOD/TFs/meta_analysis_tf_activity/sjogren/tf_gsva_sjogren_vs_healthy_whole blood_complete_results.tsv", sep = '\t')
sjogren_blood <- sjogren_blood[sjogren_blood$numStudies >= 2, ]
#sjogren_blood <- sjogren_blood[sjogren_blood$effectSizeFDR < 0.05, ]
sjogren_tissue <- read.csv("TISSUE/TFs/meta_analysis_tf_activity/sjogren/SALIVARY.GLAND_tf_gsva_sjogren_vs_healthy_tissue_complete_results.tsv", sep = '\t')
sjogren_tissue <- sjogren_tissue[sjogren_tissue$numStudies >=2, ]
#sjogren_tissue <- sjogren_tissue[sjogren_tissue$effectSizeFDR < 0.05, ]


merged <- merge(sjogren_blood, sjogren_tissue, by = 'gene')
correlation <- psych::corr.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')


#SSc - cell types
scleroderma_blood <- read.csv("BLOOD/meta_analysis_cell_proportions/scleroderma/cell_prop_scleroderma_vs_healthy_whole blood_complete_results.tsv", sep = '\t')
#scleroderma_blood <- scleroderma_blood[scleroderma_blood$numStudies >= 2, ]
#scleroderma_blood <- scleroderma_blood[scleroderma_blood$effectSizeFDR < 0.05, ]
scleroderma_tissue <- read.csv("TISSUE/meta_analysis_cell_proportions/scleroderma/SKIN_scleroderma_vs_healthy_tissue_complete_results.tsv", sep = '\t')
#scleroderma_tissue <- scleroderma_tissue[scleroderma_tissue$numStudies >=2, ]
#scleroderma_tissue <- scleroderma_tissue[scleroderma_tissue$effectSizeFDR < 0.05, ]

merged <- merge(scleroderma_blood, scleroderma_tissue, by = 'gene')
correlation <- psych::corr.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')
