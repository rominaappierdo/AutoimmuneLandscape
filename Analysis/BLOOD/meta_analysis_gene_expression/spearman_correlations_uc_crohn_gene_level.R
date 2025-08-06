setwd("Analysis/BLOOD/meta_analysis_gene_expression")

crohn <- read.csv("crohn/crohn_vs_healthy_whole blood_complete_results.tsv", sep = '\t')
crohn <- crohn[crohn$effectSizeFDR < 0.05, ]
crohn <- crohn[crohn$numStudies >= 2, ]


uc <- read.csv("uc/uc_vs_healthy__whole blood_complete_results.tsv", sep = '\t')
uc <- uc[uc$effectSizeFDR < 0.05, ]
uc <- uc[uc$numStudies >= 2, ]

merged <- merge(crohn, uc , by = 'gene')

cor(merged$effectSize.x, merged$effectSize.y)

cor.test(merged$effectSize.x, merged$effectSize.y, method = 'spearman')
