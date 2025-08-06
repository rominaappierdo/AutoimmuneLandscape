setwd('/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/lupus_kidney_enrichment_analysis')

df <- read.csv("../TISSUE/meta_analysis_gene_expression/lupus/KIDNEY_lupus_vs_healthy_tissue_complete_results.tsv", sep ='\t')
#write.table(df$gene, file = 'background_genes_LN.txt', sep = '\n', quote = F, row.names = F, col.names = F)
df <- df[df$effectSizeFDR < 0.05, ]
df <- df[df$numStudies >=2, ]

df_pos <- df[df$effectSize > 0, ]
df_neg <- df[df$effectSize < 0, ]

write.table(df_pos$gene, file = 'upreg_genes_LN.txt', sep = '\n', quote = F, row.names = F, col.names = F)
