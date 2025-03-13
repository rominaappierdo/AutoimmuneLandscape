setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/visualizeLandscape/script_to_select_significant_features_tissue")

#cytokines
cyto <- readLines("../../../signaling_molecules/curated_kines/all_kines.txt")

significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_gene_expression", recursive = T, pattern = '*complete*', full.names = T)) {
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  #t<-t[t["effectSize"]<= -1 | t["effectSize"]>= 1,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)
significant_genes <- significant_genes[significant_genes %in% cyto]

write.table(significant_genes, "significant_cyto_in_at_least_one_disease.txt",quote = F, col.names = F, row.names = F, sep = '\n')


#cell types
redundant_ct <- readLines("../../../gene_sets/redundant_cell_types.txt")
significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_cell_proportions/", recursive = T, pattern = '*complete*', full.names = T)) {
  print(file)
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  ##t<-t[t["effectSize"]<= -1 | t["effectSize"]>= 1,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)
significant_genes <- significant_genes[!significant_genes %in% redundant_ct]
write.table(significant_genes, "significant_cell_in_at_least_one_disease.txt",quote = F, col.names = F, row.names = F, sep = '\n')

#hallmark
significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_pathways_activity/hallmarks/", recursive = T, pattern = '*complete*', full.names = T)) {
  print(file)
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  #t<-t[t["effectSize"]<= -1 | t["effectSize"]>= 1,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)
write.table(significant_genes, "significant_hall_in_at_least_one_disease.txt",quote = F, col.names = F, row.names = F, sep = '\n')

#kegg
significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_pathways_activity/kegg/", recursive = T, pattern = '*complete*', full.names = T)) {
  print(file)
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  #t<-t[t["effectSize"]<= -1 | t["effectSize"]>= 1,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)
write.table(significant_genes, "significant_kegg_in_at_least_one_disease.txt",quote = F, col.names = F, row.names = F, sep = '\n')

#TF
significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/TFs/meta_analysis_tf_activity/", recursive = T, pattern = '*complete*', full.names = T)) {
  print(file)
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  #t<-t[t["effectSize"]<= -1 | t["effectSize"]>= 1,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)
write.table(significant_genes, "significant_tf_in_at_least_one_disease.txt",quote = F, col.names = F, row.names = F, sep = '\n')

#mirna
significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/miRNA/meta_analysis_mirna_activity", recursive = T, pattern = '*complete*', full.names = T)) {
  print(file)
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  #t<-t[t["effectSize"]<= -1 | t["effectSize"]>= 1,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)
write.table(significant_genes, "significant_mirna_in_at_least_one_disease.txt",quote = F, col.names = F, row.names = F, sep = '\n')


#wikipathways

significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_pathways_activity/wikipathways/", recursive = T, pattern = '*complete*', full.names = T)) {
  print(file)
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  #t<-t[t["effectSize"]<= -1 | t["effectSize"]>= 1,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)
write.table(significant_genes, "significant_wiki_in_at_least_one_disease.txt",quote = F, col.names = F, row.names = F, sep = '\n')

