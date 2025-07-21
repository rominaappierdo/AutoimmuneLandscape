setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/suppl_figs_featurization")
options(repr.plot.width=15, repr.plot.height=15)
library(dendextend)
library(ComplexHeatmap)
library(grid)

immunerelated <- c(readLines("../../gene_sets/hallmark_immune_related.txt"), readLines("../../gene_sets/kegg_immune_related.txt"), readLines("../../gene_sets/wiki_immune_related.txt"))

significant_genes<-c()
for (file in list.files("../../Analysis/TISSUE/meta_analysis_pathways_activity/merged_pathways_results/", recursive = T, pattern = '*complete*', full.names = T)) {
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  t <- t[abs(t$effectSize) > 1, ]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)



data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files("../../Analysis/TISSUE/meta_analysis_pathways_activity/merged_pathways_results", recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/|_")[[1]][13]
  tissue <- strsplit(file, "\\/|_")[[1]][12]
  disease <- paste0(disease, "_", tissue)
 
  df<-read.csv(file, sep='\t')
  df<-df[df$gene %in% significant_genes, c("effectSize", "gene")]
  df$disease<-rep(disease, nrow(df))
  data<-rbind(data, df)
}

#reshaping
data<-tidyr::spread(data, "disease", "effectSize")
rownames(data)<-data$gene
data$gene<-NULL

data <- data[rownames(data) %in% immunerelated, ]

n_studies<- c(27, 5, 7, 23, 3, 2, 6, 30)

plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))


data[data > 2]<- 2
data[data < -2]<- -2

colnames(data) <- c("Crohn's disease (Gut)", "Lupus (Kidney)", "Lupus (Skin)", "Psoriasis (Skin)", "Rheumatoid arthritis (Synovium)", 
                    "Scleroderma (Skin)", "Sjogren's syndrome (Salivary glands)", "Ulcerative colitis (Gut)")



column_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F), which = 'column')

ht <- Heatmap(as.matrix(data), 
              name = "pooled_effectSize", 
              cluster_rows = row_dend, 
              cluster_columns = column_dend, 
              na_col = 'black',
              column_names_rot = 45,
              bottom_annotation = column_ha,
              rect_gp = gpar(col = "black", lwd = 0.3),
              col = circlize::colorRamp2(c(-1, 0, 1), c("#56B4E9", "#FFFFFF", "#FF7400")),
              column_names_gp = gpar(fontsize = 13),
              row_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = T,
              heatmap_width = unit(1.8, "cm")*ncol(data),
              heatmap_height = unit(3, "cm")*ncol(data)
)
draw(ht, heatmap_legend_side = 'left')


jpeg("ht_pathways_tissue.jpeg", units="in", width=15, height=8, res=600)
draw(ht, heatmap_legend_side = 'left')
dev.off()
