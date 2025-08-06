setwd("paper_figs/suppl_figs_featurization")
options(repr.plot.width=15, repr.plot.height=15)
library(dendextend)
library(ComplexHeatmap)
library(grid)

significant_genes<-c()
for (file in list.files("../../Analysis/TISSUE/meta_analysis_cell_proportions/", recursive = T, pattern = '*complete*', full.names = T)) {
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)



data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files("../../Analysis/TISSUE/meta_analysis_cell_proportions/", recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/|_")[[1]][10]
  tissue <- strsplit(file, "\\/|_")[[1]][11]
  disease <- paste0(tissue, "_", disease)
  df<-read.csv(file, sep='\t')
  df<-df[df$gene %in% significant_genes, c("effectSize", "gene")]
  df$disease<-rep(disease, nrow(df))
  data<-rbind(data, df)
}

#reshaping
data<-tidyr::spread(data, "disease", "effectSize")
rownames(data)<-data$gene
data$gene<-NULL

redundant_cell_types <- readLines("../../gene_sets/redundant_cell_types.txt")
data <- data[!rownames(data) %in% redundant_cell_types, ]

n_studies<- c(27, 30, 5, 6, 7, 23, 2, 3)

plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))


data[data > 1]<- 1
data[data < -1]<- -1
print(colnames(data))
colnames(data) <- c("Crohn's disease (Gut)", "Ulcerative colitis (Gut)", "Lupus (Kidney)", "Sjogren's syndrome (Salivary glands)",
                    "Lupus (Skin)", "Psoriasis (Skin)", "Scleroderma (Skin)", "Rheumatoid arthritis (Synovium)")

rownames(data) <- c("Basophils", "CD14+ monocytes", "CD4+ T cells", "CD56dim NK cells",
                    "CD8+ T cells", "Eosinophils", "Hematopoietic progenitors", "Macrophages M0",
                    "Macrophages M1", "Macrophages M2", "MAST cells", "memory B cells", "Myeloid dendritic cells", "naive B cells",
                    "Neutrophils",  "Plasma cells", "Plasmacytoid dendritic cells")

column_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F, add_numbers = T, numbers_rot = 0), which = 'column')

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
              row_names_gp = gpar(fontsize = 13),
              show_heatmap_legend = T,
              heatmap_width = unit(1.7, "cm")*ncol(data),
              heatmap_height = unit(2.5, "cm")*ncol(data)
)
ht


jpeg("ht_cell_tissue.jpeg", units="in", width=10, height=10, res=600)
draw(ht)
dev.off()
