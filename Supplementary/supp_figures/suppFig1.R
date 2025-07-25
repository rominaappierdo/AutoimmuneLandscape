setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Supplementary/supp_figures")
options(repr.plot.width=15, repr.plot.height=15)
library(dendextend)
library(ComplexHeatmap)
library(grid)

significant_genes<-c()
for (file in list.files("../../Analysis/BLOOD/meta_analysis_cell_proportions/", recursive = T, pattern = '*complete*', full.names = T)) {
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)



data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files("../../Analysis/BLOOD/meta_analysis_cell_proportions/", recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/")[[1]][7]
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

n_studies<- c(8, 4, 2, 16, 4, 3, 5, 2, 3, 6)

plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))

data[data > 1]<- 1
data[data < -1]<- -1

colnames(data) <- c("Crohn's disease", "Type 1 diabetes", "Juvenile idiopathic arthritis", "Lupus", "Multiple sclerosis", 
                    "Psoriasis", "Rheumatoid arthritis", "Scleroderma", "Sjogren's syndrome", 
                    "Ulcerative colitis")

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
              row_names_gp = gpar(fontsize = 13),
              show_heatmap_legend = F,
              heatmap_width = unit(1.7, "cm")*ncol(data),
              heatmap_height = unit(2.35, "cm")*ncol(data)
)


jpeg("SuppFig1.jpeg", units="cm", width=23, height=25, res=600)
draw(ht)
dev.off()
