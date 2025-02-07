setwd('/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/TFs/meta_analysis_tf_activity/')
options(repr.plot.width=15, repr.plot.height=15)

library(ComplexHeatmap)
library(grid)

significant_genes<-c()
for (file in list.files(getwd(), recursive = T, pattern = '*complete*', full.names = T)) {
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  t<-t[t["effectSize"]<=-0.8 | t["effectSize"]>= 0.8,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)

#write.table(significant_genes, file = '../../../../visualizeLandscape/tf_significant_genes_high_es.txt', sep = '\n', quote = F, col.names = F, row.names = F)

data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files(getwd(), recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/")[[1]][10]
  df<-read.csv(file, sep='\t')
  df<-df[df$gene %in% significant_genes, c("effectSize", "gene")]
  df$disease<-rep(disease, nrow(df))
  data<-rbind(data, df)
}

#reshaping
data<-tidyr::spread(data, "disease", "effectSize")
rownames(data)<-data$gene
data$gene<-NULL
n_studies<- c(8, 4, 2, 16, 4, 3, 5, 2, 3, 6) 

#data<-t(data)

plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))

saveRDS(column_dend, file = "/mnt/data3/appierdo/PhD-project/paper/paper_figs/fig2/column_dend_tf.rds")
saveRDS(row_dend, file = "/mnt/data3/appierdo/PhD-project/paper/paper_figs/fig2/row_dend_tf.rds")

data[data > 1]<- 1
data[data < -1]<- -1

column_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F), which = 'column')
colnames(data) <- c("Crohn's disease", "Type 1 diabetes", "Juvenile idiopathic arthritis", "Lupus", "Multiple sclerosis", 
                    "Psoriasis", "Rheumatoid arthritis", "Scleroderma", "Sjogren's syndrome", 
                    "Ulcerative colitis")

ht<-Heatmap(as.matrix(data), 
            name = "pooled_effectSize", 
            cluster_rows = row_dend, 
            cluster_columns = column_dend, 
            na_col = 'black',
            column_names_rot = 45,
            bottom_annotation = column_ha,
            col = circlize::colorRamp2(c(-1, 0, 1), c("#56B4E9", "#FFFFFF", "#FF7400")),
            show_row_names = T,
            use_raster = F,
            show_heatmap_legend = F,
            heatmap_width = unit(1.4, "cm")*ncol(data),
            heatmap_height = unit(23.5, "cm"),
            row_names_gp = gpar(fontsize = 9),
            column_names_gp = gpar(fontsize = 14),
            rect_gp = gpar(col = "black", lwd = 0.2)
            
            )

ht = draw(ht)
tiff("ht_tf.tiff", units="in", width=10, height=10, res=400)
draw(ht)
dev.off()

stop()
write.table(data, file = '/mnt/data3/appierdo/PhD-project/paper/paper_figs/fig2/tf.tsv', sep = '\t', row.names = T, col.names = NA, quote = F)

#write.table(rownames(data), file = 'tf_visualizaLandscape.txt', sep = '\n', quote = F, col.names = F, row.names = F)
