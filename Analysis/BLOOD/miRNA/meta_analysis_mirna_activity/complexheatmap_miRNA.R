setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/miRNA/meta_analysis_mirna_activity/")
options(repr.plot.width=15, repr.plot.height=15)

library(ComplexHeatmap)
library(grid)
library(dendextend)

significant_genes<-c()
for (file in list.files(getwd(), recursive = T, pattern = '*complete*', full.names = T)) {
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  t<-t[t["effectSize"]<= -0.8 | t["effectSize"]>= 0.8,]
  significant_genes<-append(significant_genes, t$gene)
}
significant_genes <- unique(significant_genes)

#write.table(significant_genes, file = '../../../../visualizeLandscape/mirna_significant_genes_high_es.txt', sep = '\n', quote = F, col.names = F, row.names = F)


data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files(getwd(), recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/")[[1]][10]
  print(disease)
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

data<-na.omit(data)
plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))

saveRDS(column_dend, file = "/mnt/data3/appierdo/PhD-project/paper/paper_figs/fig2/column_dend_mirna.rds")
saveRDS(row_dend, file = "/mnt/data3/appierdo/PhD-project/paper/paper_figs/fig2/row_dend_mirna.rds")

data[data > 1]<- 1
data[data < -1]<- -1

colnames(data) <- c("Crohn's disease", "Type 1 diabetes", "Juvenile idiopathic arthritis", "Lupus", "Multiple sclerosis", 
                    "Psoriasis", "Rheumatoid arthritis", "Scleroderma", "Sjogren's syndrome", 
                    "Ulcerative colitis")
column_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F), which = 'column')

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
        heatmap_width = unit(1.4, "cm")*ncol(data),
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 14),
        rect_gp = gpar(col = "black", lwd = 0.2)
)

ht = draw(ht)
tiff("ht_mirna.tiff", units="in", width=10, height=10, res=400)
draw(ht)

stop()
write.table(data, file = '/mnt/data3/appierdo/PhD-project/paper/paper_figs/fig2/mirna.tsv', sep = '\t', row.names = T, col.names = NA, quote = F)

stop()
# 
# 
# clusters <- cutree(as.hclust(row_dend(ht)), k = 2)
# cluster1 <- names(clusters[clusters == 1]) #smaller cluster
# cluster2 <- names(clusters[clusters == 2]) #larger cluster
# 
# write.table(rownames(data), file='background.txt', quote = F, row.names = F, col.names = F, sep='\n')
# write.table(cluster1, file='cluster1.txt', quote = F, row.names = F, col.names = F, sep='\n')
# write.table(cluster2, file='cluster2.txt', quote = F, row.names = F, col.names = F, sep='\n')

# 
# 
row_dend = color_branches(row_dend, k = 2)

ht<-Heatmap(as.matrix(data),
            #name = "effectSize",
            name = "legend",
            cluster_rows = row_dend,
            cluster_columns = column_dend,
            na_col = 'black',
            column_names_rot = 45,
            #bottom_annotation = column_ha,
            col = circlize::colorRamp2(c(-1, 0, 1), c("#56B4E9", "#FFFFFF", "#FF7400")),
            show_row_names = F,
            show_column_names = F,
            use_raster = F,
            heatmap_width = unit(1.5, "cm")*ncol(data),
            column_names_gp = grid::gpar(fontsize = 14),
            row_names_gp = grid::gpar(fontsize = 9)
            )

draw(ht)


# clusters <- cutree(as.hclust(row_dend(ht)), k = 8)
# cluster_downreg <- names(clusters[clusters == 2])
# "ICAM2" %in% cluster_downreg
# write.table(cluster_downreg, file='cluster_downreg.txt', quote = F, row.names = F, col.names = F, sep='\n')
