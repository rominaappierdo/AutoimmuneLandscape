setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig3")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
options(repr.plot.width=15, repr.plot.height=15)
library(ComplexHeatmap)
library(grid)
library(dendextend)

significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("../../Analysis/TISSUE/meta_analysis_gene_expression", recursive = T, pattern = '*complete*', full.names = T)) {
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"] < 0.05,]
  #t<-t[t$heterogeneityPval > 0.05,]
  t<-t[t["effectSize"]<= -0.8 | t["effectSize"]>= 0.8,]
  significant_genes<-append(significant_genes, t$gene)
}


sign_genes_in_50_perc_diseases <- names(table(significant_genes)[table(significant_genes) >= 2])


data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files("../../Analysis/TISSUE/meta_analysis_gene_expression", recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/")[[1]][6]
  tissue <- strsplit(file, "\\/|_")[[1]][10]
  disease <- paste0(disease, "_", tissue, sep = '')
  df<-read.csv(file, sep='\t')
  df<-df[df$gene %in% sign_genes_in_50_perc_diseases, c("effectSize", "gene")]
  #df<-df[df$gene %in% significant_genes, c("effectSize", "gene")]
  
  df$disease<-rep(disease, nrow(df))
  data<-rbind(data, df)
}

#reshaping
data<-tidyr::spread(data, "disease", "effectSize")
rownames(data)<-data$gene
data$gene<-NULL

n_studies<- c(27, 5, 7, 23, 3, 2, 6, 30)
data<-na.omit(data)
#data <- as.data.frame(t(data))
#plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))
#column_dend = color_branches(column_dend, k = 2)
row_dend = color_branches(row_dend, k = 2, col = c('blue', 'red'))

# data[data > 2]<- 2
# data[data < -2]<- -2


row_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F, add_numbers = T, numbers_rot = 0), which = 'column')

colnames(data) <- c("Crohn's disease (Gut)", "Lupus (Kidney)", "Lupus (Skin)", "Psoriasis (Skin)",
                    "Rheumatoid arthritis (Synovium)", "Scleroderma (Skin)", "Sjogren's syndrome (Salivary glands)",
                    "Ulcerative colitis (Gut)")

ht<-Heatmap(as.matrix(data), 
            name = "pooled_effectSize", 
            cluster_rows = row_dend, 
            cluster_columns = column_dend, 
            na_col = 'black',
            column_names_rot = 45,
            bottom_annotation = row_ha,
            #col = circlize::colorRamp2(c(-1, 0, 1), c("#56B4E9", "#FFFFFF", "#FF7400")),
            show_row_names = F,
            show_column_names = T,
            use_raster = F,
            column_names_gp = gpar(fontsize = 11),
            #row_names_gp = gpar(fontsize = 9),
            heatmap_width = unit(0.95, "cm")*ncol(data),
            heatmap_height = unit(0.013, "cm")*nrow(data),
            rect_gp = gpar(col = "black", lwd = 0.01),
            show_heatmap_legend = T
)


ht <- draw(ht)
ht


jpeg("ht_genes.jpeg", units="in", width=8, height=10, res=600)
draw(ht)
dev.off()


# 
# clusters <- cutree(as.hclust(row_dend(ht)), k = 2)
# cluster1 <- names(clusters[clusters == 1]) 
# cluster2 <- names(clusters[clusters == 2]) 
# 
# write.table(cluster1, file='cluster1.txt', quote = F, row.names = F, col.names = F, sep='\n')
# write.table(cluster2, file='cluster2.txt', quote = F, row.names = F, col.names = F, sep='\n')