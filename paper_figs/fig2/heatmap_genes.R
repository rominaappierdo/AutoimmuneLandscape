setwd('/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig2')
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")
options(repr.plot.width=15, repr.plot.height=15)
library(ComplexHeatmap)
library(grid)
library(dendextend)


significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("../../Analysis/BLOOD/meta_analysis_gene_expression", recursive = T, pattern = '*complete*', full.names = T)) {
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"] < 0.05,]
  t<-t[t["effectSize"]<= -0.8 | t["effectSize"]>= 0.8,]
  significant_genes<-append(significant_genes, t$gene)
}

data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files("../../Analysis/BLOOD/meta_analysis_gene_expression", recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/")[[1]][6]
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
row_dend <- color_branches(row_dend, k=2, col = c('blue', 'red'))


column_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F, add_numbers = T, numbers_rot = 0), which = 'column')

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
            show_row_names = F,
            use_raster = F,
            column_names_gp = gpar(fontsize = 11),
            row_names_gp = gpar(fontsize = 8),
            heatmap_width = unit(0.8, "cm")*ncol(data),
            heatmap_height = unit(0.025, "cm")*nrow(data),
            rect_gp = gpar(col = "black", lwd = 0.03), 
            row_title = 'Genes'
)


ht = draw(ht)


jpeg("ht_genes.jpeg", units="in", width=8, height=9, res=600)
draw(ht)
dev.off()

