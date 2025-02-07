setwd('/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/TFs/meta_analysis_tf_activity/')
source("../../../../utils_functions/for_plotting/for_pheatmap.R")
options(repr.plot.width=15, repr.plot.height=15)

library(pheatmap)
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

data<-t(data)

plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))

data[data > 1]<- 1
data[data < -1]<- -1

heatmap<-pheatmap::pheatmap(as.matrix(data),
                            cluster_cols = column_dend,
                            cluster_rows = row_dend, 
                            na_col = 'black',
                            cellwidth =12
                            )
# svg('heat2_deconv.svg')
#ggsave(file="heat1.svg", plot=heatmap1, width=10, height=8)

#library(ggplot2)
#ggsave(file="heat2.svg", plot=heatmap2, width=10, height=8)
