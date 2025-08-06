setwd("Analysis/BLOOD/pipeline_landscape_disease_level_final")
library(ComplexHeatmap)
library(dendextend)

correlation_matrix <- read.csv("input_correlation_mat.tsv", sep = '\t', row.names = 1, check.names = F)
dim(correlation_matrix)

#compute correlation matrix and use it as.a distance measure
correlation_matrix<- Hmisc::rcorr(as.matrix(correlation_matrix), type = 'spearman')$r

distance_matrix <- as.dist((1 - correlation_matrix)/2)

nclusters <- 15
column_dend = hclust(distance_matrix)
row_dend = hclust(distance_matrix)
row_dend = color_branches(row_dend, k = nclusters)


# row_dend = color_branches(row_dend, k = nclusters)
# labels(row_dend) <- rep('', 700)
# plot(row_dend)

p <- ComplexHeatmap::Heatmap(as.matrix(correlation_matrix),
                             show_row_names = F,
                             cluster_rows = row_dend, 
                             cluster_columns = column_dend,
                             column_names_rot = 50,
                             row_split = nclusters,
                             col = circlize::colorRamp2(c(-0.9, 0, 0.9), c("#56B4E9", "#FFFFFF", "#FF7400")),
                             row_names_gp = grid::gpar(fontsize = 10),
                             column_names_gp = grid::gpar(fontsize = 2),
                             show_heatmap_legend = F,
                             heatmap_legend_param = list(title = 'correlation'),
                             na_col = 'black',
                             show_column_names = F,
                             heatmap_width = unit(0.25, "mm")*ncol(correlation_matrix),
                             heatmap_height = unit(0.25, "mm")*ncol(correlation_matrix)
)
p <- draw(p, padding = unit(c(2, 2, 2, 2), "mm"))


#saveRDS(row_dend, file= 'row_dend.rds')

cluster_membership <- as.data.frame(cutree(column_dend, k= nclusters))
colnames(cluster_membership) <- c('cluster')
cluster_membership$feature <- rownames(cluster_membership)
clusters<- split(cluster_membership$feature, cluster_membership$cluster)
#saveRDS(clusters, file= 'clusters.rds')

