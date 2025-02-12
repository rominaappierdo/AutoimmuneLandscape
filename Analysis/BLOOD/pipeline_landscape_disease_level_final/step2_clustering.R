setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final")
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
p<- draw(p, padding = unit(c(2, 2, 2, 2), "mm"))


saveRDS(row_dend, file= 'row_dend.rds')

cluster_membership <- as.data.frame(cutree(column_dend, k= nclusters))
colnames(cluster_membership) <- c('cluster')
cluster_membership$feature <- rownames(cluster_membership)
clusters<- split(cluster_membership$feature, cluster_membership$cluster)
saveRDS(clusters, file= 'clusters.rds')


stop()



#transpose
nclusters <- 15
column_dend = hclust(distance_matrix)
row_dend = hclust(distance_matrix)
row_dend = color_branches(row_dend, k = nclusters)


p3 <- ComplexHeatmap::Heatmap(as.matrix(correlation_matrix),
                             show_row_names = F,
                             cluster_rows = row_dend, 
                             cluster_columns = column_dend,
                             column_names_rot = 50,
                             row_split=nclusters,
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

p3<- draw(p3, padding = unit(c(2, 2, 2, 2), "mm"))


stop()
#save cluster membership
cluster_membership <- as.data.frame(cutree(column_dend, k= nclusters))
colnames(cluster_membership) <- c('cluster')
cluster_membership$feature <- rownames(cluster_membership)
head(cluster_membership)

clusters<- split(cluster_membership$feature, cluster_membership$cluster)
saveRDS(clusters, file= 'clusters.rds')




firsts<-c()
for (i in 1:length(clusters)) {
  firsts<- append(firsts, paste0(clusters[[i]][1], '_', i))
}

correlation_matrix2 <- correlation_matrix

new_colnames<- c()
for (c in colnames(correlation_matrix2)){

  if (any(grepl(c, firsts))) {
    new_colnames <- append(new_colnames, firsts[which(grepl(c, firsts))])
  } else {
    new_colnames <- append(new_colnames, "")
  }
}

colnames(correlation_matrix2) <- new_colnames
p <- ComplexHeatmap::Heatmap(as.matrix(correlation_matrix2),
                             show_row_names = F,
                             cluster_rows = row_dend, 
                             cluster_columns = column_dend,
                             column_names_rot = 50,
                             column_split=nclusters,
                             col = circlize::colorRamp2(c(-0.9, 0, 0.9), c("#56B4E9", "#FFFFFF", "#FF7400")),
                             row_names_gp = grid::gpar(fontsize = 10),
                             column_names_gp = grid::gpar(fontsize = 10),
                             show_heatmap_legend = F,
                             heatmap_legend_param = list(title = 'correlation'),
                             na_col = 'black',
                             show_column_names = T,
                             heatmap_width = unit(0.25, "mm")*ncol(correlation_matrix),
                             heatmap_height = unit(0.25, "mm")*ncol(correlation_matrix)
)
p<- draw(p, padding = unit(c(2, 2, 2, 2), "mm"))
stop()
##########plot clusters on effect size matrix
matrix <- read.csv("mat.tsv", sep = '\t', row.names = 1, check.names = F)
ComplexHeatmap::Heatmap(as.matrix(matrix),
                        cluster_columns = column_dend,
                        column_split = 16,
                        show_column_names = F)

###########

clusters[12]
