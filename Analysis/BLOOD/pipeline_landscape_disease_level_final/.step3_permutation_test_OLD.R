setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final")
library(dendextend)
set.seed(3333)

row_dend <- readRDS("row_dend.rds")
nclusters <- 15

heigths <- sort(heights_per_k.dendrogram(as.dendrogram(row_dend)), T)
height <- as.numeric(heigths[names(heigths) == nclusters])

#plot(row_dend, labels=FALSE)
plot.new()
nperm <- 100000
out <- c()
for (i in c(1:nperm)) {
  #print(i)
  correlation_matrix <- read.csv("input_correlation_mat.tsv", sep = '\t', row.names = 1, check.names = F)
  correlation_matrix<- Hmisc::rcorr(as.matrix(correlation_matrix), type = 'spearman')$r
  #shuffle here
  shuffled_vector <- sample(c(correlation_matrix))
  shuffled_mat <- matrix(shuffled_vector,nrow = nrow(correlation_matrix),ncol = ncol(correlation_matrix))
  distance_matrix <- as.dist((1 - shuffled_mat)/2)
  row_dend = hclust(distance_matrix)
  numCluster <- length(rect.hclust(row_dend, h=height, border=0))
  out <- append(out, numCluster)
}

saveRDS(out, file = "res_permutation_100000_test_clusters_number.rds")

print('Done')

