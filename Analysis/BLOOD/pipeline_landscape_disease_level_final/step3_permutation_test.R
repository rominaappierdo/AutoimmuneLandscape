setwd("Analysis/BLOOD/pipeline_landscape_disease_level_final/")
library(dendextend)
set.seed(1234)

row_dend <- readRDS("row_dend.rds") #clustering of 700 features
nclusters <- 15

heigths <- sort(heights_per_k.dendrogram(as.dendrogram(row_dend)), T) #heights of splits
height <- as.numeric(heigths[names(heigths) == nclusters]) #height of split at ncluster

plot(row_dend, leaflab = "none")
abline(h=height,col="red")

input_mat <- read.csv("input_correlation_mat.tsv", sep = '\t', row.names = 1, check.names = F) #diseases on rows and 700 features on columns 
correlation_matrix <- Hmisc::rcorr(as.matrix(input_mat), type = 'spearman')$r

nperm <- 100000
out <- c()

for (i in c(1:nperm)) {
  print(i)
  
  
  upper_idx <- which(upper.tri(correlation_matrix))
  permuted_values <- sample(correlation_matrix[upper_idx])
 
  
  # Create new permuted symmetric matrix
  corr_perm <- matrix(0, nrow = ncol(input_mat), ncol = ncol(input_mat))
  rownames(corr_perm) <- colnames(input_mat)
  colnames(corr_perm) <- colnames(input_mat)
  
  corr_perm[upper.tri(corr_perm)] <- permuted_values
  corr_perm[lower.tri(corr_perm)] <- t(corr_perm)[lower.tri(corr_perm)]
  diag(corr_perm) <- 1  # correlation with self
  
  distance_matrix <- as.dist((1 - corr_perm)/2)
  
  row_dend = hclust(distance_matrix)
  numCluster <- length(rect.hclust(row_dend, h=height, border=0))
  out <- append(out, numCluster)
  
}

pval <- (sum(out <= nclusters) + 1) / (nperm + 1)
print(pval)

write.table(pval, file = 'pval.txt', row.names = F, col.names = F, quote = F, sep = '\n')

saveRDS(out, file = "res_permutation_100000_test_clusters_number.rds")

print('Done')


