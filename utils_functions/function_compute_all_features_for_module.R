library(GSVA)
source("/mnt/data3/appierdo/PhD-project/utils_functions/deconvolution_with_immunoStates.R")

compute_features <- function(expr, module, deconvolution = T, ncore = 8) {
  counter_feature_set <- 0
  list_all_features <- list()
  
  cyto <- expr[paste0(rownames(expr), "_cyto") %in% module,]
  cyto <- scale(cyto)
  if (nrow(cyto) >0) {
    counter_feature_set <- counter_feature_set+1
    list_all_features[[counter_feature_set]] <- cyto
  }
  
  print('doing deconvolution')
  if (deconvolution == T) {
    cell_prop <- as.data.frame(compute_cell_proportions(expr))
    cell_prop <- t(cell_prop[, colnames(cell_prop) %in% module, drop=F])
    if (nrow(cell_prop) > 0) {
      counter_feature_set <- counter_feature_set+1
      list_all_features[[counter_feature_set]] <- cell_prop
    }
  }
  print('done deconvolution')
  
  print("doing hallmark")
  hall <- readRDS("/mnt/data3/appierdo/PhD-project/gene_sets/hallmarks_list_GENE_SYMBOL.RData")
  hall <- hall[names(hall) %in% module]
  if (length(hall) >0) {
    counter_feature_set <- counter_feature_set+1
    hall_df <- t(scale(t(gsva(as.matrix(expr), gset.idx.list = hall, min.sz=5, parallel.sz=ncore, method = 'ssgsea', ssgsea.norm = F, verbose = F))))
    list_all_features[[counter_feature_set]] <- hall_df
  }
  print('done hallmark')
  
  print("doing kegg")
  kegg <- readRDS("/mnt/data3/appierdo/PhD-project/gene_sets/KEGG_list_GENE_SYMBOL.RData")
  kegg <- kegg[names(kegg) %in% module]
  
  if (length(kegg) >0) {
    counter_feature_set <- counter_feature_set+1
    kegg_df <- t(scale(t(gsva(as.matrix(expr), gset.idx.list = kegg, min.sz=5, parallel.sz=ncore, method = 'ssgsea', ssgsea.norm = F, verbose = F))))
    list_all_features[[counter_feature_set]] <- kegg_df
  }
  print("done kegg")
  
  print("doing wp")
  btm <- readRDS("/mnt/data3/appierdo/PhD-project/gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData")
  btm <- btm[names(btm) %in% module]
  if (length(btm) >0) {
    counter_feature_set <- counter_feature_set+1
    btm_df <- t(scale(t(gsva(as.matrix(expr), gset.idx.list = btm, min.sz=5, parallel.sz=ncore, method = 'ssgsea', ssgsea.norm = F, verbose = F))))
    list_all_features[[counter_feature_set]] <- btm_df
  }
  print('done wp')
  print("doing tf")
  tf <- readRDS("/mnt/data3/appierdo/PhD-project/gene_sets/collecTRI_regulons_dorothea_08122023.rds")
  tf <- tf[paste0(names(tf), "_tf") %in% module]
  if (length(tf) >0) {
    counter_feature_set <- counter_feature_set+1
    tf_df <- t(scale(t(gsva(as.matrix(expr), gset.idx.list = tf, min.sz=5, parallel.sz=ncore, method = 'ssgsea', ssgsea.norm = F, verbose = F))))
    list_all_features[[counter_feature_set]] <- tf_df
  }
  print("done tf")
  
  print("doing mirna")
  mirna <- readRDS("/mnt/data3/appierdo/PhD-project/gene_sets/miRTarBase_regulons.rds")
  mirna <- mirna[names(mirna) %in% module]
  if (length(mirna) >0) {
    counter_feature_set <- counter_feature_set+1
    mirna_df <- t((scale(t(gsva(as.matrix(expr), gset.idx.list = mirna, min.sz=5, parallel.sz=ncore, method = 'ssgsea', ssgsea.norm = F, verbose = F)))))* (-1)
    list_all_features[[counter_feature_set]] <- mirna_df
  }
  print("done mirna")
  
  total_df <- as.data.frame(do.call("rbind", list_all_features))
  
  return(total_df)
}


