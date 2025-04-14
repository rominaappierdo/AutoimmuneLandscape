setwd("/mnt/data3/appierdo/PhD-project/paper/paper_figs/module_validation_figs/crohn/")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/PhD-project/utils_functions/function_compute_all_features_for_module.R")
source("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_validation/readLines_all_modules.R")


gse <- 'GSE16879'


########## set module
list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <-  c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                   "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/Transendothelial_migration", 
                   "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/Innate_immune_response", 
                   "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired_cell_proliferation", "Lipid metabolism")

mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")


########## get expression data
data <- getGEOData(gse)
expr <- data$originalData$GSE16879$expr
map_probes <- (as.data.frame(data$originalData$GSE16879$keys))
colnames(map_probes) <- 'keys'
expr <- merge(expr, map_probes, by = 'row.names')
expr$Row.names <- NULL
expr <- na.omit(expr)
expr <- mean_expression_gene_rowise(expr, colname = 'keys')
expr <- order_by_col_names(expr)

########### get metadata

meta <- GEOquery::getGEO(gse)
meta <- meta$GSE16879_series_matrix.txt.gz@phenoData@data
meta <- sort_by_col(meta, colname = 'geo_accession')



meta$tissue <- meta$`tissue:ch1`
meta$disease <- meta$`disease:ch1`
meta$treatment <- meta$`before or after first infliximab treatment:ch1`
meta$treatment <- ifelse(grepl("Not app", meta$treatment), 'control', ifelse(grepl('After', meta$treatment), 'after', 'before'))
meta$patient_id <- unlist(lapply(meta$title, function(x) strsplit(x, '_')[[1]][1]))
meta$response <- ifelse(grepl("Not app", meta$`response to infliximab:ch1`), 'control', ifelse(grepl('Yes', meta$`response to infliximab:ch1`), 'T', 'F'))


pheno <- meta[, c("geo_accession", "patient_id", "tissue", "disease", "treatment", "response")]

if (!identical(colnames(expr), pheno$geo_accession)){
  stop()
}

pca.res <- prcomp(t(expr))

fviz_pca_ind(pca.res, habillage = pheno$tissue, mean.point = F, geom = 'point', pointsize=5)
fviz_pca_ind(pca.res, habillage = pheno$disease, mean.point = F, geom = 'point',  pointsize=5)
fviz_pca_ind(pca.res, habillage = pheno$response, mean.point = F, geom = 'point',  pointsize=5)
fviz_pca_ind(pca.res, habillage = pheno$treatment, mean.point = F, geom = 'point')

#QUESTION #1: Can modules capture treatment effect? 

#we have to split the analysis in colon vs ileum
#ILEUM
pheno_ileum <- pheno[pheno$tissue == 'Ileum', ]

#get only the AFTER treatment ones
pheno_ileum <- pheno_ileum[pheno_ileum$treatment %in% c("control", "after"), ]

expr_ileum <- expr[, colnames(expr) %in% pheno_ileum$geo_accession]

table(pheno_ileum$disease)

pca.res.ileum <- prcomp(t(expr_ileum))
fviz_pca_ind(pca.res.ileum, habillage = pheno_ileum$disease, mean.point = F, geom = 'point', pointsize = 4)#, addEllipses = T)
fviz_pca_ind(pca.res.ileum, habillage = pheno_ileum$response, mean.point = F, geom = 'point', pointsize = 4, addEllipses = T)
#fviz_pca_ind(pca.res.ileum, habillage = pheno_ileum$treatment, mean.point = F, geom = 'point', pointsize = 4, addEllipses = T)



list_plot <- list()
c <- 0
j <- 0
for (module in list_modules) {
  c <- c+1
  # if (!c %in% sig_ones) {
  #   next
  # }
  j <- j+1
  
  total_df <- compute_features(expr = expr_ileum, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  mod_name <- module_names[[c]]
  
  merged <- merge(score_df, pheno_ileum, by = 'row.names')
  merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
  merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))
  
  p <- ggplot(merged, aes(x = disease, y = module, fill = response))+
    geom_boxplot() +
    theme_classic() +
    stat_compare_means(aes(group = response), 
                       method = "t.test",
                       label = "p.format") +
    ggtitle(mod_name)
  
  list_plot[[j]] <- p
  
}

do.call(grid.arrange,c(list_plot, ncol = 5))

stop()

################################################################################################################################################
#Question #2: Can modules capture differences in response at baseline to predict response to treatment?

#we have to split the analysis in colon vs ileum
#ILEUM
pheno_ileum <- pheno[pheno$tissue == 'Ileum', ]

#get only the AFTER treatment ones
pheno_ileum <- pheno_ileum[pheno_ileum$treatment %in% c("control", "before"), ]

expr_ileum <- expr[, colnames(expr) %in% pheno_ileum$geo_accession]

table(pheno_ileum$disease)

pca.res.ileum <- prcomp(t(expr_ileum))
fviz_pca_ind(pca.res.ileum, habillage = pheno_ileum$disease, mean.point = F, geom = 'point', pointsize = 4)#, addEllipses = T)
fviz_pca_ind(pca.res.ileum, habillage = pheno_ileum$response, mean.point = F, geom = 'point', pointsize = 4, addEllipses = T)
#fviz_pca_ind(pca.res.ileum, habillage = pheno_ileum$treatment, mean.point = F, geom = 'point', pointsize = 4, addEllipses = T)



list_plot <- list()
c <- 0
j <- 0
for (module in list_modules) {
  c <- c+1
  # if (!c %in% sig_ones) {
  #   next
  # }
  j <- j+1
  
  total_df <- compute_features(expr = expr_ileum, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  mod_name <- module_names[[c]]
  
  merged <- merge(score_df, pheno_ileum, by = 'row.names')
  merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
  merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))
  
  p <- ggplot(merged, aes(x = disease, y = module, fill = response))+
    geom_boxplot() +
    theme_classic() +
    stat_compare_means(aes(group = response), 
                       method = "t.test",
                       label = "p.format") +
    ggtitle(mod_name)
  
  list_plot[[j]] <- p
  
}

do.call(grid.arrange,c(list_plot, ncol = 5))


stop()




################################################################################################
#AUC Inflammation module on both after and before treatment

pheno_ileum <- pheno[pheno$tissue == 'Ileum', ]

#BEFORE

#get only the BEFORE treatment ones
pheno_ileum_before <- pheno_ileum[pheno_ileum$treatment %in% c("control", "before"), ]
expr_ileum_before <- expr[, colnames(expr) %in% pheno_ileum_before$geo_accession]


#I want to have a matrix which is 15 x number of samples
dim(expr_ileum_before)
dim(pheno_ileum_before)
identical(colnames(expr_ileum_before), pheno_ileum_before$geo_accession)

#compute modules scores
c <- 0
for (module in list_modules[13]) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_ileum_before, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  mod_name <- 'Inflammation_module'
  colnames(score_df) <- mod_name
  
  if (c == 1) {
    final_score_df <- score_df
  } else {
    final_score_df <- cbind(final_score_df, score_df)
  }
  
}

merged_before <- merge(final_score_df, pheno_ileum_before, by = 'row.names')
merged_before <- merged_before[merged_before$tissue =='Ileum',]
merged_before <- merged_before[!merged_before$response == 'control',]
rownames(merged_before) <- merged_before$Row.names
merged_before$Row.names <- NULL
merged_before$response <- factor(ifelse(merged_before$response, 'T', 'F'), levels = c('F','T'))


roc_obj <- roc(merged_before$response, merged_before$Inflammation_module)
auc_value <- auc(roc_obj)
print(auc_value)

ggroc(roc_obj) +
  ggtitle("ROC Curve - Before Treatment - Inflammation Module") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme_classic() +
  annotate("text", x = 0.98, y = 1, 
           label = paste("AUC:", round(auc_value, 3)), 
           size = 5, color = "blue")



#AFTER

#get only the AFTER treatment ones
pheno_ileum_after <- pheno_ileum[pheno_ileum$treatment %in% c("control", "after"), ]
expr_ileum_after <- expr[, colnames(expr) %in% pheno_ileum_after$geo_accession]

#I want to have a matrix which is 15 x number of samples
dim(expr_ileum_after)
dim(pheno_ileum_after)
identical(colnames(expr_ileum_after), pheno_ileum_after$geo_accession)

#compute modules scores
c <- 0
for (module in list_modules[13]) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_ileum_after, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  mod_name <- 'Inflammation_module'
  colnames(score_df) <- mod_name
  
  if (c == 1) {
    final_score_df <- score_df
  } else {
    final_score_df <- cbind(final_score_df, score_df)
  }
  
}

merged_after <- merge(final_score_df, pheno_ileum_after, by = 'row.names')
merged_after <- merged_after[merged_after$tissue =='Ileum',]
merged_after <- merged_after[!merged_after$response == 'control',]
rownames(merged_after) <- merged_after$Row.names
merged_after$Row.names <- NULL
merged_after$response <- factor(ifelse(merged_after$response, 'T', 'F'), levels = c('F','T'))


roc_obj <- roc(merged_after$response, merged_after$Inflammation_module)
auc_value <- auc(roc_obj)
print(auc_value)

ggroc(roc_obj) +
  ggtitle("ROC Curve - After Treatment - Inflammation Module") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme_classic() +
  annotate("text", x = 0.98, y = 1, 
           label = paste("AUC:", round(auc_value, 3)), 
           size = 5, color = "blue")

