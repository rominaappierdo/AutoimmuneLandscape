setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig8")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
source("../../utils_functions/function_compute_all_features_for_module.R")
source("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")


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
expr <- expr[, !colnames(expr) == 'GSM364647']

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
pheno <- pheno[!pheno$geo_accession == 'GSM364647', ]
if (!identical(colnames(expr), pheno$geo_accession)){
  stop()
}

pca.res <- prcomp(t(expr))

fviz_pca_ind(pca.res, habillage = pheno$tissue, mean.point = F, geom = 'point', pointsize=5)
fviz_pca_ind(pca.res, habillage = pheno$disease, mean.point = F, geom = 'point',  pointsize=5)
fviz_pca_ind(pca.res, habillage = pheno$response, mean.point = F, geom = 'point',  pointsize=5)
fviz_pca_ind(pca.res, habillage = pheno$treatment, mean.point = F, geom = 'point')

#QUESTION #1: Can modules capture treatment effect? 

pheno_colon <- pheno[pheno$tissue == 'Colon', ]
#get only the AFTER treatment ones
pheno_colon <- pheno_colon[pheno_colon$treatment %in% c("control", "after"), ]
expr_colon <- expr[, colnames(expr) %in% pheno_colon$geo_accession]

pca.res.after <- prcomp(t(expr_colon))
fviz_pca_ind(pca.res.after, habillage = pheno_colon$response, mean.point = F, geom = 'point', pointsize=5)


#only inflammation module
total_df <- compute_features(expr = expr_colon, module = list_modules[[13]], deconvolution = T)
total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
total_df <- na.omit(total_df)
score_df <- as.data.frame(colSums(total_df))
colnames(score_df) <- "module"
mod_name <-  module_names[[13]]

merged <- merge(score_df, pheno_colon, by = 'row.names')
merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))

p1 <- ggplot(merged, aes(x = disease, y = module, fill = response))+
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c('bisque', '#add98e', '#ff800e'), name = 'Response', labels = c('Control', 'Non responder', 'Responder')) +
  stat_compare_means(aes(group = response), 
                     method = "wilcox.test",
                     label = "p.format") +
  ggtitle(paste0(mod_name, " | ", "After-treatment samples"))+
  ylab("module score") +
  xlab("")

p1


############################baseline, predictive potential

pheno_colon <- pheno[pheno$tissue == 'Colon', ]

#get only the BEFORE treatment ones
pheno_colon <- pheno_colon[pheno_colon$treatment %in% c("control", "before"), ]
expr_colon <- expr[, colnames(expr) %in% pheno_colon$geo_accession]

pca.res.before <- prcomp(t(expr_colon))


total_df <- compute_features(expr = expr_colon, module = list_modules[[13]], deconvolution = T)
total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
total_df <- na.omit(total_df)
score_df <- as.data.frame(colSums(total_df))
colnames(score_df) <- "module"
mod_name <- module_names[[13]]

merged <- merge(score_df, pheno_colon, by = 'row.names')
merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))

p2 <- ggplot(merged, aes(x = disease, y = module, fill = response))+
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(aes(group = response), 
                     method = "wilcox.test",
                     label = "p.format") +
  ggtitle(paste0(mod_name, " | ", "Before-treatment samples")) +
  scale_fill_manual(values = c('bisque', '#add98e', '#ff800e'), name = 'Response', labels = c('Control', 'Non responder', 'Responder'))

p2




jpeg("fig8A.jpeg", width = 16,height = 14, units = "cm", res = 600)
p1
dev.off()


jpeg("fig8B.jpeg", width = 16,height = 14, units = "cm", res = 600)
p2
dev.off()


stop()

#Interestingly, module activity in responders approached levels observed in healthy controls (actually non significant)


pheno_colon <- pheno[pheno$tissue == 'Colon', ]

#get only the After treatment ones
pheno_colon <- pheno_colon[pheno_colon$treatment %in% c("control", "after"), ]
#get only controls and responders
pheno_colon <- pheno_colon[pheno_colon$response %in% c("control", "T"), ]
expr_colon <- expr[, colnames(expr) %in% pheno_colon$geo_accession]


total_df <- compute_features(expr = expr_colon, module = list_modules[[13]], deconvolution = T)
total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
total_df <- na.omit(total_df)
score_df <- as.data.frame(colSums(total_df))
colnames(score_df) <- "module"
mod_name <- module_names[[13]]

merged <- merge(score_df, pheno_colon, by = 'row.names')
merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))


comparisons_list <- list(c("Control", "UC"), c("Control", "CD"), c("UC", "CD"))


p3 <- ggplot(merged, aes(x = disease, y = module, fill = response))+
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(comparisons = comparisons_list,
                     method = "t.test",
                     label = "p.format") +
  ggtitle(paste0(mod_name, " | ", "After-treatment samples")) +
  scale_fill_manual(values = c('bisque', '#add98e', '#ff800e'), name = 'Response')

p3
