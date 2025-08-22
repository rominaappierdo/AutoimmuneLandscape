setwd("Supplementary/supp_figures")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
library(caret)
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
source("../../utils_functions/function_compute_all_features_for_module.R")
source("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")


gse <- 'GSE12251'

########## set module
list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <-  c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                   "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/Transendothelial_migration", 
                   "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/Innate_immune_response", 
                   "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired_cell_proliferation", "Lipid metabolism")

mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")


########## get expression data
data <- getGEOData(gse)
expr <- data$originalData$GSE12251$expr
map_probes <- (as.data.frame(data$originalData$GSE12251$keys))
colnames(map_probes) <- 'keys'
expr <- merge(expr, map_probes, by = 'row.names')
expr$Row.names <- NULL
expr <- na.omit(expr)
expr <- mean_expression_gene_rowise(expr, colname = 'keys')
expr <- order_by_col_names(expr)

########### get metadata

meta <- GEOquery::getGEO(gse)
meta <- meta$GSE12251_series_matrix.txt.gz@phenoData@data
meta <- sort_by_col(meta, colname = 'geo_accession')

meta$response <- meta$`WK8RSPHM:ch1`
meta$disease <- rep('UC', nrow(meta)) 
meta$tissue <- rep('Colon', nrow(meta)) 
meta$treatment <- rep('before', nrow(meta)) 



if (!identical(colnames(expr), meta$geo_accession)){
  stop('different order')
}


pca.res <- prcomp(t(expr))
fviz_pca_ind(pca.res, habillage = meta$response, mean.point = F, geom = 'point') #gene expression no difference


module <- list_modules[[13]]
total_df <- compute_features(expr = expr, module = module, deconvolution = T)
total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
total_df <- na.omit(total_df)
score_df <- as.data.frame(colSums(total_df))
mod_name <- 'Inflammation_module'
colnames(score_df) <- mod_name
merged <- merge(score_df, meta, by = 'row.names')
merged$response <- factor(merged$response, levels = c("No", "Yes"))



roc_obj <- roc(merged$response , merged$Inflammation_module, direction = ">")
auc_value <- auc(roc_obj)
auc_value


p <- ggroc(roc_obj) +
  #ggtitle("GSE12251 validation - UC - before treatment Inflamm Module - response label") +
  ggtitle("GSE12251 Validation Dataset") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme_classic() +
  annotate("text", x = 0.95, y = 1, 
           label = paste("AUC =", round(auc_value, 3)), 
           size = 5, color = "blue")


jpeg(file = 'SuppFig8A.jpeg', width = 18 , height = 15, units = 'cm', res = 600)
p 
dev.off()
