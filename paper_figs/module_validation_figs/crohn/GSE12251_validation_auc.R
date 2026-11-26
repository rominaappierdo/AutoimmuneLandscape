setwd("/mnt/data3/appierdo/PhD-project/paper/paper_figs/module_validation_figs/crohn")
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
library(caret)

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


ggroc(roc_obj) +
  ggtitle("GSE12251 validation - UC - before treatment Inflamm Module - response label") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme_classic() +
  annotate("text", x = 0.95, y = 1, 
           label = paste("AUC =", round(auc_value, 3)), 
           size = 5, color = "blue")


#increase text
ggroc(roc_obj) +
  ggtitle("GSE12251 validation - UC - before treatment Inflammation Module - response label") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme_classic(base_size = 18) +  # increase base text size
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold")
  ) +
  annotate(
    "text",
    x = 0.8, y = 1,
    label = paste("AUC =", round(auc_value, 3)),
    size = 9,       # increased annotation text size
    color = "blue"
  )

stop()


# 
# set.seed(123)
# #Convert the outcome variable to a factor with valid level names
# merged$class <- factor(merged$response, levels = c("Yes", "No"))
# 
# # Define cross-validation control with 10 folds
# train_control <- trainControl(
#   method = "cv",                # Cross-validation
#   number = 10,                  # Number of folds
#   classProbs = TRUE,            # Include class probabilities
#   summaryFunction = twoClassSummary, # Evaluate using ROC, Sensitivity, Specificity
#   savePredictions = "final"
# )
# 
# # Train a logistic regression model with 10-fold cross-validation
# cv_model <- train(
#   class ~ module,               # Formula
#   data = merged,                # Dataset
#   method = "glm",               # Logistic regression
#   family = binomial,            # Specify binomial family
#   trControl = train_control,    # Cross-validation settings
#   metric = "ROC"                # Optimize for ROC
# )
# 
# # Print the cross-validation results
# print(cv_model)
# 
# cat("Average ROC (10-fold CV):", max(cv_model$results$ROC, na.rm = TRUE), "\n")
# 
# 
# cv_model$trainingData$probabilities <- predict(cv_model$finalModel, newdata = cv_model$trainingData, type = "response")
# 
# cv_model$trainingData$.outcome <- factor(cv_model$trainingData$.outcome, levels = c('No', 'Yes'))
# roc_curve <- roc(cv_model$trainingData$.outcome, cv_model$trainingData$probabilities)
# plot(roc_curve, main = "ROC Curve")
# auc(roc_curve)  # Compute AUC
# 
# 
# ggroc(roc_curve) +
#   ggtitle(paste0("GSE12251 - ROC Curve\nInflammation module - AUC: ", round(auc(roc_curve), 3))) +
#   xlab("False Positive Rate (1 - Specificity)") +
#   ylab("True Positive Rate (Sensitivity)") +
#   theme_classic()
# 
# stop()
# 
# 
# 
# predictions <- cv_model$pred
# 
# folds <- unique(cv_model$pred$Resample)
# 
# 
# plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
#      xlab = "False Positive Rate", ylab = "True Positive Rate",
#      main = "ROC Curves for 10-fold CV")
# abline(a = 0, b = 1, lty = 2, col = "gray")
# 
# 
# 
# for (fold in folds) {
#   fold_predictions <- predictions[predictions$Resample == fold, ]
#   pred_prob <- cv_model$trainingData$probabilities[predictions$Resample == fold]
#   
#   roc_obj <- roc(fold_predictions$obs, pred_prob, levels = rev(levels(fold_predictions$obs)))
#   lines(roc_obj, col = "blue", lwd = 1.5, lty = 2)
# }
# 
# 
