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


if (!identical(colnames(expr), meta$geo_accession)){
  stop('different order')
}


pca.res <- prcomp(t(expr))
fviz_pca_ind(pca.res, habillage = meta$response, mean.point = F, geom = 'point') #gene expression no difference

list.plots<-list()
c <-0
for (module in list_modules){
  c<-c+1
  total_df <- compute_features(expr = expr, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  mod_name <- module_names[[c]]
  colnames(score_df) <- 'module'
  merged <- merge(score_df, meta, by = 'row.names')
  
  set.seed(123)
  #Convert the outcome variable to a factor with valid level names
  merged$class <- factor(merged$response, levels = c("Yes", "No"))
  
  # Define cross-validation control with 10 folds
  train_control <- trainControl(
    method = "cv",                # Cross-validation
    number = 10,                  # Number of folds
    classProbs = TRUE,            # Include class probabilities
    summaryFunction = twoClassSummary, # Evaluate using ROC, Sensitivity, Specificity
    savePredictions = "final"
  )
  
  # Train a logistic regression model with 10-fold cross-validation
  cv_model <- train(
    class ~ module,               # Formula
    data = merged,                # Dataset
    method = "glm",               # Logistic regression
    family = binomial,            # Specify binomial family
    trControl = train_control,    # Cross-validation settings
    metric = "ROC"                # Optimize for ROC
  )
  
  # Print the cross-validation results
  print(cv_model)
  
  cat("Average ROC (10-fold CV):", max(cv_model$results$ROC, na.rm = TRUE), "\n")
 
  
  cv_model$trainingData$probabilities <- predict(cv_model$finalModel, newdata = cv_model$trainingData, type = "response")
  
  cv_model$trainingData$.outcome <- factor(cv_model$trainingData$.outcome, levels = c('No', 'Yes'))
  roc_curve <- roc(cv_model$trainingData$.outcome, cv_model$trainingData$probabilities)
  plot(roc_curve, main = "ROC Curve")
  auc(roc_curve)  # Compute AUC
  
  
  p <- ggroc(roc_curve) +
    ggtitle(mod_name) +
    xlab("False Positive Rate (1 - Specificity)") +
    ylab("True Positive Rate (Sensitivity)") +
    theme_classic()
  
  list.plots[[c]] <- roc_curve
}


ggroc(list.plots) +
  ggtitle("Multiple ROC Curves") +
  theme_minimal() +
  scale_color_manual(values = rainbow(length(list.plots)),
                                     labels = paste0(names(list.plots), " (AUC = ", round(auc_values, 3), ")"))
  #labs(color = "Model") 

stop()



predictions <- cv_model$pred

folds <- unique(cv_model$pred$Resample)


plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "ROC Curves for 10-fold CV")
abline(a = 0, b = 1, lty = 2, col = "gray")



for (fold in folds) {
  fold_predictions <- predictions[predictions$Resample == fold, ]
  pred_prob <- cv_model$trainingData$probabilities[predictions$Resample == fold]
  
  roc_obj <- roc(fold_predictions$obs, pred_prob, levels = rev(levels(fold_predictions$obs)))
  lines(roc_obj, col = "blue", lwd = 1.5, lty = 2)
}


