setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/module_validation_figs/crohn")
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
#COLON
pheno_colon <- pheno[pheno$tissue == 'Colon', ]

#get only the AFTER treatment ones
pheno_colon <- pheno_colon[pheno_colon$treatment %in% c("control", "after"), ]
expr_colon <- expr[, colnames(expr) %in% pheno_colon$geo_accession]

table(pheno_colon$disease)

pca.res.colon <- prcomp(t(expr_colon))
fviz_pca_ind(pca.res.colon, habillage = pheno_colon$disease, mean.point = F, geom = 'point', pointsize = 4)#, addEllipses = T)
fviz_pca_ind(pca.res.colon, habillage = pheno_colon$response, mean.point = F, geom = 'point', pointsize = 4)#, addEllipses = T)

#fviz_pca_ind(pca.res.colon, habillage = pheno_colon$treatment, mean.point = F, geom = 'point', pointsize = 4, addEllipses = T)


#I want to compute only the significant ones for before treatment

#sig_ones <- c(6, 8, 9, 13, 15) #both significant


list_plot <- list()
c <- 0
j <- 0
for (module in list_modules) {
  c <- c+1
  j <- j+1
  
  total_df <- compute_features(expr = expr_colon, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  mod_name <- module_names[[c]]
  
  merged <- merge(score_df, pheno_colon, by = 'row.names')
  merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
  merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))
  
  p <- ggplot(merged, aes(x = disease, y = module, fill = response))+
    geom_boxplot() +
    theme_classic() +
    #scale_fill_tableau(palette = "Classic Green-Orange 12", name = 'Response')+
    scale_fill_manual(values = c('bisque', '#add98e', '#ff800e'), name = 'Response', labels = c('Control', 'Non responder', 'Responder')) +
    stat_compare_means(aes(group = response), 
                       method = "t.test",
                       label = "p.format") +
    ggtitle(mod_name) +
    ylab("module score") +
    xlab("")# +
    #theme(legend.text=element_text(size=18), 
     #     legend.title=element_text(size=18)
      #    )#legend.position = 'none')

  #legend
  # p <- p+ theme(legend.text=element_text(size=16), legend.title = element_text(size=16), legend.position="bottom")
  # 
  # jpeg(file = 'legend_colon_after_treatment.jpeg', units = 'cm', width = 20, height = 20, res = 600)
  # p
  # dev.off()
  
  list_plot[[j]] <- p
  
}

jpeg(file = 'boxplot_all_modules_colon.jpeg', units = 'cm', width = 57, height = 30, res = 200)
do.call(grid.arrange,c(list_plot, ncol = 5))
dev.off()

#saveRDS(list_plot, 'list_plot_before_treat_colon.rds')


#sig at least 1

list_plot_new <- list()
list_plot_new[[1]] <- list_plot[[13]]
list_plot_new[[2]] <- list_plot[[4]]
list_plot_new[[3]] <- list_plot[[5]]
list_plot_new[[4]] <- list_plot[[9]]
list_plot_new[[5]] <- list_plot[[11]]
list_plot_new[[6]] <- list_plot[[12]]
list_plot_new[[7]] <- list_plot[[15]]
list_plot_new[[8]] <- list_plot[[1]]
list_plot_new[[9]] <- list_plot[[2]]
list_plot_new[[10]] <- list_plot[[3]]
list_plot_new[[11]] <- list_plot[[6]]
list_plot_new[[12]] <- list_plot[[7]]
list_plot_new[[13]] <- list_plot[[8]]
list_plot_new[[14]] <- list_plot[[14]]

for (i in 1:length(list_plot_new)) {
  plot <- list_plot_new[[i]]
  
  new_plot <-  plot + 
    theme(
      axis.title.y = element_text(size = 18),
      plot.title = ggtext::element_markdown(size = 19, hjust = 0.5, face = 'bold', margin=margin(0,0,20,0)),
      axis.text.x = element_text(size = 18)
      
    )
  
  new_plot$layers[[2]]$aes_params$textsize <- 14 #not working
  
  
  list_plot_new[[i]] <- new_plot
  
}


jpeg(file = 'boxplot_sig_modules_colon_after_treatment.jpeg', units = 'cm', width = 68, height = 47, res = 600)
do.call(grid.arrange,c(list_plot_new, ncol = 4))
dev.off()

jpeg(file = 'boxplot_sig_modules_colon_after_treatment_only_six_for_phdmeeting.jpeg', units = 'cm', width = 50, height = 35, res = 600)
do.call(grid.arrange,c(list_plot_new[1:6], ncol = 3))
dev.off()


################################################################################################################################################
#Question #2: Can modules capture differences in response at baseline to predict response to treatment?

#we have to split the analysis in colon vs ileum
#COLON
pheno_colon <- pheno[pheno$tissue == 'Colon', ]

#get only the BEFORE treatment ones
pheno_colon <- pheno_colon[pheno_colon$treatment %in% c("control", "before"), ]

expr_colon <- expr[, colnames(expr) %in% pheno_colon$geo_accession]

table(pheno_colon$disease)

pca.res.colon <- prcomp(t(expr_colon))
#fviz_pca_ind(pca.res.colon, habillage = pheno_colon$disease, mean.point = F, geom = 'point', pointsize = 4)#, addEllipses = T)
#fviz_pca_ind(pca.res.colon, habillage = pheno_colon$response, mean.point = F, geom = 'point', pointsize = 4, addEllipses = T)
#fviz_pca_ind(pca.res.colon, habillage = pheno_colon$treatment, mean.point = F, geom = 'point', pointsize = 4, addEllipses = T)


#I want to compute only the significant ones for before treatment
#sig_ones <- c(6, 7, 8, 9, 11, 12, 13, 14, 15) #at least 1
#sig_ones <- c(6, 8, 9, 13, 15) #both significant


list_plot <- list()
c <- 0
j <- 0
for (module in list_modules) {
  c <- c+1
  # if (!c %in% sig_ones) {
  #   next
  # }
  j <- j+1
  
  total_df <- compute_features(expr = expr_colon, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  mod_name <- module_names[[c]]
  
  merged <- merge(score_df, pheno_colon, by = 'row.names')
  merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
  merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))
  
  p <- ggplot(merged, aes(x = disease, y = module, fill = response))+
    geom_boxplot() +
    theme_classic() +
    stat_compare_means(aes(group = response), 
                       method = "t.test",
                       label = "p.format") +
    ggtitle(mod_name) +
    scale_fill_manual(values = c('bisque', '#add98e', '#ff800e'), name = 'Response')
  
  list_plot[[j]] <- p
  
}


#saveRDS(list_plot, 'list_plot_colon_before.rds')
list_plot <- readRDS("list_plot_colon_before.rds")
list_plot_new <- list()
list_plot_new[[1]] <- list_plot[[13]]
list_plot_new[[2]] <- list_plot[[9]]
list_plot_new[[3]] <- list_plot[[11]]
list_plot_new[[4]] <- list_plot[[12]]
list_plot_new[[5]] <- list_plot[[15]]
list_plot_new[[6]] <- list_plot[[6]]
list_plot_new[[7]] <- list_plot[[7]]
list_plot_new[[8]] <- list_plot[[8]]
list_plot_new[[9]] <- list_plot[[14]]

for (i in 1:length(list_plot_new)) {
  plot <- list_plot_new[[i]]
  
  new_plot <-  plot + 
    theme(
      axis.title.y = element_text(size = 22),
      plot.title = ggtext::element_markdown(size = 19, hjust = 0.5, face = 'bold', margin=margin(0,0,20,0)),
      axis.text.x = element_text(size = 18)
      
    ) +
    xlab('')
  
  #new_plot$layers[[2]]$aes_params$textsize <- 12 #not working
  
  
  list_plot_new[[i]] <- new_plot
  
}


jpeg(file = 'boxplot_sig_modules_colon_before_treatment.jpeg', units = 'cm', width = 35, height = 20, res = 600)
do.call(grid.arrange,c(list_plot_new, ncol = 4))
dev.off()

jpeg(file = 'boxplot_sig_modules_colon_before_treatment_only4_phdmeeting.jpeg', units = 'cm', width = 55, height = 40, res = 600)
do.call(grid.arrange,c(list_plot_new[1:4], ncol = 2))
dev.off()


stop()

################################################################################################
#AUC Inflammation module on both after and before treatment

pheno_colon <- pheno[pheno$tissue == 'Colon', ]

#BEFORE

#get only the BEFORE treatment ones
pheno_colon_before <- pheno_colon[pheno_colon$treatment %in% c("control", "before"), ]
expr_colon_before <- expr[, colnames(expr) %in% pheno_colon_before$geo_accession]


#I want to have a matrix which is 15 x number of samples
dim(expr_colon_before)
dim(pheno_colon_before)
identical(colnames(expr_colon_before), pheno_colon_before$geo_accession)

#compute modules scores
c <- 0
for (module in list_modules[13]) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_colon_before, module = module, deconvolution = T)
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

merged_before <- merge(final_score_df, pheno_colon_before, by = 'row.names')
merged_before <- merged_before[merged_before$tissue =='Colon',]
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
pheno_colon_after <- pheno_colon[pheno_colon$treatment %in% c("control", "after"), ]
expr_colon_after <- expr[, colnames(expr) %in% pheno_colon_after$geo_accession]

#I want to have a matrix which is 15 x number of samples
dim(expr_colon_after)
dim(pheno_colon_after)
identical(colnames(expr_colon_after), pheno_colon_after$geo_accession)

#compute modules scores
c <- 0
for (module in list_modules[13]) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_colon_after, module = module, deconvolution = T)
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

merged_after <- merge(final_score_df, pheno_colon_after, by = 'row.names')
merged_after <- merged_after[merged_after$tissue =='Colon',]
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

stop()


#plot before and after treatment

pheno_colon <- pheno[pheno$tissue == 'Colon', ]
pheno_after_before <- pheno_colon[pheno_colon$treatment %in% c("before", "after"), ]
expr_after_before <- expr[, colnames(expr) %in% pheno_after_before$geo_accession]

identical(colnames(expr_after_before), pheno_after_before$geo_accession)

c <- 0
for (module in list_modules[13]) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_after_before, module = module, deconvolution = T)
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


merged_after_before <- merge(final_score_df, pheno_after_before, by = 'row.names')
merged_after_before$treatment <- factor(merged_after_before$treatment, levels = c('before', 'after'))

line_plot <- ggplot(merged_after_before, aes(y = Inflammation_module, x = treatment, group = patient_id, color = response))+
  geom_point(size = 4) +
  theme_classic()+
  geom_line() +
  ggtitle("Colon") +
  scale_color_tableau(palette = "Nuriel Stone", name = 'Response') +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)
        )+
  xlab("Treatment")+
  ylab("Inflammation Module")

jpeg(file = 'line_plot_inflammation_module_colon.jpeg', units = 'cm', width = 20, height = 20, res = 300)
line_plot
dev.off()


merged_after_before$treatment <- varhandle::unfactor(merged_after_before$treatment)
merged_after_before$treatment[merged_after_before$treatment == 'before'] <- 'Before'
merged_after_before$treatment[merged_after_before$treatment == 'after'] <- 'After'
merged_after_before$treatment <- factor(merged_after_before$treatment, levels = c('Before', 'After'))

ggplot(merged_after_before, aes(y = Inflammation_module, x = disease, fill = response)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~treatment)+
  theme(strip.text.x = element_text(size = 14, face = 'bold'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_fill_manual(values = c("#add98e", "#ff800e"), name = 'Response')+
  xlab("")+
  ylab("Inflammation module") +
  stat_compare_means(aes(group = response), 
                     method = "t.test",
                     label = "p.format")


b <- ggplot(merged_after_before, aes(y = Inflammation_module, x = treatment, fill = response)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~disease)+
  theme(strip.text.x = element_text(size = 14, face = 'bold'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_fill_manual(values = c("#add98e", "#ff800e"), name = 'Response')+
  xlab("")+
  ylab("Inflammation module") +
  stat_compare_means(aes(group = response), 
                     method = "t.test",
                     label = "p.format")

jpeg(file = 'boxplot_before_after_treat_comparison_colon_facet_disease.jpeg', units = 'cm', width = 25, height = 20, res = 300)
b
dev.off()


c <- ggplot(merged_after_before, aes(y = Inflammation_module, x = disease, fill = response)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~treatment)+
  theme(strip.text.x = element_text(size = 14, face = 'bold'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_fill_manual(values = c("#add98e", "#ff800e"), name = 'Response')+
  xlab("")+
  ylab("Inflammation module") +
  stat_compare_means(aes(group = response), 
                     method = "t.test",
                     label = "p.format")

jpeg(file = 'boxplot_before_after_treat_comparison_colon_facet_treatment.jpeg', units = 'cm', width = 25, height = 20, res = 300)
c
dev.off()


# # Caricare i pacchetti necessari
# library(randomForest)
# library(caret)
# 
# # Simulare il dataset di addestramento
# set.seed(123)
# train_data <- merged[, c(1:5, 50)]
# 
# # Define cross-validation parameters
# control <- trainControl(
#   method = "cv",        # Cross-validation
#   number = 5,           # Number of folds
#   verboseIter = TRUE,   # Show progress
#   classProbs = TRUE,    # Compute class probabilities
#   summaryFunction = twoClassSummary # Use AUC as the metric
# )
# 
# # Train Random Forest model with cross-validation
# 
# rf_model_cv <- train(
#   label ~ ., 
#   data = train_data,
#   method = "rf",        # Random Forest
#   metric = "ROC",       # Use AUC for performance
#   trControl = control,
#   tuneLength = 3        # Number of hyperparameter combinations to try
# )
# 
# # Print model results
# print(rf_model_cv)
# 
# # Extract the best model
# best_model <- rf_model_cv$finalModel
# # Creare il modello di Random Forest
# #rf_model <- randomForest(label ~ ., data = train_data, importance = TRUE, ntree = 100)
# 
# saveRDS(best_model, file='rf_model_cv.rds')
# 
# stop()
#  
# library(caret)
# #tidymodels
# set.seed(123)
# #Convert the outcome variable to a factor with valid level names
# merged$class <- as.factor(ifelse(grepl("T", merged$response), "Class1", "Class2"))
# 
# # Define cross-validation control with 10 folds
# train_control <- trainControl(
#   method = "cv",                # Cross-validation
#   number = 10,                  # Number of folds
#   classProbs = TRUE,            # Include class probabilities
#   summaryFunction = twoClassSummary # Evaluate using ROC, Sensitivity, Specificity
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
# # Evaluate performance
# cat("Average ROC (10-fold CV):", max(cv_model$results$ROC, na.rm = TRUE), "\n")
# 
# 
# 
# 
# 
# # ggplot(merged, aes(x = disease, y = module, fill = treatment))+
# #   geom_boxplot() +
# #   theme_classic()
