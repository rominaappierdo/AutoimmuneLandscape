setwd("paper_figs/module_validation_figs/psoriasis_GSE117239")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/function_compute_all_features_for_module.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files//readLines_all_modules.R")


gse <- 'GSE117239'


########## set module
list_modules <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
module_names <-  c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment", 
                   "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/Transendothelial_migration", 
                   "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/Innate_immune_response", 
                   "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired_cell_proliferation", "Lipid metabolism")

mods <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12", "mod13", "mod14", "mod15")


########## get expression data
data <- getGEOData(gse)
expr <- data$originalData$GSE117239$expr
map_probes <- (as.data.frame(data$originalData$GSE117239$keys))
colnames(map_probes) <- 'keys'
expr <- merge(expr, map_probes, by = 'row.names')
expr$Row.names <- NULL
expr <- na.omit(expr)
expr <- mean_expression_gene_rowise(expr, colname = 'keys')
expr <- order_by_col_names(expr)

########### get metadata

meta <- GEOquery::getGEO(gse)
meta <- meta$GSE117239_series_matrix.txt.gz@phenoData@data
meta <- sort_by_col(meta, colname = 'geo_accession')



meta$tissue <- meta$`tissue:ch1` #non-lesional/lesional
meta$disease <- meta$`condition:ch1` #psoriasis
meta$timepoint <- meta$`timepoint:ch1` #baseline week1 week12
meta$treatment <- meta$`treatment:ch1`
meta$patient_id <- unlist(lapply(meta$characteristics_ch1, function(x) gsub("subjectid: ", "", x)[[1]][1]))
meta$response <- meta$`pasi75 response at week 12:ch1`
meta$pasi <- meta$`pasi:ch1`
meta$pasi[meta$pasi == 'NA'] <- 'baseline'

pheno <- meta[, c("geo_accession", "patient_id", "tissue", "disease", "treatment", "response", "timepoint", "pasi")]

if (!identical(colnames(expr), pheno$geo_accession)){
  stop()
}

length(unique(meta$patient_id)) #89 unique patients

pca.res <- prcomp(t(expr))


fviz_pca_ind(pca.res, habillage = pheno$tissue, mean.point = F, geom = 'point')
fviz_pca_ind(pca.res, habillage = pheno$response, mean.point = F, geom = 'point')
fviz_pca_ind(pca.res, habillage = pheno$treatment, mean.point = F, geom = 'point')
fviz_pca_ind(pca.res, habillage = pheno$timepoint, mean.point = F, geom = 'point')


#get lesional samples
pheno_lesional <- pheno[pheno$tissue == 'Lesional', ]
expr_lesional <- expr[, colnames(expr) %in% pheno_lesional$geo_accession]
pca.res.lesional <- prcomp(t(expr_lesional))

fviz_pca_ind(pca.res.lesional, habillage = pheno_lesional$tissue, mean.point = F, geom = 'point')
fviz_pca_ind(pca.res.lesional, habillage = pheno_lesional$timepoint, mean.point = F, geom = 'point')
fviz_pca_ind(pca.res.lesional, habillage = pheno_lesional$treatment, mean.point = F, geom = 'point')


#separate by drug in lesional
pheno_lesional_etanercept <- pheno_lesional[grepl("etaner",pheno_lesional$treatment, ignore.case = T), ]
expr_lesional_etanercept <- expr[, colnames(expr) %in% pheno_lesional_etanercept$geo_accession]
pca.res.lesional_etanercept <- prcomp(t(expr_lesional_etanercept))

fviz_pca_ind(pca.res.lesional_etanercept, habillage = pheno_lesional_etanercept$timepoint, mean.point = F, geom = 'point')
fviz_pca_ind(pca.res.lesional_etanercept, habillage = pheno_lesional_etanercept$response, mean.point = F, geom = 'point')



pheno_lesional_ustekinumab <- pheno_lesional[grepl("uste",pheno_lesional$treatment, ignore.case = T), ]
expr_lesional_ustekinumab <- expr[, colnames(expr) %in% pheno_lesional_ustekinumab$geo_accession]
pca.res.lesional_ustekinumab <- prcomp(t(expr_lesional_ustekinumab))

fviz_pca_ind(pca.res.lesional_ustekinumab, habillage = pheno_lesional_ustekinumab$timepoint, mean.point = F, geom = 'point')
fviz_pca_ind(pca.res.lesional_ustekinumab, habillage = pheno_lesional_ustekinumab$response, mean.point = F, geom = 'point')



#check differences in drugs for inflamm module

list_plot <- list()
c <- 0
for (module in list_modules) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_lesional, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  mod_name <- module_names[[c]]
  
  merged <- merge(score_df, pheno_lesional, by = 'row.names')
  #merged$disease <- factor(merged$disease, levels = c('Control', "CD", "UC"))
  #merged$treatment <- factor(merged$treatment, levels = c('control', "before", "after"))
  
  timepoint_comparisons <- list(
    c("Baseline", "Week 1"),
    c("Baseline", "Week 12"),
    c("Week 1", "Week 12")
  )
  
  
   p <-   ggplot(merged, aes(x = timepoint, y = module, fill = response)) +
   geom_boxplot() +
   stat_compare_means(
     #comparisons = timepoint_comparisons, # Pairwise comparisons
     #method = "anova",                  # Statistical test
     label = "p.format",                 # Display significance stars
     aes(group = response)              # Group by timepoint for comparison
   ) +
   facet_wrap(~treatment) + 
   theme_classic() +
   ggtitle('Inflammation module - Lesional skin pso')
   
   
   
  list_plot[[c]] <- p
  
  # pdf(file = paste0('module_', c, '.pdf'), width = 15, height = 10)
  # print(p)
  # dev.off()
}


pheno_corr <- pheno[!pheno$pasi == 'baseline', ]
expr_corr <- expr[, colnames(expr) %in% pheno_corr$geo_accession]


module <- list_modules[[13]]
total_df <- compute_features(expr = expr_corr, module = module, deconvolution = T)
total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
total_df <- na.omit(total_df)
score_df <- as.data.frame(colSums(total_df))
mod_name <- 'Inflammation_module'
colnames(score_df) <- mod_name
merged <- merge(score_df, pheno_corr, by = 'row.names')
merged$response <- factor(merged$response, levels = c("NR", "R"))
merged$pasi <- as.numeric(merged$pasi)

cor(as.numeric(merged$Inflammation_module), as.numeric(merged$pasi), method = 'pearson')

pl1 <- ggplot(merged, aes(y = pasi, x = timepoint, fill = response)) +
  geom_boxplot()+
  stat_compare_means(method = 't.test') +
  theme_classic() +
  #scale_fill_tableau(palette = 'Superfishel Stone')
  scale_fill_tableau(palette = 'Nuriel Stone', name = 'Response') +
  ylab("PASI")+
  xlab("") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))+
  ggtitle("PASI score at baseline vs week 12")

pl1

# jpeg(file = "PASI_score_at_baseline_vs_week_12.jpeg", units = 'cm', res = 300, width = 25, height = 20)
# pl1
# dev.off()

pl2 <- ggplot(merged, aes(y = pasi, x = timepoint, fill = response)) +
  geom_boxplot()+
  stat_compare_means(method = 't.test', size = 4.5) +
  theme_classic() +
  #scale_fill_tableau(palette = 'Superfishel Stone')
  scale_fill_tableau(palette = 'Nuriel Stone', name = 'Response') +
  ylab("PASI")+
  xlab("") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        strip.text.x = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+
 # ggtitle("PASI score at baseline vs week 12")+
  facet_wrap(~treatment)


pl2



stat.test <- merged %>%
  group_by(treatment, timepoint) %>%
  t_test(pasi ~ response) %>%
  mutate(label = formatC(p, format = "e", digits = 2)) %>%  # scientific, 2 decimals
  add_xy_position(x = "timepoint", dodge = 0.8)


# Plot
pl2.2 <- ggplot(merged, aes(x = timepoint, y = pasi, fill = response)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  stat_pvalue_manual(
    stat.test,
    label = "label",
    hide.ns = FALSE,               # <- KEEP non-significant brackets
    bracket.size = 0.4,
    tip.length = 0.01,
    inherit.aes = FALSE
  ) +
  theme_classic() +
  scale_fill_tableau(palette = "Nuriel Stone", name = "Response") +
  ylab("PASI") + xlab("") +
  theme(
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 12),
    axis.text.x  = element_text(size = 14),
    strip.text.x = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  facet_wrap(~ treatment)

pl2.2


jpeg(file = "PASI_score_at_baseline_vs_week_12_facet_treatment.jpeg", units = 'cm', res = 600, width = 30, height = 15)
pl2.2
dev.off()


#significant reductions in PASI scores were observed in responders treated with ustekinumab (both 45 mg and 90 mg) and etanercept  
#pl3
new_merged <- merged
new_merged <- new_merged[new_merged$response == 'R', ]

pl3 <- ggplot(new_merged, aes(y = pasi, x = timepoint, fill = response)) +
  geom_boxplot()+
  stat_compare_means(method = 't.test', size = 4.5) +
  theme_classic() +
  #scale_fill_tableau(palette = 'Superfishel Stone')
  scale_fill_tableau(palette = 'Nuriel Stone', name = 'Response') +
  ylab("PASI")+
  xlab("") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        strip.text.x = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+
  # ggtitle("PASI score at baseline vs week 12")+
  facet_wrap(~treatment)

pl3


new_merged <- merged
new_merged2 <- new_merged[new_merged$response == 'NR', ]

pl4 <- ggplot(new_merged2, aes(y = pasi, x = timepoint, fill = response)) +
  geom_boxplot()+
  stat_compare_means(method = 't.test', size = 4.5) +
  theme_classic() +
  #scale_fill_tableau(palette = 'Superfishel Stone')
  scale_fill_tableau(palette = 'Nuriel Stone', name = 'Response') +
  ylab("PASI")+
  xlab("") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        strip.text.x = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+
  # ggtitle("PASI score at baseline vs week 12")+
  facet_wrap(~treatment)

pl4


#correlation inflammation module
corr_plot <- ggplot(merged, aes(x = Inflammation_module, y = pasi)) +
  geom_point() +
  theme_classic()+
  geom_smooth(method = 'lm') +
  stat_cor(method = "spearman", label.x = -18, label.y = 58, size = 6, cor.coef.name = "ρ" )+
  facet_wrap(~treatment)+
  theme(
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 14),    # Increase axis text size
    legend.title = element_text(size = 14), # Increase legend title size
    legend.text = element_text(size = 14),  # Increase legend text size
    strip.text = element_text(size = 14, face = 'bold')    # Increase facet label size
  ) +
  ylab('PASI') +
  xlab('Inflammation module')

corr_plot
stop()
# jpeg(file = "corr_plot.jpeg", units = 'cm', res = 600, width = 30, height = 15)
# corr_plot
# dev.off()


#CORRELATIONS ALL MODULES

pheno_corr <- pheno[!pheno$pasi == 'baseline', ]
expr_corr <- expr[, colnames(expr) %in% pheno_corr$geo_accession]

z <- 0
list_plots2 <- list()
for (module in list_modules) {
  z <- z+1
  total_df <- compute_features(expr = expr_corr, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  mod_name <- 'module'
  colnames(score_df) <- mod_name
  merged <- merge(score_df, pheno_corr, by = 'row.names')
  merged$response <- factor(merged$response, levels = c("NR", "R"))
  merged$pasi <- as.numeric(merged$pasi)
  
  x <- ggplot(merged, aes(x = module, y = pasi)) +
    geom_point() +
    theme_classic()+
    geom_smooth(method = 'lm') +
    stat_cor(method = "spearman", size = 5)+
    facet_wrap(~treatment)+
    theme(
      axis.title = element_text(size = 14),   # Increase axis title size
      axis.text = element_text(size = 14),    # Increase axis text size
      legend.title = element_text(size = 14), # Increase legend title size
      legend.text = element_text(size = 14),  # Increase legend text size
      strip.text = element_text(size = 14, face = 'bold')    # Increase facet label size
    ) +
    ggtitle(module_names[[z]]) +
    ylab("PASI")+
    xlab("Module score")
  
  list_plots2[[z]] <- x
 
}


do.call(grid.arrange,c(list_plots2, ncol = 3))







stop()


#select for week12 timepoint lesional

pheno_week12 <- pheno_lesional[pheno_lesional$timepoint == 'Week 12', ]
expr_week12_lesional <- expr[, colnames(expr) %in% pheno_week12$geo_accession]



list_plot <- list()
c <- 0
for (module in list_modules[13]) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_week12_lesional, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  mod_name <- module_names[[c]]
  
  merged <- merge(score_df, pheno_week12, by = 'row.names')
  
  p <- ggplot(merged, aes(x = treatment, y = module, fill = response))+
    geom_boxplot() +
    theme_classic()+
    stat_compare_means(aes(group = response)
                      )+
    ggtitle('Inflammation module - Lesional skin pso week 12 - responder vs non responder')
  
  p
  
  list_plot[[c]] <- p

  
}





#select for baseline timepoint lesional

pheno_baseline <- pheno_lesional[pheno_lesional$timepoint == 'Baseline', ]
expr_baseline_lesional <- expr[, colnames(expr) %in% pheno_baseline$geo_accession]

list_plot <- list()
c <- 0
for (module in list_modules) {
  c <- c+1
  
  total_df <- compute_features(expr = expr_baseline_lesional, module = module, deconvolution = T)
  total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable
  
  total_df <- na.omit(total_df)
  score_df <- as.data.frame(colSums(total_df))
  colnames(score_df) <- "module"
  mod_name <- module_names[[c]]
  
  merged <- merge(score_df, pheno_baseline, by = 'row.names')
  
  p <- ggplot(merged, aes(x = disease, y = module, fill = response))+
    geom_boxplot() +
    theme_classic()+
    stat_compare_means(
      aes(group = response),
      method = 't.test'
    ) +
    ggtitle(mod_name)
  
  p
  
  list_plot[[c]] <- p
  
  
}

do.call(grid.arrange,c(list_plot, ncol = 5))

