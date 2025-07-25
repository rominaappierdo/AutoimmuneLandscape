setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Supplementary/supp_figures")
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

#QUESTION #1: Can modules capture treatment effect? 

#we have to split the analysis in colon vs ileum
#COLON
pheno_colon <- pheno[pheno$tissue == 'Colon', ]

#get only the AFTER treatment ones
pheno_colon <- pheno_colon[pheno_colon$treatment %in% c("control", "after"), ]
expr_colon <- expr[, colnames(expr) %in% pheno_colon$geo_accession]

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
    xlab("") +
    theme(legend.text=element_text(size=18),
      legend.title=element_text(size=18),
      legend.position = 'none')
  
  list_plot[[j]] <- p
  
  #save legend
  # l <- p+ theme(legend.text=element_text(size=16), legend.title = element_text(size=16), legend.position="bottom")
  # jpeg(file = 'legend_colon_after_treatment.jpeg', units = 'cm', width = 20, height = 20, res = 600)
  # l
  # dev.off()
  
  
}




jpeg(file = 'boxplot_all_modules_colon.jpeg', units = 'cm', width = 57, height = 30, res = 200)
do.call(grid.arrange,c(list_plot, ncol = 5))
dev.off()

saveRDS(list_plot, 'list_plots_after_treat_colon.rds')


#significant modules at least in 1 comparison

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


jpeg(file = 'SuppFig7A.jpeg', units = 'cm', width = 68, height = 47, res = 600)
do.call(grid.arrange,c(list_plot_new, ncol = 4))
dev.off()

