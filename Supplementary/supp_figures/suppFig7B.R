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

#Question #2: Can modules capture differences in response at baseline to predict response to treatment?

#we have to split the analysis in colon vs ileum
#COLON
pheno_colon <- pheno[pheno$tissue == 'Colon', ]

#get only the BEFORE treatment ones
pheno_colon <- pheno_colon[pheno_colon$treatment %in% c("control", "before"), ]
expr_colon <- expr[, colnames(expr) %in% pheno_colon$geo_accession]


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
    scale_fill_manual(values = c('bisque', '#add98e', '#ff800e'), name = 'Response')+
    theme(legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          legend.position = 'none')
  
  list_plot[[j]] <- p
  
}


saveRDS(list_plot, 'list_plot_colon_before.rds')
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


jpeg(file = 'SuppFig7B.jpeg', units = 'cm', width = 60, height = 38, res = 600)
do.call(grid.arrange,c(list_plot_new, ncol = 4))
dev.off()

