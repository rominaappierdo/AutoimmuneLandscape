setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/modules_applications/perturbation_experiment_ifn/")
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

gse<-'GSE78193'

########## set module
module_name <- "mod5"
module <- readLines("../../BLOOD/pipeline_landscape_disease_level_final/modules_files/mod5/mod5.txt")
module_interpretation <- 'IFN-g module'

########## get expression data
data <- getGEOData(gse)
expr <- data$originalData$GSE78193_GPL6848$expr
map_probes <- (as.data.frame(data$originalData$GSE78193_GPL6848$keys))
colnames(map_probes) <- 'keys'
expr <- merge(expr, map_probes, by = 'row.names')
expr$Row.names <- NULL
expr <- na.omit(expr)
expr <- mean_expression_gene_rowise(expr, colname = 'keys')
expr <- order_by_col_names(expr)


########## get metadata
meta <- fread("/mnt/data3/appierdo/PhD-project/Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv")
meta <- meta[meta$series_id==gse]
meta <- meta[meta$cell_type=='cultured sample', ]
meta <- sort_by_col(meta, colname = 'geo_accession')

if (!identical(meta$geo_accession, colnames(expr))) {
  stop("different samples from expression and metadata")
}

########## compute features 
total_df <- compute_features(expr = expr, module = module, deconvolution = T)
total_df <- t(scale(t(total_df))) #scale feature across sample, make samples comparable

rownames(total_df) <- gsub("KEGG_|HALLMARK_|WP_", "", rownames(total_df))

#compute module scores
score_df <- as.data.frame(colSums(total_df))
colnames(score_df) <- "module"
score_df$Condition <- paste0(meta$timepoint, "_", meta$treatment)

score_df$Condition <- factor(score_df$Condition, levels = c("t0_unstimulated", "t0_IFN-gamma",
                                                            "24h_unstimulated", "24h_IFN-gamma",
                                                            "48h_unstimulated", "48h_IFN-gamma"))

score_df$timepoint <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][1]}))
score_df$stimulation <- unlist(lapply(as.vector(score_df$Condition), function(x) {strsplit(x, '_')[[1]][2]}))
score_df$facet = factor(score_df$timepoint, levels = c("t0", "24h", "48h"))


tableau_colors <- ggthemes::tableau_color_pal("Tableau 10")(10)  # Get the full palette of 10 colors
# Subset and reverse the colors (e.g., use the first 5 colors and reverse them)
custom_colors <- tableau_colors[6:10]

my_comparisons <- list(c("IFN-gamma", "unstimulated"))

score_df$stimulation <- factor(score_df$stimulation, levels = c("unstimulated", "IFN-gamma"))
box <- ggplot(score_df, aes(x = stimulation, y = module)) +
  geom_boxplot(aes(fill = stimulation)) +
  stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test") +
  ggtitle("IFN module")+
  ylab("module score")+
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(
          size = 14, color = "white", face = 'bold')
  ) +
  #scale_fill_tableau(direction = -1) +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~ facet)

box



g <- ggplot_gtable(ggplot_build(box))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#377eb8","#e41a1b",  "#4eaf4a")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- fills[k]
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  
  k <- k+1
}
grid.draw(g)