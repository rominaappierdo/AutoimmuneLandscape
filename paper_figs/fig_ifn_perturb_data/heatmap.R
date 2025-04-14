setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig_ifn_perturb_data")
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
module <- readLines("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/mod5/mod5.txt")
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


palette_for_anno <- list(
  Timepoint = c('t0' = '#377eb8', '24h' = '#e41a1b', '48h' = '#4eaf4a'),
  Stimulation = c("IFN-gamma" = "#B07AA1", "unstimulated" = "#EDC948")
)

anno <- HeatmapAnnotation(
  Timepoint = meta$timepoint,
  Stimulation = meta$treatment,
  col = palette_for_anno
)

rownames(total_df) <- gsub("KEGG_|HALLMARK_|WP_", "", rownames(total_df))

p <- Heatmap(as.matrix(total_df),
             bottom_annotation = anno,
             column_names_rot = 45,
             column_split = factor(meta$timepoint, levels = c('t0', '24h', '48h')) ,
             column_names_gp = grid::gpar(fontsize = 8),
             row_names_gp = grid::gpar(fontsize = 11),
             heatmap_height = unit(0.45, "cm")*nrow(total_df),
             heatmap_width = unit(0.7, "cm")*ncol(total_df),
             heatmap_legend_param = list(
               title = "Zscores",
               labels_gp = grid::gpar(fontsize = 40)  # Change legend labels size
             ),
             cluster_column_slices = F,
             #col = circlize::colorRamp2(c(-1, 0, 1), c("orange1", "#FFFFFF", "purple")),
             show_heatmap_legend = F,
             show_column_names = F ,
             column_title = paste0(gse, " | ", module_interpretation),
             rect_gp = gpar(col = "black", lwd = 0.3),
             column_gap = unit(3, "mm"),
             column_title_gp = gpar(fontsize = 18, fontface = "bold")
             
             
             
             
)
p<-draw(p, annotation_legend_side = "left", heatmap_legend_side = "bottom")



png(filename = 'ht.png', res = 600, width = 10000, height = 7000)
p
dev.off()
