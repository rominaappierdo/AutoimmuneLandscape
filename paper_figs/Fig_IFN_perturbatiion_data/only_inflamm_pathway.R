setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/Fig_IFN_perturbatiion_data")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(GSVA)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")
#source("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/module_validation_figs/psoriasis_GSE117239/Untitled.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/function_compute_all_features_for_module.R")

gse<-'GSE78193'

hallmark <- readRDS("/mnt/data3/appierdo/AutoimmuneLandscape/gene_sets/hallmarks_list_GENE_SYMBOL.RData")
names(hallmark)[grepl("inflamm", names(hallmark), ignore.case = T)]

infl <- hallmark['HALLMARK_INFLAMMATORY_RESPONSE']





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

hall_df <- t(scale(t(gsva(
  as.matrix(expr), 
  gset.idx.list = infl, 
  min.sz=5, 
  parallel.sz=12, 
  method = 'ssgsea', 
  ssgsea.norm = F, 
  verbose = F))))


hall_df <- as.data.frame(t(hall_df))
colnames(hall_df) <- 'ssGSEA score'



hall_df$timepoint <- factor(meta$timepoint, levels = c('t0', '24h', '48h'))
hall_df$treatment <- factor(meta$treatment, levels = c("unstimulated", "IFN-gamma"))

tableau_colors <- ggthemes::tableau_color_pal("Tableau 10")(10)  # Get the full palette of 10 colors
# Subset and reverse the colors (e.g., use the first 5 colors and reverse them)
custom_colors <- tableau_colors[6:10]

box <- ggplot(hall_df, aes(x = treatment, y = `ssGSEA score`, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~timepoint) +
  theme_classic() +
  stat_compare_means() +
  scale_fill_manual(values = custom_colors) +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(
          size = 14, color = "white", face = 'bold')
  ) +
  ggtitle('HALLMARK_INFLAMMATORY_RESPONSE')

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

stop()
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
               title_gp = grid::gpar(fontsize = 40, fontface = "bold"),  # title size
               labels_gp = grid::gpar(fontsize = 30) 
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

