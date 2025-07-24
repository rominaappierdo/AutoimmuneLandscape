setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/visualizeLandscape")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
library(MetaIntegrator)
library(ComplexHeatmap)

immunoStates <- immunoStatesMatrix
redundant_cell_types <- readLines("../../gene_sets/redundant_cell_types.txt")
cells_mat <- data.frame(matrix(nrow = length(colnames(immunoStates)), ncol = 0, dimnames = list(colnames(immunoStates))))
source.dir <- "../TISSUE/meta_analysis_cell_proportions/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  file_name <- strsplit(file, '/')[[1]][6]
  disease <- strsplit(file_name, '_')[[1]][2]
  tissue <- strsplit(file_name, '_')[[1]][1]
  disease <- paste0(disease, "_", tissue, sep = '')
  cell_prop <- read.csv(file, sep = '\t')
  cell_prop <- cell_prop[!cell_prop$gene %in% redundant_cell_types, ]
  cell_prop <- cell_prop[, c("gene", "effectSize")]
  cell_prop <- cbind(cell_prop, 'disease' = rep(disease,nrow(cell_prop)))
  cell_prop <- reshape(cell_prop, idvar = "gene", timevar = "disease", direction = "wide")
  rownames(cell_prop) <- cell_prop$gene
  cell_prop$gene <- NULL
  colnames(cell_prop) <- disease
  cells_mat <- merge(cells_mat, cell_prop, by = "row.names", all = T)
  rownames(cells_mat) <- cells_mat$Row.names
  cells_mat$Row.names <- NULL
}

#cell_prop <- na.omit(cell_prop)
cells_mat <- cells_mat[rownames(cells_mat) %in% readLines("script_to_select_significant_features_tissue/significant_cell_in_at_least_one_disease.txt"), ]


hallmarks_list <- names(readRDS("../../gene_sets/hallmarks_list_GENE_SYMBOL.RData"))
hall_mat <- data.frame(matrix(nrow = length(hallmarks_list), ncol = 0, dimnames = list(hallmarks_list)))
source.dir <- "../TISSUE/meta_analysis_pathways_activity/hallmarks/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  file_name <- strsplit(file, '/')[[1]][7]
  disease <- strsplit(file_name, '_')[[1]][3]
  tissue <- strsplit(file_name, '_')[[1]][1]
  disease <- paste0(disease, "_", tissue, sep = '')
  hall <- read.csv(file, sep = '\t')
  hall <- hall[, c("gene", "effectSize")]
  hall <- cbind(hall, 'disease' = rep(disease,nrow(hall)))
  hall <- reshape(hall, idvar = "gene", timevar = "disease", direction = "wide")
  rownames(hall) <- hall$gene
  hall$gene <- NULL
  colnames(hall) <- disease
  hall_mat <- merge(hall_mat, hall, by = "row.names", all = T)
  rownames(hall_mat) <- hall_mat$Row.names
  hall_mat$Row.names <- NULL
}

hall_mat <- na.omit(hall_mat)
hall_mat <- hall_mat[rownames(hall_mat) %in% readLines("script_to_select_significant_features_tissue/significant_hall_in_at_least_one_disease.txt"), ]


kegg_list <- names(readRDS("../../gene_sets/KEGG_list_GENE_SYMBOL.RData"))
kegg_mat <- data.frame(matrix(nrow = length(kegg_list), ncol = 0, dimnames = list(kegg_list)))
source.dir <- "../TISSUE//meta_analysis_pathways_activity/kegg/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  file_name <- strsplit(file, '/')[[1]][8]
  disease <- strsplit(file_name, '_')[[1]][3]
  tissue <- strsplit(file_name, '_')[[1]][1]
  disease <- paste0(disease, "_", tissue, sep = '')
  kegg <- read.csv(file, sep = '\t')
  kegg <- kegg[, c("gene", "effectSize")]
  kegg <- cbind(kegg, 'disease' = rep(disease,nrow(kegg)))
  kegg <- reshape(kegg, idvar = "gene", timevar = "disease", direction = "wide")
  rownames(kegg) <- kegg$gene
  kegg$gene <- NULL
  colnames(kegg) <- disease
  kegg_mat <- merge(kegg_mat, kegg, by = "row.names", all = T)
  rownames(kegg_mat) <- kegg_mat$Row.names
  kegg_mat$Row.names <- NULL
}
kegg_mat <- na.omit(kegg_mat)
kegg_mat <- kegg_mat[rownames(kegg_mat) %in% readLines("script_to_select_significant_features_tissue/significant_kegg_in_at_least_one_disease.txt"), ]
#rownames(kegg_mat) <- unlist(lapply(rownames(kegg_mat), function(x) gsub("KEGG_", "", x)))



tf_list <- names(readRDS("../../gene_sets/collecTRI_regulons_dorothea_08122023.rds"))
tf_mat <- data.frame(matrix(nrow = length(tf_list), ncol = 0, dimnames = list(tf_list)))
source.dir <- "../TISSUE/TFs/meta_analysis_tf_activity/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  file_name <- strsplit(file, '/')[[1]][7]
  disease <- strsplit(file_name, '_')[[1]][4]
  tissue <- strsplit(file_name, '_')[[1]][1]
  disease <- paste0(disease, "_", tissue, sep = '')
  tfs <- read.csv(file, sep = '\t')
  tfs <- tfs[, c("gene", "effectSize")]
  tfs <- cbind(tfs, 'disease' = rep(disease,nrow(tfs)))
  tfs <- reshape(tfs, idvar = "gene", timevar = "disease", direction = "wide")
  rownames(tfs) <- tfs$gene
  tfs$gene <- NULL
  colnames(tfs) <- disease
  tf_mat <- merge(tf_mat, tfs, by = "row.names", all = T)
  rownames(tf_mat) <- tf_mat$Row.names
  tf_mat$Row.names <- NULL
}
tf_mat <- tf_mat[rownames(tf_mat) %in% readLines("script_to_select_significant_features_tissue/significant_tf_in_at_least_one_disease.txt"), ]
rownames(tf_mat) <- unlist(lapply(rownames(tf_mat), function(x) paste0(x, '_tf')))
tf_mat <- na.omit(tf_mat)
colnames(tf_mat)[colnames(tf_mat) == 'gsva_SYNOVIUM'] <- 'ra_SYNOVIUM'

mirna_list <- names(readRDS("../../gene_sets/miRTarBase_regulons.rds"))
mirna_mat <- data.frame(matrix(nrow = length(mirna_list), ncol = 0, dimnames = list(mirna_list)))
source.dir <- "../TISSUE/miRNA/meta_analysis_mirna_activity/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  print(file)
  file_name <- strsplit(file, '/')[[1]][7]
  disease <- strsplit(file_name, '_')[[1]][4]
  tissue <- strsplit(file_name, '_')[[1]][1]
  disease <- paste0(disease, "_", tissue, sep = '')
  mirnas <- read.csv(file, sep = '\t')
  mirnas <- mirnas[, c("gene", "effectSize")]
  mirnas <- cbind(mirnas, 'disease' = rep(disease,nrow(mirnas)))
  mirnas <- reshape(mirnas, idvar = "gene", timevar = "disease", direction = "wide")
  rownames(mirnas) <- mirnas$gene
  mirnas$gene <- NULL
  colnames(mirnas) <- disease
  mirna_mat <- merge(mirna_mat, mirnas, by = "row.names", all = T)
  rownames(mirna_mat) <- mirna_mat$Row.names
  mirna_mat$Row.names <- NULL
}

mirna_mat <- mirna_mat[rownames(mirna_mat) %in% readLines("script_to_select_significant_features_tissue/significant_mirna_in_at_least_one_disease.txt"), ]
mirna_mat <- na.omit(mirna_mat)
colnames(mirna_mat)[colnames(mirna_mat) == 'vs_SYNOVIUM'] <- 'ra_SYNOVIUM'

cyto_list <- readLines("../../signaling_molecules/curated_kines/all_kines.txt")
cyto_mat <- data.frame(matrix(nrow = length(cyto_list), ncol = 0, dimnames = list(cyto_list)))

source.dir <- "../TISSUE/meta_analysis_gene_expression/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  file_name <- strsplit(file, '/')[[1]][6]
  disease <- strsplit(file_name, '_')[[1]][2]
  tissue <- strsplit(file_name, '_')[[1]][1]
  disease <- paste0(disease, "_", tissue, sep = '')
  cyto <- read.csv(file, sep = '\t')
  cyto <- cyto[cyto$gene %in% cyto_list,]
  cyto <- cyto[, c("gene", "effectSize")]
  cyto <- cbind(cyto, 'disease' = rep(disease,nrow(cyto)))
  cyto <- reshape(cyto, idvar = "gene", timevar = "disease", direction = "wide")
  rownames(cyto) <- cyto$gene
  cyto$gene <- NULL
  colnames(cyto) <- disease
  cyto_mat <- merge(cyto_mat, cyto, by = "row.names", all = T)
  rownames(cyto_mat) <- cyto_mat$Row.names
  cyto_mat$Row.names <- NULL
}

cyto_mat <- na.omit(cyto_mat)
cyto_mat <- cyto_mat[rownames(cyto_mat) %in% readLines("script_to_select_significant_features_tissue/significant_cyto_in_at_least_one_disease.txt"), ]
rownames(cyto_mat) <- unlist(lapply(rownames(cyto_mat), function(x) paste0(x, '_cyto')))


common_columns <- Reduce(intersect, list(colnames(cells_mat),colnames(hall_mat),colnames(kegg_mat), colnames(cyto_mat), colnames(mirna_mat), colnames(tf_mat)))

cells_mat <- cells_mat[, colnames(cells_mat) %in% common_columns]
hall_mat <- hall_mat[, colnames(hall_mat) %in% common_columns]
kegg_mat <- kegg_mat[, colnames(kegg_mat) %in% common_columns]
cyto_mat <- cyto_mat[, colnames(cyto_mat) %in% common_columns]
mirna_mat <- mirna_mat[, colnames(mirna_mat) %in% common_columns]
tf_mat <- tf_mat[, colnames(tf_mat) %in% common_columns]

colnames(cells_mat)
colnames(kegg_mat)
colnames(cyto_mat)
colnames(hall_mat)
colnames(tf_mat)
colnames(mirna_mat)

df <- do.call("rbind", list((cells_mat), (hall_mat), (kegg_mat), (cyto_mat), (mirna_mat), (tf_mat)))
column_split <- c(rep('Cell\nproportions', nrow(cells_mat)), rep('Pathways activity', nrow(hall_mat) + nrow(kegg_mat)), rep('Cytokines\nexpression', nrow(cyto_mat)), rep('miRNAs activity', nrow(mirna_mat)), rep('TFs activity', nrow(tf_mat)))

df <- t((df))
rownames(df) <- c("Crohn's disease (Gut)", 'Lupus (Kidney)', 'Lupus (Skin)', 'Psoriasis (Skin)', "Rheumatoid arthritis (Synovium)",  'Scleroderma (Skin)', "Sjogren's syndrome (Salivary glands)", 'Ulcerative colitis (Gut)')


ht <- Heatmap(
  as.matrix(df),
  column_split = column_split,
  name = 'pooled effect size',
  column_names_gp = gpar(fontsize = 2),
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  #bottom_annotation = column_ha,
  #column_gap = unit(2, "mm"),
  column_title_gp = gpar(col = c("darkgreen","darkviolet", "orange", "steelblue", "red"), font = rep(2, 5)),
  heatmap_height = unit(10, "mm")*nrow(df),
  heatmap_width = unit(0.4, "mm")*ncol(df),
  na_col = 'black',
  show_column_names = F,
  heatmap_legend_param = list(
    legend_direction = "horizontal")
)

ht <- draw(ht, heatmap_legend_side = "bottom")



jpeg("ht_lanadscape_tissue.jpeg", units="in", width=18, height=5, res=600)
draw(ht)
dev.off()


