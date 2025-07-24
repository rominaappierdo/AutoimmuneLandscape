setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/visualizeLandscape")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
library(MetaIntegrator)
library(ComplexHeatmap)


nDiseases <- 10
cols <- c("crohn", "diabetes", "jia", "lupus", "multiple_sclerosis", "psoriasis", "ra", "scleroderma", "sjogren", "uc")
immunoStates <- immunoStatesMatrix
redundant_cell_types <- readLines("../../gene_sets/redundant_cell_types.txt")
cells_mat <- data.frame(matrix(nrow = length(colnames(immunoStates)), ncol = 0, dimnames = list(colnames(immunoStates))))

source.dir <- "../BLOOD/meta_analysis_cell_proportions/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  disease <- strsplit(file, '/')[[1]][5]
  
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
cells_mat <- na.omit(cells_mat)
cells_mat <- cells_mat[rownames(cells_mat) %in% readLines("script_to_select_significant_features_blood/significant_cell_in_at_least_one_disease.txt"), ]



hallmarks_list <- names(readRDS("../../gene_sets/hallmarks_list_GENE_SYMBOL.RData"))
hall_mat <- data.frame(matrix(nrow = length(hallmarks_list), ncol = 0, dimnames = list(hallmarks_list)))
source.dir <- "../BLOOD/meta_analysis_pathways_activity/hallmarks/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  disease <- strsplit(file, '/')[[1]][6]
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
hall_mat <- hall_mat[rownames(hall_mat) %in% readLines("script_to_select_significant_features_blood/significant_hall_in_at_least_one_disease.txt"), ]
hall_mat <- hall_mat[rownames(hall_mat) %in% readLines("../../gene_sets/hallmark_immune_related.txt"), ]


kegg_list <- names(readRDS("../../gene_sets/KEGG_list_GENE_SYMBOL.RData"))
kegg_mat <- data.frame(matrix(nrow = length(kegg_list), ncol = 0, dimnames = list(kegg_list)))
source.dir <- "../BLOOD/meta_analysis_pathways_activity/kegg/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  disease <- strsplit(file, '/')[[1]][6]
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
kegg_mat <- kegg_mat[rownames(kegg_mat) %in% readLines("script_to_select_significant_features_blood/significant_kegg_in_at_least_one_disease.txt"), ]
#rownames(kegg_mat) <- unlist(lapply(rownames(kegg_mat), function(x) gsub("KEGG_", "", x)))
kegg_mat <- na.omit(kegg_mat)
kegg_mat <- kegg_mat[rownames(kegg_mat) %in% readLines("../../gene_sets/kegg_immune_related.txt"), ]



wiki_list <- names(readRDS("../../gene_sets/WIKIPATHWAYS_list_GENE_SYMBOL.RData"))
wiki_mat <- data.frame(matrix(nrow = length(wiki_list), ncol = 0, dimnames = list(wiki_list)))
source.dir <- "../BLOOD/meta_analysis_pathways_activity/wikipathways/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  disease <- strsplit(file, '/')[[1]][6]
  wiki <- read.csv(file, sep = '\t')
  wiki <- wiki[, c("gene", "effectSize")]
  wiki <- cbind(wiki, 'disease' = rep(disease,nrow(wiki)))
  wiki <- reshape(wiki, idvar = "gene", timevar = "disease", direction = "wide")
  rownames(wiki) <- wiki$gene
  wiki$gene <- NULL
  colnames(wiki) <- disease
  wiki_mat <- merge(wiki_mat, wiki, by = "row.names", all = T)
  rownames(wiki_mat) <- wiki_mat$Row.names
  wiki_mat$Row.names <- NULL
}
wiki_mat <- na.omit(wiki_mat)
wiki_mat <- wiki_mat[rownames(wiki_mat) %in% readLines("script_to_select_significant_features_blood/significant_wiki_in_at_least_one_disease.txt"), ]
wiki_mat <- wiki_mat[rownames(wiki_mat) %in% readLines("../../gene_sets/wiki_immune_related.txt"), ]


tf_list <- names(readRDS("../../gene_sets/collecTRI_regulons_dorothea_08122023.rds"))
tf_mat <- data.frame(matrix(nrow = length(tf_list), ncol = 0, dimnames = list(tf_list)))
source.dir <- "../BLOOD/TFs/meta_analysis_tf_activity/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  disease <- strsplit(file, '/')[[1]][6]
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
tf_mat <- tf_mat[rownames(tf_mat) %in% readLines("script_to_select_significant_features_blood/significant_tf_in_at_least_one_disease.txt"), ]
rownames(tf_mat) <- unlist(lapply(rownames(tf_mat), function(x) paste0(x, '_tf')))
tf_mat <- na.omit(tf_mat)


mirna_list <- names(readRDS("../../gene_sets/miRTarBase_regulons.rds"))
mirna_mat <- data.frame(matrix(nrow = length(mirna_list), ncol = 0, dimnames = list(mirna_list)))
source.dir <- "../BLOOD/miRNA/meta_analysis_mirna_activity/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  disease <- strsplit(file, '/')[[1]][6]
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
mirna_mat <- mirna_mat[rownames(mirna_mat) %in% readLines("script_to_select_significant_features_blood/significant_mirna_in_at_least_one_disease.txt"), ]
mirna_mat <- na.omit(mirna_mat)



cyto_list <- readLines("../../signaling_molecules/curated_kines/all_kines.txt")
cyto_mat <- data.frame(matrix(nrow = length(cyto_list), ncol = 0, dimnames = list(cyto_list)))
source.dir <- "../BLOOD/meta_analysis_gene_expression/"
for (file in list.files(source.dir, pattern = '.*complete_res.*', recursive = T, full.names = T)) {
  disease <- strsplit(file, '/')[[1]][5]
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
cyto_mat <- cyto_mat[rownames(cyto_mat) %in% readLines("script_to_select_significant_features_blood/significant_cyto_in_at_least_one_disease.txt"), ]
rownames(cyto_mat) <- unlist(lapply(rownames(cyto_mat), function(x) paste0(x, '_cyto')))
cyto_mat <- na.omit(cyto_mat)



df <- do.call("rbind", list((cells_mat), (hall_mat), (kegg_mat), (cyto_mat), (mirna_mat), (tf_mat), (wiki_mat)))


#column_split <- c(rep('cells', nrow(cells_mat)), rep('HALLMARK\npathways', nrow(hall_mat)), rep('KEGG\npathways', nrow(kegg_mat)), rep('cytokines', nrow(cyto_mat)), rep('miRNAs', nrow(mirna_mat)), rep('TFs', nrow(tf_mat)))

column_split <- c(rep('Cell\nproportions', nrow(cells_mat)), rep('Pathways activity', nrow(hall_mat)+nrow(kegg_mat)+nrow(wiki_mat)), rep('Cytokines\nexpression', nrow(cyto_mat)), rep('miRNAs\nactivity', nrow(mirna_mat)), rep('TFs activity', nrow(tf_mat)))


df <- rename_single_column(df, 'diabetes', 'type1_diabetes')
df <- t((df))

rownames(df) <- c("Crohn's disease", "Type 1 diabetes", "Juvenile idiopathic arthritis", "Lupus", "Multiple sclerosis", 
                  "Psoriasis", "Rheumatoid arthritis", "Scleroderma", "Sjogren's syndrome", 
                  "Ulcerative colitis")



ht <- Heatmap(
  as.matrix(df),
  column_split = column_split,
  name = 'pooled effect size',
  #cluster_rows = F, 
  #cluster_columns = F,
  column_names_gp = gpar(fontsize = 2),
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  #bottom_annotation = column_ha,
  #column_gap = unit(2, "mm"),
  heatmap_height = unit(6.4, "mm")*nrow(df),
  heatmap_width = unit(0.5, "mm")*ncol(df),
  na_col = 'black',
  show_column_names = F,
  rect_gp = gpar(col = "black", lwd = 0.03),
  show_heatmap_legend = T,
  column_title_gp = gpar(col = c("steelblue","orange", "darkviolet", "darkgreen",  "red"), font = rep(2, 5)),
  heatmap_legend_param = list(
    legend_direction = "horizontal")
)

ht <- draw(ht, heatmap_legend_side = "bottom")




jpeg("ht_lanadscape_blood.jpeg", units="in", width=15, height=5, res=600)
draw(ht)
dev.off()



