setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Supplementary/supp_figures")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
options(repr.plot.width=15, repr.plot.height=15)
library(ComplexHeatmap)
library(grid)
library(dendextend)

redundant <- readLines("../../gene_sets/redundant_cell_types.txt")

significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0
for (file in list.files("../../Analysis/TISSUE/meta_analysis_cell_proportions/", recursive = T, pattern = '*complete*', full.names = T)) {
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t<-t[t["effectSizeFDR"]<0.05,]
  significant_genes<-append(significant_genes, t$gene)
}

tail(sort(table(significant_genes)), 100)


data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in list.files("../../Analysis/TISSUE/meta_analysis_cell_proportions/", recursive = T, pattern = '*complete*', full.names = T)) {
  disease<-strsplit(file, "\\/")[[1]][7]
  tissue<-strsplit(file, "\\/|_")[[1]][11]
  disease <- paste0(disease, "_", tissue, sep = '')
  df<-read.csv(file, sep='\t')
  df<-df[df$gene %in% significant_genes, c("effectSize", "gene")]
  
  df$disease<-rep(disease, nrow(df))
  data<-rbind(data, df)
}

table(data$gene)

#reshaping
data<-tidyr::spread(data, "disease", "effectSize")
rownames(data)<-data$gene
data$gene<-NULL

data <- data[!rownames(data) %in% redundant, ]
n_studies<- c(27, 5, 7, 23, 3, 2, 6, 30)


plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))
#row_dend = color_branches(row_dend, k = 2)

data[data > 2]<- 2
data[data < -2]<- -2


column_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F, add_numbers = T, numbers_rot = 0), which = 'column')

colnames(data) <- c("Crohn's disease (Gut)", "Lupus (Kidney)", "Lupus (Skin)", "Psoriasis (Skin)",
                    "Rheumatoid arthritis (Synovium)", "Scleroderma (Skin)", "Sjogren's syndrome (Salivary glands)",
                    "Ulcerative colitis (Gut)")

anno_lineage <- c("myeloid", "lymphoid", "lymphoid", "lymphoid", "myeloid", "myeloid", "myeloid", "myeloid", "myeloid",
                  "myeloid", "lymphoid", "myeloid", "lymphoid", "lymphoid")

palette_for_anno <- list(lineage = c('myeloid' = 'black', 'lymphoid' = 'lightgrey'))

anno <- rowAnnotation(lineage = anno_lineage, col = palette_for_anno, show_annotation_name = FALSE , simple_anno_size = unit(0.3, "cm"))

rownames(data) <- c("Basophils", "CD14+ monocytes", "CD4+ alpha beta T cells","CD56dim NK cells",
                    "CD8+ alpha beta T cells", "Eosinophils", "Hematopoietic progenitors", "Macrophages M0", "Macrophages M1", "Macrophages M2","MAST cells",
                    "Memory B cells","Myeloid dendritic cells", "Naive B cells", "Neutrophils", "Plasma cells", "Plasmacytoid dendritic cells")


ht<-Heatmap(as.matrix(data), 
            name = "pooled ES", 
            cluster_rows = row_dend, 
            cluster_columns = column_dend, 
            na_col = 'black',
            column_names_rot = 45,
            bottom_annotation = column_ha,
            col = circlize::colorRamp2(c(-1, 0, 1), c("#56B4E9", "#FFFFFF", "#FF7400")),
            show_row_names = T,
            use_raster = F,
            column_names_gp = gpar(fontsize = 11),
            row_names_gp = gpar(fontsize = 10),
            heatmap_width = unit(1.7, "cm")*ncol(data),
            heatmap_height = unit(1.3, "cm")*nrow(data),
            rect_gp = gpar(col = "black", lwd = 0.2),
            #right_annotation = anno
)


ht = draw(ht)

jpeg("SuppFig6.jpeg", width = 22, height = 23, units = 'cm', res = 600)
draw(ht)
dev.off()

