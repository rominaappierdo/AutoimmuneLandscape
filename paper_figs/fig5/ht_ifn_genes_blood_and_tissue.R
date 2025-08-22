#setwd("paper_figs/fig5")
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
library(ggplot2)
library(viridis)
library(ComplexHeatmap)

ifn_genes <- readRDS("/mnt/data3/appierdo/PhD-project/gene_sets/hallmarks_list_GENE_SYMBOL.RData")
ifn_genes <- ifn_genes[grepl("interf", names(ifn_genes), ignore.case = T)]

ifna <- ifn_genes[[1]]
ifng <- ifn_genes[[2]]

ifn_genes <- unique(unname(unlist(ifn_genes)))


blood_files <- list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/meta_analysis_gene_expression/", pattern = 'complete', full.names = T, recursive = T)

b <- 0
for (file in blood_files) {
  b <- b+1
  disease <- strsplit(file, "/")[[1]][10]
  df <- read.csv(file, sep ='\t')
  df <- df[df$gene %in% ifn_genes, c("gene", "effectSize")]
  colnames(df) <- c('gene', disease)
  if (b ==1 ) {
    merged <- df
  } else {
    merged <- merge(merged, df, by = 'gene', all =T)
  }
  
}
rownames(merged) <- merged$gene
merged$gene <- NULL

merged_blood <- merged

blood_df <- wide_to_long(merged, name_row_column = 'genes', name_column_column = 'diseases', name_value_column = 'effectSize')

p1 <- ggplot(blood_df, aes(y = diseases, x = genes)) +
   geom_point(aes(size = effectSize, color = effectSize)) +
   theme_classic()

tissue_files <- list.files("/mnt/data3/appierdo/AutoimmuneLandscape//Analysis/TISSUE/meta_analysis_gene_expression/", pattern = 'complete', full.names = T, recursive = T)

t <- 0
for (file in tissue_files) {
  t <- t+1
  disease <- strsplit(file, "/")[[1]][11]
  tissue <- strsplit(file, "/|_")[[1]][15]
 
  disease <- paste0(disease ,"_", tissue)
  
  df <- read.csv(file, sep ='\t')
  #df <- df[df$numStudies >=2, ]
  #df <- df[df$effectSizeFDR < 0.05, ]
  df <- df[df$gene %in% ifn_genes, c("gene", "effectSize")]
  colnames(df) <- c('gene', disease)
  if (t ==1 ) {
    merged <- df
  } else {
    merged <- merge(merged, df, by = 'gene', all =T)
  }
  
}
rownames(merged) <- merged$gene
merged$gene <- NULL

merged_tissue <- merged

merged <- cbind(merged_blood, merged_tissue)

split <- ifelse(colnames(merged) %in% colnames(merged_blood), "Blood", "Tissue")

merged[merged > 2]<- 2
merged[merged < -2]<- -2
merged <- as.data.frame(t(merged))


rownames(merged)[rownames(merged) == 'crohn'] <- "Crohn's disease"
rownames(merged)[rownames(merged) =='diabetes'] <- "Type I diabetes"
rownames(merged)[rownames(merged) =='ra'] <- "Rheumatoid arthritis"
rownames(merged)[rownames(merged) =='psoriasis'] <- "Psoriasis"
rownames(merged)[rownames(merged) =='sjogren'] <- "Sjogren's syndrome"
rownames(merged)[rownames(merged) =='uc'] <- "Ulcerative colitis"
rownames(merged)[rownames(merged) =='multiple_sclerosis'] <- "Multiple sclerosis"
rownames(merged)[rownames(merged) =='scleroderma'] <- "Scleroderma"
rownames(merged)[rownames(merged) =='jia'] <- "Juvenile idiopathic arthritis"
rownames(merged)[rownames(merged) =='lupus'] <- "Lupus"

rownames(merged)[rownames(merged) =='crohn_GUT'] <- "Crohn's disease (Gut)"
rownames(merged)[rownames(merged) =='ra_SYNOVIUM'] <- "Rheumatoid arthritis (Synovium)"
rownames(merged)[rownames(merged) =='psoriasis_SKIN'] <- "Psoriasis (Skin)"
rownames(merged)[rownames(merged) =='sjogren_SALIVARY.GLAND'] <- "Sjogren's syndrome (Salivary glands)"
rownames(merged)[rownames(merged) =='uc_GUT'] <- "Ulcerative colitis (Gut)"
rownames(merged)[rownames(merged) =='scleroderma_SKIN'] <- "Scleroderma (Skin)"
rownames(merged)[rownames(merged) =='lupus_SKIN'] <- "Lupus (Skin)"
rownames(merged)[rownames(merged) =='lupus_KIDNEY'] <- "Lupus (Kidney)"



ht <- Heatmap(as.matrix(merged),
        row_split = split,
        col = circlize::colorRamp2(c(-2, 0, 2), c('blue', "#FFFFFF", "red")),
        row_gap = unit(3.5, "mm"),
        na_col = 'black',
        heatmap_height = unit(0.6, "cm")*nrow(merged),
        heatmap_width = unit(0.15, "cm")*ncol(merged),
        rect_gp = gpar(col = "black", lwd = 0.1),
        show_column_names = F,cluster_row_slices = F,
        heatmap_legend_param = list(title = "Effect Size"),
        column_title = "genes",
        column_title_side = "bottom"
        
)


mean(unlist(merged['Lupus', ]))
mean(unlist(merged['Lupus (Kidney)', ]))
mean(unlist(merged['Lupus (Skin)', ]))

# jpeg("fig5_ht.jpeg", width = 38,height = 18, units = "cm", res = 600)
# draw(ht)
# dev.off()
