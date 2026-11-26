# ============================================================
# AutoimmuneLandscape
# LESION ANALYSIS - STEP 3 PLOTS
# Affected vs Unaffected tissue module activity
# Clustered heatmaps: FDR and nominal p-values
# ============================================================

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. LOAD RESULTS
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat_tissue/"

output_dir <- file.path(
  base_dir,
  "module_clinical_associations",
  "lesion",
  "step3_lesion_associations"
)

results <- as.data.table(readRDS(
  file.path(
    output_dir,
    "lesion_module_associations_by_GSE_disease.rds"
  )
))

module_cols <- paste0("mod", 1:15)

stopifnot(all(c(
  "series_id",
  "disease",
  "module",
  "beta_Affected_vs_Unaffected",
  "pvalue",
  "FDR"
) %in% names(results)))

if (nrow(results) == 0) {
  stop("No results available for plotting.")
}

# ============================================================
# 2. MODULE NAMES
# ============================================================

module_names <- c(
  mod1  = "Th1 response",
  mod2  = "Monocyte/macrophage-NK crosstalk",
  mod3  = "T cell activation and recruitment",
  mod4  = "Anti-inflammatory response",
  mod5  = "IFN response",
  mod6  = "IL6 trans-signaling / Transendothelial migration",
  mod7  = "Lymphoid cell development",
  mod8  = "Endothelial dysfunction",
  mod9  = "TLR / NF-kB / Innate response",
  mod10 = "Lipid metabolism",
  mod11 = "Inflammation / Wound healing",
  mod12 = "B cell response",
  mod13 = "Impaired cell proliferation",
  mod14 = "Integrins",
  mod15 = "ECM remodeling"
)

results[, group := paste(series_id, disease, sep = " | ")]

# ============================================================
# 3. BUILD MATRICES
# ============================================================

make_matrix <- function(value_col) {
  
  mat <- dcast(
    results,
    module ~ group,
    value.var = value_col
  )
  
  mat <- as.data.frame(mat)
  
  rownames(mat) <- mat$module
  
  mat <- as.matrix(
    mat[, -1, drop = FALSE]
  )
  
  storage.mode(mat) <- "double"
  
  # Keep biological module order
  mat <- mat[
    module_cols[module_cols %in% rownames(mat)],
    ,
    drop = FALSE
  ]
  
  return(mat)
}

beta_mat <- make_matrix("beta_Affected_vs_Unaffected")
p_mat    <- make_matrix("pvalue")
fdr_mat  <- make_matrix("FDR")

# Ensure identical matrix dimensions and ordering

stopifnot(
  identical(rownames(beta_mat), rownames(p_mat)),
  identical(rownames(beta_mat), rownames(fdr_mat)),
  identical(colnames(beta_mat), colnames(p_mat)),
  identical(colnames(beta_mat), colnames(fdr_mat))
)

# Replace module IDs with biological names

rownames(beta_mat) <- module_names[rownames(beta_mat)]
rownames(p_mat)    <- module_names[rownames(p_mat)]
rownames(fdr_mat)  <- module_names[rownames(fdr_mat)]

# ============================================================
# 4. HIERARCHICAL CLUSTERING
# ============================================================

# Use the same clustering for both heatmaps.
# Missing values are replaced with zero ONLY for clustering.

cluster_mat <- beta_mat

cluster_mat[!is.finite(cluster_mat)] <- 0

row_clust <- if (nrow(cluster_mat) >= 2) {
  hclust(dist(cluster_mat), method = "complete")
} else {
  FALSE
}

col_clust <- if (ncol(cluster_mat) >= 2) {
  hclust(dist(t(cluster_mat)), method = "complete")
} else {
  FALSE
}

# ============================================================
# 5. COLORS
# ============================================================

# Blue = higher in unaffected
# White = no difference
# Red = higher in affected

lim <- max(abs(beta_mat), na.rm = TRUE)

if (!is.finite(lim) || lim == 0) {
  lim <- 1
}

col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B6FB6", "white", "#C94B4B")
)

# ============================================================
# 6. HEATMAP FUNCTION
# ============================================================

make_heatmap <- function(pvalues, title) {
  
  Heatmap(
    beta_mat,
    
    name = "Affected - Unaffected",
    
    col = col_fun,
    
    na_col = "#E5E5E5",
    
    cluster_rows = row_clust,
    cluster_columns = col_clust,
    
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    
    column_names_rot = 60,
    column_names_gp = gpar(fontsize = 7),
    
    column_title = title,
    
    rect_gp = gpar(
      col = "white",
      lwd = 0.3
    ),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      
      if (
        is.finite(pvalues[i, j]) &&
        pvalues[i, j] < 0.05
      ) {
        
        grid.text(
          "*",
          x,
          y,
          gp = gpar(
            fontsize = 10,
            fontface = "bold"
          )
        )
        
      }
    },
    
    heatmap_legend_param = list(
      title = "Affected - Unaffected",
      at = c(-lim, 0, lim),
      labels = round(c(-lim, 0, lim), 2)
    )
  )
}

# ============================================================
# 7. FDR HEATMAP
# ============================================================

ht_fdr <- make_heatmap(
  fdr_mat,
  "Lesion-associated module activity | FDR < 0.05"
)

pdf(
  file.path(
    output_dir,
    "lesion_heatmap_clustered_FDR.pdf"
  ),
  width = 20,
  height = 9
)

draw(ht_fdr)
dev.off()

# ============================================================
# 8. NOMINAL P-VALUE HEATMAP
# ============================================================

ht_nominal <- make_heatmap(
  p_mat,
  "Lesion-associated module activity | Nominal p < 0.05"
)

pdf(
  file.path(
    output_dir,
    "lesion_heatmap_clustered_nominal.pdf"
  ),
  width = 20,
  height = 9
)

draw(ht_nominal)
dev.off()

# ============================================================
# 9. DISPLAY
# ============================================================

draw(ht_fdr)
draw(ht_nominal)

cat("\n====================================\n")
cat("LESION ANALYSIS - STEP 3 PLOTS COMPLETED\n")
cat("====================================\n")

cat("Modules:", nrow(beta_mat), "\n")
cat("GSE x disease groups:", ncol(beta_mat), "\n")

cat("\nPlots saved in:\n", output_dir, "\n")
cat("====================================\n")