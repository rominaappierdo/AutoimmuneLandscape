library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. LOAD RESULTS
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat"

output_dir <- file.path(
  base_dir,
  "module_clinical_associations",
  "sex",
  "step3_sex_associations"
)

results <- as.data.table(readRDS(
  file.path(output_dir, "sex_module_associations_by_GSE_disease.rds")
))

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
  
  dcast(
    results,
    module ~ group,
    value.var = value_col
  ) |>
    as.data.frame() |>
    (\(x) {
      rownames(x) <- x$module
      as.matrix(x[, -1, drop = FALSE])
    })()
}

beta_mat <- make_matrix("beta_F_vs_M")
p_mat    <- make_matrix("pvalue")
fdr_mat  <- make_matrix("FDR")

# Biological module labels
rownames(beta_mat) <- module_names[rownames(beta_mat)]
rownames(p_mat)    <- module_names[rownames(p_mat)]
rownames(fdr_mat)  <- module_names[rownames(fdr_mat)]

# ============================================================
# 4. HIERARCHICAL CLUSTERING
# ============================================================

# Use the same clustering for both heatmaps
# Missing values, if any, are replaced only for clustering

cluster_mat <- beta_mat
cluster_mat[!is.finite(cluster_mat)] <- 0

row_clust <- hclust(dist(cluster_mat), method = "complete")
col_clust <- hclust(dist(t(cluster_mat)), method = "complete")

# ============================================================
# 5. COLORS
# ============================================================

lim <- max(abs(beta_mat), na.rm = TRUE)

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
    name = "Female - Male",
    col = col_fun,
    
    cluster_rows = row_clust,
    cluster_columns = col_clust,
    
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    
    column_names_rot = 60,
    column_names_gp = gpar(fontsize = 7),
    
    column_title = title,
    
    rect_gp = gpar(col = "white", lwd = 0.3),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (is.finite(pvalues[i, j]) &&
          pvalues[i, j] < 0.05) {
        
        grid.text(
          "*", x, y,
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      }
    },
    
    heatmap_legend_param = list(
      title = "Female - Male"
    )
  )
}

# ============================================================
# 7. FDR HEATMAP
# ============================================================

ht_fdr <- make_heatmap(
  fdr_mat,
  "Sex-associated module activity | FDR < 0.05"
)

pdf(
  file.path(output_dir, "sex_heatmap_clustered_FDR.pdf"),
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
  "Sex-associated module activity | Nominal p < 0.05"
)

pdf(
  file.path(output_dir, "sex_heatmap_clustered_nominal.pdf"),
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