
# ============================================================
# AutoimmuneLandscape - Module scores quality control
# ============================================================

library(data.table)
library(ggplot2)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat_tissue/")

scores_dir <- file.path(getwd(), "module_scores_final")
qc_dir <- file.path(getwd(), "module_scores_QC")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. LOAD COVERAGE QC
# ------------------------------------------------------------

qc <- fread(file.path(scores_dir, "module_score_coverage_QC.tsv"))

# Summary across datasets
coverage_summary <- qc[, .(
  n_datasets = .N,
  mean_coverage = mean(coverage, na.rm = TRUE),
  min_coverage = min(coverage, na.rm = TRUE),
  n_datasets_below_50 = sum(coverage < 0.5),
  n_datasets_below_80 = sum(coverage < 0.8),
  n_datasets_with_scores = sum(n_scored_GSM > 0)
), by = module]

coverage_summary[, module := factor(
  module, levels = paste0("mod", 1:15)
)]

print(coverage_summary)

fwrite(
  coverage_summary,
  file.path(qc_dir, "coverage_summary.tsv"),
  sep = "\t"
)

# ------------------------------------------------------------
# 2. COVERAGE HEATMAP
# ------------------------------------------------------------

qc[, module := factor(
  module, levels = paste0("mod", 1:15)
)]

p_coverage <- ggplot(
  qc,
  aes(x = module, y = dataset_file, fill = coverage)
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    limits = c(0, 1)
  ) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Module feature coverage across datasets",
    x = "Module",
    y = "Dataset",
    fill = "Coverage"
  ) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(qc_dir, "module_coverage_heatmap.pdf"),
  p_coverage,
  width = 10,
  height = 17
)

# ------------------------------------------------------------
# 3. LOAD ALL MODULE SCORES
# ------------------------------------------------------------

score_files <- list.files(
  scores_dir,
  pattern = "_module_scores\\.rds$",
  full.names = TRUE
)

stopifnot(length(score_files) > 0)

scores_list <- lapply(score_files, function(f) {
  
  x <- as.matrix(readRDS(f))
  
  dt <- as.data.table(x, keep.rownames = "GSM")
  
  dt[, dataset_file := sub(
    "_module_scores\\.rds$",
    "_features.rds",
    basename(f)
  )]
  
  dt
})

scores_all <- rbindlist(scores_list, fill = TRUE)

module_cols <- paste0("mod", 1:15)

# Save combined scores
fwrite(
  scores_all,
  file.path(qc_dir, "all_module_scores.tsv"),
  sep = "\t"
)

cat("Total rows:", nrow(scores_all), "\n")
cat("Unique GSM:", uniqueN(scores_all$GSM), "\n")
cat("Datasets:", uniqueN(scores_all$dataset_file), "\n")

# ------------------------------------------------------------
# 4. SCORE DISTRIBUTIONS
# ------------------------------------------------------------

scores_long <- melt(
  scores_all,
  id.vars = c("GSM", "dataset_file"),
  measure.vars = module_cols,
  variable.name = "module",
  value.name = "score"
)

scores_long[, module := factor(
  module, levels = module_cols
)]

# Summary statistics
score_summary <- scores_long[, .(
  n_valid = sum(is.finite(score)),
  mean = mean(score, na.rm = TRUE),
  sd = sd(score, na.rm = TRUE),
  median = median(score, na.rm = TRUE)
), by = .(dataset_file, module)]

fwrite(
  score_summary,
  file.path(qc_dir, "score_distribution_summary.tsv"),
  sep = "\t"
)

# Plot each module across datasets
for (mod in module_cols) {
  
  d <- scores_long[
    module == mod & is.finite(score)
  ]
  
  if (nrow(d) == 0) next
  
  p <- ggplot(
    d,
    aes(x = reorder(dataset_file, score, median),
        y = score)
  ) +
    geom_boxplot(outlier.size = 0.3) +
    coord_flip() +
    theme_classic(base_size = 10) +
    labs(
      title = paste("Score distribution:", mod),
      x = "Dataset",
      y = "Module score"
    ) +
    theme(axis.text.y = element_text(size = 5))
  
  ggsave(
    file.path(qc_dir, paste0(mod, "_distribution.pdf")),
    p,
    width = 9,
    height = 16
  )
}

# ------------------------------------------------------------
# 5. PCA EXPLORATORY ANALYSIS
# ------------------------------------------------------------

# Keep complete cases for the first exploratory PCA
pca_data <- scores_all[
  complete.cases(scores_all[, ..module_cols])
]

cat("Samples retained for PCA:", nrow(pca_data), "\n")

if (nrow(pca_data) >= 3) {
  
  x <- as.matrix(pca_data[, ..module_cols])
  
  # Scale modules globally for PCA
  nonconstant <- apply(x, 2, sd) > 0
  x <- x[, nonconstant, drop = FALSE]
  
  if (ncol(x) >= 2) {
    
    pca <- prcomp(x, center = TRUE, scale. = TRUE)
    
    pca_plot <- data.table(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      dataset_file = pca_data$dataset_file
    )
    
    variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    
    p <- ggplot(
      pca_plot,
      aes(x = PC1, y = PC2, color = dataset_file)
    ) +
      geom_point(alpha = 0.5, size = 0.8) +
      theme_classic() +
      guides(color = "none") +
      labs(
        title = "PCA of module scores",
        x = paste0("PC1 (", round(variance[1], 1), "%)"),
        y = paste0("PC2 (", round(variance[2], 1), "%)")
      )
    
    ggsave(
      file.path(qc_dir, "module_scores_PCA.pdf"),
      p,
      width = 8,
      height = 6
    )
    
    fwrite(
      pca_plot,
      file.path(qc_dir, "module_scores_PCA_coordinates.tsv"),
      sep = "\t"
    )
  }
}

cat("\nQC completed. Results saved in:", qc_dir, "\n")






library(ComplexHeatmap)
library(circlize)

module_cols <- paste0("mod", 1:15)

cor_mat <- cor(
  as.matrix(scores_all[, ..module_cols]),
  use = "pairwise.complete.obs",
  method = "pearson"
)

pdf(
  file.path(qc_dir, "module_correlation_heatmap.pdf"),
  width = 9,
  height = 8
)

Heatmap(
  cor_mat,
  name = "Pearson r",
  col = colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  ),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(
      sprintf("%.2f", cor_mat[i, j]),
      x, y,
      gp = grid::gpar(fontsize = 8)
    )
  }
)

dev.off()







library(data.table)
library(ggplot2)

module_pairs <- as.data.table(
  t(combn(names(modules), 2))
)

setnames(module_pairs, c("mod1", "mod2"))

module_pairs[, `:=`(
  n_shared = mapply(
    function(a, b) length(intersect(modules[[a]], modules[[b]])),
    mod1, mod2
  ),
  jaccard = mapply(
    function(a, b) {
      length(intersect(modules[[a]], modules[[b]])) /
        length(union(modules[[a]], modules[[b]]))
    },
    mod1, mod2
  )
)]

module_pairs[
  order(-jaccard),
  .(mod1, mod2, n_shared, jaccard)
][1:20]




# ============================================================
# 8. MODULE CORRELATION HEATMAP WITH BIOLOGICAL NAMES
# ============================================================

library(ComplexHeatmap)
library(circlize)

module_names <- c(
  mod1  = "Th1 response",
  mod2  = "Monocyte/macrophage - NK crosstalk",
  mod3  = "T cell activation and recruitment",
  mod4  = "Anti-inflammatory response",
  mod5  = "IFN response",
  mod6  = "IL6 trans-signaling / Transendothelial migration",
  mod7  = "Lymphoid cell development",
  mod8  = "Endothelial dysfunction",
  mod9  = "TLR / NF-kB / Innate immunity",
  mod10 = "Lipid metabolism",
  mod11 = "Inflammation / Wound healing",
  mod12 = "B cell response",
  mod13 = "Integrins",
  mod14 = "Impaired cell proliferation",
  mod15 = "ECM remodeling"
)

# Copy original correlation matrix
cor_mat_named <- cor_mat

# Replace module IDs with biological names
rownames(cor_mat_named) <- module_names[rownames(cor_mat)]
colnames(cor_mat_named) <- module_names[colnames(cor_mat)]

# Save heatmap
pdf(
  file.path(qc_dir, "module_correlation_heatmap_named.pdf"),
  width = 15,
  height = 11
)

Heatmap(
  cor_mat_named,
  name = "Pearson r",
  col = colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  ),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 10),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(
      sprintf("%.2f", cor_mat_named[i, j]),
      x, y,
      gp = grid::gpar(fontsize = 8)
    )
  }
)

dev.off()

