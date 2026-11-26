
# ============================================================
# AutoimmuneLandscape
# STEP 4 - EXPLORATORY PCA OF MODULE SCORES
# ============================================================

library(data.table)
library(ggplot2)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat")

input_dir <- file.path(getwd(), "endotype_analysis")
output_dir <- file.path(input_dir, "PCA")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)

# ============================================================
# 1. LOAD DATA
# ============================================================

data_all <- as.data.table(readRDS(
  file.path(input_dir, "all_samples_with_metadata.rds")
))

# Keep only samples with all 15 finite scores
complete <- rowSums(
  is.finite(as.matrix(data_all[, ..module_cols]))
) == length(module_cols)

data_pca <- copy(data_all[complete])

cat("Total samples:", nrow(data_all), "\n")
cat("Samples retained:", nrow(data_pca), "\n")
cat("Samples excluded:", sum(!complete), "\n")

# ============================================================
# 2. EXCLUDE NON-RELEVANT CONDITIONS
# ============================================================

# These conditions are not included in the main
# autoimmune / immune-mediated disease analysis.

exclude_diseases <- c(
  "staphylococcus infection",
  "streptococcus infection",
  "down syndrome",
  "prediabetic",
  "fever without source",
  "irritable bowel syndrome",
  "diverticular disease",
  "other"
)

data_pca <- data_pca[
  !disease %in% exclude_diseases
]

cat("Samples after disease filtering:", nrow(data_pca), "\n")

# ============================================================
# 3. PCA
# ============================================================

# Scores were already standardized within each dataset.
# We do not standardize again across the pooled samples.
#
# prcomp(center = TRUE, scale. = FALSE) centers each
# module globally, preserving its relative variance.

X <- as.matrix(data_pca[, ..module_cols])

pca <- prcomp(
  X,
  center = TRUE,
  scale. = FALSE
)

variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

cat("\nVariance explained:\n")
print(round(variance[1:5], 2))

# Add PCA coordinates to metadata
data_pca[, `:=`(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3]
)]

# ============================================================
# 4. PCA BY DISEASE
# ============================================================

p_disease <- ggplot(
  data_pca,
  aes(x = PC1, y = PC2, color = disease)
) +
  geom_point(alpha = 0.45, size = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCA - Disease",
    x = paste0("PC1 (", round(variance[1], 1), "%)"),
    y = paste0("PC2 (", round(variance[2], 1), "%)"),
    color = "Disease"
  )

ggsave(
  file.path(output_dir, "PCA_by_disease.pdf"),
  p_disease,
  width = 12,
  height = 8
)

# ============================================================
# 5. PCA BY DATASET
# ============================================================

p_dataset <- ggplot(
  data_pca,
  aes(x = PC1, y = PC2, color = series_id)
) +
  geom_point(alpha = 0.5, size = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCA - Dataset",
    x = paste0("PC1 (", round(variance[1], 1), "%)"),
    y = paste0("PC2 (", round(variance[2], 1), "%)"),
    color = "GSE"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(output_dir, "PCA_by_dataset.pdf"),
  p_dataset,
  width = 10,
  height = 8
)

# ============================================================
# 6. PCA BY HEALTHY / DISEASE
# ============================================================

p_condition <- ggplot(
  data_pca,
  aes(x = PC1, y = PC2, color = condition)
) +
  geom_point(alpha = 0.4, size = 0.9) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCA - Healthy vs Disease",
    x = paste0("PC1 (", round(variance[1], 1), "%)"),
    y = paste0("PC2 (", round(variance[2], 1), "%)")
  )

ggsave(
  file.path(output_dir, "PCA_by_condition.pdf"),
  p_condition,
  width = 9,
  height = 7
)

# ============================================================
# 7. DISEASE CENTROIDS
# ============================================================

# Calculate mean PC coordinates for each disease.
# Show only groups with at least 20 samples.

centroids <- data_pca[
  , .(
    n_samples = .N,
    PC1 = mean(PC1),
    PC2 = mean(PC2)
  ),
  by = disease
][n_samples >= 20]

p_centroids <- ggplot(
  centroids,
  aes(x = PC1, y = PC2, label = disease)
) +
  geom_point(aes(size = n_samples), alpha = 0.7) +
  geom_text(
    check_overlap = TRUE,
    hjust = -0.1,
    size = 3
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Disease centroids in PCA space",
    x = paste0("PC1 (", round(variance[1], 1), "%)"),
    y = paste0("PC2 (", round(variance[2], 1), "%)"),
    size = "Samples"
  )

ggsave(
  file.path(output_dir, "PCA_disease_centroids.pdf"),
  p_centroids,
  width = 12,
  height = 9
)

# ============================================================
# 8. MODULE LOADINGS
# ============================================================

# Identify which modules contribute to each PC.

loadings <- as.data.table(
  pca$rotation[, 1:3, drop = FALSE],
  keep.rownames = "module"
)

fwrite(
  loadings,
  file.path(output_dir, "PCA_module_loadings.tsv"),
  sep = "\t"
)

fwrite(
  data_pca[, c(
    "series_id", "geo_accession",
    "disease", "condition",
    "PC1", "PC2", "PC3"
  ), with = FALSE],
  file.path(output_dir, "PCA_coordinates.tsv"),
  sep = "\t"
)

saveRDS(
  data_pca,
  file.path(output_dir, "PCA_samples.rds")
)

cat("\nPCA completed!\n")
cat("Output:", output_dir, "\n")

# Display plots
print(p_disease)
print(p_dataset)
print(p_condition)
print(p_centroids)
