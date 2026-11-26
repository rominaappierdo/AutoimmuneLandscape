
# ============================================================
# AutoimmuneLandscape
# STEP 5B - RESIDUAL STRUCTURE AFTER PC1 REMOVAL
# ============================================================

library(data.table)
library(ggplot2)
library(cluster)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat_tissue/")

input_dir <- file.path(getwd(), "endotype_analysis", "clustering")
output_dir <- file.path(input_dir, "PC1_residual_analysis")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)
set.seed(123)


# ============================================================
# 1. LOAD PATIENTS AND ORIGINAL MODULE SCORES
# ============================================================

patients <- as.data.table(readRDS(
  file.path(input_dir, "patients_with_clusters.rds")
))

stopifnot(all(c(
  module_cols, "series_id", "geo_accession",
  "disease", "cluster"
) %in% names(patients)))

X <- as.matrix(patients[, ..module_cols])
storage.mode(X) <- "double"

stopifnot(all(is.finite(X)))

# Preserve the original clustering for comparison.
patients[, original_cluster := as.character(cluster)]

cat("Patients:", nrow(patients), "\n")
cat("Datasets:", uniqueN(patients$series_id), "\n")


# ============================================================
# 2. REMOVE THE FIRST PRINCIPAL COMPONENT
# ============================================================

# PCA on original module scores.
# No additional scaling.

pca_original <- prcomp(
  X,
  center = TRUE,
  scale. = FALSE
)

original_variance <- 100 *
  pca_original$sdev^2 / sum(pca_original$sdev^2)

cat("\nOriginal PC1 variance:",
    round(original_variance[1], 2), "%\n")

# Reconstruct the component of each sample explained
# by PC1 and subtract it from the centered matrix.
#
# X_centered = PC1 contribution + residual matrix

X_centered <- sweep(
  X, 2, pca_original$center, "-"
)

pc1_contribution <- tcrossprod(
  pca_original$x[, 1],
  pca_original$rotation[, 1]
)

X_residual <- X_centered - pc1_contribution

colnames(X_residual) <- module_cols

stopifnot(all(is.finite(X_residual)))

# Numerical check: residuals should be orthogonal to PC1.
orthogonality <- max(abs(
  crossprod(pca_original$x[, 1], X_residual)
))

cat("PC1 residual orthogonality check:",
    signif(orthogonality, 4), "\n")

# Save residual module scores separately.
saveRDS(
  X_residual,
  file.path(output_dir, "module_scores_without_PC1.rds")
)


# ============================================================
# 3. PCA OF THE RESIDUAL MATRIX
# ============================================================

# The new first residual PC corresponds to original PC2,
# up to numerical precision and arbitrary axis sign.
#
# Residual PCA reveals the remaining variation but
# does not create new information.

pca_residual <- prcomp(
  X_residual,
  center = TRUE,
  scale. = FALSE
)

residual_variance <- 100 *
  pca_residual$sdev^2 / sum(pca_residual$sdev^2)

patients[, `:=`(
  residual_PC1 = pca_residual$x[, 1],
  residual_PC2 = pca_residual$x[, 2]
)]

cat("\nResidual PCA variance explained:\n")
print(round(residual_variance[1:5], 2))


# ============================================================
# 4. VISUALIZE RESIDUAL PCA
# ============================================================

p_disease <- ggplot(
  patients,
  aes(
    x = residual_PC1,
    y = residual_PC2,
    color = disease
  )
) +
  geom_point(alpha = 0.45, size = 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Residual PCA by disease",
    x = paste0(
      "Residual PC1 (",
      round(residual_variance[1], 1), "%)"
    ),
    y = paste0(
      "Residual PC2 (",
      round(residual_variance[2], 1), "%)"
    )
  )

p_dataset <- ggplot(
  patients,
  aes(
    x = residual_PC1,
    y = residual_PC2,
    color = series_id
  )
) +
  geom_point(alpha = 0.45, size = 0.8) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Residual PCA by dataset",
    x = paste0(
      "Residual PC1 (",
      round(residual_variance[1], 1), "%)"
    ),
    y = paste0(
      "Residual PC2 (",
      round(residual_variance[2], 1), "%)"
    )
  )

ggsave(
  file.path(output_dir, "residual_PCA_by_disease.pdf"),
  p_disease, width = 9, height = 7
)

ggsave(
  file.path(output_dir, "residual_PCA_by_dataset.pdf"),
  p_dataset, width = 9, height = 7
)


# ============================================================
# 5. EVALUATE K-MEANS ON RESIDUAL SCORES
# ============================================================

# K-means is applied to all 15 residual dimensions.
# PC1 has been removed, but the other dimensions
# retain their original relative variances.
#
# Use a fixed subset for silhouette evaluation.

k_values <- 2:8

idx <- sample(
  seq_len(nrow(X_residual)),
  min(2000L, nrow(X_residual))
)

distances <- dist(
  X_residual[idx, , drop = FALSE]
)

results_list <- vector("list", length(k_values))
fits <- vector("list", length(k_values))

names(fits) <- paste0("K", k_values)

for (i in seq_along(k_values)) {
  
  k <- k_values[i]
  
  cat("Evaluating residual K =", k, "\n")
  
  fit <- kmeans(
    X_residual,
    centers = k,
    nstart = 50,
    iter.max = 200
  )
  
  if (fit$ifault != 0) {
    warning("K-means issue at K = ", k)
  }
  
  fits[[i]] <- fit
  
  if (uniqueN(fit$cluster[idx]) >= 2) {
    sil <- silhouette(
      fit$cluster[idx],
      distances
    )
    mean_sil <- mean(sil[, "sil_width"])
  } else {
    mean_sil <- NA_real_
  }
  
  results_list[[i]] <- data.table(
    k = k,
    silhouette = mean_sil,
    explained_ss = fit$betweenss / fit$totss,
    min_cluster_size = min(tabulate(fit$cluster, nbins = k)),
    max_cluster_size = max(tabulate(fit$cluster, nbins = k)),
    ifault = fit$ifault
  )
}

results <- rbindlist(results_list)

print(results)

fwrite(
  results,
  file.path(output_dir, "residual_clustering_k_evaluation.tsv"),
  sep = "\t"
)


# ============================================================
# 6. SELECT PRELIMINARY K
# ============================================================

valid_results <- results[
  is.finite(silhouette) & ifault == 0
]

if (nrow(valid_results) == 0) {
  stop("No valid residual clustering solution")
}

best_k <- valid_results[
  which.max(silhouette), k
]

fit_final <- fits[[paste0("K", best_k)]]

patients[, residual_cluster := factor(
  paste0("R", fit_final$cluster),
  levels = paste0("R", seq_len(best_k))
)]

cat("\nBest residual K:", best_k, "\n")


# ============================================================
# 7. CLUSTER PLOTS
# ============================================================

p_sil <- ggplot(
  results,
  aes(x = k, y = silhouette)
) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = k_values) +
  theme_classic(base_size = 12) +
  labs(
    title = "Residual clustering: silhouette",
    x = "Number of clusters",
    y = "Mean silhouette"
  )

p_clusters <- ggplot(
  patients,
  aes(
    x = residual_PC1,
    y = residual_PC2,
    color = residual_cluster
  )
) +
  geom_point(alpha = 0.45, size = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0(
      "Residual clustering (K = ", best_k, ")"
    ),
    x = "Residual PC1",
    y = "Residual PC2",
    color = "Cluster"
  )

ggsave(
  file.path(output_dir, "residual_silhouette_by_k.pdf"),
  p_sil, width = 7, height = 5
)

ggsave(
  file.path(output_dir, "residual_PCA_clusters.pdf"),
  p_clusters, width = 9, height = 7
)


# ============================================================
# 8. COMPARE ORIGINAL AND RESIDUAL CLUSTERS
# ============================================================

cluster_comparison <- patients[
  , .N,
  by = .(original_cluster, residual_cluster)
]

cluster_comparison[
  , percentage :=
    100 * N / sum(N),
  by = original_cluster
]

print(cluster_comparison)

fwrite(
  cluster_comparison,
  file.path(output_dir, "original_vs_residual_clusters.tsv"),
  sep = "\t"
)

# Check whether residual clusters simply reproduce
# the original high/low PC1 grouping.


# ============================================================
# 9. RESIDUAL CLUSTER COMPOSITION
# ============================================================

cluster_sizes <- patients[
  , .N,
  by = residual_cluster
][order(residual_cluster)]

cluster_disease <- patients[
  , .N,
  by = .(residual_cluster, disease)
]

cluster_disease[
  , percentage_within_disease :=
    100 * N / sum(N),
  by = disease
]

cluster_dataset <- patients[
  , .N,
  by = .(residual_cluster, series_id)
]

cluster_dataset[
  , percentage_within_dataset :=
    100 * N / sum(N),
  by = series_id
]

fwrite(
  cluster_sizes,
  file.path(output_dir, "residual_cluster_sizes.tsv"),
  sep = "\t"
)

fwrite(
  cluster_disease,
  file.path(output_dir, "residual_cluster_disease.tsv"),
  sep = "\t"
)

fwrite(
  cluster_dataset,
  file.path(output_dir, "residual_cluster_dataset.tsv"),
  sep = "\t"
)


# ============================================================
# 10. MODULE PROFILES OF RESIDUAL CLUSTERS
# ============================================================

# IMPORTANT:
# Characterize residual clusters using ORIGINAL module
# scores, not only residual values, to retain
# biological interpretability.

original_module_profiles <- patients[
  , lapply(.SD, mean),
  by = residual_cluster,
  .SDcols = module_cols
]

fwrite(
  original_module_profiles,
  file.path(output_dir, "residual_clusters_original_module_profiles.tsv"),
  sep = "\t"
)

# Also inspect residual module profiles.
residual_dt <- as.data.table(X_residual)
residual_dt[, residual_cluster := patients$residual_cluster]

residual_module_profiles <- residual_dt[
  , lapply(.SD, mean),
  by = residual_cluster,
  .SDcols = module_cols
]

fwrite(
  residual_module_profiles,
  file.path(output_dir, "residual_clusters_residual_module_profiles.tsv"),
  sep = "\t"
)

profiles_long <- melt(
  residual_module_profiles,
  id.vars = "residual_cluster",
  variable.name = "module",
  value.name = "mean_score"
)

p_heatmap <- ggplot(
  profiles_long,
  aes(
    x = module,
    y = residual_cluster,
    fill = mean_score
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Residual module profiles by cluster",
    x = "Module",
    y = "Residual cluster",
    fill = "Mean residual"
  )

ggsave(
  file.path(output_dir, "residual_cluster_module_heatmap.pdf"),
  p_heatmap,
  width = 10,
  height = max(4, 0.6 * best_k + 2)
)


# ============================================================
# 11. SAVE RESULTS
# ============================================================

saveRDS(
  patients,
  file.path(output_dir, "patients_residual_clusters.rds")
)

saveRDS(
  fit_final,
  file.path(output_dir, "residual_kmeans_model.rds")
)

saveRDS(
  pca_original,
  file.path(output_dir, "original_PCA_model.rds")
)

saveRDS(
  pca_residual,
  file.path(output_dir, "residual_PCA_model.rds")
)

fwrite(
  patients[, .(
    series_id,
    geo_accession,
    disease,
    original_cluster,
    residual_cluster,
    residual_PC1,
    residual_PC2
  )],
  file.path(output_dir, "residual_cluster_assignments.tsv"),
  sep = "\t"
)


# ============================================================
# 12. FINAL OUTPUT
# ============================================================

cat("\n====================================\n")
cat("STEP 5B COMPLETED\n")
cat("Patients:", nrow(patients), "\n")
cat("Original PC1 variance:",
    round(original_variance[1], 2), "%\n")
cat("Best residual K:", best_k, "\n")
cat("Output:", output_dir, "\n")
cat("====================================\n")

print(results)
print(cluster_sizes)
print(cluster_comparison)

print(p_sil)
print(p_clusters)
print(p_heatmap)
