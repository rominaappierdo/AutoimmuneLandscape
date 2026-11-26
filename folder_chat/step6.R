
# ============================================================
# AutoimmuneLandscape
# STEP 5A - EXPLORATORY ENDOTYPE CLUSTERING
# ============================================================

library(data.table)
library(ggplot2)
library(cluster)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat")

input_dir <- file.path(getwd(), "endotype_analysis")
output_dir <- file.path(input_dir, "clustering")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)

set.seed(123)


# ============================================================
# 1. LOAD DATA
# ============================================================

data_all <- as.data.table(readRDS(
  file.path(input_dir, "all_samples_with_metadata.rds")
))

stopifnot(all(c(
  module_cols,
  "series_id",
  "geo_accession",
  "disease",
  "condition"
) %in% names(data_all)))

# Exclude conditions unrelated to the main
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

patients <- copy(data_all[
  condition == "Disease" &
    !is.na(disease) &
    !disease %in% exclude_diseases
])

cat("Patients before score filtering:", nrow(patients), "\n")


# ============================================================
# 2. REMOVE INCOMPLETE MODULE SCORE PROFILES
# ============================================================

# Retain patients with 15 finite module scores.
# Exclude NA, NaN and Inf values.

X_initial <- as.matrix(patients[, ..module_cols])

storage.mode(X_initial) <- "double"

complete_idx <- rowSums(
  is.finite(X_initial)
) == length(module_cols)

excluded <- patients[!complete_idx]

patients <- patients[complete_idx]

cat("Patients retained:", nrow(patients), "\n")
cat("Patients excluded:", nrow(excluded), "\n")
cat("Datasets:", uniqueN(patients$series_id), "\n")
cat("Disease labels:", uniqueN(patients$disease), "\n")

fwrite(
  excluded[, .N, by = .(series_id, disease)],
  file.path(output_dir, "excluded_incomplete_samples.tsv"),
  sep = "\t"
)

if (nrow(patients) < 10) {
  stop("Too few patients for clustering")
}


# ============================================================
# 3. PREPARE MODULE SCORE MATRIX
# ============================================================

X <- as.matrix(patients[, ..module_cols])

storage.mode(X) <- "double"

stopifnot(
  ncol(X) == 15,
  all(is.finite(X))
)

# No additional scaling:
# module scores were generated from features
# standardized within their original datasets.
#
# Important: module variances may still differ,
# influencing Euclidean distances and K-means.

module_variance <- data.table(
  module = module_cols,
  mean = colMeans(X),
  sd = apply(X, 2, sd),
  variance = apply(X, 2, var)
)

print(module_variance)

fwrite(
  module_variance,
  file.path(output_dir, "module_variance.tsv"),
  sep = "\t"
)

if (any(module_variance$variance == 0)) {
  stop("One or more modules have zero variance")
}


# ============================================================
# 4. PCA FOR VISUALIZATION
# ============================================================

pca <- prcomp(
  X,
  center = TRUE,
  scale. = FALSE
)

variance_explained <- 100 *
  pca$sdev^2 / sum(pca$sdev^2)

patients[, `:=`(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)]

cat("\nPCA variance explained:\n")
print(round(variance_explained[1:5], 2))


# ============================================================
# 5. EVALUATE K = 2:8
# ============================================================

# Fit K-means on all patients.
#
# Evaluate silhouette on the SAME random subset
# of up to 2000 patients for every K.
#
# This reduces computation and makes comparisons
# between K values more consistent.

k_values <- 2:8

if (nrow(X) <= max(k_values)) {
  stop("Not enough patients for the requested K range")
}

idx <- sample(
  seq_len(nrow(X)),
  size = min(2000L, nrow(X))
)

distances <- dist(
  X[idx, , drop = FALSE],
  method = "euclidean"
)

results_list <- vector("list", length(k_values))
fits <- vector("list", length(k_values))

names(fits) <- paste0("K", k_values)

for (i in seq_along(k_values)) {
  
  k <- k_values[i]
  
  cat("\nEvaluating K =", k, "\n")
  
  fit <- kmeans(
    X,
    centers = k,
    nstart = 50,
    iter.max = 200,
    algorithm = "Hartigan-Wong"
  )
  
  if (fit$ifault != 0) {
    warning(
      "K-means convergence/numerical issue at K = ",
      k,
      " (ifault = ",
      fit$ifault,
      ")"
    )
  }
  
  fits[[i]] <- fit
  
  # Silhouette of the fixed evaluation subset,
  # using cluster assignments from the full-data fit.
  #
  # A cluster may occasionally be absent from the
  # subset; record NA rather than failing.
  
  sampled_clusters <- fit$cluster[idx]
  
  if (length(unique(sampled_clusters)) >= 2) {
    
    sil <- silhouette(
      sampled_clusters,
      distances
    )
    
    mean_sil <- mean(sil[, "sil_width"])
    
  } else {
    
    mean_sil <- NA_real_
    
  }
  
  results_list[[i]] <- data.table(
    k = k,
    silhouette = mean_sil,
    withinss = fit$tot.withinss,
    betweenss = fit$betweenss,
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
  file.path(output_dir, "clustering_k_evaluation.tsv"),
  sep = "\t"
)


# ============================================================
# 6. PLOT SILHOUETTE AND EXPLAINED VARIANCE
# ============================================================

p_sil <- ggplot(
  results,
  aes(x = k, y = silhouette)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = k_values) +
  theme_classic(base_size = 12) +
  labs(
    title = "K-means: silhouette by number of clusters",
    x = "Number of clusters (K)",
    y = "Mean silhouette width"
  )

p_ss <- ggplot(
  results,
  aes(x = k, y = explained_ss)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = k_values) +
  theme_classic(base_size = 12) +
  labs(
    title = "Variance explained by K-means",
    x = "Number of clusters (K)",
    y = "Between-cluster / total sum of squares"
  )

ggsave(
  file.path(output_dir, "silhouette_by_k.pdf"),
  p_sil,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_dir, "explained_variance_by_k.pdf"),
  p_ss,
  width = 7,
  height = 5
)


# ============================================================
# 7. SELECT PRELIMINARY BEST K
# ============================================================

# Select K with the highest mean silhouette.
#
# This is an exploratory choice, not yet a
# biologically validated number of endot ypes.

valid_results <- results[
  is.finite(silhouette) & ifault == 0
]

if (nrow(valid_results) == 0) {
  stop("No valid clustering solution available")
}

best_k <- valid_results[
  which.max(silhouette),
  k
]

cat("\nBest K by silhouette:", best_k, "\n")

# Reuse the already evaluated fit.
# This ensures that the selected clustering
# corresponds exactly to the evaluated solution.

fit_final <- fits[[paste0("K", best_k)]]

patients[, cluster := factor(
  paste0("C", fit_final$cluster),
  levels = paste0("C", seq_len(best_k))
)]


# ============================================================
# 8. PCA COLORED BY CLUSTER
# ============================================================

p_clusters <- ggplot(
  patients,
  aes(
    x = PC1,
    y = PC2,
    color = cluster
  )
) +
  geom_point(
    alpha = 0.5,
    size = 0.8
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0(
      "Exploratory endotype clustering (K = ",
      best_k,
      ")"
    ),
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    color = "Cluster"
  )

ggsave(
  file.path(output_dir, "PCA_clusters.pdf"),
  p_clusters,
  width = 9,
  height = 7
)


# ============================================================
# 9. CLUSTER SIZES
# ============================================================

cluster_sizes <- patients[
  , .(
    n_samples = .N,
    percentage = 100 * .N / nrow(patients)
  ),
  by = cluster
][order(cluster)]

print(cluster_sizes)

fwrite(
  cluster_sizes,
  file.path(output_dir, "cluster_sizes.tsv"),
  sep = "\t"
)


# ============================================================
# 10. CLUSTER COMPOSITION BY DISEASE
# ============================================================

cluster_disease <- patients[
  , .N,
  by = .(cluster, disease)
]

# Percentage of each disease assigned to each cluster.
# Denominator: all patients with that disease.

cluster_disease[
  , percentage_within_disease :=
    100 * N / sum(N),
  by = disease
]

# Percentage of each cluster represented by each disease.
# Denominator: all patients in that cluster.

cluster_disease[
  , percentage_within_cluster :=
    100 * N / sum(N),
  by = cluster
]

fwrite(
  cluster_disease,
  file.path(output_dir, "cluster_disease.tsv"),
  sep = "\t"
)


# ============================================================
# 11. CLUSTER COMPOSITION BY DATASET
# ============================================================

cluster_dataset <- patients[
  , .N,
  by = .(cluster, series_id)
]

cluster_dataset[
  , percentage_within_dataset :=
    100 * N / sum(N),
  by = series_id
]

cluster_dataset[
  , percentage_within_cluster :=
    100 * N / sum(N),
  by = cluster
]

fwrite(
  cluster_dataset,
  file.path(output_dir, "cluster_dataset.tsv"),
  sep = "\t"
)


# ============================================================
# 12. MODULE PROFILES BY CLUSTER
# ============================================================

# Mean activity of each module within each cluster.

module_profiles <- patients[
  , lapply(.SD, mean),
  by = cluster,
  .SDcols = module_cols
]

print(module_profiles)

fwrite(
  module_profiles,
  file.path(output_dir, "cluster_module_profiles.tsv"),
  sep = "\t"
)

# Long-format table for plotting module profiles.

module_profiles_long <- melt(
  module_profiles,
  id.vars = "cluster",
  variable.name = "module",
  value.name = "mean_score"
)

p_modules <- ggplot(
  module_profiles_long,
  aes(
    x = module,
    y = cluster,
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
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean module scores by cluster",
    x = "Module",
    y = "Cluster",
    fill = "Mean score"
  )

ggsave(
  file.path(output_dir, "cluster_module_heatmap.pdf"),
  p_modules,
  width = 10,
  height = max(4, 0.6 * best_k + 2)
)


# ============================================================
# 13. SAVE RESULTS
# ============================================================

saveRDS(
  patients,
  file.path(output_dir, "patients_with_clusters.rds")
)

saveRDS(
  fit_final,
  file.path(output_dir, "kmeans_model.rds")
)

saveRDS(
  fits,
  file.path(output_dir, "kmeans_all_K_models.rds")
)

fwrite(
  patients[, c(
    "series_id",
    "geo_accession",
    "disease",
    "cluster",
    module_cols
  ), with = FALSE],
  file.path(output_dir, "patients_cluster_assignments.tsv"),
  sep = "\t"
)


# ============================================================
# 14. FINAL OUTPUT
# ============================================================

cat("\n====================================\n")
cat("STEP 5A COMPLETED\n")
cat("Patients:", nrow(patients), "\n")
cat("Datasets:", uniqueN(patients$series_id), "\n")
cat("Best K:", best_k, "\n")
cat("Output:", output_dir, "\n")
cat("====================================\n")

print(results)
print(cluster_sizes)
print(p_sil)
print(p_ss)
print(p_clusters)
print(p_modules)









# Correlation between PC1 and module scores
pc1_correlations <- sapply(
  module_cols,
  function(m) cor(patients[[m]], patients$PC1)
)

print(round(sort(pc1_correlations), 3))

# PC1 distribution by cluster
ggplot(patients, aes(x = PC1, fill = cluster)) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  labs(
    title = "PC1 distribution by cluster",
    x = "PC1",
    y = "Density"
  )

# PC1 vs PC2 colored by cluster
ggplot(patients, aes(PC1, PC2, color = cluster)) +
  geom_point(alpha = 0.35, size = 0.7) +
  theme_classic()







# Correlation between PC1 and each module
pc1_cor <- data.table(
  module = module_cols,
  correlation = sapply(
    module_cols,
    function(m) cor(patients[[m]], patients$PC1)
  )
)

print(pc1_cor[order(-abs(correlation))])

# Variance explained by PC1
round(variance_explained[1], 2)




# Correlation of each module with global PC1
# separately within each dataset

pc1_by_gse <- patients[
  , lapply(.SD, function(x) {
    if (.N < 20 || sd(x) == 0 || sd(PC1) == 0)
      return(NA_real_)
    
    cor(x, PC1)
  }),
  by = series_id,
  .SDcols = module_cols
]

# Number of datasets showing the same direction
# as the overall correlation

direction_summary <- data.table(
  module = module_cols,
  n_datasets = sapply(
    module_cols,
    function(m) sum(!is.na(pc1_by_gse[[m]]))
  ),
  n_negative = sapply(
    module_cols,
    function(m) sum(pc1_by_gse[[m]] < 0, na.rm = TRUE)
  ),
  n_positive = sapply(
    module_cols,
    function(m) sum(pc1_by_gse[[m]] > 0, na.rm = TRUE)
  )
)

print(direction_summary)



