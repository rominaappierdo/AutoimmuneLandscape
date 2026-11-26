# ============================================================
# AutoimmuneLandscape
# ENDOTYPE ANALYSIS BY DATASET
# STEP 3 - K-means clustering and silhouette analysis
# Whole blood - Patients only
# ============================================================
setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/endotype_analysis_by disease")

library(data.table)
library(cluster)
library(ggplot2)

set.seed(123)

# ============================================================
# 1. DIRECTORIES
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/"

analysis_dir <- file.path(
  base_dir,
  "endotype_analysis_by disease"
)

input_dir <- file.path(
  analysis_dir,
  "step1_data_preparation"
)

step2_dir <- file.path(
  analysis_dir,
  "step2_module_correlations"
)

output_dir <- file.path(
  analysis_dir,
  "step3_kmeans_clustering"
)

plot_dir <- file.path(output_dir, "silhouette_plots")
cluster_dir <- file.path(output_dir, "cluster_assignments")
model_dir <- file.path(output_dir, "kmeans_models")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)

k_values <- 2:8

nstart <- 50

max_silhouette_samples <- 2000


# ============================================================
# 2. LOAD DATA
# ============================================================

patients <- as.data.table(readRDS(
  file.path(input_dir, "eligible_patients.rds")
))

# Use exactly the same datasets analyzed in Step 2

eligible_gse <- readRDS(
  file.path(step2_dir, "analyzed_GSEs.rds")
)

patients <- patients[
  series_id %in% eligible_gse
]

cat("Samples:", nrow(patients), "\n")
cat("Datasets:", length(eligible_gse), "\n")

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "disease",
  module_cols
) %in% names(patients)))


# ============================================================
# 3. INITIALIZE RESULTS
# ============================================================

silhouette_results <- list()

best_k_results <- list()

cluster_assignments <- list()

cluster_size_results <- list()


# ============================================================
# 4. K-MEANS WITHIN EACH GSE
# ============================================================

for (gse in eligible_gse) {
  
  cat("\n====================================\n")
  cat("Processing:", gse, "\n")
  cat("====================================\n")
  
  gse_data <- copy(
    patients[series_id == gse]
  )
  
  X <- as.matrix(gse_data[, ..module_cols])
  
  storage.mode(X) <- "double"
  
  stopifnot(all(is.finite(X)))
  
  # No additional scaling
  # Use original module scores
  
  n <- nrow(X)
  
  cat("Samples:", n, "\n")
  
  cat(
    "Diseases:",
    paste(unique(gse_data$disease), collapse = " | "),
    "\n"
  )
  
  # Use the same samples for silhouette evaluation
  # across all tested values of K
  
  idx <- if (n > max_silhouette_samples) {
    
    sort(sample(
      seq_len(n),
      max_silhouette_samples
    ))
    
  } else {
    
    seq_len(n)
    
  }
  
  # Compute distances only once per dataset
  
  distance_matrix <- dist(
    X[idx, , drop = FALSE],
    method = "euclidean"
  )
  
  gse_results <- list()
  
  gse_models <- list()
  
  
  # ==========================================================
  # 5. TEST DIFFERENT VALUES OF K
  # ==========================================================
  
  for (k in k_values) {
    
    if (k >= n || k > nrow(unique(X))) {
      next
    }
    
    cat("Testing K =", k, "\n")
    
    fit <- tryCatch(
      
      kmeans(
        X,
        centers = k,
        nstart = nstart,
        iter.max = 100
      ),
      
      error = function(e) {
        
        cat("K-means failed:", conditionMessage(e), "\n")
        
        return(NULL)
        
      }
      
    )
    
    if (is.null(fit)) {
      next
    }
    
    if (fit$ifault != 0) {
      cat("K-means convergence issue. Skipping K =", k, "\n")
      next
    }
    
    sampled_clusters <- fit$cluster[idx]
    
    # All clusters must be represented in the
    # silhouette evaluation subset
    
    if (uniqueN(sampled_clusters) != k) {
      next
    }
    
    sil <- silhouette(
      sampled_clusters,
      distance_matrix
    )
    
    sil_values <- sil[, "sil_width"]
    
    cluster_sizes <- table(fit$cluster)
    
    gse_results[[as.character(k)]] <- data.table(
      
      series_id = gse,
      
      k = k,
      
      n_samples = n,
      
      n_silhouette_samples = length(idx),
      
      silhouette = mean(sil_values),
      
      median_silhouette = median(sil_values),
      
      proportion_negative_silhouette =
        mean(sil_values < 0),
      
      tot_withinss = fit$tot.withinss,
      
      between_ss_ratio =
        fit$betweenss / fit$totss,
      
      min_cluster_size = min(cluster_sizes),
      
      max_cluster_size = max(cluster_sizes),
      
      min_cluster_proportion =
        min(cluster_sizes) / n
      
    )
    
    gse_models[[as.character(k)]] <- fit
    
  }
  
  
  # ==========================================================
  # 6. SELECT OPTIMAL K
  # ==========================================================
  
  gse_results <- rbindlist(gse_results)
  
  if (nrow(gse_results) == 0) {
    
    cat("No valid clustering solutions for", gse, "\n")
    
    next
    
  }
  
  valid_results <- gse_results[
    is.finite(silhouette)
  ]
  
  if (nrow(valid_results) == 0) {
    next
  }
  
  best_result <- valid_results[
    which.max(silhouette)
  ]
  
  best_k <- best_result$k
  
  best_fit <- gse_models[[as.character(best_k)]]
  
  cat("Optimal K:", best_k, "\n")
  
  cat(
    "Best silhouette:",
    round(best_result$silhouette, 3),
    "\n"
  )
  
  print(table(best_fit$cluster))
  
  
  # ==========================================================
  # 7. SAVE CLUSTER ASSIGNMENTS
  # ==========================================================
  
  gse_data[
    , cluster := best_fit$cluster
  ]
  
  gse_data[
    , optimal_k := best_k
  ]
  
  gse_data[
    , mean_silhouette := best_result$silhouette
  ]
  
  cluster_assignments[[gse]] <- gse_data
  
  saveRDS(
    gse_data,
    file.path(
      cluster_dir,
      paste0(gse, "_cluster_assignments.rds")
    )
  )
  
  fwrite(
    gse_data,
    file.path(
      cluster_dir,
      paste0(gse, "_cluster_assignments.tsv")
    ),
    sep = "\t"
  )
  
  saveRDS(
    best_fit,
    file.path(
      model_dir,
      paste0(gse, "_best_kmeans_model.rds")
    )
  )
  
  
  # ==========================================================
  # 8. CLUSTER SIZE SUMMARY
  # ==========================================================
  
  cluster_sizes <- gse_data[
    , .(
      n_samples = .N,
      proportion = .N / n
    ),
    by = cluster
  ]
  
  cluster_sizes[
    , series_id := gse
  ]
  
  cluster_sizes[
    , optimal_k := best_k
  ]
  
  cluster_size_results[[gse]] <- cluster_sizes
  
  
  # ==========================================================
  # 9. STORE RESULTS
  # ==========================================================
  
  silhouette_results[[gse]] <- gse_results
  
  best_result[
    , n_diseases := uniqueN(gse_data$disease)
  ]
  
  best_result[
    , diseases := paste(
      sort(unique(gse_data$disease)),
      collapse = " | "
    )
  ]
  
  best_k_results[[gse]] <- best_result
  
  
  # ==========================================================
  # 10. SILHOUETTE PLOT FOR EACH GSE
  # ==========================================================
  
  p <- ggplot(
    gse_results,
    aes(x = k, y = silhouette)
  ) +
    
    geom_line(
      linewidth = 0.7,
      color = "grey40"
    ) +
    
    geom_point(
      size = 3,
      color = "#3B6FB6"
    ) +
    
    geom_point(
      data = best_result,
      color = "#C94B4B",
      size = 4
    ) +
    
    geom_hline(
      yintercept = 0.255,
      linetype = "dashed",
      color = "grey50"
    ) +
    
    scale_x_continuous(
      breaks = k_values
    ) +
    
    theme_classic() +
    
    labs(
      
      title = paste0(
        gse,
        " | n = ", n
      ),
      
      subtitle = paste0(
        "Optimal K = ", best_k,
        " | Silhouette = ",
        round(best_result$silhouette, 3)
      ),
      
      x = "Number of clusters (K)",
      
      y = "Mean silhouette score"
      
    )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(gse, "_silhouette.pdf")
    ),
    p,
    width = 6,
    height = 4.5
  )
  
}


# ============================================================
# 11. COMBINE RESULTS
# ============================================================

silhouette_results <- rbindlist(silhouette_results)

best_k_results <- rbindlist(best_k_results)

cluster_assignments <- rbindlist(
  cluster_assignments,
  fill = TRUE
)

cluster_size_results <- rbindlist(
  cluster_size_results
)

if (nrow(best_k_results) == 0) {
  stop("No valid clustering results.")
}

setorder(best_k_results, -silhouette)

print(best_k_results, nrows = Inf)


# ============================================================
# 12. SAVE RESULTS
# ============================================================

fwrite(
  silhouette_results,
  file.path(output_dir, "silhouette_results_all_K.tsv"),
  sep = "\t"
)

fwrite(
  best_k_results,
  file.path(output_dir, "optimal_K_by_GSE.tsv"),
  sep = "\t"
)

fwrite(
  cluster_size_results,
  file.path(output_dir, "cluster_sizes_by_GSE.tsv"),
  sep = "\t"
)

saveRDS(
  silhouette_results,
  file.path(output_dir, "silhouette_results_all_K.rds")
)

saveRDS(
  best_k_results,
  file.path(output_dir, "optimal_K_by_GSE.rds")
)

saveRDS(
  cluster_assignments,
  file.path(output_dir, "all_cluster_assignments.rds")
)

saveRDS(
  cluster_size_results,
  file.path(output_dir, "cluster_sizes_by_GSE.rds")
)


# ============================================================
# 13. DISTRIBUTION OF OPTIMAL K
# ============================================================

k_summary <- best_k_results[
  , .(
    n_datasets = .N
  ),
  by = k
]

k_summary <- merge(
  data.table(k = k_values),
  k_summary,
  by = "k",
  all.x = TRUE
)

k_summary[
  is.na(n_datasets),
  n_datasets := 0L
]

k_summary[
  , percentage :=
    100 * n_datasets / nrow(best_k_results)
]

print(k_summary)

fwrite(
  k_summary,
  file.path(output_dir, "optimal_K_distribution.tsv"),
  sep = "\t"
)

p_k <- ggplot(
  k_summary,
  aes(x = factor(k), y = n_datasets)
) +
  
  geom_col(
    fill = "#3B6FB6",
    width = 0.7
  ) +
  
  geom_text(
    aes(label = n_datasets),
    vjust = -0.5
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  
  theme_classic() +
  
  labs(
    title = "Distribution of optimal K across datasets",
    x = "Optimal number of clusters",
    y = "Number of datasets"
  )

ggsave(
  file.path(output_dir, "optimal_K_distribution.pdf"),
  p_k,
  width = 7,
  height = 5
)


# ============================================================
# 14. BEST SILHOUETTE SCORE BY GSE
# ============================================================

p_best <- ggplot(
  best_k_results,
  aes(
    x = reorder(series_id, silhouette),
    y = silhouette,
    fill = factor(k)
  )
) +
  
  geom_col(width = 0.75) +
  
  geom_hline(
    yintercept = 0.255,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    color = "grey40"
  ) +
  
  coord_flip() +
  
  theme_classic() +
  
  labs(
    
    title = "Best silhouette score within each dataset",
    
    subtitle = "Dashed line: global analysis (0.255)",
    
    x = NULL,
    
    y = "Mean silhouette score",
    
    fill = "Optimal K"
    
  )

ggsave(
  file.path(output_dir, "best_silhouette_by_GSE.pdf"),
  p_best,
  width = 10,
  height = 13
)


# ============================================================
# 15. SILHOUETTE DISTRIBUTION ACROSS K
# ============================================================

p_box <- ggplot(
  silhouette_results,
  aes(x = factor(k), y = silhouette)
) +
  
  geom_boxplot(
    fill = "#B5C9E2",
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.12,
    alpha = 0.5,
    size = 1.5
  ) +
  
  geom_hline(
    yintercept = 0.255,
    linetype = "dashed"
  ) +
  
  theme_classic() +
  
  labs(
    title = "Silhouette scores across datasets",
    x = "Number of clusters (K)",
    y = "Mean silhouette score"
  )

ggsave(
  file.path(output_dir, "silhouette_distribution_by_K.pdf"),
  p_box,
  width = 8,
  height = 5
)


# ============================================================
# 16. ALL SILHOUETTE CURVES
# ============================================================

p_curves <- ggplot(
  silhouette_results,
  aes(
    x = k,
    y = silhouette,
    group = series_id
  )
) +
  
  geom_line(
    alpha = 0.25,
    color = "grey40"
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    linewidth = 1.2,
    color = "#C94B4B"
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "point",
    size = 3,
    color = "#C94B4B"
  ) +
  
  scale_x_continuous(
    breaks = k_values
  ) +
  
  theme_classic() +
  
  labs(
    title = "Silhouette profiles across datasets",
    subtitle = "Red line: median silhouette across datasets",
    x = "Number of clusters (K)",
    y = "Mean silhouette score"
  )

ggsave(
  file.path(output_dir, "all_silhouette_curves.pdf"),
  p_curves,
  width = 8,
  height = 5
)


# ============================================================
# 17. SILHOUETTE VS SAMPLE SIZE
# ============================================================

p_size <- ggplot(
  best_k_results,
  aes(
    x = n_samples,
    y = silhouette,
    color = factor(k)
  )
) +
  
  geom_point(size = 3) +
  
  geom_hline(
    yintercept = 0.255,
    linetype = "dashed"
  ) +
  
  scale_x_log10() +
  
  theme_classic() +
  
  labs(
    title = "Cluster separation and dataset sample size",
    x = "Number of samples (log scale)",
    y = "Best silhouette score",
    color = "Optimal K"
  )

ggsave(
  file.path(output_dir, "silhouette_vs_sample_size.pdf"),
  p_size,
  width = 7,
  height = 5
)


# ============================================================
# 18. FINAL SUMMARY
# ============================================================

cat("\n====================================\n")
cat("ENDOTYPE ANALYSIS BY DATASET\n")
cat("STEP 3 COMPLETED\n")
cat("====================================\n")

cat(
  "Datasets analyzed:",
  nrow(best_k_results), "\n"
)

cat(
  "Median best silhouette:",
  median(best_k_results$silhouette), "\n"
)

cat(
  "Minimum best silhouette:",
  min(best_k_results$silhouette), "\n"
)

cat(
  "Maximum best silhouette:",
  max(best_k_results$silhouette), "\n"
)

cat(
  "Datasets with optimal K = 2:",
  sum(best_k_results$k == 2), "\n"
)

cat(
  "Datasets with silhouette > 0.255:",
  sum(best_k_results$silhouette > 0.255), "\n"
)

cat(
  "Datasets with silhouette >= 0.5:",
  sum(best_k_results$silhouette >= 0.5), "\n"
)

cat(
  "Datasets with silhouette >= 0.7:",
  sum(best_k_results$silhouette >= 0.7), "\n"
)

cat(
  "Datasets with a cluster smaller than 5 samples:",
  sum(best_k_results$min_cluster_size < 5), "\n"
)

cat("\nOptimal K distribution:\n")

print(k_summary)

cat("\nOutput directory:\n", output_dir, "\n")

cat("====================================\n")