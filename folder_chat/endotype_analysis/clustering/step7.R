
library(data.table)
library(ggplot2)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/")

input_dir <- file.path(getwd(), "endotype_analysis", "clustering")
output_dir <- file.path(input_dir, "LODO_validation")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

patients <- as.data.table(readRDS(
  file.path(input_dir, "patients_with_clusters.rds")
))

module_cols <- paste0("mod", 1:15)

X <- as.matrix(patients[, ..module_cols])
storage.mode(X) <- "double"

stopifnot(all(is.finite(X)))

# Global reference PC1
global_pca <- prcomp(X, center = TRUE, scale. = FALSE)
global_loading <- global_pca$rotation[, 1]

gse_ids <- unique(patients$series_id)

results <- list()
projections <- list()

for (gse in gse_ids) {
  
  train_idx <- patients$series_id != gse
  test_idx <- patients$series_id == gse
  
  if (sum(train_idx) < 20 || sum(test_idx) < 20) next
  
  X_train <- X[train_idx, , drop = FALSE]
  X_test <- X[test_idx, , drop = FALSE]
  
  # Train PCA without the held-out dataset
  pca_train <- prcomp(
    X_train,
    center = TRUE,
    scale. = FALSE
  )
  
  loading <- pca_train$rotation[, 1]
  
  # Align the arbitrary PC1 sign to the global reference
  if (sum(loading * global_loading) < 0) {
    loading <- -loading
  }
  
  # Project held-out samples using training centering
  pc1_test <- as.vector(
    sweep(X_test, 2, pca_train$center, "-") %*% loading
  )
  
  # Correlations between projected PC1 and module scores
  correlations <- sapply(
    module_cols,
    function(m) {
      x <- X_test[, m]
      if (sd(x) == 0 || sd(pc1_test) == 0) {
        return(NA_real_)
      }
      cor(x, pc1_test)
    }
  )
  
  results[[gse]] <- data.table(
    series_id = gse,
    n_samples = sum(test_idx),
    module = module_cols,
    correlation = as.numeric(correlations)
  )
  
  projections[[gse]] <- data.table(
    series_id = gse,
    geo_accession = patients$geo_accession[test_idx],
    PC1_LODO = pc1_test
  )
}

lodo_results <- rbindlist(results)
lodo_projections <- rbindlist(projections)

fwrite(
  lodo_results,
  file.path(output_dir, "LODO_module_correlations.tsv"),
  sep = "\t"
)

fwrite(
  lodo_projections,
  file.path(output_dir, "LODO_PC1_projections.tsv"),
  sep = "\t"
)

# Summary across held-out datasets
direction_summary <- lodo_results[
  , .(
    n_datasets = sum(!is.na(correlation)),
    n_negative = sum(correlation < 0, na.rm = TRUE),
    n_positive = sum(correlation > 0, na.rm = TRUE),
    median_correlation = median(correlation, na.rm = TRUE)
  ),
  by = module
]

setorder(direction_summary, module)

print(direction_summary)

fwrite(
  direction_summary,
  file.path(output_dir, "LODO_direction_summary.tsv"),
  sep = "\t"
)

# Heatmap of correlations
p <- ggplot(
  lodo_results,
  aes(x = module, y = series_id, fill = correlation)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Leave-one-dataset-out PC1 validation",
    x = "Module",
    y = "Held-out GSE",
    fill = "Pearson r"
  )

ggsave(
  file.path(output_dir, "LODO_module_correlation_heatmap.pdf"),
  p,
  width = 10,
  height = max(8, 0.22 * uniqueN(lodo_results$series_id))
)

print(p)
