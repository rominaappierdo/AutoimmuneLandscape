setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/endotype_analysis_by disease")

# ============================================================
# AutoimmuneLandscape
# ENDOTYPE ANALYSIS BY DATASET
# STEP 2 - Module correlations within each GSE
# Whole blood - Patients only
# ============================================================

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)

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

output_dir <- file.path(
  analysis_dir,
  "step2_module_correlations"
)

heatmap_dir <- file.path(output_dir, "heatmaps")
matrix_dir <- file.path(output_dir, "correlation_matrices")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)


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


# ============================================================
# 3. LOAD DATA
# ============================================================

patients <- as.data.table(readRDS(
  file.path(input_dir, "eligible_patients.rds")
))

gse_summary <- as.data.table(readRDS(
  file.path(input_dir, "GSE_summary.rds")
))

# Exclude datasets with repeated subjects
# or incomplete subject identification

eligible_gse <- gse_summary[
  eligible == TRUE &
    !is.na(n_repeated_subjects) &
    n_repeated_subjects == 0 &
    n_samples_with_subject_id == n_samples,
  series_id
]

patients <- patients[
  series_id %in% eligible_gse
]

cat("Samples:", nrow(patients), "\n")
cat("Datasets:", length(eligible_gse), "\n")

if (length(eligible_gse) == 0) {
  stop("No eligible datasets.")
}

fwrite(
  gse_summary[!series_id %in% eligible_gse],
  file.path(output_dir, "excluded_datasets.tsv"),
  sep = "\t"
)

saveRDS(
  eligible_gse,
  file.path(output_dir, "analyzed_GSEs.rds")
)


# ============================================================
# 4. HEATMAP COLORS
# ============================================================

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#3B6FB6", "white", "#C94B4B")
)


# ============================================================
# 5. HEATMAP FUNCTION
# ============================================================

plot_correlation_heatmap <- function(
    mat,
    title,
    output_file,
    cluster = TRUE,
    legend_title = "Pearson r"
) {
  
  plot_mat <- mat
  
  rownames(plot_mat) <- module_names[rownames(mat)]
  colnames(plot_mat) <- module_names[colnames(mat)]
  
  # Same clustering for rows and columns
  
  clust <- if (cluster) {
    hclust(dist(plot_mat), method = "complete")
  } else {
    FALSE
  }
  
  ht <- Heatmap(
    plot_mat,
    
    name = legend_title,
    col = col_fun,
    
    cluster_rows = clust,
    cluster_columns = clust,
    
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 60,
    
    column_title = title,
    
    rect_gp = gpar(col = "white", lwd = 0.5),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      
      grid.text(
        sprintf("%.2f", plot_mat[i, j]),
        x, y,
        gp = gpar(
          fontsize = 8,
          col = ifelse(
            abs(plot_mat[i, j]) > 0.7,
            "white",
            "black"
          )
        )
      )
      
    }
  )
  
  pdf(output_file, width = 13, height = 11)
  draw(ht)
  dev.off()
}


# ============================================================
# 6. CORRELATION ANALYSIS WITHIN EACH GSE
# ============================================================

cor_matrices <- list()
cor_results <- list()
dataset_results <- list()

for (gse in eligible_gse) {
  
  cat("\nProcessing:", gse, "\n")
  
  gse_data <- patients[series_id == gse]
  
  X <- as.matrix(gse_data[, ..module_cols])
  storage.mode(X) <- "double"
  
  # Pearson correlations between module scores
  
  cor_mat <- cor(
    X,
    method = "pearson"
  )
  
  stopifnot(all(is.finite(cor_mat)))
  
  cor_matrices[[gse]] <- cor_mat
  
  # Save individual correlation matrix
  
  saveRDS(
    cor_mat,
    file.path(matrix_dir, paste0(gse, "_correlation.rds"))
  )
  
  # Extract unique module pairs
  # Exclude diagonal and duplicated correlations
  
  idx <- which(
    upper.tri(cor_mat),
    arr.ind = TRUE
  )
  
  cor_values <- cor_mat[idx]
  
  cor_results[[gse]] <- data.table(
    series_id = gse,
    module1 = rownames(cor_mat)[idx[, 1]],
    module2 = colnames(cor_mat)[idx[, 2]],
    correlation = cor_values,
    n_samples = nrow(gse_data)
  )
  
  # Summary statistics for each dataset
  
  dataset_results[[gse]] <- data.table(
    series_id = gse,
    n_samples = nrow(gse_data),
    n_diseases = uniqueN(gse_data$disease),
    
    diseases = paste(
      sort(unique(gse_data$disease)),
      collapse = " | "
    ),
    
    mean_correlation = mean(cor_values),
    median_correlation = median(cor_values),
    
    mean_absolute_correlation = mean(abs(cor_values)),
    
    proportion_positive = mean(cor_values > 0),
    
    proportion_strong_positive = mean(cor_values >= 0.7),
    
    proportion_strong_negative = mean(cor_values <= -0.7),
    
    proportion_strong_absolute = mean(abs(cor_values) >= 0.7)
  )
  
  # Generate individual heatmap
  
  plot_correlation_heatmap(
    mat = cor_mat,
    
    title = paste0(
      gse,
      " | n = ", nrow(gse_data),
      "\n",
      paste(unique(gse_data$disease), collapse = " / ")
    ),
    
    output_file = file.path(
      heatmap_dir,
      paste0(gse, "_module_correlations.pdf")
    )
  )
}


# ============================================================
# 7. SAVE ALL CORRELATION RESULTS
# ============================================================

cor_results <- rbindlist(cor_results)

dataset_results <- rbindlist(dataset_results)

setorder(dataset_results, -mean_absolute_correlation)

fwrite(
  cor_results,
  file.path(output_dir, "module_correlations_by_GSE.tsv"),
  sep = "\t"
)

fwrite(
  dataset_results,
  file.path(output_dir, "correlation_summary_by_GSE.tsv"),
  sep = "\t"
)

saveRDS(
  cor_matrices,
  file.path(output_dir, "all_correlation_matrices.rds")
)

saveRDS(
  cor_results,
  file.path(output_dir, "module_correlations_by_GSE.rds")
)

print(dataset_results, nrows = Inf)


# ============================================================
# 8. CORRELATION SUMMARY ACROSS DATASETS
# ============================================================

# Each dataset contributes equally to the summary.

pair_summary <- cor_results[
  , .(
    n_datasets = .N,
    
    mean_r = mean(correlation),
    median_r = median(correlation),
    
    sd_r = sd(correlation),
    IQR_r = IQR(correlation),
    
    min_r = min(correlation),
    max_r = max(correlation),
    
    proportion_positive = mean(correlation > 0),
    proportion_negative = mean(correlation < 0),
    
    proportion_strong_positive = mean(correlation >= 0.7),
    
    proportion_strong_negative = mean(correlation <= -0.7)
  ),
  by = .(module1, module2)
]

pair_summary[
  , module1_name := module_names[module1]
]

pair_summary[
  , module2_name := module_names[module2]
]

setorder(pair_summary, -median_r)

fwrite(
  pair_summary,
  file.path(output_dir, "module_pair_summary.tsv"),
  sep = "\t"
)

cat("\nMost positively correlated module pairs:\n")

print(head(pair_summary, 15))


# ============================================================
# 9. BUILD CONSENSUS CORRELATION MATRIX
# ============================================================

# Median Pearson correlation across datasets

median_mat <- matrix(
  NA_real_,
  nrow = length(module_cols),
  ncol = length(module_cols),
  dimnames = list(module_cols, module_cols)
)

diag(median_mat) <- 1

for (i in seq_len(nrow(pair_summary))) {
  
  m1 <- pair_summary$module1[i]
  m2 <- pair_summary$module2[i]
  
  r <- pair_summary$median_r[i]
  
  median_mat[m1, m2] <- r
  median_mat[m2, m1] <- r
}

saveRDS(
  median_mat,
  file.path(output_dir, "median_correlation_matrix.rds")
)

fwrite(
  as.data.table(median_mat, keep.rownames = "module"),
  file.path(output_dir, "median_correlation_matrix.tsv"),
  sep = "\t"
)


# ============================================================
# 10. CONSENSUS CORRELATION HEATMAP
# ============================================================

plot_correlation_heatmap(
  
  mat = median_mat,
  
  title = paste0(
    "Median module correlations across ",
    length(eligible_gse),
    " datasets"
  ),
  
  output_file = file.path(
    output_dir,
    "median_module_correlation_heatmap.pdf"
  ),
  
  legend_title = "Median Pearson r"
)


# ============================================================
# 11. CORRELATION REPRODUCIBILITY ACROSS GSEs
# ============================================================

# Each column is a GSE.
# Each row is one of the 105 unique module pairs.

cor_results[
  , module_pair := paste(module1, module2, sep = " | ")
]

wide <- dcast(
  cor_results,
  module_pair ~ series_id,
  value.var = "correlation"
)

wide_df <- as.data.frame(wide)

rownames(wide_df) <- wide_df$module_pair

cor_profile_mat <- as.matrix(
  wide_df[, -1, drop = FALSE]
)

# Correlate the 105 correlation coefficients
# between every pair of datasets

gse_similarity <- cor(
  cor_profile_mat,
  method = "pearson"
)

saveRDS(
  gse_similarity,
  file.path(output_dir, "GSE_correlation_similarity.rds")
)

fwrite(
  as.data.table(gse_similarity, keep.rownames = "series_id"),
  file.path(output_dir, "GSE_correlation_similarity.tsv"),
  sep = "\t"
)

if (ncol(gse_similarity) >= 2) {
  
  gse_clust <- hclust(
    dist(gse_similarity),
    method = "complete"
  )
  
  ht_similarity <- Heatmap(
    
    gse_similarity,
    
    name = "Pearson r",
    col = col_fun,
    
    cluster_rows = gse_clust,
    cluster_columns = gse_clust,
    
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 60,
    
    column_title = "Similarity of module correlation profiles"
  )
  
  pdf(
    file.path(output_dir, "GSE_correlation_similarity_heatmap.pdf"),
    width = 13,
    height = 12
  )
  
  draw(ht_similarity)
  dev.off()
}


# ============================================================
# 12. DISTRIBUTION OF CORRELATIONS BY DATASET
# ============================================================

dataset_order <- dataset_results$series_id

cor_results[
  , series_id := factor(
    series_id,
    levels = rev(dataset_order)
  )
]

p_distribution <- ggplot(
  cor_results,
  aes(x = series_id, y = correlation)
) +
  
  geom_boxplot(
    outlier.size = 0.5,
    fill = "#B5C9E2"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(-1, 1)
  ) +
  
  theme_classic() +
  
  labs(
    x = NULL,
    y = "Pearson correlation",
    title = "Distribution of module correlations within each dataset"
  )

ggsave(
  file.path(output_dir, "correlation_distribution_by_GSE.pdf"),
  p_distribution,
  width = 10,
  height = 13
)


# ============================================================
# 13. MEAN ABSOLUTE CORRELATION BY DATASET
# ============================================================

p_mean <- ggplot(
  dataset_results,
  aes(
    x = reorder(series_id, mean_absolute_correlation),
    y = mean_absolute_correlation
  )
) +
  
  geom_col(fill = "#3B6FB6") +
  
  coord_flip() +
  
  theme_classic() +
  
  labs(
    x = NULL,
    y = "Mean absolute Pearson correlation",
    title = "Overall module correlation strength by dataset"
  )

ggsave(
  file.path(output_dir, "mean_absolute_correlation_by_GSE.pdf"),
  p_mean,
  width = 9,
  height = 13
)


# ============================================================
# 14. FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("ENDOTYPE ANALYSIS BY DATASET\n")
cat("STEP 2 COMPLETED\n")
cat("========================================\n")

cat("Datasets analyzed:", length(eligible_gse), "\n")

cat("Samples analyzed:", nrow(patients), "\n")

cat("Module pairs per dataset:",
    choose(length(module_cols), 2), "\n")

cat("Total correlations:", nrow(cor_results), "\n")

cat(
  "Median correlation across all dataset-module pairs:",
  median(cor_results$correlation),
  "\n"
)

cat(
  "Median similarity between dataset correlation profiles:",
  median(gse_similarity[upper.tri(gse_similarity)]),
  "\n"
)

cat("\nOutput directory:\n", output_dir, "\n")

cat("========================================\n")