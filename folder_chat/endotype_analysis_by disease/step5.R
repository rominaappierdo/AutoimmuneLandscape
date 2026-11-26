setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/endotype_analysis_by disease")

# ============================================================
# AutoimmuneLandscape
# ENDOTYPE ANALYSIS BY DATASET
# STEP 5 - Hierarchical clustering of patients
# Whole blood - Patients only
# ============================================================

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)

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
  "step3_kmeans_clustering"
)

output_dir <- file.path(
  analysis_dir,
  "step5_hierarchical_clustering"
)

heatmap_dir <- file.path(
  output_dir,
  "heatmaps"
)

model_dir <- file.path(
  output_dir,
  "hierarchical_models"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)

# Clustering parameters

patient_distance <- "euclidean"
patient_method <- "ward.D2"

module_distance <- "euclidean"
module_method <- "complete"

# Use original module scores
scale_modules <- FALSE


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
  file.path(
    input_dir,
    "all_cluster_assignments.rds"
  )
))

best_k_results <- as.data.table(readRDS(
  file.path(
    input_dir,
    "optimal_K_by_GSE.rds"
  )
))

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "condition",
  "disease",
  "cluster",
  module_cols
) %in% names(patients)))

# Additional safety check: disease patients only

patients <- patients[
  condition == "Disease"
]

eligible_gse <- sort(
  unique(patients$series_id)
)

cat("Samples:", nrow(patients), "\n")
cat("Datasets:", length(eligible_gse), "\n")

stopifnot(length(eligible_gse) > 0)


# ============================================================
# 4. COLOR PALETTES
# ============================================================

# Consistent disease colors across all datasets

disease_levels <- sort(
  unique(as.character(patients$disease))
)

disease_colors <- setNames(
  grDevices::hcl.colors(
    length(disease_levels),
    palette = "Dark 3"
  ),
  disease_levels
)

# K-means cluster colors

cluster_levels <- sort(
  unique(as.character(patients$cluster))
)

cluster_colors <- setNames(
  grDevices::hcl.colors(
    length(cluster_levels),
    palette = "Set 2"
  ),
  cluster_levels
)

sex_colors <- c(
  F = "#D88DA8",
  M = "#729CC8",
  Unknown = "grey85"
)


# ============================================================
# 5. COMMON HEATMAP COLOR SCALE
# ============================================================

# Use the same color scale across all datasets.
# Extreme values are clipped for visualization only.

all_scores <- as.matrix(
  patients[, ..module_cols]
)

storage.mode(all_scores) <- "double"

stopifnot(all(is.finite(all_scores)))

if (scale_modules) {
  
  # Common scale for within-GSE standardized scores
  
  lim <- 3
  
} else {
  
  # Symmetric color scale based on the 99th percentile
  # of absolute module scores across all samples
  
  lim <- as.numeric(
    quantile(
      abs(all_scores),
      probs = 0.99,
      na.rm = TRUE
    )
  )
  
}

if (!is.finite(lim) || lim == 0) {
  lim <- 1
}

col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B6FB6", "white", "#C94B4B")
)

cat("Heatmap color limits:", -lim, lim, "\n")


# ============================================================
# 6. INITIALIZE RESULTS
# ============================================================

hierarchical_models <- list()

hierarchical_assignments <- list()

dataset_summary <- list()


# ============================================================
# 7. HIERARCHICAL CLUSTERING WITHIN EACH GSE
# ============================================================

for (gse in eligible_gse) {
  
  cat("\n====================================\n")
  cat("Processing:", gse, "\n")
  cat("====================================\n")
  
  gse_data <- copy(
    patients[series_id == gse]
  )
  
  stopifnot(
    !anyDuplicated(gse_data$geo_accession)
  )
  
  n <- nrow(gse_data)
  
  gse_info <- best_k_results[
    series_id == gse
  ]
  
  best_k <- gse_info$k[1]
  
  best_silhouette <- gse_info$silhouette[1]
  
  cat("Samples:", n, "\n")
  cat("Optimal K-means K:", best_k, "\n")
  cat("K-means silhouette:", best_silhouette, "\n")
  
  
  # ==========================================================
  # 8. BUILD MODULE SCORE MATRIX
  # ==========================================================
  
  # Rows = patients
  # Columns = modules
  
  X <- as.matrix(
    gse_data[, ..module_cols]
  )
  
  storage.mode(X) <- "double"
  
  rownames(X) <- gse_data$geo_accession
  
  stopifnot(all(is.finite(X)))
  
  # Optional standardization within each GSE
  
  if (scale_modules) {
    
    X <- scale(X)
    
    stopifnot(all(is.finite(X)))
    
  }
  
  # Transpose for the heatmap
  # Rows = modules
  # Columns = patients
  
  mat <- t(X)
  
  rownames(mat) <- module_names[
    rownames(mat)
  ]
  
  
  # ==========================================================
  # 9. HIERARCHICAL CLUSTERING OF PATIENTS
  # ==========================================================
  
  patient_dist <- dist(
    X,
    method = patient_distance
  )
  
  patient_hclust <- hclust(
    patient_dist,
    method = patient_method
  )
  
  
  # ==========================================================
  # 10. HIERARCHICAL CLUSTERING OF MODULES
  # ==========================================================
  
  module_dist <- dist(
    mat,
    method = module_distance
  )
  
  module_hclust <- hclust(
    module_dist,
    method = module_method
  )
  
  
  # ==========================================================
  # 11. OPTIONAL CUT OF THE PATIENT DENDROGRAM
  # ==========================================================
  
  # Cut the hierarchical tree using the optimal K
  # previously identified by K-means.
  #
  # This is used only to compare the two methods.
  # The hierarchical tree itself is not constrained by K.
  
  hc_clusters <- cutree(
    patient_hclust,
    k = best_k
  )
  
  gse_data[
    , hierarchical_cluster :=
      unname(hc_clusters[geo_accession])
  ]
  
  gse_data[
    , hierarchical_order :=
      match(
        geo_accession,
        patient_hclust$labels[
          patient_hclust$order
        ]
      )
  ]
  
  hierarchical_assignments[[gse]] <- gse_data
  
  
  # ==========================================================
  # 12. BUILD SAMPLE ANNOTATIONS
  # ==========================================================
  
  annotation_data <- data.frame(
    
    Disease = as.character(
      gse_data$disease
    ),
    
    Kmeans = as.character(
      gse_data$cluster
    ),
    
    Hierarchical = as.character(
      gse_data$hierarchical_cluster
    ),
    
    row.names = gse_data$geo_accession
    
  )
  
  annotation_colors <- list(
    
    Disease = disease_colors[
      unique(annotation_data$Disease)
    ],
    
    Kmeans = cluster_colors[
      unique(annotation_data$Kmeans)
    ],
    
    Hierarchical = cluster_colors[
      unique(annotation_data$Hierarchical)
    ]
    
  )
  
  
  # Add sex if available
  
  if ("gender" %in% names(gse_data)) {
    
    sex <- as.character(gse_data$gender)
    
    sex[
      is.na(sex) |
        !sex %in% c("F", "M")
    ] <- "Unknown"
    
    if (any(sex %in% c("F", "M"))) {
      
      annotation_data$Sex <- sex
      
      annotation_colors$Sex <- sex_colors[
        unique(sex)
      ]
      
    }
    
  }
  
  
  # Add age if available and informative
  
  if ("age" %in% names(gse_data)) {
    
    age <- suppressWarnings(
      as.numeric(
        as.character(gse_data$age)
      )
    )
    
    age[!is.finite(age)] <- NA_real_
    
    if (
      sum(!is.na(age)) >= 2 &&
      length(unique(age[!is.na(age)])) > 1
    ) {
      
      annotation_data$Age <- age
      
      age_range <- range(
        age,
        na.rm = TRUE
      )
      
      annotation_colors$Age <- colorRamp2(
        age_range,
        c("#F7E7A9", "#80501A")
      )
      
    }
    
  }
  
  top_ha <- HeatmapAnnotation(
    
    df = annotation_data,
    
    col = annotation_colors,
    
    na_col = "grey85",
    
    annotation_name_side = "left",
    
    annotation_name_gp = gpar(
      fontsize = 9
    ),
    
    simple_anno_size = unit(3, "mm")
    
  )
  
  
  # ==========================================================
  # 13. HEATMAP TITLE
  # ==========================================================
  
  disease_title <- paste(
    sort(unique(gse_data$disease)),
    collapse = " / "
  )
  
  heatmap_title <- paste0(
    
    gse,
    
    " | n = ", n,
    
    "\n",
    
    disease_title,
    
    "\n",
    
    "K-means K = ", best_k,
    
    " | Silhouette = ",
    
    round(best_silhouette, 3)
    
  )
  
  
  # ==========================================================
  # 14. GENERATE HIERARCHICAL CLUSTERING HEATMAP
  # ==========================================================
  
  ht <- Heatmap(
    
    mat,
    
    name = "Module score",
    
    col = col_fun,
    
    top_annotation = top_ha,
    
    cluster_rows = module_hclust,
    
    cluster_columns = patient_hclust,
    
    # Preserve the original dendrogram ordering
    
    row_dend_reorder = FALSE,
    
    column_dend_reorder = FALSE,
    
    show_row_dend = TRUE,
    
    show_column_dend = TRUE,
    
    row_dend_width = unit(20, "mm"),
    
    column_dend_height = unit(25, "mm"),
    
    row_names_side = "left",
    
    row_names_gp = gpar(
      fontsize = 10
    ),
    
    # Show sample names only in small datasets
    
    show_column_names = n <= 60,
    
    column_names_gp = gpar(
      fontsize = 5
    ),
    
    column_names_rot = 90,
    
    column_title = heatmap_title,
    
    column_title_gp = gpar(
      fontsize = 12,
      fontface = "bold"
    ),
    
    rect_gp = gpar(
      col = NA
    ),
    
    use_raster = n > 200,
    
    heatmap_legend_param = list(
      
      title = "Module score",
      
      at = c(-lim, 0, lim),
      
      labels = c(
        round(-lim, 2),
        0,
        round(lim, 2)
      )
      
    )
    
  )
  
  
  # ==========================================================
  # 15. SAVE HEATMAP
  # ==========================================================
  
  pdf(
    
    file.path(
      heatmap_dir,
      paste0(gse, "_hierarchical_heatmap.pdf")
    ),
    
    width = 16,
    
    height = 8
    
  )
  
  draw(
    
    ht,
    
    heatmap_legend_side = "right",
    
    annotation_legend_side = "right",
    
    merge_legends = TRUE
    
  )
  
  dev.off()
  
  
  # ==========================================================
  # 16. SAVE SECOND HEATMAP WITH CLUSTER SEPARATION
  # ==========================================================
  
  # Same hierarchical clustering, but with visible
  # gaps between branches defined by cutree().
  
  ht_split <- Heatmap(
    
    mat,
    
    name = "Module score",
    
    col = col_fun,
    
    top_annotation = top_ha,
    
    cluster_rows = module_hclust,
    
    cluster_columns = patient_hclust,
    
    row_dend_reorder = FALSE,
    
    column_dend_reorder = FALSE,
    
    column_split = best_k,
    
    column_gap = unit(3, "mm"),
    
    show_column_names = FALSE,
    
    show_row_dend = TRUE,
    
    show_column_dend = TRUE,
    
    row_names_side = "left",
    
    row_names_gp = gpar(
      fontsize = 10
    ),
    
    column_title = "Hierarchical clusters",
    
    rect_gp = gpar(col = NA),
    
    use_raster = n > 200
    
  )
  
  pdf(
    
    file.path(
      heatmap_dir,
      paste0(gse, "_hierarchical_heatmap_split.pdf")
    ),
    
    width = 16,
    
    height = 8
    
  )
  
  draw(
    
    ht_split,
    
    column_title = heatmap_title,
    
    heatmap_legend_side = "right",
    
    annotation_legend_side = "right",
    
    merge_legends = TRUE
    
  )
  
  dev.off()
  
  
  # ==========================================================
  # 17. SAVE HIERARCHICAL CLUSTERING MODELS
  # ==========================================================
  
  model <- list(
    
    series_id = gse,
    
    patient_hclust = patient_hclust,
    
    module_hclust = module_hclust,
    
    patient_distance = patient_distance,
    
    patient_method = patient_method,
    
    module_distance = module_distance,
    
    module_method = module_method,
    
    scaled_modules = scale_modules,
    
    k = best_k,
    
    hierarchical_clusters = hc_clusters
    
  )
  
  hierarchical_models[[gse]] <- model
  
  saveRDS(
    
    model,
    
    file.path(
      model_dir,
      paste0(gse, "_hierarchical_model.rds")
    )
    
  )
  
  
  # ==========================================================
  # 18. DATASET SUMMARY
  # ==========================================================
  
  hc_sizes <- table(
    gse_data$hierarchical_cluster
  )
  
  dataset_summary[[gse]] <- data.table(
    
    series_id = gse,
    
    n_samples = n,
    
    diseases = disease_title,
    
    k = best_k,
    
    kmeans_silhouette = best_silhouette,
    
    min_hierarchical_cluster_size = min(hc_sizes),
    
    max_hierarchical_cluster_size = max(hc_sizes),
    
    hierarchical_cluster_sizes = paste(
      hc_sizes,
      collapse = " / "
    )
    
  )
  
  cat("Hierarchical cluster sizes:\n")
  
  print(hc_sizes)
  
}


# ============================================================
# 19. COMBINE RESULTS
# ============================================================

hierarchical_assignments <- rbindlist(
  hierarchical_assignments,
  fill = TRUE
)

dataset_summary <- rbindlist(
  dataset_summary
)

setorder(
  dataset_summary,
  -kmeans_silhouette
)

setorder(
  hierarchical_assignments,
  series_id,
  hierarchical_order
)


# ============================================================
# 20. SAVE RESULTS
# ============================================================

saveRDS(
  
  hierarchical_assignments,
  
  file.path(
    output_dir,
    "all_hierarchical_cluster_assignments.rds"
  )
  
)

saveRDS(
  
  hierarchical_models,
  
  file.path(
    output_dir,
    "all_hierarchical_models.rds"
  )
  
)

fwrite(
  
  dataset_summary,
  
  file.path(
    output_dir,
    "hierarchical_clustering_summary.tsv"
  ),
  
  sep = "\t"
  
)

fwrite(
  
  hierarchical_assignments,
  
  file.path(
    output_dir,
    "all_hierarchical_cluster_assignments.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 21. FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("ENDOTYPE ANALYSIS BY DATASET\n")
cat("STEP 5 COMPLETED\n")
cat("========================================\n")

cat(
  "Datasets analyzed:",
  nrow(dataset_summary),
  "\n"
)

cat(
  "Samples analyzed:",
  nrow(hierarchical_assignments),
  "\n"
)

cat(
  "Heatmaps generated:",
  2 * nrow(dataset_summary),
  "\n"
)

cat(
  "Patient clustering:",
  patient_distance,
  "+",
  patient_method,
  "\n"
)

cat(
  "Module clustering:",
  module_distance,
  "+",
  module_method,
  "\n"
)

cat(
  "Module standardization:",
  scale_modules,
  "\n"
)

print(
  dataset_summary,
  nrows = Inf
)

cat("\nOutput directory:\n", output_dir, "\n")

cat("========================================\n")