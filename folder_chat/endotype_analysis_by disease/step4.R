setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/endotype_analysis_by disease")

# ============================================================
# AutoimmuneLandscape
# ENDOTYPE ANALYSIS BY DATASET
# STEP 4 - Cluster characterization
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
  "step3_kmeans_clustering"
)

output_dir <- file.path(
  analysis_dir,
  "step4_cluster_characterization"
)

individual_plot_dir <- file.path(
  output_dir,
  "individual_GSE_plots"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  individual_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

module_cols <- paste0("mod", 1:15)

min_cluster_size <- 5


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
# 3. LOAD STEP 3 RESULTS
# ============================================================

patients <- as.data.table(readRDS(
  file.path(input_dir, "all_cluster_assignments.rds")
))

best_k_results <- as.data.table(readRDS(
  file.path(input_dir, "optimal_K_by_GSE.rds")
))

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "condition",
  "disease",
  "cluster",
  module_cols
) %in% names(patients)))

stopifnot(all(patients$condition == "Disease"))

cat("Original datasets:", uniqueN(patients$series_id), "\n")
cat("Original samples:", nrow(patients), "\n")


# ============================================================
# 4. SELECT COMPARABLE CLUSTERING SOLUTIONS
# ============================================================

# Include only datasets with:
# - Optimal K = 2
# - At least 5 samples in each cluster

eligible_gse <- best_k_results[
  k == 2 &
    min_cluster_size >= 5,
  series_id
]

excluded_gse <- best_k_results[
  !series_id %in% eligible_gse
]

fwrite(
  excluded_gse,
  file.path(output_dir, "excluded_GSEs.tsv"),
  sep = "\t"
)

patients <- copy(
  patients[series_id %in% eligible_gse]
)

cat("Eligible datasets:", length(eligible_gse), "\n")
cat("Eligible samples:", nrow(patients), "\n")

if (length(eligible_gse) < 2) {
  stop("Insufficient eligible datasets.")
}


# ============================================================
# 5. CALCULATE OVERALL MODULE ACTIVITY
# ============================================================

# Average activity across the 15 original module scores.
# Used only to orient cluster labels consistently.

patients[
  , overall_activity := rowMeans(
    as.matrix(.SD)
  ),
  .SDcols = module_cols
]


# ============================================================
# 6. ALIGN CLUSTER LABELS ACROSS DATASETS
# ============================================================

cluster_activity <- patients[
  , .(
    mean_activity = mean(overall_activity),
    n_samples = .N
  ),
  by = .(series_id, cluster)
]

# Order clusters by increasing overall activity

setorder(
  cluster_activity,
  series_id,
  mean_activity,
  cluster
)

cluster_activity[
  , cluster_label := c("Low", "High"),
  by = series_id
]

# Flag datasets with identical overall cluster means

cluster_activity[
  , orientation_tie := diff(mean_activity) == 0,
  by = series_id
]

patients <- merge(
  patients,
  cluster_activity[
    , .(series_id, cluster, cluster_label)
  ],
  by = c("series_id", "cluster")
)

patients[
  , cluster_label := factor(
    cluster_label,
    levels = c("Low", "High")
  )
]

fwrite(
  cluster_activity,
  file.path(output_dir, "cluster_label_mapping.tsv"),
  sep = "\t"
)

saveRDS(
  patients,
  file.path(output_dir, "aligned_cluster_assignments.rds")
)

print(cluster_activity, nrows = Inf)


# ============================================================
# 7. CONVERT TO LONG FORMAT
# ============================================================

long_data <- melt(
  patients,
  
  id.vars = c(
    "series_id",
    "geo_accession",
    "disease",
    "cluster",
    "cluster_label"
  ),
  
  measure.vars = module_cols,
  
  variable.name = "module",
  
  value.name = "score",
  
  variable.factor = FALSE
)

long_data[
  , module_name := module_names[module]
]


# ============================================================
# 8. MEAN MODULE ACTIVITY PER CLUSTER
# ============================================================

cluster_profiles <- long_data[
  , .(
    mean_score = mean(score),
    median_score = median(score),
    sd_score = sd(score),
    n_samples = .N
  ),
  by = .(
    series_id,
    cluster_label,
    module,
    module_name
  )
]

fwrite(
  cluster_profiles,
  file.path(output_dir, "cluster_module_profiles.tsv"),
  sep = "\t"
)

saveRDS(
  cluster_profiles,
  file.path(output_dir, "cluster_module_profiles.rds")
)


# ============================================================
# 9. CALCULATE HIGH - LOW DIFFERENCES
# ============================================================

module_differences <- dcast(
  cluster_profiles,
  series_id + module + module_name ~ cluster_label,
  value.var = "mean_score"
)

module_differences[
  , delta := High - Low
]

module_differences <- merge(
  module_differences,
  best_k_results[
    , .(
      series_id,
      silhouette,
      n_samples,
      n_diseases,
      diseases
    )
  ],
  by = "series_id"
)

fwrite(
  module_differences,
  file.path(output_dir, "module_differences_by_GSE.tsv"),
  sep = "\t"
)

saveRDS(
  module_differences,
  file.path(output_dir, "module_differences_by_GSE.rds")
)


# ============================================================
# 10. BUILD MODULE DIFFERENCE MATRIX
# ============================================================

wide <- dcast(
  module_differences,
  module ~ series_id,
  value.var = "delta"
)

wide <- as.data.frame(wide)

rownames(wide) <- wide$module

delta_mat <- as.matrix(
  wide[, -1, drop = FALSE]
)

delta_mat <- delta_mat[
  module_cols,
  ,
  drop = FALSE
]

rownames(delta_mat) <- module_names[
  rownames(delta_mat)
]

stopifnot(all(is.finite(delta_mat)))

saveRDS(
  delta_mat,
  file.path(output_dir, "module_difference_matrix.rds")
)


# ============================================================
# 11. HEATMAP COLORS
# ============================================================

lim <- max(abs(delta_mat))

if (lim == 0) {
  lim <- 1
}

col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B6FB6", "white", "#C94B4B")
)


# ============================================================
# 12. ANNOTATIONS
# ============================================================

gse_info <- best_k_results[
  match(colnames(delta_mat), series_id)
]

top_ha <- HeatmapAnnotation(
  
  Silhouette = gse_info$silhouette,
  
  Samples = anno_barplot(
    gse_info$n_samples,
    gp = gpar(fill = "grey60", col = NA),
    border = FALSE
  ),
  
  Diagnoses = ifelse(
    gse_info$n_diseases == 1,
    "Single",
    "Multiple"
  ),
  
  col = list(
    
    Silhouette = colorRamp2(
      c(0, 0.25, 0.5),
      c("white", "#F2C66D", "#A65E18")
    ),
    
    Diagnoses = c(
      Single = "#91B8A4",
      Multiple = "#9C85BB"
    )
    
  ),
  
  annotation_name_gp = gpar(fontsize = 9)
)


# ============================================================
# 13. HEATMAP OF MODULE DIFFERENCES
# ============================================================

row_clust <- hclust(
  dist(delta_mat),
  method = "complete"
)

col_clust <- hclust(
  dist(t(delta_mat)),
  method = "complete"
)

ht_delta <- Heatmap(
  
  delta_mat,
  
  name = "High - Low",
  
  col = col_fun,
  
  top_annotation = top_ha,
  
  cluster_rows = row_clust,
  
  cluster_columns = col_clust,
  
  row_names_side = "left",
  
  row_names_gp = gpar(fontsize = 9),
  
  column_names_gp = gpar(fontsize = 8),
  
  column_names_rot = 60,
  
  rect_gp = gpar(col = "white", lwd = 0.3),
  
  column_title = paste0(
    "Module activity differences across ",
    length(eligible_gse),
    " datasets"
  )
  
)

pdf(
  file.path(output_dir, "module_differences_heatmap.pdf"),
  width = 18,
  height = 9
)

draw(ht_delta)

dev.off()


# ============================================================
# 14. HEATMAP ORDERED BY SILHOUETTE
# ============================================================

gse_order <- best_k_results[
  series_id %in% eligible_gse
][order(-silhouette), series_id]

ht_ordered <- Heatmap(
  
  delta_mat[, gse_order, drop = FALSE],
  
  name = "High - Low",
  
  col = col_fun,
  
  cluster_rows = row_clust,
  
  cluster_columns = FALSE,
  
  row_names_side = "left",
  
  row_names_gp = gpar(fontsize = 9),
  
  column_names_gp = gpar(fontsize = 8),
  
  column_names_rot = 60,
  
  column_title =
    "Module differences | Datasets ordered by silhouette",
  
  rect_gp = gpar(col = "white", lwd = 0.3)
  
)

pdf(
  file.path(
    output_dir,
    "module_differences_ordered_by_silhouette.pdf"
  ),
  width = 18,
  height = 8
)

draw(ht_ordered)

dev.off()


# ============================================================
# 15. HEATMAP OF MEAN CLUSTER PROFILES
# ============================================================

cluster_profiles[
  , profile_id := paste(
    series_id,
    cluster_label,
    sep = " | "
  )
]

profile_wide <- dcast(
  cluster_profiles,
  module ~ profile_id,
  value.var = "mean_score"
)

profile_wide <- as.data.frame(profile_wide)

rownames(profile_wide) <- profile_wide$module

profile_mat <- as.matrix(
  profile_wide[, -1, drop = FALSE]
)

profile_mat <- profile_mat[
  module_cols,
  ,
  drop = FALSE
]

rownames(profile_mat) <- module_names[
  rownames(profile_mat)
]

# Keep Low and High clusters next to each other

profile_order <- unlist(
  lapply(gse_order, function(gse) {
    paste(gse, c("Low", "High"), sep = " | ")
  })
)

profile_mat <- profile_mat[
  ,
  profile_order,
  drop = FALSE
]

profile_lim <- max(abs(profile_mat))

if (profile_lim == 0) {
  profile_lim <- 1
}

profile_col_fun <- colorRamp2(
  c(-profile_lim, 0, profile_lim),
  c("#3B6FB6", "white", "#C94B4B")
)

profile_ha <- HeatmapAnnotation(
  
  Cluster = rep(
    c("Low", "High"),
    times = length(gse_order)
  ),
  
  col = list(
    Cluster = c(
      Low = "#3B6FB6",
      High = "#C94B4B"
    )
  )
)

ht_profiles <- Heatmap(
  
  profile_mat,
  
  name = "Mean score",
  
  col = profile_col_fun,
  
  top_annotation = profile_ha,
  
  cluster_rows = row_clust,
  
  cluster_columns = FALSE,
  
  column_split = factor(
    rep(gse_order, each = 2),
    levels = gse_order
  ),
  
  cluster_column_slices = FALSE,
  
  column_gap = unit(2, "mm"),
  
  show_column_names = FALSE,
  
  column_title_gp = gpar(fontsize = 7),
  
  row_names_side = "left",
  
  row_names_gp = gpar(fontsize = 9),
  
  rect_gp = gpar(col = "white", lwd = 0.3)
  
)

pdf(
  file.path(output_dir, "mean_cluster_profiles_heatmap.pdf"),
  width = 24,
  height = 8
)

draw(ht_profiles)

dev.off()


# ============================================================
# 16. SUMMARY OF MODULE DIFFERENCES ACROSS GSEs
# ============================================================

# Each dataset contributes equally.

module_summary <- module_differences[
  , .(
    
    n_datasets = .N,
    
    mean_delta = mean(delta),
    
    median_delta = median(delta),
    
    Q1 = quantile(delta, 0.25),
    
    Q3 = quantile(delta, 0.75),
    
    sd_delta = sd(delta),
    
    median_absolute_delta = median(abs(delta)),
    
    n_higher_high = sum(delta > 0),
    
    n_higher_low = sum(delta < 0),
    
    proportion_positive = mean(delta > 0),
    
    proportion_negative = mean(delta < 0)
    
  ),
  by = .(module, module_name)
]

setorder(module_summary, -median_delta)

fwrite(
  module_summary,
  file.path(output_dir, "module_difference_summary.tsv"),
  sep = "\t"
)

print(module_summary, nrows = Inf)


# ============================================================
# 17. MEDIAN MODULE DIFFERENCES - BARPLOT
# ============================================================

p_median <- ggplot(
  module_summary,
  
  aes(
    x = reorder(module_name, median_delta),
    y = median_delta,
    fill = median_delta > 0
  )
) +
  
  geom_col(width = 0.7) +
  
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 0.2
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  scale_fill_manual(
    values = c(
      "TRUE" = "#C94B4B",
      "FALSE" = "#3B6FB6"
    ),
    guide = "none"
  ) +
  
  coord_flip() +
  
  theme_classic() +
  
  labs(
    title = "Median module differences across datasets",
    subtitle = "Error bars represent the interquartile range",
    x = NULL,
    y = "Median difference (High - Low)"
  )

ggsave(
  file.path(output_dir, "median_module_differences.pdf"),
  p_median,
  width = 10,
  height = 7
)


# ============================================================
# 18. DISTRIBUTION OF MODULE DIFFERENCES
# ============================================================

module_differences[
  , module_name := factor(
    module_name,
    levels = rev(module_summary$module_name)
  )
]

p_distribution <- ggplot(
  module_differences,
  
  aes(
    x = module_name,
    y = delta
  )
) +
  
  geom_boxplot(
    fill = "grey90",
    outlier.shape = NA,
    width = 0.6
  ) +
  
  geom_jitter(
    aes(color = silhouette),
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  
  scale_color_gradient(
    low = "#B5C9E2",
    high = "#173B70"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  coord_flip() +
  
  theme_classic() +
  
  labs(
    title = "Distribution of module differences across GSEs",
    x = NULL,
    y = "Mean score difference (High - Low)",
    color = "Silhouette"
  )

ggsave(
  file.path(output_dir, "module_difference_distribution.pdf"),
  p_distribution,
  width = 11,
  height = 8
)


# ============================================================
# 19. DIRECTION CONSISTENCY ACROSS DATASETS
# ============================================================

direction_data <- melt(
  
  module_summary,
  
  id.vars = c("module", "module_name"),
  
  measure.vars = c(
    "proportion_positive",
    "proportion_negative"
  ),
  
  variable.name = "direction",
  
  value.name = "proportion"
)

direction_data[
  , direction := fifelse(
    direction == "proportion_positive",
    "Higher in High cluster",
    "Higher in Low cluster"
  )
]

direction_data[
  , module_name := factor(
    module_name,
    levels = rev(module_summary$module_name)
  )
]

p_direction <- ggplot(
  
  direction_data,
  
  aes(
    x = module_name,
    y = proportion,
    fill = direction
  )
) +
  
  geom_col(width = 0.75) +
  
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed"
  ) +
  
  scale_fill_manual(
    values = c(
      "Higher in High cluster" = "#C94B4B",
      "Higher in Low cluster" = "#3B6FB6"
    )
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  
  coord_flip() +
  
  theme_classic() +
  
  theme(legend.position = "bottom") +
  
  labs(
    title = "Consistency of module activity differences",
    x = NULL,
    y = "Proportion of datasets",
    fill = NULL
  )

ggsave(
  file.path(output_dir, "module_direction_consistency.pdf"),
  p_direction,
  width = 10,
  height = 7
)


# ============================================================
# 20. SIMILARITY OF CLUSTER DIFFERENCES ACROSS GSEs
# ============================================================

# Compare the 15-dimensional High - Low profiles
# between each pair of datasets.

# Exclude constant profiles, for which Pearson r
# would be undefined.

valid_profiles <- apply(
  delta_mat,
  2,
  sd
) > 0

similarity_mat <- cor(
  delta_mat[, valid_profiles, drop = FALSE],
  method = "pearson"
)

saveRDS(
  similarity_mat,
  file.path(output_dir, "GSE_cluster_similarity.rds")
)

fwrite(
  as.data.table(
    similarity_mat,
    keep.rownames = "series_id"
  ),
  file.path(output_dir, "GSE_cluster_similarity.tsv"),
  sep = "\t"
)

similarity_col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#3B6FB6", "white", "#C94B4B")
)

if (ncol(similarity_mat) >= 2) {
  
  similarity_clust <- hclust(
    as.dist(1 - similarity_mat),
    method = "complete"
  )
  
  ht_similarity <- Heatmap(
    
    similarity_mat,
    
    name = "Pearson r",
    
    col = similarity_col_fun,
    
    cluster_rows = similarity_clust,
    
    cluster_columns = similarity_clust,
    
    row_names_gp = gpar(fontsize = 8),
    
    column_names_gp = gpar(fontsize = 8),
    
    column_names_rot = 60,
    
    column_title =
      "Similarity of cluster difference profiles across GSEs"
    
  )
  
  pdf(
    file.path(output_dir, "GSE_cluster_similarity_heatmap.pdf"),
    width = 13,
    height = 12
  )
  
  draw(ht_similarity)
  
  dev.off()
  
}


# ============================================================
# 21. CONSENSUS CLUSTER DIFFERENCE PROFILE
# ============================================================

# Median High - Low difference for each module.

consensus_profile <- apply(
  delta_mat,
  1,
  median
)

# Leave-one-GSE-out consensus:
# Each GSE is compared with the median profile
# calculated from all OTHER datasets.

consensus_similarity <- rbindlist(
  lapply(colnames(delta_mat), function(gse) {
    
    reference_profile <- apply(
      delta_mat[
        ,
        colnames(delta_mat) != gse,
        drop = FALSE
      ],
      1,
      median
    )
    
    observed_profile <- delta_mat[, gse]
    
    r <- if (
      sd(reference_profile) > 0 &&
      sd(observed_profile) > 0
    ) {
      cor(
        observed_profile,
        reference_profile,
        method = "pearson"
      )
    } else {
      NA_real_
    }
    
    data.table(
      series_id = gse,
      correlation_with_consensus = r
    )
    
  })
)

consensus_similarity <- merge(
  consensus_similarity,
  best_k_results[
    , .(
      series_id,
      silhouette,
      n_samples,
      diseases
    )
  ],
  by = "series_id"
)

setorder(
  consensus_similarity,
  -correlation_with_consensus
)

fwrite(
  consensus_similarity,
  file.path(output_dir, "consensus_similarity_by_GSE.tsv"),
  sep = "\t"
)


# ============================================================
# 22. CONSENSUS SIMILARITY BARPLOT
# ============================================================

p_consensus <- ggplot(
  
  consensus_similarity,
  
  aes(
    x = reorder(
      series_id,
      correlation_with_consensus
    ),
    
    y = correlation_with_consensus,
    
    fill = correlation_with_consensus > 0
  )
) +
  
  geom_col(width = 0.75) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  scale_fill_manual(
    values = c(
      "TRUE" = "#C94B4B",
      "FALSE" = "#3B6FB6"
    ),
    guide = "none"
  ) +
  
  coord_flip() +
  
  theme_classic() +
  
  labs(
    title = "Similarity to the consensus cluster profile",
    subtitle = "Leave-one-GSE-out consensus",
    x = NULL,
    y = "Pearson correlation"
  )

ggsave(
  file.path(output_dir, "consensus_similarity_by_GSE.pdf"),
  p_consensus,
  width = 8,
  height = 12
)


# ============================================================
# 23. SILHOUETTE VS CONSENSUS SIMILARITY
# ============================================================

p_silhouette <- ggplot(
  
  consensus_similarity,
  
  aes(
    x = silhouette,
    y = correlation_with_consensus
  )
  
) +
  
  geom_point(
    aes(size = n_samples),
    color = "#3B6FB6",
    alpha = 0.75
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 0.255,
    linetype = "dashed"
  ) +
  
  theme_classic() +
  
  labs(
    title = "Cluster separation and molecular reproducibility",
    x = "Mean silhouette score",
    y = "Correlation with consensus profile",
    size = "Samples"
  )

ggsave(
  file.path(output_dir, "silhouette_vs_consensus_similarity.pdf"),
  p_silhouette,
  width = 8,
  height = 6
)


# ============================================================
# 24. INDIVIDUAL DATASET PLOTS
# ============================================================

for (gse in eligible_gse) {
  
  gse_long <- copy(
    long_data[series_id == gse]
  )
  
  gse_long[
    , module_name := factor(
      module_name,
      levels = unname(module_names)
    )
  ]
  
  gse_silhouette <- best_k_results[
    series_id == gse,
    silhouette
  ]
  
  # Distribution of individual module scores
  
  p_modules <- ggplot(
    
    gse_long,
    
    aes(
      x = cluster_label,
      y = score,
      fill = cluster_label
    )
    
  ) +
    
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      width = 0.12,
      size = 0.7,
      alpha = 0.3
    ) +
    
    facet_wrap(
      ~ module_name,
      scales = "free_y",
      ncol = 5
    ) +
    
    scale_fill_manual(
      values = c(
        Low = "#3B6FB6",
        High = "#C94B4B"
      )
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    ) +
    
    labs(
      title = paste0(gse, " | Module activity by cluster"),
      subtitle = paste0(
        "Silhouette = ",
        round(gse_silhouette, 3)
      ),
      x = NULL,
      y = "Module score"
    )
  
  ggsave(
    file.path(
      individual_plot_dir,
      paste0(gse, "_module_boxplots.pdf")
    ),
    p_modules,
    width = 15,
    height = 9
  )
  
  
  # Overall activity distribution
  
  gse_data <- patients[series_id == gse]
  
  p_activity <- ggplot(
    
    gse_data,
    
    aes(
      x = overall_activity,
      fill = cluster_label
    )
    
  ) +
    
    geom_histogram(
      bins = 25,
      position = "identity",
      alpha = 0.6
    ) +
    
    scale_fill_manual(
      values = c(
        Low = "#3B6FB6",
        High = "#C94B4B"
      )
    ) +
    
    theme_classic() +
    
    labs(
      title = paste0(gse, " | Overall module activity"),
      subtitle = paste0(
        "Silhouette = ",
        round(gse_silhouette, 3)
      ),
      x = "Mean activity across 15 modules",
      y = "Number of samples",
      fill = "Cluster"
    )
  
  ggsave(
    file.path(
      individual_plot_dir,
      paste0(gse, "_overall_activity.pdf")
    ),
    p_activity,
    width = 7,
    height = 5
  )
  
}


# ============================================================
# 25. DISEASE COMPOSITION OF CLUSTERS
# ============================================================

multi_disease_gse <- patients[
  , .(
    n_diseases = uniqueN(disease)
  ),
  by = series_id
][n_diseases > 1, series_id]

disease_composition <- patients[
  series_id %in% multi_disease_gse,
  
  .(
    n_samples = .N
  ),
  
  by = .(
    series_id,
    cluster_label,
    disease
  )
]

disease_composition[
  , proportion := n_samples / sum(n_samples),
  by = .(series_id, cluster_label)
]

fwrite(
  disease_composition,
  file.path(output_dir, "disease_composition_by_cluster.tsv"),
  sep = "\t"
)


# ============================================================
# 26. DISEASE COMPOSITION PLOTS
# ============================================================

for (gse in multi_disease_gse) {
  
  plot_data <- disease_composition[
    series_id == gse
  ]
  
  p_disease <- ggplot(
    
    plot_data,
    
    aes(
      x = cluster_label,
      y = proportion,
      fill = disease
    )
    
  ) +
    
    geom_col(
      width = 0.65,
      color = "white"
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format()
    ) +
    
    theme_classic() +
    
    labs(
      title = paste0(
        gse,
        " | Disease composition by cluster"
      ),
      x = "Cluster",
      y = "Proportion of samples",
      fill = "Disease"
    )
  
  ggsave(
    file.path(
      individual_plot_dir,
      paste0(gse, "_disease_composition.pdf")
    ),
    p_disease,
    width = 10,
    height = 5
  )
  
}


# ============================================================
# 27. ASSOCIATION BETWEEN CLUSTER AND DISEASE
# ============================================================

# Exploratory association test for datasets
# containing multiple diagnoses.

set.seed(123)

disease_associations <- rbindlist(
  
  lapply(multi_disease_gse, function(gse) {
    
    gse_data <- patients[series_id == gse]
    
    tab <- table(
      gse_data$cluster_label,
      gse_data$disease
    )
    
    chi <- suppressWarnings(
      chisq.test(tab, correct = FALSE)
    )
    
    # Monte Carlo Fisher test for sparse tables
    
    if (any(chi$expected < 5)) {
      
      test <- fisher.test(
        tab,
        simulate.p.value = TRUE,
        B = 10000
      )
      
      test_name <- "Fisher Monte Carlo"
      
    } else {
      
      test <- chi
      
      test_name <- "Chi-squared"
      
    }
    
    # Cramer's V: strength of association
    
    cramer_v <- sqrt(
      as.numeric(chi$statistic) /
        (
          sum(tab) *
            min(nrow(tab) - 1, ncol(tab) - 1)
        )
    )
    
    data.table(
      series_id = gse,
      n_samples = nrow(gse_data),
      n_diseases = uniqueN(gse_data$disease),
      test = test_name,
      cramer_v = cramer_v,
      pvalue = test$p.value
    )
    
  })
  
)

if (nrow(disease_associations) > 0) {
  
  disease_associations[
    , FDR := p.adjust(
      pvalue,
      method = "BH"
    )
  ]
  
  setorder(disease_associations, FDR)
  
  fwrite(
    disease_associations,
    file.path(output_dir, "cluster_disease_associations.tsv"),
    sep = "\t"
  )
  
  print(disease_associations)
  
}


# ============================================================
# 28. DISEASE ASSOCIATION BARPLOT
# ============================================================

if (nrow(disease_associations) > 0) {
  
  disease_associations[
    , significant := FDR < 0.05
  ]
  
  p_association <- ggplot(
    
    disease_associations,
    
    aes(
      x = reorder(series_id, cramer_v),
      y = cramer_v,
      fill = significant
    )
    
  ) +
    
    geom_col(width = 0.7) +
    
    scale_fill_manual(
      values = c(
        "TRUE" = "#C94B4B",
        "FALSE" = "grey70"
      )
    ) +
    
    coord_flip() +
    
    theme_classic() +
    
    labs(
      title = "Association between cluster and diagnosis",
      x = NULL,
      y = "Cramer's V",
      fill = "FDR < 0.05"
    )
  
  ggsave(
    file.path(output_dir, "cluster_disease_association.pdf"),
    p_association,
    width = 8,
    height = 6
  )
  
}


# ============================================================
# 29. FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("ENDOTYPE ANALYSIS BY DATASET\n")
cat("STEP 4 COMPLETED\n")
cat("========================================\n")

cat(
  "Datasets analyzed:",
  length(eligible_gse), "\n"
)

cat(
  "Samples analyzed:",
  nrow(patients), "\n"
)

cat(
  "Datasets containing multiple diseases:",
  length(multi_disease_gse), "\n"
)

cat(
  "Median pairwise cluster profile similarity:",
  median(
    similarity_mat[upper.tri(similarity_mat)],
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median similarity to leave-one-GSE-out consensus:",
  median(
    consensus_similarity$correlation_with_consensus,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Datasets positively correlated with consensus:",
  sum(
    consensus_similarity$correlation_with_consensus > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat("\nModule difference summary:\n")

print(module_summary, nrows = Inf)

cat("\nConsensus similarity:\n")

print(consensus_similarity, nrows = Inf)

cat("\nDisease associations:\n")

print(disease_associations)

cat("\nOutput directory:\n", output_dir, "\n")

cat("========================================\n")