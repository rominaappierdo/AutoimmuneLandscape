
# ============================================================
# AutoimmuneLandscape
# STEP 3 - METADATA AND MODULE SCORES INTEGRATION
# ============================================================

library(data.table)
library(ggplot2)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat")

scores_dir <- file.path(getwd(), "module_scores_final")
output_dir <- file.path(getwd(), "endotype_analysis")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)


# ============================================================
# 1. LOAD METADATA
# ============================================================

# complete_metadata() must be available in the R session.
# Each row represents a metadata annotation, which may
# not necessarily correspond to a unique biological sample.

meta <- as.data.table(complete_metadata())

stopifnot(all(
  c("series_id", "geo_accession", "disease") %in% names(meta)
))

meta[, `:=`(
  series_id = trimws(as.character(series_id)),
  geo_accession = trimws(as.character(geo_accession)),
  disease = trimws(as.character(disease))
)]

# Remove records without valid study/sample identifiers.
meta <- meta[
  !is.na(series_id) & nzchar(series_id) &
    !is.na(geo_accession) & nzchar(geo_accession)
]

cat("Metadata records:", nrow(meta), "\n")
cat("Unique GSM:", uniqueN(meta$geo_accession), "\n")


# ============================================================
# 2. CLASSIFY HEALTHY / DISEASE
# ============================================================

# Classification based on the disease annotation.
# Explicit healthy labels are recognized.
# Generic "control" labels are accepted only when they
# are standalone or clearly describe healthy controls.
#
# All other nonempty labels are provisionally classified
# as Disease, including non-autoimmune conditions.

classify_condition <- function(x) {
  
  x <- tolower(trimws(as.character(x)))
  
  result <- rep(NA_character_, length(x))
  valid <- !is.na(x) & nzchar(x)
  
  healthy <- valid & (
    x %in% c(
      "healthy",
      "control",
      "healthy control",
      "healthy controls",
      "normal",
      "normal control",
      "normal controls"
    ) |
      grepl("\\bhealthy\\b", x) |
      grepl("^non[- ]ibd controls?(,|$)", x)
  )
  
  result[healthy] <- "Healthy"
  result[valid & !healthy] <- "Disease"
  
  return(result)
}

meta[, condition := classify_condition(disease)]

# Save all observed clinical labels for manual inspection.
disease_labels <- meta[
  , .N, by = .(disease, condition)
][order(condition, -N)]

fwrite(
  disease_labels,
  file.path(output_dir, "disease_labels_QC.tsv"),
  sep = "\t"
)

cat("\nInitial condition distribution:\n")
print(meta[, .N, by = condition])


# ============================================================
# 3. RESOLVE DUPLICATE METADATA
# ============================================================

# GSE = dataset
# GSM = sample
#
# A sample may have multiple annotations.
# Different labels do not necessarily imply a conflict.
#
# Example:
# "healthy" and "non-IBD control, healthy"
# both indicate Healthy.
#
# However, Disease vs Healthy for the same GSM
# is treated as ambiguous.

metadata_check <- meta[, .(
  n_records = .N,
  n_conditions = uniqueN(condition, na.rm = FALSE),
  n_disease_labels = uniqueN(disease, na.rm = FALSE),
  disease_original_labels = paste(
    sort(unique(disease[!is.na(disease)])),
    collapse = " | "
  )
), by = .(series_id, geo_accession)]

ambiguous <- metadata_check[n_conditions > 1]

fwrite(
  ambiguous,
  file.path(output_dir, "ambiguous_metadata.tsv"),
  sep = "\t"
)

cat("\nAmbiguous samples:", nrow(ambiguous), "\n")

# Exclude samples with contradictory Healthy/Disease labels.
if (nrow(ambiguous) > 0) {
  meta <- meta[
    !ambiguous,
    on = .(series_id, geo_accession)
  ]
}

# Remove samples without an interpretable condition.
meta <- meta[!is.na(condition)]

# Preserve all original disease labels for each sample.
meta[, disease_original_labels := paste(
  sort(unique(disease)),
  collapse = " | "
), by = .(series_id, geo_accession)]

# Standardize healthy labels.
meta[condition == "Healthy", disease := "healthy"]

# Check whether any diseased sample still has
# multiple distinct disease annotations.
disease_conflicts <- meta[
  condition == "Disease",
  .(
    n_labels = uniqueN(disease),
    labels = paste(sort(unique(disease)), collapse = " | ")
  ),
  by = .(series_id, geo_accession)
][n_labels > 1]

fwrite(
  disease_conflicts,
  file.path(output_dir, "disease_label_conflicts.tsv"),
  sep = "\t"
)

# Do not select a disease label arbitrarily.
# Keep one row per GSM only if the disease annotation
# is unambiguous after standardization.
if (nrow(disease_conflicts) > 0) {
  meta <- meta[
    !disease_conflicts,
    on = .(series_id, geo_accession)
  ]
}

# For remaining duplicated annotations, retain one row.
# Their Healthy/Disease status and disease label agree.
setorder(meta, series_id, geo_accession)

meta_clean <- unique(
  meta,
  by = c("series_id", "geo_accession")
)

stopifnot(!anyDuplicated(
  meta_clean[, .(series_id, geo_accession)]
))

cat("Unique metadata samples:", nrow(meta_clean), "\n")


# ============================================================
# 4. LOAD MODULE SCORES
# ============================================================

score_files <- list.files(
  scores_dir,
  pattern = "_module_scores\\.rds$",
  full.names = TRUE
)

if (length(score_files) == 0) {
  stop("No module score files found")
}

scores_list <- lapply(score_files, function(f) {
  
  x <- as.matrix(readRDS(f))
  
  if (!all(module_cols %in% colnames(x))) {
    stop("Missing module columns in: ", basename(f))
  }
  
  if (is.null(rownames(x))) {
    stop("Missing GSM identifiers in: ", basename(f))
  }
  
  dt <- as.data.table(
    x[, module_cols, drop = FALSE],
    keep.rownames = "geo_accession"
  )
  
  dt[, source_file := basename(f)]
  
  # Extract GSE from the source filename.
  dt[, series_id := sub(
    ".*(GSE[0-9]+).*",
    "\\1",
    source_file
  )]
  
  if (any(!grepl("^GSE[0-9]+$", dt$series_id))) {
    stop("Invalid GSE identifier: ", basename(f))
  }
  
  return(dt)
})

scores <- rbindlist(scores_list, fill = TRUE)

cat("\nScore records:", nrow(scores), "\n")
cat("Source files:", length(score_files), "\n")


# ============================================================
# 5. CHECK DUPLICATE SCORES
# ============================================================

# The same GSE + GSM may appear in multiple source files.
# First, identify all duplicate records.

scores[, n_valid_modules := rowSums(
  is.finite(as.matrix(.SD))
), .SDcols = module_cols]

score_duplicates <- scores[
  , .(
    n_records = .N,
    n_sources = uniqueN(source_file)
  ),
  by = .(series_id, geo_accession)
][n_records > 1]

fwrite(
  score_duplicates,
  file.path(output_dir, "duplicate_score_records.tsv"),
  sep = "\t"
)

cat("\nDuplicated scored samples:",
    nrow(score_duplicates), "\n")

# Compare duplicate score vectors.
# If two finite values disagree, mark the sample
# as discordant rather than silently selecting one.

score_conflicts <- scores[
  , {
    x <- as.matrix(.SD)
    
    inconsistent <- any(vapply(
      seq_len(ncol(x)),
      function(j) {
        v <- x[, j]
        v <- v[is.finite(v)]
        
        length(v) > 1 &&
          (max(v) - min(v) > 1e-8)
      },
      logical(1)
    ))
    
    .(score_conflict = inconsistent)
  },
  by = .(series_id, geo_accession),
  .SDcols = module_cols
][score_conflict == TRUE]

fwrite(
  score_conflicts,
  file.path(output_dir, "discordant_score_records.tsv"),
  sep = "\t"
)

cat("Samples with discordant scores:",
    nrow(score_conflicts), "\n")

# Exclude discordant score records for now.
if (nrow(score_conflicts) > 0) {
  scores <- scores[
    !score_conflicts,
    on = .(series_id, geo_accession)
  ]
}

# For compatible duplicate records, keep the one
# containing the highest number of valid module scores.
setorder(
  scores,
  series_id,
  geo_accession,
  -n_valid_modules,
  source_file
)

scores_clean <- unique(
  scores,
  by = c("series_id", "geo_accession")
)

cat("Unique scored samples:", nrow(scores_clean), "\n")


# ============================================================
# 6. MERGE SCORES AND METADATA
# ============================================================

# Inner join: retain samples present in both sources.

data_all <- merge(
  scores_clean,
  meta_clean,
  by = c("series_id", "geo_accession"),
  all = FALSE
)

stopifnot(!anyDuplicated(
  data_all[, .(series_id, geo_accession)]
))

cat("\nMerged samples:", nrow(data_all), "\n")

cat("\nCondition distribution:\n")
print(data_all[, .N, by = condition])


# ============================================================
# 7. SEPARATE PATIENTS AND CONTROLS
# ============================================================

patients <- data_all[condition == "Disease"]
controls <- data_all[condition == "Healthy"]

# Keep complete score profiles for downstream clustering.
patients_complete <- patients[
  complete.cases(patients[, ..module_cols])
]

controls_complete <- controls[
  complete.cases(controls[, ..module_cols])
]

cat("\nPatients:", nrow(patients), "\n")
cat("Patients with all 15 scores:",
    nrow(patients_complete), "\n")
cat("Healthy controls:", nrow(controls), "\n")
cat("Healthy controls with all 15 scores:",
    nrow(controls_complete), "\n")


# ============================================================
# 8. SAVE RESULTS
# ============================================================

saveRDS(
  data_all,
  file.path(output_dir, "all_samples_with_metadata.rds")
)

saveRDS(
  patients,
  file.path(output_dir, "patients.rds")
)

saveRDS(
  patients_complete,
  file.path(output_dir, "patients_complete.rds")
)

saveRDS(
  controls,
  file.path(output_dir, "healthy_controls.rds")
)

saveRDS(
  controls_complete,
  file.path(output_dir, "healthy_controls_complete.rds")
)

fwrite(
  data_all[, c(
    "series_id",
    "geo_accession",
    "disease",
    "disease_original_labels",
    "condition",
    module_cols
  ), with = FALSE],
  file.path(output_dir, "all_samples_module_scores.tsv"),
  sep = "\t"
)


# ============================================================
# 9. DISEASE SUMMARY
# ============================================================

summary_disease <- data_all[, .(
  n_samples = .N,
  n_complete = sum(complete.cases(.SD)),
  n_datasets = uniqueN(series_id)
), by = .(disease, condition),
.SDcols = module_cols]

setorder(summary_disease, -n_samples)

fwrite(
  summary_disease,
  file.path(output_dir, "disease_summary.tsv"),
  sep = "\t"
)

print(summary_disease)


# ============================================================
# 10. DISEASE DISTRIBUTION PLOT
# ============================================================

p <- ggplot(
  summary_disease,
  aes(
    x = reorder(disease, n_samples),
    y = n_samples,
    fill = condition
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    title = "Samples by disease",
    x = NULL,
    y = "Number of samples",
    fill = "Condition"
  )

ggsave(
  file.path(output_dir, "disease_distribution.pdf"),
  p,
  width = 10,
  height = max(7, 0.3 * nrow(summary_disease))
)

cat("\nINTEGRATION COMPLETED\n")
cat("Output directory:", output_dir, "\n")




# ============================================================
# 3. HARMONIZE DISEASE LABELS AND RESOLVE DUPLICATES
# ============================================================

# Preserve original annotations before harmonization
meta[, disease_original := disease]

# Harmonize known synonymous labels
meta[, disease := tolower(trimws(disease))]

meta[, disease := gsub(
  "rheumatoid arthiritis",
  "rheumatoid arthritis",
  disease,
  fixed = TRUE
)]

meta[disease == "palmoplantar pustular psoriasis",
     disease := "psoriasis (palmoplantar pustular)"]

meta[disease == "diffuse cutaneous systemic sclerosis",
     disease := "systemic sclerosis (diffuse cutaneous)"]

meta[disease == "dermatomyositis with muscular perifascicular atrophy",
     disease := "dermatomyositis with perifascicular atrophy"]

meta[disease == "antiphospholipid syndrome",
     disease := "primary antiphospholipid syndrome"]

meta[disease == "primary antiphospholipid syndrome, systemic lupus erythematosus",
     disease := "systemic lupus erythematosus, antiphospholipid syndrome"]

# Preserve all original labels for each GSM
meta[, disease_original_labels := paste(
  sort(unique(disease_original[!is.na(disease_original)])),
  collapse = " | "
), by = .(series_id, geo_accession)]

# Standardize healthy labels
meta[condition == "Healthy", disease := "healthy"]

# Check remaining conflicts after harmonization
conflicts <- meta[, .(
  n_conditions = uniqueN(condition, na.rm = FALSE),
  n_diseases = uniqueN(disease, na.rm = FALSE),
  labels = paste(sort(unique(disease)), collapse = " | ")
), by = .(series_id, geo_accession)]

print(conflicts[n_conditions > 1 | n_diseases > 1])

# Stop if unresolved conflicts remain
if (nrow(conflicts[n_conditions > 1 | n_diseases > 1]) > 0) {
  stop("Unresolved metadata conflicts")
}

# Now duplicate annotations are consistent
meta_clean <- unique(
  meta,
  by = c("series_id", "geo_accession")
)

cat("Unique metadata samples:", nrow(meta_clean), "\n")













# Numero di campioni
data_all[, .N, by = condition]

# Numero di campioni con tutti i 15 score
data_all[, .(
  n_samples = .N,
  n_complete = sum(complete.cases(.SD)),
  n_GSE = uniqueN(series_id)
), by = condition, .SDcols = module_cols]

# Distribuzione per diagnosi
summary_disease

# Verifica che ogni GSM compaia una sola volta
data_all[, .N, by = .(series_id, geo_accession)][N > 1]





# Missing module scores by disease
missing_by_disease <- data_all[, lapply(.SD, function(x) {
  sum(!is.finite(x))
}), by = disease, .SDcols = module_cols]

print(
  missing_by_disease[
    disease %in% c(
      "sjogren's syndrome",
      "multiple sclerosis (relapsing-remitting)"
    )
  ]
)

# Missing scores by dataset for Sjögren
data_all[
  disease == "sjogren's syndrome",
  c(
    list(n_samples = .N),
    lapply(.SD, function(x) sum(!is.finite(x)))
  ),
  by = series_id,
  .SDcols = module_cols
]








# Keep only samples with 15 finite module scores
complete_scores <- rowSums(
  is.finite(as.matrix(data_all[, ..module_cols]))
) == length(module_cols)

data_complete <- data_all[complete_scores]

# QC: excluded samples by dataset and disease
excluded_samples <- data_all[!complete_scores]

excluded_samples[
  , .N,
  by = .(series_id, disease)
][order(-N)]

# Final dataset for downstream analysis
saveRDS(
  data_complete,
  file.path(output_dir, "all_samples_complete.rds")
)

cat("Samples retained:", nrow(data_complete), "\n")
cat("Samples excluded:", nrow(excluded_samples), "\n")

