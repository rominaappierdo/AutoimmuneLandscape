# ============================================================
# AutoimmuneLandscape
# LESION ANALYSIS - STEP 1
# Data preparation and quality control
# Affected vs Unaffected tissues
# ============================================================

library(data.table)

# ============================================================
# 1. DIRECTORIES
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat_tissue/"

input_file <- file.path(
  base_dir,
  "endotype_analysis",
  "all_samples_with_metadata.rds"
)

output_dir <- file.path(
  base_dir,
  "module_clinical_associations",
  "lesion"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)


# ============================================================
# 2. LOAD DATA
# ============================================================

data <- as.data.table(readRDS(input_file))

cat("Original samples:", nrow(data), "\n")
cat("Original datasets:", uniqueN(data$series_id), "\n")

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "disease",
  "sample_source",
  module_cols
) %in% names(data)))


# ============================================================
# 3. EXCLUDE HEALTHY CONTROLS
# ============================================================

cat("\nDisease categories before filtering:\n")
print(data[, .N, by = disease][order(-N)], nrows = Inf)

healthy_labels <- c(
  "healthy",
  "healthy control",
  "healthy controls",
  "control",
  "controls",
  "normal"
)

data[, disease_clean := tolower(trimws(
  as.character(disease)
))]

data <- data[
  !is.na(disease_clean) &
    disease_clean != "" &
    !disease_clean %in% healthy_labels
]

cat("\nSamples after excluding healthy controls:",
    nrow(data), "\n")

cat("Datasets retained:", uniqueN(data$series_id), "\n")


# ============================================================
# 4. CLASSIFY AFFECTED VS UNAFFECTED TISSUES
# ============================================================

data[, sample_source_clean := tolower(trimws(
  as.character(sample_source)
))]

# Unaffected tissue patterns
unaffected_pattern <- paste(
  "unaffected",
  "unaffect",
  "non[-_ ]?lesion",
  "non[-_ ]?inflam",
  "uninflam",
  "non[-_ ]?affect",
  "non[-_ ]?involv",
  sep = "|"
)

# Affected tissue patterns
affected_pattern <- paste(
  "affected",
  "lesion",
  "inflam",
  "involved",
  sep = "|"
)

# Check unaffected patterns first:
# "unaffected" contains "affected"
# "non-lesional" contains "lesion"

data[, lesion_status := fifelse(
  grepl(unaffected_pattern, sample_source_clean),
  "Unaffected",
  fifelse(
    grepl(affected_pattern, sample_source_clean),
    "Affected",
    NA_character_
  )
)]


# ============================================================
# 5. INSPECT SAMPLE CLASSIFICATION
# ============================================================

classification_check <- data[
  , .N,
  by = .(
    disease,
    sample_source,
    lesion_status
  )
][order(disease, lesion_status, -N)]

print(classification_check, nrows = Inf)

fwrite(
  classification_check,
  file.path(output_dir, "sample_source_classification.tsv"),
  sep = "\t"
)

cat("\nClassification summary:\n")

print(
  data[
    , .N,
    by = lesion_status
  ][order(lesion_status)]
)

# IMPORTANT:
# Inspect sample_source_classification.tsv to verify that
# the assigned labels are biologically appropriate.
#
# Samples without an explicit affected/unaffected label
# remain NA and are excluded from the comparison.


# ============================================================
# 6. FILTER SAMPLES
# ============================================================

lesion_data <- copy(data[
  lesion_status %in% c("Affected", "Unaffected")
])

cat("\nSamples with lesion status:", nrow(lesion_data), "\n")

# Keep samples with all 15 module scores available

X <- as.matrix(lesion_data[, ..module_cols])
storage.mode(X) <- "double"

complete_idx <- rowSums(is.finite(X)) == length(module_cols)

lesion_data <- lesion_data[complete_idx]

cat("Samples with complete module scores:",
    nrow(lesion_data), "\n")

cat("Datasets retained:",
    uniqueN(lesion_data$series_id), "\n")


# ============================================================
# 7. OVERALL LESION STATUS DISTRIBUTION
# ============================================================

lesion_summary <- lesion_data[
  , .(
    n_samples = .N,
    percentage = 100 * .N / nrow(lesion_data)
  ),
  by = lesion_status
][order(lesion_status)]

print(lesion_summary)

fwrite(
  lesion_summary,
  file.path(output_dir, "lesion_summary.tsv"),
  sep = "\t"
)


# ============================================================
# 8. LESION DISTRIBUTION BY DATASET
# ============================================================

lesion_by_gse <- lesion_data[
  , .(
    n_affected = sum(lesion_status == "Affected"),
    n_unaffected = sum(lesion_status == "Unaffected"),
    n_total = .N
  ),
  by = series_id
]

lesion_by_gse[
  , both_conditions :=
    n_affected > 0 & n_unaffected > 0
]

setorder(lesion_by_gse, -n_total)

print(lesion_by_gse, nrows = Inf)

fwrite(
  lesion_by_gse,
  file.path(output_dir, "lesion_by_GSE.tsv"),
  sep = "\t"
)


# ============================================================
# 9. LESION DISTRIBUTION BY DISEASE
# ============================================================

lesion_by_disease <- lesion_data[
  , .(
    n_affected = sum(lesion_status == "Affected"),
    n_unaffected = sum(lesion_status == "Unaffected"),
    n_total = .N
  ),
  by = disease
]

setorder(lesion_by_disease, -n_total)

print(lesion_by_disease, nrows = Inf)

fwrite(
  lesion_by_disease,
  file.path(output_dir, "lesion_by_disease.tsv"),
  sep = "\t"
)


# ============================================================
# 10. LESION DISTRIBUTION BY GSE x DISEASE
# ============================================================

lesion_by_group <- lesion_data[
  , .(
    n_affected = sum(lesion_status == "Affected"),
    n_unaffected = sum(lesion_status == "Unaffected"),
    n_total = .N
  ),
  by = .(series_id, disease)
]

lesion_by_group[
  , both_conditions :=
    n_affected > 0 & n_unaffected > 0
]

setorder(lesion_by_group, -n_total)

print(lesion_by_group, nrows = Inf)

fwrite(
  lesion_by_group,
  file.path(output_dir, "lesion_by_GSE_disease.tsv"),
  sep = "\t"
)


# ============================================================
# 11. IDENTIFY GROUPS SUITABLE FOR COMPARISON
# ============================================================

min_samples <- 5

eligible_groups <- lesion_by_group[
  n_affected >= min_samples &
    n_unaffected >= min_samples
]

cat("\nEligible GSE x disease groups:",
    nrow(eligible_groups), "\n")

print(eligible_groups, nrows = Inf)

fwrite(
  eligible_groups,
  file.path(output_dir, "eligible_GSE_disease_groups.tsv"),
  sep = "\t"
)


# ============================================================
# 12. CHECK POTENTIALLY PAIRED SAMPLES
# ============================================================

if ("subject_id" %in% names(lesion_data)) {
  
  paired_check <- lesion_data[
    !is.na(subject_id) &
      trimws(as.character(subject_id)) != "",
    .(
      n_affected = sum(lesion_status == "Affected"),
      n_unaffected = sum(lesion_status == "Unaffected")
    ),
    by = .(series_id, disease, subject_id)
  ]
  
  paired_check[
    , has_both :=
      n_affected > 0 & n_unaffected > 0
  ]
  
  cat("\nPotentially paired subjects:",
      sum(paired_check$has_both), "\n")
  
  print(
    paired_check[has_both == TRUE],
    nrows = Inf
  )
  
  fwrite(
    paired_check,
    file.path(output_dir, "paired_subjects_check.tsv"),
    sep = "\t"
  )
  
} else {
  
  cat("\nNo subject_id column available.",
      "Pairing cannot be assessed.\n")
  
}


# ============================================================
# 13. SAVE FILTERED DATA
# ============================================================

saveRDS(
  lesion_data,
  file.path(output_dir, "lesion_data.rds")
)


# ============================================================
# 14. FINAL SUMMARY
# ============================================================

cat("\n====================================\n")
cat("LESION ANALYSIS - STEP 1 COMPLETED\n")
cat("====================================\n")

cat("Samples:", nrow(lesion_data), "\n")

cat(
  "Affected:",
  sum(lesion_data$lesion_status == "Affected"),
  "\n"
)

cat(
  "Unaffected:",
  sum(lesion_data$lesion_status == "Unaffected"),
  "\n"
)

cat(
  "Datasets:",
  uniqueN(lesion_data$series_id),
  "\n"
)

cat(
  "Datasets with both conditions:",
  sum(lesion_by_gse$both_conditions),
  "\n"
)

cat(
  "GSE x disease groups with both conditions:",
  sum(lesion_by_group$both_conditions),
  "\n"
)

cat(
  "Eligible groups (>= 5 per condition):",
  nrow(eligible_groups),
  "\n"
)

cat("Output:", output_dir, "\n")
cat("====================================\n")