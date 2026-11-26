# ============================================================
# AutoimmuneLandscape
# LESION ANALYSIS - STEP 2
# Affected vs Unaffected tissue distribution
# by dataset and clinical condition
# ============================================================

library(data.table)

# ============================================================
# 1. DIRECTORIES
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat_tissue/"

input_dir <- file.path(
  base_dir,
  "module_clinical_associations",
  "lesion"
)

output_dir <- input_dir

# ============================================================
# 2. LOAD DATA
# ============================================================

lesion_data <- as.data.table(readRDS(
  file.path(input_dir, "lesion_data.rds")
))

module_cols <- paste0("mod", 1:15)

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "lesion_status",
  "condition",
  "disease",
  module_cols
) %in% names(lesion_data)))

cat("Samples:", nrow(lesion_data), "\n")
cat("Datasets:", uniqueN(lesion_data$series_id), "\n")

# ============================================================
# 3. LESION DISTRIBUTION BY GSE AND CONDITION
# ============================================================

lesion_by_gse_condition <- lesion_data[
  , .(
    n_affected = sum(lesion_status == "Affected"),
    n_unaffected = sum(lesion_status == "Unaffected"),
    n_total = .N
  ),
  by = .(series_id, condition)
]

lesion_by_gse_condition[
  , both_conditions := n_affected > 0 & n_unaffected > 0
]

setorder(lesion_by_gse_condition, series_id, condition)

print(lesion_by_gse_condition, nrows = Inf)

fwrite(
  lesion_by_gse_condition,
  file.path(output_dir, "lesion_by_GSE_condition.tsv"),
  sep = "\t"
)

# ============================================================
# 4. LESION DISTRIBUTION BY GSE AND DISEASE
# ============================================================

lesion_by_gse_disease <- lesion_data[
  , .(
    n_affected = sum(lesion_status == "Affected"),
    n_unaffected = sum(lesion_status == "Unaffected"),
    n_total = .N
  ),
  by = .(series_id, disease)
]

lesion_by_gse_disease[
  , both_conditions := n_affected > 0 & n_unaffected > 0
]

setorder(lesion_by_gse_disease, series_id, disease)

print(lesion_by_gse_disease, nrows = Inf)

fwrite(
  lesion_by_gse_disease,
  file.path(output_dir, "lesion_by_GSE_disease.tsv"),
  sep = "\t"
)

# ============================================================
# 5. IDENTIFY GROUPS ELIGIBLE FOR COMPARISON
# ============================================================

min_samples <- 5

eligible_groups <- lesion_by_gse_disease[
  n_affected >= min_samples &
    n_unaffected >= min_samples
]

cat(
  "\nEligible GSE-disease groups (>= 5 per condition):",
  nrow(eligible_groups),
  "\n"
)

print(eligible_groups, nrows = Inf)

fwrite(
  eligible_groups,
  file.path(output_dir, "eligible_GSE_disease_groups.tsv"),
  sep = "\t"
)

# ============================================================
# 6. SUMMARY
# ============================================================

cat("\n====================================\n")
cat("LESION ANALYSIS - STEP 2 COMPLETED\n")
cat("====================================\n")

cat(
  "GSE-condition groups with both lesion statuses:",
  sum(lesion_by_gse_condition$both_conditions),
  "\n"
)

cat(
  "GSE-disease groups with both lesion statuses:",
  sum(lesion_by_gse_disease$both_conditions),
  "\n"
)

cat(
  "Eligible GSE-disease groups (>= 5 per status):",
  nrow(eligible_groups),
  "\n"
)

cat("Output:", output_dir, "\n")
cat("====================================\n")