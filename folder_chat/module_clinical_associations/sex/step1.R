# ============================================================
# AutoimmuneLandscape
# SEX ANALYSIS - STEP 1
# Data preparation and quality control
# ============================================================

library(data.table)

# ============================================================
# 1. DIRECTORIES
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat"

input_file <- file.path(
  base_dir,
  "endotype_analysis",
  "all_samples_with_metadata.rds"
)

output_dir <- file.path(
  base_dir,
  "module_clinical_associations",
  "sex"
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
  "gender",
  "disease",
  "condition",
  module_cols
) %in% names(data)))


# ============================================================
# 3. FILTER SAMPLES
# ============================================================

# Keep only samples with known sex (F/M)

sex_data <- copy(data[
  gender %in% c("F", "M")
])

cat("\nSamples with known sex:", nrow(sex_data), "\n")

# Keep samples with all 15 module scores available

X <- as.matrix(sex_data[, ..module_cols])
storage.mode(X) <- "double"

complete_idx <- rowSums(is.finite(X)) == length(module_cols)

sex_data <- sex_data[complete_idx]

cat("Samples with complete module scores:", nrow(sex_data), "\n")
cat("Datasets retained:", uniqueN(sex_data$series_id), "\n")


# ============================================================
# 4. SEX DISTRIBUTION
# ============================================================

sex_summary <- sex_data[
  , .(
    n_samples = .N,
    percentage = 100 * .N / nrow(sex_data)
  ),
  by = gender
][order(gender)]

print(sex_summary)

fwrite(
  sex_summary,
  file.path(output_dir, "sex_summary.tsv"),
  sep = "\t"
)


# ============================================================
# 5. SEX DISTRIBUTION BY DATASET
# ============================================================

sex_by_gse <- sex_data[
  , .(
    n_female = sum(gender == "F"),
    n_male = sum(gender == "M"),
    n_total = .N
  ),
  by = series_id
]

sex_by_gse[
  , `:=`(
    percentage_female = 100 * n_female / n_total,
    percentage_male = 100 * n_male / n_total,
    both_sexes = n_female > 0 & n_male > 0
  )
]

setorder(sex_by_gse, -n_total)

print(sex_by_gse)

fwrite(
  sex_by_gse,
  file.path(output_dir, "sex_by_GSE.tsv"),
  sep = "\t"
)


# ============================================================
# 6. SEX DISTRIBUTION BY DISEASE
# ============================================================

sex_by_disease <- sex_data[
  , .(
    n_female = sum(gender == "F"),
    n_male = sum(gender == "M"),
    n_total = .N
  ),
  by = disease
]

setorder(sex_by_disease, -n_total)

fwrite(
  sex_by_disease,
  file.path(output_dir, "sex_by_disease.tsv"),
  sep = "\t"
)


# ============================================================
# 7. SEX DISTRIBUTION BY CONDITION
# ============================================================

sex_by_condition <- sex_data[
  , .N,
  by = .(condition, gender)
][order(condition, gender)]

print(sex_by_condition)

fwrite(
  sex_by_condition,
  file.path(output_dir, "sex_by_condition.tsv"),
  sep = "\t"
)


# ============================================================
# 8. SAVE FILTERED DATA
# ============================================================

saveRDS(
  sex_data,
  file.path(output_dir, "sex_data.rds")
)


# ============================================================
# 9. FINAL SUMMARY
# ============================================================

cat("\n====================================\n")
cat("SEX ANALYSIS - STEP 1 COMPLETED\n")
cat("====================================\n")

cat("Samples:", nrow(sex_data), "\n")
cat("Females:", sum(sex_data$gender == "F"), "\n")
cat("Males:", sum(sex_data$gender == "M"), "\n")

cat("Datasets:", uniqueN(sex_data$series_id), "\n")

cat(
  "Datasets with both sexes:",
  sum(sex_by_gse$both_sexes),
  "\n"
)

cat("Output:", output_dir, "\n")
cat("====================================\n")