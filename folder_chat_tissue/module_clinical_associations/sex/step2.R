# ============================================================
# AutoimmuneLandscape
# SEX ANALYSIS - STEP 2
# Sex distribution by dataset and clinical condition
# ============================================================

library(data.table)

# ============================================================
# 1. DIRECTORIES
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat_tissue/"

input_dir <- file.path(
  base_dir,
  "module_clinical_associations",
  "sex"
)

output_dir <- input_dir

# ============================================================
# 2. LOAD DATA
# ============================================================

sex_data <- as.data.table(readRDS(
  file.path(input_dir, "sex_data.rds")
))

module_cols <- paste0("mod", 1:15)

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "gender",
  "condition",
  "disease",
  module_cols
) %in% names(sex_data)))

cat("Samples:", nrow(sex_data), "\n")
cat("Datasets:", uniqueN(sex_data$series_id), "\n")

# ============================================================
# 3. SEX DISTRIBUTION BY GSE AND CONDITION
# ============================================================

sex_by_gse_condition <- sex_data[
  , .(
    n_female = sum(gender == "F"),
    n_male = sum(gender == "M"),
    n_total = .N
  ),
  by = .(series_id, condition)
]

sex_by_gse_condition[
  , both_sexes := n_female > 0 & n_male > 0
]

setorder(sex_by_gse_condition, series_id, condition)

print(sex_by_gse_condition)

fwrite(
  sex_by_gse_condition,
  file.path(output_dir, "sex_by_GSE_condition.tsv"),
  sep = "\t"
)

# ============================================================
# 4. SEX DISTRIBUTION BY GSE AND DISEASE
# ============================================================

sex_by_gse_disease <- sex_data[
  , .(
    n_female = sum(gender == "F"),
    n_male = sum(gender == "M"),
    n_total = .N
  ),
  by = .(series_id, disease)
]

sex_by_gse_disease[
  , both_sexes := n_female > 0 & n_male > 0
]

setorder(sex_by_gse_disease, series_id, disease)

fwrite(
  sex_by_gse_disease,
  file.path(output_dir, "sex_by_GSE_disease.tsv"),
  sep = "\t"
)

# ============================================================
# 5. SUMMARY
# ============================================================

cat("\n====================================\n")
cat("SEX ANALYSIS - STEP 2 COMPLETED\n")
cat("====================================\n")

cat(
  "GSE-condition groups with both sexes:",
  sum(sex_by_gse_condition$both_sexes),
  "\n"
)

cat(
  "GSE-disease groups with both sexes:",
  sum(sex_by_gse_disease$both_sexes),
  "\n"
)

cat("Output:", output_dir, "\n")
cat("====================================\n")