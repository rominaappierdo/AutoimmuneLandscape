# ============================================================
# AutoimmuneLandscape
# SEX ANALYSIS - STEP 3
# Association between sex and module activity
# Patients only - Stratified by GSE and disease
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

output_dir <- file.path(input_dir, "step3_sex_associations")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. LOAD DATA
# ============================================================

sex_data <- as.data.table(readRDS(
  file.path(input_dir, "sex_data.rds")
))

module_cols <- paste0("mod", 1:15)

stopifnot(all(c(
  "series_id", "geo_accession", "gender",
  "condition", "disease", module_cols
) %in% names(sex_data)))

# ============================================================
# 3. FILTER PATIENTS
# ============================================================

patients <- copy(
  sex_data[
    condition == "Disease" &
      gender %in% c("F", "M") &
      !is.na(disease) &
      nzchar(trimws(disease))
  ]
)

# Require complete module scores
X <- as.matrix(patients[, ..module_cols])
storage.mode(X) <- "double"

patients <- patients[
  rowSums(is.finite(X)) == length(module_cols)
]

patients[, gender := factor(gender, levels = c("M", "F"))]

cat("Patients:", nrow(patients), "\n")
cat("GSE:", uniqueN(patients$series_id), "\n")

# ============================================================
# 4. IDENTIFY ELIGIBLE GSE x DISEASE GROUPS
# ============================================================

min_per_sex <- 5

group_summary <- patients[
  , .(
    n_female = sum(gender == "F"),
    n_male = sum(gender == "M"),
    n_total = .N
  ),
  by = .(series_id, disease)
]

group_summary[
  , eligible := n_female >= min_per_sex &
    n_male >= min_per_sex
]

setorder(group_summary, -n_total)

fwrite(
  group_summary,
  file.path(output_dir, "group_eligibility.tsv"),
  sep = "\t"
)

eligible_groups <- group_summary[
  eligible == TRUE,
  .(series_id, disease)
]

patients_eligible <- merge(
  patients,
  eligible_groups,
  by = c("series_id", "disease")
)

cat("Eligible groups:", nrow(eligible_groups), "\n")
cat("Eligible patients:", nrow(patients_eligible), "\n")

# ============================================================
# 5. LINEAR MODELS WITHIN GSE x DISEASE
# ============================================================

results <- patients_eligible[
  , {
    group_data <- .SD
    
    rbindlist(lapply(module_cols, function(mod) {
      
      df <- data.frame(
        score = as.numeric(group_data[[mod]]),
        gender = group_data$gender
      )
      
      fit <- lm(score ~ gender, data = df)
      
      coef_table <- summary(fit)$coefficients
      
      beta <- unname(coef_table["genderF", "Estimate"])
      se <- unname(coef_table["genderF", "Std. Error"])
      p <- unname(coef_table["genderF", "Pr(>|t|)"])
      
      ci <- confint(fit, "genderF", level = 0.95)
      
      data.table(
        module = mod,
        n_female = sum(df$gender == "F"),
        n_male = sum(df$gender == "M"),
        mean_female = mean(df$score[df$gender == "F"]),
        mean_male = mean(df$score[df$gender == "M"]),
        beta_F_vs_M = beta,
        se = se,
        CI_lower = ci[1],
        CI_upper = ci[2],
        pvalue = p
      )
    }))
  },
  by = .(series_id, disease)
]

# ============================================================
# 6. MULTIPLE TESTING CORRECTION
# ============================================================

results[
  , FDR := p.adjust(pvalue, method = "BH"),
  by = .(series_id, disease)
]

results[
  , direction := fifelse(
    beta_F_vs_M > 0,
    "Higher in females",
    fifelse(
      beta_F_vs_M < 0,
      "Higher in males",
      "No difference"
    )
  )
]

setorder(results, series_id, disease, module)

# ============================================================
# 7. SAVE RESULTS
# ============================================================

fwrite(
  results,
  file.path(output_dir, "sex_module_associations_by_GSE_disease.tsv"),
  sep = "\t"
)

saveRDS(
  results,
  file.path(output_dir, "sex_module_associations_by_GSE_disease.rds")
)

# ============================================================
# 8. SUMMARY
# ============================================================

significant_results <- results[
  !is.na(FDR) & FDR < 0.05
]

significant_summary <- significant_results[
  , .(
    n_significant = .N,
    n_higher_females = sum(beta_F_vs_M > 0),
    n_higher_males = sum(beta_F_vs_M < 0)
  ),
  by = module
]

all_modules <- data.table(module = module_cols)

significant_summary <- merge(
  all_modules,
  significant_summary,
  by = "module",
  all.x = TRUE
)

for (col in c(
  "n_significant",
  "n_higher_females",
  "n_higher_males"
)) {
  set(
    significant_summary,
    which(is.na(significant_summary[[col]])),
    col,
    0L
  )
}

setorder(significant_summary, -n_significant)

fwrite(
  significant_summary,
  file.path(output_dir, "significant_associations_summary.tsv"),
  sep = "\t"
)

cat("\n====================================\n")
cat("SEX ANALYSIS - STEP 3 COMPLETED\n")
cat("====================================\n")

cat("Eligible GSE x disease groups:", nrow(eligible_groups), "\n")
cat("Eligible patients:", nrow(patients_eligible), "\n")
cat("Models fitted:", nrow(results), "\n")
cat("Significant associations (FDR < 0.05):",
    nrow(significant_results), "\n")

cat("\nSignificant associations by module:\n")
print(significant_summary)

cat("\nResults saved in:\n", output_dir, "\n")
cat("====================================\n")