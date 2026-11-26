# ============================================================
# AutoimmuneLandscape
# LESION ANALYSIS - STEP 3
# Association between lesion status and module activity
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
  "lesion"
)

output_dir <- file.path(
  input_dir,
  "step3_lesion_associations"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)

min_per_status <- 5


# ============================================================
# 2. LOAD DATA
# ============================================================

lesion_data <- as.data.table(readRDS(
  file.path(input_dir, "lesion_data.rds")
))

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "lesion_status",
  "disease",
  module_cols
) %in% names(lesion_data)))

cat("Samples:", nrow(lesion_data), "\n")
cat("Datasets:", uniqueN(lesion_data$series_id), "\n")


# ============================================================
# 3. FILTER PATIENTS
# ============================================================

# Healthy controls were excluded in Step 1.
# Apply the exclusion again as a safeguard.

healthy_labels <- c(
  "healthy",
  "healthy control",
  "healthy controls",
  "control",
  "controls",
  "normal"
)

patients <- copy(lesion_data[
  lesion_status %in% c("Affected", "Unaffected") &
    !is.na(disease) &
    nzchar(trimws(as.character(disease))) &
    !tolower(trimws(as.character(disease))) %in%
    healthy_labels
])

# Require complete module scores

X <- as.matrix(patients[, ..module_cols])
storage.mode(X) <- "double"

patients <- patients[
  rowSums(is.finite(X)) == length(module_cols)
]

# Unaffected is the reference category

patients[
  , lesion_status := factor(
    lesion_status,
    levels = c("Unaffected", "Affected")
  )
]

cat("Patient tissue samples:", nrow(patients), "\n")
cat("GSE:", uniqueN(patients$series_id), "\n")


# ============================================================
# 4. IDENTIFY ELIGIBLE GSE x DISEASE GROUPS
# ============================================================

group_summary <- patients[
  , .(
    n_affected = sum(lesion_status == "Affected"),
    n_unaffected = sum(lesion_status == "Unaffected"),
    n_total = .N
  ),
  by = .(series_id, disease)
]

group_summary[
  , eligible := n_affected >= min_per_status &
    n_unaffected >= min_per_status
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
cat("Eligible samples:", nrow(patients_eligible), "\n")


# ============================================================
# 5. CHECK PAIRED SAMPLE STRUCTURE
# ============================================================

# A valid subject_id must identify the same patient
# across affected and unaffected tissue samples.

if (!"subject_id" %in% names(patients_eligible)) {
  patients_eligible[, subject_id := NA_character_]
}

patients_eligible[
  , subject_id := trimws(as.character(subject_id))
]

patients_eligible[
  subject_id %in% c("", "NA", "N/A", "unknown", "Unknown"),
  subject_id := NA_character_
]

pairing_summary <- patients_eligible[
  , {
    group_data <- .SD
    
    valid_ids <- group_data[!is.na(subject_id)]
    
    subject_counts <- valid_ids[
      , .(
        n_affected = sum(lesion_status == "Affected"),
        n_unaffected = sum(lesion_status == "Unaffected")
      ),
      by = subject_id
    ]
    
    n_paired <- sum(
      subject_counts$n_affected > 0 &
        subject_counts$n_unaffected > 0
    )
    
    # Strict pairing:
    # Every sample must have a subject ID;
    # each subject must have exactly one sample
    # from each tissue condition.
    
    fully_paired <- (
      nrow(group_data) > 0 &&
        all(!is.na(group_data$subject_id)) &&
        nrow(subject_counts) > 0 &&
        all(subject_counts$n_affected == 1) &&
        all(subject_counts$n_unaffected == 1)
    )
    
    .(
      n_samples = nrow(group_data),
      n_subjects_with_id = nrow(subject_counts),
      n_paired_subjects = n_paired,
      fully_paired = fully_paired,
      pairing_status = if (fully_paired) {
        "Fully paired"
      } else if (n_paired > 0) {
        "Partially paired / repeated samples"
      } else {
        "No confirmed pairing"
      }
    )
  },
  by = .(series_id, disease)
]

print(pairing_summary, nrows = Inf)

fwrite(
  pairing_summary,
  file.path(output_dir, "pairing_summary.tsv"),
  sep = "\t"
)


# ============================================================
# 6. LINEAR MODELS WITHIN GSE x DISEASE
# ============================================================

# Fully paired groups:
#   Paired analysis using within-subject differences.
#
# Groups without confirmed pairing:
#   Ordinary linear regression.
#
# Partially paired groups:
#   Flagged and not fitted automatically because
#   an unpaired model would ignore dependence.

results <- patients_eligible[
  , {
    group_data <- .SD
    
    pairing_info <- pairing_summary[
      series_id == .BY$series_id &
        disease == .BY$disease
    ]
    
    pairing_status <- pairing_info$pairing_status[1]
    
    if (pairing_status == "Partially paired / repeated samples") {
      
      NULL
      
    } else {
      
      rbindlist(lapply(module_cols, function(mod) {
        
        df <- data.table(
          score = as.numeric(group_data[[mod]]),
          lesion_status = group_data$lesion_status,
          subject_id = group_data$subject_id
        )
        
        n_affected <- sum(df$lesion_status == "Affected")
        n_unaffected <- sum(df$lesion_status == "Unaffected")
        
        mean_affected <- mean(
          df$score[df$lesion_status == "Affected"]
        )
        
        mean_unaffected <- mean(
          df$score[df$lesion_status == "Unaffected"]
        )
        
        if (pairing_status == "Fully paired") {
          
          # Paired samples: calculate affected-unaffected
          # differences within each patient.
          
          wide <- dcast(
            df,
            subject_id ~ lesion_status,
            value.var = "score"
          )
          
          wide[
            , difference := Affected - Unaffected
          ]
          
          fit <- lm(difference ~ 1, data = wide)
          
          coef_name <- "(Intercept)"
          model_type <- "Paired"
          
        } else {
          
          # Independent samples
          
          fit <- lm(
            score ~ lesion_status,
            data = df
          )
          
          coef_name <- "lesion_statusAffected"
          model_type <- "Unpaired"
        }
        
        coef_table <- summary(fit)$coefficients
        
        beta <- unname(
          coef_table[coef_name, "Estimate"]
        )
        
        se <- unname(
          coef_table[coef_name, "Std. Error"]
        )
        
        p <- unname(
          coef_table[coef_name, "Pr(>|t|)"]
        )
        
        ci <- confint(
          fit,
          coef_name,
          level = 0.95
        )
        
        data.table(
          module = mod,
          model_type = model_type,
          n_affected = n_affected,
          n_unaffected = n_unaffected,
          n_paired_subjects =
            pairing_info$n_paired_subjects[1],
          mean_affected = mean_affected,
          mean_unaffected = mean_unaffected,
          beta_Affected_vs_Unaffected = beta,
          se = se,
          CI_lower = ci[1],
          CI_upper = ci[2],
          pvalue = p
        )
        
      }))
      
    }
  },
  by = .(series_id, disease)
]


# ============================================================
# 7. MULTIPLE TESTING CORRECTION
# ============================================================

if (nrow(results) > 0) {
  
  results[
    , FDR := p.adjust(pvalue, method = "BH"),
    by = .(series_id, disease)
  ]
  
  results[
    , direction := fifelse(
      beta_Affected_vs_Unaffected > 0,
      "Higher in affected",
      fifelse(
        beta_Affected_vs_Unaffected < 0,
        "Higher in unaffected",
        "No difference"
      )
    )
  ]
  
  setorder(results, series_id, disease, module)
}


# ============================================================
# 8. SAVE RESULTS
# ============================================================

fwrite(
  results,
  file.path(
    output_dir,
    "lesion_module_associations_by_GSE_disease.tsv"
  ),
  sep = "\t"
)

saveRDS(
  results,
  file.path(
    output_dir,
    "lesion_module_associations_by_GSE_disease.rds"
  )
)


# ============================================================
# 9. SIGNIFICANT ASSOCIATIONS SUMMARY
# ============================================================

significant_results <- results[
  !is.na(FDR) & FDR < 0.05
]

significant_summary <- significant_results[
  , .(
    n_significant = .N,
    n_higher_affected =
      sum(beta_Affected_vs_Unaffected > 0),
    n_higher_unaffected =
      sum(beta_Affected_vs_Unaffected < 0)
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
  "n_higher_affected",
  "n_higher_unaffected"
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
  file.path(
    output_dir,
    "significant_associations_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 10. SUMMARY
# ============================================================

cat("\n====================================\n")
cat("LESION ANALYSIS - STEP 3 COMPLETED\n")
cat("====================================\n")

cat(
  "Eligible GSE x disease groups:",
  nrow(eligible_groups), "\n"
)

cat(
  "Eligible samples:",
  nrow(patients_eligible), "\n"
)

cat(
  "Fully paired groups:",
  sum(pairing_summary$fully_paired), "\n"
)

cat(
  "Groups excluded due to partial pairing:",
  sum(
    pairing_summary$pairing_status ==
      "Partially paired / repeated samples"
  ),
  "\n"
)

cat("Models fitted:", nrow(results), "\n")

cat(
  "Significant associations (FDR < 0.05):",
  nrow(significant_results), "\n"
)

cat("\nSignificant associations by module:\n")
print(significant_summary)

cat("\nResults saved in:\n", output_dir, "\n")
cat("====================================\n")