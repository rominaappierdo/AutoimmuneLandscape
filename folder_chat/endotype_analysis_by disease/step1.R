# ============================================================
# AutoimmuneLandscape
# ENDOTYPE ANALYSIS BY DATASET
# STEP 1 - Data preparation and quality control
# Whole blood - Patients only
# ============================================================

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/endotype_analysis_by disease")
library(data.table)

# ============================================================
# 1. DIRECTORIES
# ============================================================

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/"

# Previous analysis: READ ONLY

input_file <- file.path(
  base_dir,
  "endotype_analysis",
  "all_samples_with_metadata.rds"
)

# New analysis directory

output_dir <- file.path(
  base_dir,
  "endotype_analysis_by disease"
)

step_dir <- file.path(
  output_dir,
  "step1_data_preparation"
)

gse_dir <- file.path(
  step_dir,
  "datasets"
)

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)

module_cols <- paste0("mod", 1:15)

min_samples <- 30


# ============================================================
# 2. LOAD DATA
# ============================================================

data_all <- as.data.table(readRDS(input_file))

cat("Original samples:", nrow(data_all), "\n")
cat("Original GSEs:", uniqueN(data_all$series_id), "\n")

stopifnot(all(c(
  "series_id",
  "geo_accession",
  "disease",
  "condition",
  module_cols
) %in% names(data_all)))

cat("\nOriginal disease distribution:\n")

print(
  data_all[, .N, by = disease][order(-N)],
  nrows = Inf
)

cat("\nOriginal condition distribution:\n")

print(data_all[, .N, by = condition])


# ============================================================
# 3. CHECK TISSUE
# ============================================================

# The input file was generated in the previous
# whole-blood analysis.

if ("sample_source" %in% names(data_all)) {
  
  cat("\nSample sources:\n")
  
  print(
    data_all[, .N, by = sample_source][order(-N)],
    nrows = Inf
  )
  
  stopifnot(
    all(
      tolower(trimws(as.character(data_all$sample_source))) ==
        "whole blood"
    )
  )
  
}


# ============================================================
# 4. FILTER PATIENTS
# ============================================================

# Use the same disease exclusion criteria
# as in the previous endotype analysis.

exclude_diseases <- c(
  "staphylococcus infection",
  "streptococcus infection",
  "down syndrome",
  "prediabetic",
  "fever without source",
  "irritable bowel syndrome",
  "diverticular disease",
  "other"
)

patients <- copy(data_all[
  condition == "Disease" &
    !is.na(disease) &
    nzchar(trimws(as.character(disease))) &
    !disease %in% exclude_diseases &
    !is.na(series_id) &
    nzchar(trimws(as.character(series_id)))
])

cat("\nDisease samples:", nrow(patients), "\n")

cat(
  "GSEs after patient selection:",
  uniqueN(patients$series_id),
  "\n"
)


# ============================================================
# 5. REMOVE INCOMPLETE MODULE SCORE PROFILES
# ============================================================

X <- as.matrix(patients[, ..module_cols])

storage.mode(X) <- "double"

complete_idx <- rowSums(
  is.finite(X)
) == length(module_cols)

excluded <- patients[!complete_idx]

patients <- patients[complete_idx]

cat("\nSamples with incomplete module scores:",
    nrow(excluded), "\n")

cat("Samples with complete module scores:",
    nrow(patients), "\n")

cat("GSEs retained:",
    uniqueN(patients$series_id), "\n")

fwrite(
  excluded[, .N, by = .(series_id, disease)],
  file.path(step_dir, "excluded_incomplete_samples.tsv"),
  sep = "\t"
)

if (nrow(patients) == 0) {
  stop("No patients retained after filtering.")
}


# ============================================================
# 6. CHECK DUPLICATED SAMPLES
# ============================================================

duplicate_samples <- patients[
  , .N,
  by = .(series_id, geo_accession)
][N > 1]

cat("\nDuplicated GSE-GSM combinations:",
    nrow(duplicate_samples), "\n")

fwrite(
  duplicate_samples,
  file.path(step_dir, "duplicate_samples.tsv"),
  sep = "\t"
)

if (nrow(duplicate_samples) > 0) {
  stop("Duplicated samples detected. Inspect before proceeding.")
}


# ============================================================
# 7. SAMPLE DISTRIBUTION BY DISEASE
# ============================================================

disease_summary <- patients[
  , .(
    n_samples = .N,
    n_datasets = uniqueN(series_id)
  ),
  by = disease
][order(-n_samples)]

cat("\nDisease summary:\n")

print(disease_summary, nrows = Inf)

fwrite(
  disease_summary,
  file.path(step_dir, "disease_summary.tsv"),
  sep = "\t"
)


# ============================================================
# 8. SAMPLE DISTRIBUTION BY GSE AND DISEASE
# ============================================================

gse_disease_summary <- patients[
  , .(
    n_samples = .N
  ),
  by = .(series_id, disease)
]

setorder(
  gse_disease_summary,
  series_id,
  -n_samples
)

cat("\nGSE x disease summary:\n")

print(gse_disease_summary, nrows = Inf)

fwrite(
  gse_disease_summary,
  file.path(step_dir, "GSE_disease_summary.tsv"),
  sep = "\t"
)


# ============================================================
# 9. DATASET SUMMARY
# ============================================================

gse_summary <- patients[
  , .(
    n_samples = .N,
    n_diseases = uniqueN(disease),
    diseases = paste(
      sort(unique(disease)),
      collapse = " | "
    )
  ),
  by = series_id
]

setorder(gse_summary, -n_samples)

cat("\nDataset summary:\n")

print(gse_summary, nrows = Inf)


# ============================================================
# 10. CHECK MODULE VARIANCE WITHIN EACH GSE
# ============================================================

module_variance <- patients[
  , {
    
    rbindlist(lapply(module_cols, function(mod) {
      
      scores <- as.numeric(.SD[[mod]])
      
      data.table(
        module = mod,
        mean = mean(scores),
        sd = sd(scores),
        variance = var(scores),
        min_score = min(scores),
        max_score = max(scores)
      )
      
    }))
    
  },
  by = series_id,
  .SDcols = module_cols
]

fwrite(
  module_variance,
  file.path(step_dir, "module_variance_by_GSE.tsv"),
  sep = "\t"
)

# Identify constant or undefined module variance

variance_summary <- module_variance[
  , .(
    n_zero_variance = sum(
      is.finite(variance) & variance == 0
    ),
    n_invalid_variance = sum(
      !is.finite(variance)
    ),
    all_modules_variable = all(
      is.finite(variance) & variance > 0
    )
  ),
  by = series_id
]

gse_summary <- merge(
  gse_summary,
  variance_summary,
  by = "series_id",
  all.x = TRUE
)


# ============================================================
# 11. CHECK REPEATED PATIENTS
# ============================================================

if ("subject_id" %in% names(patients)) {
  
  subject_data <- copy(patients)
  
  subject_data[
    , subject_id := trimws(as.character(subject_id))
  ]
  
  subject_data <- subject_data[
    !is.na(subject_id) &
      !tolower(subject_id) %in%
      c("", "na", "n/a", "unknown")
  ]
  
  subject_counts <- subject_data[
    , .(
      n_samples = .N
    ),
    by = .(series_id, subject_id)
  ]
  
  repeated_subjects <- subject_counts[
    n_samples > 1
  ]
  
  cat("\nSubjects with multiple samples:",
      nrow(repeated_subjects), "\n")
  
  fwrite(
    repeated_subjects,
    file.path(step_dir, "repeated_subjects.tsv"),
    sep = "\t"
  )
  
  subject_summary <- subject_counts[
    , .(
      n_known_subjects = .N,
      n_samples_with_subject_id = sum(n_samples),
      n_repeated_subjects = sum(n_samples > 1)
    ),
    by = series_id
  ]
  
  gse_summary <- merge(
    gse_summary,
    subject_summary,
    by = "series_id",
    all.x = TRUE
  )
  
  # Do not automatically remove repeated samples.
  # Their origin must be checked before clustering.
  
}


# ============================================================
# 12. IDENTIFY ELIGIBLE DATASETS
# ============================================================

# Eligibility criteria:
# - At least 30 samples
# - Complete module score profiles
# - All 15 modules have non-zero variance

gse_summary[
  , eligible :=
    n_samples >= min_samples &
    all_modules_variable
]

gse_summary[
  , exclusion_reason := fifelse(
    n_samples < min_samples,
    "Insufficient samples",
    fifelse(
      !all_modules_variable,
      "Constant or invalid module variance",
      "Eligible"
    )
  )
]

setorder(gse_summary, -n_samples)

eligible_gse <- gse_summary[
  eligible == TRUE,
  series_id
]

excluded_gse <- gse_summary[
  eligible == FALSE
]

cat("\nEligible datasets:", length(eligible_gse), "\n")
cat("Excluded datasets:", nrow(excluded_gse), "\n")

cat("\nEligible dataset summary:\n")

print(
  gse_summary[eligible == TRUE],
  nrows = Inf
)

fwrite(
  gse_summary,
  file.path(step_dir, "GSE_eligibility.tsv"),
  sep = "\t"
)

fwrite(
  excluded_gse,
  file.path(step_dir, "excluded_GSEs.tsv"),
  sep = "\t"
)


# ============================================================
# 13. PREPARE ELIGIBLE DATA
# ============================================================

patients_eligible <- copy(patients[
  series_id %in% eligible_gse
])

cat("\nSamples in eligible datasets:",
    nrow(patients_eligible), "\n")

cat("Diseases in eligible datasets:",
    uniqueN(patients_eligible$disease), "\n")


# ============================================================
# 14. SAVE ALL FILTERED DATA
# ============================================================

saveRDS(
  patients,
  file.path(step_dir, "all_patients_complete.rds")
)

saveRDS(
  patients_eligible,
  file.path(step_dir, "eligible_patients.rds")
)

saveRDS(
  eligible_gse,
  file.path(step_dir, "eligible_GSEs.rds")
)

saveRDS(
  gse_summary,
  file.path(step_dir, "GSE_summary.rds")
)


# ============================================================
# 15. SAVE ONE RDS FILE PER ELIGIBLE GSE
# ============================================================

for (gse in eligible_gse) {
  
  gse_data <- copy(
    patients_eligible[series_id == gse]
  )
  
  saveRDS(
    gse_data,
    file.path(gse_dir, paste0(gse, ".rds"))
  )
  
  cat(
    "Saved:",
    gse,
    "| Samples:", nrow(gse_data),
    "| Diseases:", uniqueN(gse_data$disease),
    "\n"
  )
  
}


# ============================================================
# 16. FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("ENDOTYPE ANALYSIS BY DATASET\n")
cat("STEP 1 COMPLETED\n")
cat("========================================\n")

cat("Original samples:", nrow(data_all), "\n")

cat("Disease samples with complete scores:",
    nrow(patients), "\n")

cat("Total GSEs:",
    uniqueN(patients$series_id), "\n")

cat("Minimum samples per GSE:",
    min_samples, "\n")

cat("Eligible GSEs:",
    length(eligible_gse), "\n")

cat("Excluded GSEs:",
    nrow(excluded_gse), "\n")

cat("Samples in eligible GSEs:",
    nrow(patients_eligible), "\n")

cat("Diseases in eligible GSEs:",
    uniqueN(patients_eligible$disease), "\n")

cat("\nOutput directory:\n", step_dir, "\n")

cat("========================================\n")