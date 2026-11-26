
# ============================================================
# AutoimmuneLandscape
# STEP 2 - COMPUTE MODULE SCORES FROM FEATURE MATRICES
# ============================================================

library(data.table)

# 1. DIRECTORIES ------------------------------------------------

base_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat"

modules_dir <- paste0(
  "/mnt/data3/appierdo/AutoimmuneLandscape/",
  "Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files"
)

input_dir <- file.path(base_dir, "module_scores")
output_dir <- file.path(base_dir, "module_scores_final")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# 2. LOAD THE ORIGINAL 15 MODULES -------------------------------

modules <- setNames(
  lapply(1:15, function(i) {
    readLines(
      file.path(modules_dir, paste0("mod", i),
                paste0("mod", i, ".txt"))
    )
  }),
  paste0("mod", 1:15)
)

modules <- lapply(modules, function(x) {
  unique(trimws(x[nzchar(trimws(x))]))
})

print(sapply(modules, length))


# 3. LOAD FEATURE FILE INDEX ------------------------------------

feature_files <- list.files(
  input_dir,
  pattern = "_features\\.rds$",
  full.names = TRUE
)

if (length(feature_files) == 0)
  stop("No feature files found!")

cat("Found", length(feature_files), "feature files\n")


# 4. FUNCTION TO COMPUTE MODULE SCORES --------------------------

compute_module_scores <- function(features, modules,
                                  min_coverage = 0.5) {
  
  features <- as.matrix(features)
  storage.mode(features) <- "numeric"
  
  if (is.null(rownames(features)) ||
      is.null(colnames(features))) {
    stop("Features must have GSM rows and feature columns")
  }
  
  if (anyDuplicated(colnames(features))) {
    stop("Duplicated feature names")
  }
  
  # Standardize each feature across GSM
  features_z <- scale(features)
  features_z[!is.finite(features_z)] <- NA_real_
  
  scores <- matrix(
    NA_real_,
    nrow = nrow(features_z),
    ncol = length(modules),
    dimnames = list(rownames(features_z), names(modules))
  )
  
  qc <- data.table(
    module = names(modules),
    n_expected = lengths(modules),
    n_available = 0L,
    n_valid_features = 0L,
    coverage = 0,
    n_scored_GSM = 0L
  )
  
  for (mod in names(modules)) {
    
    expected <- modules[[mod]]
    
    available <- intersect(expected, colnames(features_z))
    
    # Exclude features without any valid Z-score
    valid <- available[
      colSums(is.finite(
        features_z[, available, drop = FALSE]
      )) > 0
    ]
    
    n_expected <- length(expected)
    n_available <- length(available)
    n_valid <- length(valid)
    
    coverage <- if (n_expected > 0) {
      n_valid / n_expected
    } else {
      0
    }
    
    qc[module == mod, `:=`(
      n_available = as.integer(length(available)),
      n_valid_features = as.integer(length(valid)),
      coverage = as.numeric(length(valid) / length(expected))
    )]
    
    if (n_valid == 0 || coverage < min_coverage)
      next
    
    x <- features_z[, valid, drop = FALSE]
    
    # Require at least 50% of expected module features
    # to be valid in each individual GSM
    n_per_gsm <- rowSums(is.finite(x))
    min_required <- ceiling(min_coverage * n_expected)
    
    score <- rowMeans(x, na.rm = TRUE)
    score[n_per_gsm < min_required] <- NA_real_
    score[is.nan(score)] <- NA_real_
    
    scores[, mod] <- score
    
    qc[module == mod,
       n_scored_GSM := sum(!is.na(score))]
  }
  
  list(scores = scores, qc = qc)
}


# 5. PROCESS ALL DATASETS ---------------------------------------

qc_list <- vector("list", length(feature_files))

for (i in seq_along(feature_files)) {
  
  f <- feature_files[i]
  
  cat(
    "\n[", i, "/", length(feature_files), "]",
    basename(f), "\n"
  )
  
  # GSM x features
  features <- readRDS(f)
  
  result <- compute_module_scores(
    features = features,
    modules = modules,
    min_coverage = 0.5
  )
  
  # GSM x 15 modules
  scores <- result$scores
  
  # Save scores
  output_name <- sub(
    "_features\\.rds$",
    "_module_scores.rds",
    basename(f)
  )
  
  saveRDS(
    scores,
    file.path(output_dir, output_name)
  )
  
  # Collect QC
  qc <- copy(result$qc)
  qc[, dataset_file := basename(f)]
  qc[, n_GSM := nrow(scores)]
  
  qc_list[[i]] <- qc
  
  cat(
    "Saved:", nrow(scores), "GSM x",
    ncol(scores), "modules\n"
  )
  
  rm(features, result, scores, qc)
  gc(verbose = FALSE)
}


# 6. SAVE QUALITY CONTROL ---------------------------------------

qc_all <- rbindlist(qc_list, fill = TRUE)

fwrite(
  qc_all,
  file.path(output_dir, "module_score_coverage_QC.tsv"),
  sep = "\t"
)


# 7. SUMMARY ----------------------------------------------------

cat("\nCompleted!\n")
cat("Datasets processed:", length(feature_files), "\n")
cat("Output directory:", output_dir, "\n")

print(
  qc_all[, .(
    mean_coverage = round(mean(coverage), 3),
    datasets_with_scores = sum(n_scored_GSM > 0)
  ), by = module]
)


