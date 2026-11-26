
library(data.table)

# ============================================================
# 0. SETUP
# ============================================================

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat_tissue/")

source("/mnt/data3/appierdo/PhD-project/utils_functions/function_df_all_metadata.R")
source("/mnt/data3/appierdo/PhD-project/utils_functions/function_compute_all_features_for_module.R")
source("/mnt/data3/appierdo/AutoimmuneLandscape/folder_chat/compute_features2.R")

# ============================================================
# 1. LOAD THE 15 PREDEFINED MODULES
# ============================================================

# Load the original modules without recalculating them
modules_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files"

modules <- setNames(
  lapply(1:15, function(i) {
    readLines(file.path(
      modules_dir,
      paste0("mod", i),
      paste0("mod", i, ".txt")
    ))
  }),
  paste0("mod", 1:15)
)

# Check number of features in each module
sapply(modules, length)


# ============================================================
# 2. LOAD AND FILTER METADATA
# ============================================================

# Load complete metadata
meta <- as.data.table(complete_metadata())

# Exclude samples associated with specific cell types
sub_meta <- copy(meta[cell_type == ""])

# Define accepted baseline / untreated timepoints
valid_timepoints <- c(
  "", " ", "baseline", "pre", "NA",
  "rep1", "unstimulated", "untreated",
  "visit 1", "T0", "t0", "biopsy 1",
  "GPL96", "visit day 0", "week 0",
  "at-diagnosis (for progressor)",
  "day 0", "at diagnosis",
  "untreated, rep 1", "baseline, rep1",
  "early"
)

# Retain samples with accepted timepoints
sub_meta <- sub_meta[
  timepoint %in% valid_timepoints | is.na(timepoint)
]


# ============================================================
# 3. IDENTIFY AVAILABLE EXPRESSION DATASETS
# ============================================================

data_dir <- "/mnt/data3/appierdo/AutoimmuneLandscape/Data_Collection/downloaded_datasets"

# Find microarray and RNA-seq expression files recursively
# Microarray: GSEXXXXX_gs.tsv or GSEXXXXX_GPLXXXX_gs.tsv
# RNA-seq:    GSEXXXXX.rds
expr_files <- list.files(
  data_dir,
  pattern = "(\\.rds$|_gs\\.tsv$)",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

# Extract GSE identifiers from filenames
file_gse <- toupper(
  sub("^(GSE[0-9]+).*", "\\1", basename(expr_files))
)

# Keep only files with valid GSE identifiers
valid_files <- grepl("^GSE[0-9]+$", file_gse)

expr_files <- expr_files[valid_files]
file_gse <- file_gse[valid_files]

# Identify unique GSE with available expression data
gse_available <- sort(unique(file_gse))

cat("GSE with expression files:", length(gse_available), "\n")


# ============================================================
# 4. MATCH METADATA WITH AVAILABLE DATASETS
# ============================================================

# Identify GSE present in both metadata and expression files
gse_with_metadata <- intersect(
  gse_available,
  unique(sub_meta$series_id)
)

# Retain only datasets with expression files
sub_meta <- sub_meta[
  series_id %in% gse_with_metadata
]

cat("GSE with metadata:", uniqueN(sub_meta$series_id), "\n")
cat("Samples with metadata:", nrow(sub_meta), "\n")


# ============================================================
# 5. SELECT TISSUE SAMPLES
# ============================================================

sub_meta <- sub_meta[sample_type == "tissue"]

cat("Tissue samples:", nrow(sub_meta), "\n")
cat("Tissue GSE:", uniqueN(sub_meta$series_id), "\n")


# ============================================================
# 6. REMOVE EXACT DUPLICATE METADATA ROWS
# ============================================================

# Remove rows identical across all columns
sub_meta <- unique(sub_meta)


# ============================================================
# 7. REMOVE DUPLICATES DIFFERING ONLY IN GSM
# ============================================================

# Some samples have identical metadata but different GSM IDs.
# Randomly retain one GSM for each identical metadata profile.

set.seed(123)

# All columns except geo_accession
cols <- setdiff(names(sub_meta), "geo_accession")

# Randomize row order
sub_meta <- sub_meta[sample(.N)]

# Keep one GSM per identical metadata profile
sub_meta <- sub_meta[
  !duplicated(sub_meta[, ..cols])
]

cat("Samples after metadata deduplication:",
    nrow(sub_meta), "\n")


# ============================================================
# 8. KEEP ONE SAMPLE PER PATIENT
# ============================================================

# For longitudinal samples, retain the sample with the
# youngest recorded age.
#
# Exception: GSE99967 has non-unique subject identifiers,
# so all its GSM IDs are preserved.

meta_99967 <- copy(sub_meta[series_id == "GSE99967"])
meta_other <- copy(sub_meta[series_id != "GSE99967"])

# Convert age to numeric years
# Handles values such as "70 months", "31.4", "45"
meta_other[, age_num := suppressWarnings(
  as.numeric(
    trimws(sub(" months", "", age, fixed = TRUE))
  )
)]

# Convert ages expressed in months to years
meta_other[
  grepl("months", age, fixed = TRUE),
  age_num := age_num / 12
]

# Sort by dataset, patient, age, and GSM
# Missing or non-numeric ages are placed last
setorder(
  meta_other,
  series_id,
  subject_id,
  age_num,
  geo_accession,
  na.last = TRUE
)

# Keep the youngest sample for each patient
# Preserve all samples with empty subject_id
meta_other <- meta_other[
  subject_id == "" |
    !duplicated(meta_other[, .(series_id, subject_id)])
]

# Remove temporary age column
meta_other[, age_num := NULL]

# Recombine all samples
sub_meta <- rbindlist(
  list(meta_other, meta_99967),
  use.names = TRUE
)


# ============================================================
# 9. CHECK FINAL METADATA
# ============================================================

cat("\nFINAL METADATA\n")
cat("Samples:", nrow(sub_meta), "\n")
cat("GSE:", uniqueN(sub_meta$series_id), "\n")
cat("Disease annotations:", uniqueN(sub_meta$disease), "\n")

# Verify that no patient has multiple GSMs,
# excluding GSE99967
remaining_duplicates <- sub_meta[
  series_id != "GSE99967" & subject_id != "",
  .(n_GSM = uniqueN(geo_accession)),
  by = .(series_id, subject_id)
][n_GSM > 1]

cat("Patients with multiple GSM:",
    nrow(remaining_duplicates), "\n")

# Disease distribution
sort(table(sub_meta$disease), decreasing = TRUE)

# Dataset distribution
sort(table(sub_meta$series_id), decreasing = TRUE)


# ============================================================
# 10. IDENTIFY EXPRESSION FILES FOR SELECTED GSE
# ============================================================

# Match files using extracted GSE identifiers rather than
# filename prefixes. This handles both:
#   GSEXXXXX.rds
#   GSEXXXXX_gs.tsv
#   GSEXXXXX_GPLXXXX_gs.tsv

selected_gse <- unique(sub_meta$series_id)
expr_files_selected <- expr_files[file_gse %in% selected_gse]

# Store the GSE corresponding to each selected file
expr_file_gse <- file_gse[file_gse %in% selected_gse]

# Check number of files and represented datasets
cat("\nEXPRESSION FILES\n")
cat("Expression files:", length(expr_files_selected), "\n")
cat("GSE represented:", uniqueN(expr_file_gse), "\n")

# Check whether any GSE is missing an expression file
gse_missing <- setdiff(
  selected_gse,
  unique(expr_file_gse)
)

cat("GSE missing expression files:",
    length(gse_missing), "\n")

print(gse_missing)

# Identify GSE represented by multiple expression files
# (different GPL platforms or copies in disease directories)
multiple_files <- sort(
  table(expr_file_gse),
  decreasing = TRUE
)

print(multiple_files[multiple_files > 1])


# ============================================================
# 11. CREATE EXPRESSION FILE INDEX
# ============================================================

# Keep the complete file path associated with each GSE.
# This avoids losing information about different platforms
# or copies of the same dataset in different directories.

expr_file_index <- data.table(
  series_id = expr_file_gse,
  file_path = expr_files_selected,
  technology = ifelse(
    grepl("\\.rds$", expr_files_selected, ignore.case = TRUE),
    "RNAseq",
    "microarray"
  )
)

setorder(expr_file_index, series_id, file_path)

# Inspect file index
print(expr_file_index)

# Final consistency check
stopifnot(
  setequal(unique(expr_file_index$series_id),
           unique(sub_meta$series_id))
)


# ============================================================
# 12. INSPECT EXPRESSION FILE STRUCTURE
# ============================================================

# Select one microarray file
micro_file <- expr_file_index[
  technology == "microarray", file_path[1]
]

# Select one RNA-seq file
rna_file <- expr_file_index[
  technology == "RNAseq", file_path[1]
]

# MICROARRAY
cat("\nMICROARRAY FILE:\n", micro_file, "\n")

micro_test <- fread(micro_file, nrows = 5)

print(dim(micro_test))
print(names(micro_test)[1:min(10, ncol(micro_test))])
print(micro_test[, 1:min(5, ncol(micro_test)), with = FALSE])


# RNA-SEQ
cat("\nRNA-SEQ FILE:\n", rna_file, "\n")

rna_test <- readRDS(rna_file)

print(class(rna_test))
print(names(rna_test))

if ("expr" %in% names(rna_test)) {
  print(dim(rna_test$expr))
  print(head(colnames(rna_test$expr), 10))
  print(head(rownames(rna_test$expr), 10))
}



#######################################################################################################################################
# ============================================================
# 14. SUBSET EXPRESSION DATA TO SELECTED GSM
# ============================================================




# ============================================================
# 15. COMPUTE AND SAVE MODULE FEATURES FOR EACH DATASET
# ============================================================

# Output directory
output_dir <- file.path(getwd(), "module_scores")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# All features belonging to the 15 predefined modules
all_module_features <- unique(unlist(modules, use.names = FALSE))

# Process each expression file separately
for (i in seq_len(nrow(expr_file_index))) {
  
  gse <- expr_file_index$series_id[i]
  
  if (gse %in% c("GSE77387")) { #no corrispondenza gene symbols
    next
  }
  
  f <- expr_file_index$file_path[i]
  tech <- expr_file_index$technology[i]
  
  cat("\nProcessing", i, "/", nrow(expr_file_index), ":", gse, "\n")
  
  # Select GSM from filtered metadata
  gsm_selected <- unique(
    sub_meta[series_id == gse, geo_accession]
  )
  
  # Load expression matrix
  if (tech == "microarray") {
    
    expr <- fread(f)
    expr <- as.data.frame(expr)
    
    rownames(expr) <- expr[[1]]
    expr <- as.matrix(expr[, -1, drop = FALSE])
    
  } else {
    
    obj <- readRDS(f)
    expr <- obj$expr
    rm(obj)
    
  }
  
  # Keep only GSM present in sub_meta
  gsm_keep <- intersect(colnames(expr), gsm_selected)
  
  if (length(gsm_keep) < 2) {
    cat("Skipping: fewer than 2 matching GSM\n")
    rm(expr)
    next
  }
  
  expr <- expr[, gsm_keep, drop = FALSE]
  
  # Compute features for all 15 modules
  feature_scores <- compute_features2(
    expr = expr,
    module = all_module_features,
    deconvolution = TRUE,
    ncore = 8
  )
  
  # Create unique filename from source path
  # Avoid overwriting files belonging to different
  # platforms or disease directories
  relative_path <- substring(f, nchar(data_dir) + 2)
  file_id <- gsub("[/\\\\]", "__", relative_path)
  file_id <- sub("\\.(rds|tsv)$", "", file_id, ignore.case = TRUE)
  
  # Save GSM x features matrix
  saveRDS(
    feature_scores,
    file.path(output_dir, paste0(file_id, "_features.rds"))
  )
  
  cat(
    "Saved:",
    nrow(feature_scores), "GSM x",
    ncol(feature_scores), "features\n"
  )
  
  # Free memory
  rm(expr, feature_scores)
  gc(verbose = FALSE)
}


print("DONE")
