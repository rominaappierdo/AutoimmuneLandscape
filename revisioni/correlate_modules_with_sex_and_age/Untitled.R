
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/revisioni/correlate_modules_with_sex_and_age")

source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")

# ============================================================
# 1. METADATA
# ============================================================

meta <- complete_metadata()

sub_meta <- meta %>%
  filter(
    gender %in% c("M", "F"),
    !is.na(age),
    age != "",
    !grepl("month|M", age),
    is.na(cell_type) | cell_type == "",
    is.na(timepoint) | timepoint %in% c(
      "baseline", "pre", "NA", "rep1",
      "unstimulated", "untreated", "visit 1",
      "T0", "t0", "biopsy 1", "GPL96",
      "visit day 0", "week 0",
      "at-diagnosis (for progressor)",
      "day 0", "at diagnosis",
      "untreated, rep 1", "baseline, rep1", "early"
    )
  ) %>%
  mutate(age = suppressWarnings(as.numeric(age))) %>%
  filter(!is.na(age))

# ============================================================
# 2. FIND ALL GSE FOLDERS
# ============================================================

base_dir <- file.path(getwd(), "modules_df")

gse_dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

gse_dirs <- gse_dirs[
  grepl("^GSE[0-9]+$", basename(gse_dirs))
]

cat("Found", length(gse_dirs), "GSE folders\n")

# ============================================================
# 3. CORRELATIONS
# ============================================================

cor_results <- data.frame()

for (gse_dir in gse_dirs) {
  
  gse <- basename(gse_dir)
  cat("\nProcessing:", gse, "\n")
  
  files <- list.files(
    gse_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
  
  for (file in files) {
    
    module_name <- basename(file)
    module_name <- sub("\\.tsv$", "", module_name)
    module_name <- sub("^GSE[0-9]+_", "", module_name)
    
    module_table <- read.delim(
      file,
      row.names = 1,
      check.names = FALSE
    )
    
    # Match samples between metadata and module table
    gse_meta <- sub_meta %>%
      filter(
        series_id == gse,
        geo_accession %in% colnames(module_table)
      )
    
    # Ensure unique sample IDs
    if (anyDuplicated(gse_meta$geo_accession)) {
      warning("Duplicated samples in metadata: ", gse)
      next
    }
    
    if (nrow(gse_meta) < 3) next
    
    # Match expression columns to metadata order
    module_table <- module_table[
      , gse_meta$geo_accession, drop = FALSE
    ]
    
    stopifnot(
      identical(
        colnames(module_table),
        gse_meta$geo_accession
      )
    )
    
    # Module score
    gse_meta$module_score <- colSums(
      module_table,
      na.rm = FALSE
    )
    
    # Remove missing values
    valid <- complete.cases(
      gse_meta[, c("age", "module_score")]
    )
    
    x <- gse_meta$age[valid]
    y <- gse_meta$module_score[valid]
    
    if (
      length(x) < 3 ||
      sd(x) == 0 ||
      sd(y) == 0
    ) next
    
    # Pearson correlation
    test <- cor.test(x, y, method = "pearson")
    
    cor_results <- rbind(
      cor_results,
      data.frame(
        GSE = gse,
        Module = module_name,
        correlation = unname(test$estimate),
        pvalue = test$p.value,
        n = length(x)
      )
    )
  }
}

# ============================================================
# 4. MULTIPLE TESTING CORRECTION
# ============================================================

if (nrow(cor_results) == 0) {
  stop("No valid correlations found.")
}

cor_results <- cor_results %>%
  mutate(FDR = p.adjust(pvalue, method = "BH"))

write.csv(
  cor_results,
  "modules_age_correlations.csv",
  row.names = FALSE
)

# ============================================================
# 5. BUILD MATRICES
# ============================================================

cor_mat <- cor_results %>%
  select(Module, GSE, correlation) %>%
  pivot_wider(
    names_from = GSE,
    values_from = correlation
  ) %>%
  column_to_rownames("Module") %>%
  as.matrix()

fdr_mat <- cor_results %>%
  select(Module, GSE, FDR) %>%
  pivot_wider(
    names_from = GSE,
    values_from = FDR
  ) %>%
  column_to_rownames("Module") %>%
  as.matrix()

# Ensure identical dimensions and order
fdr_mat <- fdr_mat[
  rownames(cor_mat),
  colnames(cor_mat),
  drop = FALSE
]

# Order GSEs
cor_mat <- cor_mat[
  , order(colnames(cor_mat)), drop = FALSE
]

fdr_mat <- fdr_mat[
  rownames(cor_mat),
  colnames(cor_mat),
  drop = FALSE
]

# ============================================================
# 6. COMPLEXHEATMAP
# ============================================================

ht <- Heatmap(
  cor_mat,
  
  name = "Pearson r",
  
  col = colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  ),
  
  na_col = "grey90",
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    r <- cor_mat[i, j]
    fdr <- fdr_mat[i, j]
    
    if (!is.na(r)) {
      
      label <- sprintf("%.2f", r)
      
      if (!is.na(fdr) && fdr < 0.05) {
        label <- paste0(label, "*")
      }
      
      grid.text(
        label,
        x, y,
        gp = gpar(
          fontsize = 8,
          fontface = "bold",
          col = ifelse(abs(r) > 0.7, "white", "black")
        )
      )
    }
  },
  
  column_title = "Correlation between module scores and age",
  
  heatmap_legend_param = list(
    title = "Pearson r",
    at = c(-1, -0.5, 0, 0.5, 1)
  )
)

# ============================================================
# 7. SAVE
# ============================================================

pdf(
  "modules_age_correlation_ComplexHeatmap.pdf",
  width = max(10, ncol(cor_mat) * 0.8 + 5),
  height = 9
)

draw(ht)
dev.off()

draw(ht)

cat("\nAnalysis completed!\n")
cat("Datasets:", ncol(cor_mat), "\n")
cat("Modules:", nrow(cor_mat), "\n")
cat("Significant correlations (FDR < 0.05):",
    sum(cor_results$FDR < 0.05, na.rm = TRUE), "\n")
