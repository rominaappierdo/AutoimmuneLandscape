
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
      "untreated, rep 1", "baseline, rep1", "early", "", " "
    )
  ) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),
    Group = case_when(
      grepl("health|control", disease, ignore.case = TRUE) ~ "Healthy",
      !is.na(disease) & trimws(disease) != "" ~ "Disease",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age), !is.na(Group))

# ============================================================
# 2. FIND GSE DIRECTORIES
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
# 3. CORRELATIONS BY GROUP
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
    
    gse_meta <- sub_meta %>%
      filter(
        series_id == gse,
        geo_accession %in% colnames(module_table)
      )
    
    if (nrow(gse_meta) == 0) next
    
    if (anyDuplicated(gse_meta$geo_accession)) {
      warning("Duplicated samples in ", gse, ": ", module_name)
      next
    }
    
    # Match sample order
    module_table <- module_table[
      , gse_meta$geo_accession, drop = FALSE
    ]
    
    stopifnot(
      identical(
        colnames(module_table),
        gse_meta$geo_accession
      )
    )
    
    # Compute module scores
    gse_meta$module_score <- colSums(module_table)
    
    # Separate Healthy and Disease
    for (group in c("Healthy", "Disease")) {
      
      group_meta <- gse_meta %>%
        filter(
          Group == group,
          !is.na(age),
          !is.na(module_score)
        )
      
      if (nrow(group_meta) < 3) next
      
      x <- group_meta$age
      y <- group_meta$module_score
      
      if (sd(x) == 0 || sd(y) == 0) next
      
      test <- cor.test(
        x, y,
        method = "pearson"
      )
      
      cor_results <- rbind(
        cor_results,
        data.frame(
          GSE = gse,
          Module = module_name,
          Group = group,
          correlation = unname(test$estimate),
          pvalue = test$p.value,
          n = length(x)
        )
      )
    }
  }
}

if (nrow(cor_results) == 0) {
  stop("No valid correlations found.")
}

# ============================================================
# 4. FDR CORRECTION
# ============================================================

# Adjust separately within Healthy and Disease groups
cor_results <- cor_results %>%
  group_by(Group) %>%
  mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
  ungroup()

write.csv(
  cor_results,
  "modules_age_correlations_Healthy_vs_Disease.csv",
  row.names = FALSE
)

# ============================================================
# 5. COMPLEXHEATMAP FUNCTION
# ============================================================

plot_correlation_heatmap <- function(results, group) {
  
  df <- results %>%
    filter(Group == group)
  
  if (nrow(df) == 0) {
    message("No correlations for: ", group)
    return(invisible(NULL))
  }
  
  # Correlation matrix
  cor_mat <- df %>%
    select(Module, GSE, correlation) %>%
    pivot_wider(
      names_from = GSE,
      values_from = correlation
    ) %>%
    column_to_rownames("Module") %>%
    as.matrix()
  
  # FDR matrix
  fdr_mat <- df %>%
    select(Module, GSE, FDR) %>%
    pivot_wider(
      names_from = GSE,
      values_from = FDR
    ) %>%
    column_to_rownames("Module") %>%
    as.matrix()
  
  # Match dimensions and order
  cor_mat <- cor_mat[
    order(rownames(cor_mat)),
    order(colnames(cor_mat)),
    drop = FALSE
  ]
  
  fdr_mat <- fdr_mat[
    rownames(cor_mat),
    colnames(cor_mat),
    drop = FALSE
  ]
  
  ht <- Heatmap(
    cor_mat,
    
    name = "Pearson r",
    
    col = colorRamp2(
      c(-1, 0, 1),
      c("blue", "white", "red")
    ),
    
    na_col = "grey90",
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
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
          label, x, y,
          gp = gpar(
            fontsize = 8,
            fontface = "bold",
            col = ifelse(abs(r) > 0.7, "white", "black")
          )
        )
      }
    },
    
    column_title = paste0(
      group,
      ": Correlation between module scores and age"
    ),
    
    heatmap_legend_param = list(
      title = "Pearson r",
      at = c(-1, -0.5, 0, 0.5, 1)
    )
  )
  
  # Save PDF
  pdf(
    paste0("modules_age_correlation_", group, ".pdf"),
    width = max(9, ncol(cor_mat) * 0.8 + 5),
    height = max(7, nrow(cor_mat) * 0.35 + 3)
  )
  
  draw(ht)
  dev.off()
  
  # Display
  draw(ht)
  
  invisible(ht)
}

# ============================================================
# 6. GENERATE TWO HEATMAPS
# ============================================================

ht_healthy <- plot_correlation_heatmap(
  cor_results, "Healthy"
)

ht_disease <- plot_correlation_heatmap(
  cor_results, "Disease"
)

# ============================================================
# 7. SUMMARY
# ============================================================

print(
  cor_results %>%
    group_by(Group) %>%
    summarise(
      N_correlations = n(),
      N_datasets = n_distinct(GSE),
      N_significant_FDR = sum(FDR < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
)


