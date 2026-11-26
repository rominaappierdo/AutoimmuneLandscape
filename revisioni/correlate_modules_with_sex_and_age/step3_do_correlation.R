setwd("/mnt/data3/appierdo/AutoimmuneLandscape/revisioni/correlate_modules_with_sex_and_age")
source("/mnt/data3/appierdo/AutoimmuneLandscape/utils_functions/functions_to_process_RNAseq_Data.R")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

gses <- c("GSE9006", "GSE17590")

meta <- complete_metadata()

sub_meta <- meta[meta$gender %in% c("M", "F"), ]
sub_meta <- sub_meta[!is.na(sub_meta$age) & sub_meta$age != "", ]
sub_meta <- sub_meta[!grepl("month|M", sub_meta$age), ]
sub_meta$age <- suppressWarnings(as.numeric(sub_meta$age))
sub_meta <- sub_meta[!is.na(sub_meta$age), ]
sub_meta <- sub_meta[sub_meta$cell_type == "", ]

sub_meta <- sub_meta[
  sub_meta$timepoint %in% c(
    "baseline", "pre", NA, "NA", "rep1",
    "unstimulated", "untreated", "visit 1",
    "T0", "t0", "biopsy 1", "GPL96",
    "visit day 0", "week 0",
    "at-diagnosis (for progressor)",
    "day 0", "at diagnosis",
    "untreated, rep 1", "baseline, rep1", "early"
  ),
]

cor_results <- data.frame()

for (gse in gses) {
  
  files_modules_gse <- list.files(
    ".", pattern = gse, recursive = TRUE
  )
  
  for (file in files_modules_gse) {
    
    # Reset metadata for each module
    gse_meta <- sub_meta[sub_meta$series_id == gse, ]
    
    module_table <- read.csv(
      file, sep = "\t", row.names = 1, check.names = FALSE
    )
    
    gse_meta <- gse_meta[
      gse_meta$geo_accession %in% colnames(module_table),
    ]
    
    gse_meta <- sort_by_col(gse_meta, "geo_accession")
    module_table <- module_table[
      , gse_meta$geo_accession, drop = FALSE
    ]
    
    if (!identical(colnames(module_table), gse_meta$geo_accession)) {
      stop("Different samples")
    }
    
    # Compute module scores
    gse_meta$module <- colSums(module_table)
    
    # Pearson correlation
    valid <- complete.cases(gse_meta[, c("age", "module")])
    
    if (sum(valid) < 3 ||
        sd(gse_meta$age[valid]) == 0 ||
        sd(gse_meta$module[valid]) == 0) next
    
    cor_test <- cor.test(
      gse_meta$age[valid],
      gse_meta$module[valid],
      method = "pearson"
    )
    
    # Extract module name
    module_name <- basename(file)
    module_name <- sub("\\.(csv|tsv|txt)$", "", module_name)
    module_name <- sub("^GSE[0-9]+_", "", module_name)
    
    # Store correlation
    cor_results <- rbind(
      cor_results,
      data.frame(
        GSE = gse,
        Module = module_name,
        correlation = unname(cor_test$estimate),
        pvalue = cor_test$p.value,
        n = sum(valid)
      )
    )
    
    # Scatterplot
    p <- ggplot(gse_meta, aes(x = age, y = module)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      stat_cor(method = "pearson") +
      theme_classic() +
      ggtitle(paste(gse, module_name, sep = " | "))
    
    print(p)
  }
}

# ----------------------------
# COMPLEXHEATMAP
# ----------------------------

cor_results <- cor_results %>%
  mutate(FDR = p.adjust(pvalue, method = "BH"))

# Correlation matrix
cor_mat <- cor_results %>%
  select(Module, GSE, correlation) %>%
  pivot_wider(names_from = GSE, values_from = correlation) %>%
  tibble::column_to_rownames("Module") %>%
  as.matrix()

# FDR matrix
fdr_mat <- cor_results %>%
  select(Module, GSE, FDR) %>%
  pivot_wider(names_from = GSE, values_from = FDR) %>%
  tibble::column_to_rownames("Module") %>%
  as.matrix()

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
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 1),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    r <- cor_mat[i, j]
    fdr <- fdr_mat[i, j]
    
    if (!is.na(r)) {
      
      # Correlation coefficient + significance
      label <- sprintf("%.2f", r)
      
      if (!is.na(fdr) && fdr < 0.05) {
        label <- paste0(label, "*")
      }
      
      grid.text(
        label, x, y,
        gp = gpar(
          fontsize = 10,
          fontface = "bold",
          col = ifelse(abs(r) > 0.7, "white", "black")
        )
      )
    }
  },
  
  column_title = "Correlation between module scores and age"
)

draw(ht)



write.csv(cor_results, "modules_age_correlations.csv", row.names = FALSE)