library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("/mnt/data3/appierdo/AutoimmuneLandscape/revisioni/correlate_modules_with_sex_and_age")

# Load correlation results
cor_results <- read.csv("modules_age_correlations.csv")

# Keep only correlations with unadjusted p < 0.05
cor_results$correlation[
  is.na(cor_results$pvalue) | cor_results$pvalue >= 0.05
] <- NA_real_

# Build correlation matrix
cor_mat <- cor_results %>%
  select(Module, GSE, correlation) %>%
  pivot_wider(names_from = GSE, values_from = correlation) %>%
  column_to_rownames("Module") %>%
  as.matrix()

# Remove modules and GSEs without significant correlations
cor_mat <- cor_mat[
  rowSums(!is.na(cor_mat)) > 0,
  colSums(!is.na(cor_mat)) > 0,
  drop = FALSE
]

if (nrow(cor_mat) == 0 || ncol(cor_mat) == 0) {
  stop("No significant correlations (p < 0.05).")
}

# ComplexHeatmap
ht <- Heatmap(
  cor_mat,
  name = "Pearson r",
  
  col = colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")
  ),
  
  na_col = "white",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  
  rect_gp = gpar(col = "grey90", lwd = 0.5),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    r <- cor_mat[i, j]
    
    if (!is.na(r)) {
      grid.text(
        sprintf("%.2f", r),
        x, y,
        gp = gpar(
          fontsize = 9,
          fontface = "bold",
          col = ifelse(abs(r) > 0.7, "white", "black")
        )
      )
    }
  },
  
  column_title = "Age correlations (p < 0.05, unadjusted)"
)

# Display
draw(ht)

# Save PDF
pdf(
  "modules_age_correlations_pvalue_significant.pdf",
  width = max(8, ncol(cor_mat) * 0.8 + 4),
  height = max(6, nrow(cor_mat) * 0.35 + 2)
)

draw(ht)
dev.off()