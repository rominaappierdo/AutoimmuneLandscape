setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Supplementary/supp_figures")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
source("../../utils_functions/function_compute_all_features_for_module.R")
source("../../Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_files/readLines_all_modules.R")
set.seed(1234)

t0 <- readRDS("density_t0.rds")

p0 <- ggplot(df, aes(x = density)) +
  geom_density(fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 18) +  # increase all font sizes
  geom_vline(xintercept = g_real, color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Density Plot - Z-Score Comparison - t0",
    x = "Hedges' g",
    y = "Density"
  ) +
  annotate(
    "text",
    x = g_real - 0.2,  # move slightly left of the red line
    y = Inf, hjust = 1.1, vjust = 1.3,  # adjust placement near top
    label = paste(
      "p =", signif(p_value, digits = 3),
      "\nz =", signif(z_score, digits = 3)
    ),
    color = "red", size = 6
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )

p0
