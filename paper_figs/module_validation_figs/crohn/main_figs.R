setwd("paper_figs/module_validation_figs/crohn")
library(MetaIntegrator)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(pROC)
library(ggthemes)
library(ggpubr)
library(factoextra)
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
source("/mnt/data3/appierdo/PhD-project/utils_functions/function_compute_all_features_for_module.R")
source("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/pipeline_landscape_disease_level_final/modules_validation/readLines_all_modules.R")


gse <- 'GSE16879'

list_plot_after <- readRDS("list_plot_before_treat_colon.rds")
do.call(grid.arrange,c(list_plot_after[13], ncol = 1))

################################################################################################################################################
#Question #2: Can modules capture differences in response at baseline to predict response to treatment?

list_plot_before <- readRDS("list_plot_colon_before.rds")
do.call(grid.arrange,c(list_plot_before[13], ncol = 1))

do.call(grid.arrange,c(c(list_plot_after[13], list_plot_before[13]), ncol = 1))


