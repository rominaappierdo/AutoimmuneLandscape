setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig4")
library(gridExtra)
library(grid)
library(scales)
library(ggplot2)
final_list <- list()

files <- c("pie_chart_genes.rds",
           "pie_chart_cyto.rds",
           "pie_chart_cells.rds",
           "pie_chart_all_pathways.rds",
           "pie_chart_mirna.rds",
           "pie_chart_tf.rds")

for (file in files) {
  l <- readRDS(file)
  final_list <- append(final_list, l)
}


do.call(grid.arrange,c(final_list, ncol = 10))
#grid.lines(x = c(0, 1), y = c(0.835, 0.835), gp = gpar(col = "black", lwd = 2))
