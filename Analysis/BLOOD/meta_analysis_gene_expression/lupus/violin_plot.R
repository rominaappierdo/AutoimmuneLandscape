setwd("/mnt/data3/appierdo/PhD-project/Analysis/BLOOD/meta_analysis_gene_expression/lupus")
library(MetaIntegrator)
library(ggplot2)
library(ggpubr)
library(viridis)
library(gridExtra)
data <- qs::qread("lupus_vs_healthy_whole blood_metaintegrator_object.qs")
dir <- "/mnt/data3/appierdo/PhD-project/Data_Collection/downloaded_datasets/lupus/"
meta <- read.csv("/mnt/data3/appierdo/PhD-project/Data_Collection/raw_data/step6_starting_manage_data/annotated_metadata_8.1.2023/lupus_metadata.tsv", sep = '\t')
meta[meta==""] <- NA
list_plot <- list()
timepoint <- c('baseline', 'pre', NA,"NA", 'rep1', 'unstimulated','untreated', "visit 1",'T0', 't0', "biopsy 1", "GPL96", "visit day 0", "week 0", "at-diagnosis (for progressor)","day 0" , "at diagnosis" , "untreated, rep 1", "baseline, rep1", "early")
c <- 0
for (gse in colnames(data$metaAnalysis$datasetEffectSizes)) {
  files <- list.files(dir, recursive = T, pattern = gse, full.names = T) 
  file <- files[grepl("_gs\\.tsv|\\.rds", files)]
  sub_meta <- meta[meta$series_id == gse,]
  sub_meta <- sub_meta[grepl('heal|lupus', sub_meta$disease),]
  sub_meta[grepl("lupus", sub_meta$disease),]$disease <- 'lupus'
  sub_meta <- sub_meta[grepl('blood', sub_meta$sample_type),]
  sub_meta <- sub_meta[sub_meta$timepoint %in% timepoint,]
  
  if (grepl("rds", file)) {
    df <- readRDS(file)$expr
  } else {
    df <- read.csv(file, sep = '\t', row.names = 1)
  }
  df <- df[, colnames(df) %in% sub_meta$geo_accession]
  df <- df[rownames(df) == 'IRF9',, drop = F]
  if (nrow(df) == 0) {
    next
  } else {
    c <- c+1
    print(c)
  }
  df <- as.data.frame(t(df))
  df$class <- sub_meta$disease
  
  plot <- ggplot(df, aes(x = class, fill = class, y = IRF9)) +
    geom_violin() +
    stat_compare_means() +
    theme_classic()+
    ggtitle(gse) +
    theme(plot.title = element_text(hjust = 0.5, size = 8), legend.position='bottom')+
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE)+
    ylab("IRF9 expression") +
    xlab('') +
    guides(shape = guide_legend(override.aes = list(size = 0.5)))

  list_plot[[c]] <- plot

}

do.call(grid.arrange,c(list_plot, ncol = 5))

