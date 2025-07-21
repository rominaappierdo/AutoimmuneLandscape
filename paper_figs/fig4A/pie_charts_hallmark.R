setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig4")
library(ggplot2)
library(gridExtra)

df_blood <- data.frame()
blood_files_hall <- list.files("../../../Analysis/BLOOD/meta_analysis_pathways_activity/hallmarks/", pattern = 'complete_res', recursive = T, full.names = T)

blood_list <- list()
dis_list_blood <- c()
b <- 0
for (file in blood_files_hall) {
  b<- b+1
  disease <- strsplit(file, "/")[[1]][9]
  df <- read.csv(file, sep ='\t')
  df <- df[df$numStudies>=2,]
  df <- df[df$effectSizeFDR < 0.05,]
  #df <- df[abs(df$effectSize) > 0.8,]
  blood_list[[b]] <- df$gene
  dis_list_blood <- append(dis_list_blood, disease)
}
names(blood_list) <- dis_list_blood


tissue_files_hall <- list.files("../../../Analysis/TISSUE/meta_analysis_pathways_activity/hallmarks/", pattern = 'complete_res', recursive = T, full.names = T)

tissue_list <- list()
dis_list_tissue <- c()
t <- 0
for (file in tissue_files_hall) {
  t <- t+1
  disease <- strsplit(file, "/")[[1]][9]
  tissue <- strsplit(file, "/|_")[[1]][13]
  disease <- paste0(disease, "_", tissue, sep = '')
  df <- read.csv(file, sep ='\t')
  df <- df[df$numStudies>=2,]
  df <- df[df$effectSizeFDR < 0.05,]
  
  #df <- df[abs(df$effectSize) > 0.8,]
  tissue_list[[t]] <- df$gene
  dis_list_tissue <- append(dis_list_tissue, disease)
  
}
names(tissue_list) <- dis_list_tissue

custom_colors <- c("Blood" = "#ff8d8d", "Tissue" = "#cbe7f7", "Blood and Tissue" = "#a5292a")
list_plot_genes <- list()
c <- 0
for (name in names(blood_list)){
  c <- c+1
  print(name)
  
  n_blood <- length(blood_list[[name]])
  
  if (any(grepl(name, names(tissue_list)))) {
    n_tissue <- length(tissue_list[grepl(name, names(tissue_list))][[1]])
    n_intersection <- length(intersect(blood_list[[name]], tissue_list[grepl(name, names(tissue_list))][[1]]))
    
    tot <- length(unique(c(blood_list[[name]], tissue_list[grepl(name, names(tissue_list))][[1]])))
    
    data <- data.frame(
      category = c("Blood", "Tissue", "Blood and Tissue"),
      count = c(as.numeric(n_blood), as.numeric(n_tissue), as.numeric(n_intersection))
      
    )
    
    data$percentage <- data$count / sum(data$count) * 100
    
    
    # Add the hole in the center
    p <- ggplot(data, aes(x = 2, y = percentage, fill = category)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      theme_void() +
      theme(legend.position = "right") +
      scale_fill_manual(values = custom_colors)+
      annotate("text", x = 0, y = 0, label = tot, size = 6, color = "black") +
      theme(legend.position = "none")
    
  } else {
    tot <- length(unique(blood_list[[name]]))
    data <- data.frame(
      category = c("Blood"),
      count = c(as.numeric(n_blood))
    )
    data$percentage <- data$count / sum(data$count) * 100
    
    # Add the hole in the center
    p <- ggplot(data, aes(x = 2, y = percentage, fill = category)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      theme_void() +
      theme(legend.position = "right") +
      scale_fill_manual(values = custom_colors)+
      annotate("text", x = 0, y = 0, label = tot, size = 6, color = "black") +
      theme(legend.position = "none")
  }
  list_plot_genes[[c]] <- p
}


do.call(grid.arrange,c(list_plot_genes, ncol = 10))

saveRDS(list_plot_genes, "pie_chart_hallmark.rds")