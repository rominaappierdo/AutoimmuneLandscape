#setwd("Analysis")
setwd("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/")
blood_files <- list.files("BLOOD/meta_analysis_cell_proportions/", recursive = T, pattern = 'complet', full.names = T)

dis <- c()
ess <- c()
for (file in blood_files) {
  disease <- strsplit(file, "\\/")[[1]][4]
  
  df <- read.csv(file, sep = '\t')
  print(head(df))
  
  es <- df[df$gene == 'plasmacytoid_dendritic_cell', ]$effectSize
  
  dis <- append(dis, disease)
  ess<- append(ess, es)
  
  
}


plot(density(ess))
mean(ess)


tissue_files <- list.files("TISSUE/meta_analysis_cell_proportions", recursive = T, pattern = 'complet', full.names = T)

dis2 <- c()
ess2 <- c()
for (file in tissue_files) {
  disease <- strsplit(file, "\\/|_")[[1]][7]
  t <- strsplit(file, "\\/|_")[[1]][8]
  
  disease <- paste0(disease, "_", t)
 
  
  df <- read.csv(file, sep = '\t')
  print(head(df))
  
  es <- df[df$gene == 'plasmacytoid_dendritic_cell', ]$effectSize
  
  dis2 <- append(dis2, disease)
  ess2<- append(ess2, es)
  
  
}
mean(ess2)

plot(density(ess2))


