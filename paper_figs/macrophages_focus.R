setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs")
library(gridExtra)

#BLOOD

files <- list.files("../Analysis/BLOOD/meta_analysis_cell_proportions/", recursive = T, pattern = 'complet', full.names = T)

ess <- c()
for (file in files) {
  print(file)
  df <- read.csv(file, sep = '\t')
  es <- df[df$gene == 'macrophage_m1', ]$effectSize
  print(es)
  ess <- append(ess, es)
}

df_m1_blood <- data.frame(ess)
rownames(df_m1_blood) <- c('CD', 'T1D', 'JIA', 'SLE', 'MS', 'Ps', 'RA', 'SSc', 'SS', 'UC')
colnames(df_m1_blood) <- c('macrophage_M1_blood')


ess <- c()
for (file in files) {
  print(file)
  df <- read.csv(file, sep = '\t')
  es <- df[df$gene == 'macrophage_m2', ]$effectSize
  print(es)
  ess <- append(ess, es)
}

df_m2_blood <- data.frame(ess)
rownames(df_m2_blood) <- c('CD', 'T1D', 'JIA', 'SLE', 'MS', 'Ps', 'RA', 'SSc', 'SS', 'UC')
colnames(df_m2_blood) <- c('macrophage_M2_blood')


ess <- c()
for (file in files) {
  print(file)
  df <- read.csv(file, sep = '\t')
  es <- df[df$gene == 'macrophage_m0', ]$effectSize
  print(es)
  ess <- append(ess, es)
}

df_m0_blood <- data.frame(ess)
rownames(df_m0_blood) <- c('CD', 'T1D', 'JIA', 'SLE', 'MS', 'Ps', 'RA', 'SSc', 'SS', 'UC')
colnames(df_m0_blood) <- c('macrophage_M0_blood')


#TISSUE

files <- list.files("../Analysis/TISSUE/meta_analysis_cell_proportions/", recursive = T, pattern = 'complet', full.names = T)

ess <- c()
for (file in files) {
  df <- read.csv(file, sep = '\t')
  es <- df[df$gene == 'macrophage_m1', ]$effectSize
  print(es)
  ess <- append(ess, es)
}

p2 <- plot(density(ess))




