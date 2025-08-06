setwd("Analysis")

blood_files <- list.files("BLOOD/meta_analysis_cell_proportions/", recursive = T, pattern = 'complet', full.names = T)

dis <- c()
ess <- c()
for (file in blood_files) {
  disease <- strsplit(file, "\\/")[[1]][4]
  
  df <- read.csv(file, sep = '\t')
  print(head(df))
  
  es <- df[df$gene == 'macrophage_m1', ]$effectSize
  
  dis <- append(dis, disease)
  ess<- append(ess, es)
  
  
}


plot(density(ess))



tissue_files <- list.files("TISSUE/meta_analysis_cell_proportions", recursive = T, pattern = 'complet', full.names = T)

dis2 <- c()
ess2 <- c()
for (file in tissue_files) {
  disease <- strsplit(file, "\\/|_")[[1]][7]
  t <- strsplit(file, "\\/|_")[[1]][8]
  
  disease <- paste0(disease, "_", t)
 
  
  df <- read.csv(file, sep = '\t')
  print(head(df))
  
  es <- df[df$gene == 'macrophage_m1', ]$effectSize
  
  dis2 <- append(dis2, disease)
  ess2<- append(ess2, es)
  
  
}


plot(density(ess2))






library(ggplot2)
library(dplyr)

# Blood data
ess <- c(0.057289639, -0.110411887, 0.437787601, 1.228594825, 0.009869739, 
         0.428317246, 0.130409241, 0.580672465, 0.965366830, 0.266855842)
dis <- c("crohn", "diabetes", "jia", "lupus", "multiple_sclerosis", 
         "psoriasis", "ra", "scleroderma", "sjogren", "uc")

# Tissue data
ess2 <- c(0.7737318, 0.4498384, 2.5827669, 0.1405023, 0.4653219, 
          0.6458281, 0.6599441, 0.3074475)
dis2 <- c("crohn_GUT", "lupus_KIDNEY", "lupus_SKIN", "psoriasis_SKIN", 
          "ra_SYNOVIUM", "scleroderma_SKIN", "sjogren_SALIVARY.GLAND", "uc_GUT")

# Combine into one dataframe
df <- data.frame(
  disease = c(dis, dis2),
  effect_size = c(ess, ess2),
  source = c(rep("Blood", length(ess)), rep("Tissue", length(ess2)))
)

# Plot density
ggplot(df, aes(x = effect_size, fill = source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Effect Size Distribution by Sample Type",
       x = "Effect Size",
       y = "Density",
       fill = "Source") +
  theme_classic() +
  xlim(-1, 3)



library(ggplot2)
median(ess)
median(ess2)
# Create numeric vectors again (double-checking)
ess <- c(0.057289639, -0.110411887, 0.437787601, 1.228594825, 0.009869739, 
         0.428317246, 0.130409241, 0.580672465, 0.965366830, 0.266855842)
dis <- c("crohn", "diabetes", "jia", "lupus", "multiple_sclerosis", 
         "psoriasis", "ra", "scleroderma", "sjogren", "uc")

ess2 <- c(0.7737318, 0.4498384, 2.5827669, 0.1405023, 0.4653219, 
          0.6458281, 0.6599441, 0.3074475)
dis2 <- c("crohn_GUT", "lupus_KIDNEY", "lupus_SKIN", "psoriasis_SKIN", 
          "ra_SYNOVIUM", "scleroderma_SKIN", "sjogren_SALIVARY.GLAND", "uc_GUT")

# Combine and check structure
df <- data.frame(
  disease = c(dis, dis2),
  effect_size = as.numeric(c(ess, ess2)),
  source = factor(c(rep("Blood", length(ess)), rep("Tissue", length(ess2))))
)

# Ensure no NA values
df <- na.omit(df)

# Plot histogram
ggplot(df, aes(x = effect_size, fill = source)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 10, color = "black") +
  labs(title = "Effect Size Distribution in Blood and Tissue",
       x = "Effect Size",
       y = "Count",
       fill = "Source") +
  theme_minimal(base_size = 14)
