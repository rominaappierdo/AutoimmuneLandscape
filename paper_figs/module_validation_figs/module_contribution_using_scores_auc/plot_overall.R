setwd("paper_figs/module_validation_figs/module_contribution_using_scores_auc")
source("../../../utils_functions/functions_to_process_RNAseq_Data.R")
library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyr)
library(viridis)
library(paletteer)
library(ggtext)
files <- list.files(".", recursive = T, pattern = 'scores.rds', full.names = T)

module_names <- c("Th1 response",  "Monocyte/macrophage_NK_crosstalk", "T_cell_activation_and_recruitment",
                  "Anti-inflammatory response", "IFN response", "IL6_trans_signaling/transendothelial_migration",
                  "Lymphoid cell development", "Endothelial dysfunction", "TLR/NFKB response/innate_immune_response",
                  "B cell","Integrins", "ECM remodelling", "Inflammation/Wound healing", "Impaired cell proliferation", "Lipid metabolism")

diseases <- c()
list_plot <- list()
c <- 0
for (file in files) {
  c <- c+1
  disease <- strsplit(file, "/")[[1]][2]
  disease <- toupper(disease)
  disease <- gsub("_BLOOD", " (blood)", disease)
  disease <- gsub("_TISSUE", " (tissue)", disease)
  
  df <- readRDS(file)
  df$disease <- rep(disease, nrow(df))
  
  
  if (c == 1){
    final <- df
  } else {
    final <- rbind(final, df)
  }
  
}

ggplot(final, aes(x = disease, y = module, size = mean_aucs, color = mean_aucs)) +
  geom_point()+
  scale_size(range = c(2,10)) +
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.5)



long_blood <- final[grepl('blood', final$disease),]
long_blood <- sort_by_col(long_blood, 'disease')

long_blood$disease <- gsub("CROHN \\(blood\\)", "Crohn's disease", long_blood$disease)
long_blood$disease <- gsub("JIA \\(blood\\)", "Juvenile idiopathic arthritis", long_blood$disease)
long_blood$disease <- gsub("UC \\(blood\\)", "Ulcerative colitis", long_blood$disease)
long_blood$disease <- gsub("SJOGREN \\(blood\\)", "Sjogren's syndrome", long_blood$disease)
long_blood$disease <- gsub("SCLERODERMA \\(blood\\)", "Scleroderma", long_blood$disease)
long_blood$disease <- gsub("RA \\(blood\\)", "Rheumatoid arthritis", long_blood$disease)
long_blood$disease <- gsub("PSORIASIS \\(blood\\)", "Psoriasis", long_blood$disease)
long_blood$disease <- gsub("MULTIPLE.SCLEROSIS \\(blood\\)", "Multiple sclerosis", long_blood$disease)
long_blood$disease <- gsub("LUPUS \\(blood\\)", "Lupus", long_blood$disease)
long_blood$disease <- gsub("DIABETES \\(blood\\)", "Type I diabetes", long_blood$disease)


ggplot(long_blood, aes(x = disease, y = module, size = mean_aucs, color = mean_aucs)) +
  geom_point()+
  scale_size(range = c(2,10)) +
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.5)+
  ggtitle("BLOOD")


long_tissue <- final[grepl('tissue', final$disease),]
long_tissue <- sort_by_col(long_tissue, 'disease')

long_tissue$disease <- gsub("CROHN \\(tissue\\)", "Crohn's disease (Gut)", long_tissue$disease)
long_tissue$disease <- gsub("UC \\(tissue\\)", "Ulcerative colitis (Gut)", long_tissue$disease)
long_tissue$disease <- gsub("SJOGREN \\(tissue\\)", "Sjogren's syndrome (Salivary glands)", long_tissue$disease)
long_tissue$disease <- gsub("SCLERODERMA \\(tissue\\)", "Scleroderma (Skin)", long_tissue$disease)
long_tissue$disease <- gsub("RA \\(tissue\\)", "Rheumatoid arthritis (Synovium)", long_tissue$disease)
long_tissue$disease <- gsub("PSORIASIS \\(tissue\\)", "Psoriasis (Skin)", long_tissue$disease)
long_tissue$disease <- gsub("LUPUS \\(tissue\\)_SKIN", "Lupus (Skin)", long_tissue$disease)
long_tissue$disease <- gsub("LUPUS \\(tissue\\)_KIDNEY", "Lupus (Kidney)", long_tissue$disease)


ggplot(long_tissue, aes(x = disease, y = module,size = mean_aucs, color = mean_aucs)) +
  geom_point()+
  scale_size(range = c(2,10)) +
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.5) +
  ggtitle("TISSUE")



#to use facet

long_blood$source <- rep('Blood', nrow(long_blood))
long_tissue$source <- rep('Tissue', nrow(long_tissue))

l <- rbind(long_blood, long_tissue)
colnames(l) <- c("module", "integrated_pvals", "Mean AUCs", "disease", "source")

ggplot(l, aes(x = disease, y = module,size = `Mean AUCs`, color = `Mean AUCs`)) +
  geom_point()+
  scale_size(range = c(2,10)) +
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0.8, "cm", data = NULL),
        strip.text = element_text(
          size = 14, color = "black", face = 'bold'))+
  #scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.5) +
  facet_grid(cols = vars(source), scales = 'free') +
  xlab("") +
  ylab("") +
  scale_color_gradientn(colours = paletteer_c("ggthemes::Orange-Blue Diverging", 30), limits = c(0, 1))


#remake plot adding simpson index

# simpson_index <- function(aucs) {
#   total <- sum(aucs)
#   proportions <- aucs / total
#   D <- 1 - sum(proportions^2)
#   return(D)
# }

standard_dev <- function(aucs) {
  D <- sd(aucs)
  return(D)
}



files <- list.files(".", recursive = T, pattern = '.*simp.*.rds', full.names = T)

b <- 0
for (file in files) {
  b <- b+1
  f <- readRDS(file)
  #cvs <- as.data.frame(sapply(f, simpson_index))
  cvs <- as.data.frame(sapply(f, standard_dev))
  
  colnames(cvs) <- 'st_dev'
  disease <- strsplit(file, "/")[[1]][2]
  cvs$disease <- rep(disease, nrow(cvs))
  cvs$module <- rownames(cvs)
  
  if (b == 1){
    final_cvs <- cvs
  } else {
    final_cvs <- rbind(final_cvs, cvs)
  }
  
}

unique(final_cvs$disease)

final_cvs[final_cvs$disease == 'crohn_blood', ]$disease <- "Crohn's disease"
final_cvs[final_cvs$disease == 'lupus_blood', ]$disease <- "Lupus"
final_cvs[final_cvs$disease == 'crohn_tissue', ]$disease <- "Crohn's disease (Gut)"
final_cvs[final_cvs$disease == 'diabetes_blood', ]$disease <- "Type I diabetes"
final_cvs[final_cvs$disease == 'lupus_tissue_skin', ]$disease <- "Lupus (Skin)"
final_cvs[final_cvs$disease == 'lupus_tissue_kidney', ]$disease <- "Lupus (Kidney)"
final_cvs[final_cvs$disease == 'multiple.sclerosis_blood', ]$disease <- "Multiple sclerosis"
final_cvs[final_cvs$disease == 'psoriasis_blood', ]$disease <- "Psoriasis"
final_cvs[final_cvs$disease == 'psoriasis_tissue', ]$disease <- "Psoriasis (Skin)"
final_cvs[final_cvs$disease == 'ra_blood', ]$disease <- "Rheumatoid arthritis"
final_cvs[final_cvs$disease == 'ra_tissue', ]$disease <- "Rheumatoid arthritis (Synovium)"
final_cvs[final_cvs$disease == 'scleroderma_blood', ]$disease <- "Scleroderma"
final_cvs[final_cvs$disease == 'scleroderma_tissue', ]$disease <- "Scleroderma (Skin)"
final_cvs[final_cvs$disease == 'sjogren_blood', ]$disease <- "Sjogren's syndrome"
final_cvs[final_cvs$disease == 'sjogren_tissue', ]$disease <- "Sjogren's syndrome (Salivary glands)"
final_cvs[final_cvs$disease == 'uc_blood', ]$disease <- "Ulcerative colitis"
final_cvs[final_cvs$disease == 'uc_tissue', ]$disease <- "Ulcerative colitis (Gut)"
final_cvs[final_cvs$disease == 'jia_blood', ]$disease <- "Juvenile idiopathic arthritis"



final_merged <- merge(l, final_cvs, by = c("disease", "module"))

final_merged$`Mean AUCs`[final_merged$integrated_pvals > 0.05] <- NA



unique(final_merged$module)

ordered_vector <- c("IFN response", "Lymphoid cell development", "Endothelial dysfunction" , "TLR/NFKB response/innate_immune_response",
                    "Lipid metabolism",  "Th1 response" , "Inflammation/Wound healing", "Monocyte/macrophage_NK_crosstalk", "B cell" ,
                    "Impaired cell proliferation", "Integrins", "IL6_trans_signaling/transendothelial_migration", "T_cell_activation_and_recruitment",  
                    "Anti-inflammatory response", "ECM remodelling")                             

final_merged$module <- factor(final_merged$module, levels = rev(ordered_vector))

p <- ggplot(final_merged, aes(x = disease, y = module,size = 1/st_dev, color = `Mean AUCs`)) +
  geom_point()+
  scale_y_discrete(
    position = "right"
  )  +
  theme_classic()+
  theme(
    axis.text.y.right = element_text(hjust = 0),  # Adjust the position of y-axis text
    legend.position = "left",  # Move legends to the left
    axis.text.y = element_markdown()
  ) +# Move y-axis to the right side
  scale_size(range = c(6,10)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0.8, "cm", data = NULL),
        strip.text = element_text(
          size = 14, color = "black", face = 'bold'))+
  facet_grid(cols = vars(source), scales = 'free') +
  xlab("") +
  ylab("") +
  scale_color_gradientn(colours = paletteer_c("ggthemes::Orange-Blue Diverging", 30), limits = c(0, 1), na.value="white") 




p
#with colors
p2 <- ggplot(final_merged, aes(x = disease, y = module,size = 1/st_dev, color = `Mean AUCs`)) +
  geom_point()+
  scale_y_discrete(
    #position = "right",
    labels = c(
      "TLR/NFKB response/innate_immune_response" = "<span style='color:#ffc27f; font-size:18px;'><strong>TLR/NFKB_response/innate_immune_response</strong></span>",
      "Th1 response" = "<span style='color:#99d68c; font-size:18px;'><strong>Th1_response</strong></span>",
      "T_cell_activation_and_recruitment" = "<span style='color:#78706e; font-size:18px;'><strong>T_cell_activation_and_recruitment</strong></span>",
      "Monocyte/macrophage_NK_crosstalk" = "<span style='color:#f1ce64; font-size:18px;'><strong>Monocyte/macrophage_NK_crosstalk</strong></span>",
      "Lymphoid cell development" = "<span style='color:#a2cdeb; font-size:18px;'><strong>Lymphoid cell development</strong></span>",
      "Lipid metabolism" = "<span style='color:#58a14e; font-size:18px;'><strong>Lipid metabolism</strong></span>",
      "Integrins" = "<span style='color:#e1575a; font-size:18px;'><strong>Integrins</strong></span>",
      "Inflammation/Wound healing" = "<span style='color:#b6992d; font-size:18px;'><strong>Inflammation/Wound healing</strong></span>",
      "Impaired cell proliferation" = "<span style='color:#86bcb6; font-size:18px;'><strong>Impaired cell proliferation</strong></span>",
      "IL6_trans_signaling/transendothelial_migration" = "<span style='color:#ff9e9a; font-size:18px;'><strong>IL6_trans_signaling/transendothelial_migration</strong></span>",
      "IFN response" = "<span style='color:#4e78a7; font-size:18px;'><strong>IFN_response</strong></span>",
      "Endothelial dysfunction" = "<span style='color:#f39940; font-size:18px;'><strong>Endothelial dysfunction</strong></span>",
      "ECM remodelling" = "<span style='color:#d37295; font-size:18px;'><strong>ECM remodelling</strong></span>",
      "B cell" = "<span style='color:#489894; font-size:18px;'><strong>B cell</strong></span>",
      "Anti-inflammatory response" = "<span style='color:#b9afac; font-size:18px;'><strong>Anti-inflammatory response</strong></span>"
    )
    
  )  +
  theme_classic()+
  theme(
    #axis.text.y.right = element_text(hjust = 0),  # Adjust the position of y-axis text
    #legend.position = "left",  # Move legends to the left
    axis.text.y = element_markdown()
  ) +# Move y-axis to the right side
  scale_size(range = c(6,10)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0.8, "cm", data = NULL),
        strip.text = element_text(
          size = 14, color = "black", face = 'bold'))+
  facet_grid(cols = vars(source), scales = 'free') +
  xlab("") +
  ylab("") +
  scale_color_gradientn(colours = paletteer_c("ggthemes::Orange-Blue Diverging", 30), limits = c(0, 1), na.value="white") 

p


jpeg("fig6B.jpeg", units = 'in', width = 14, height = 8, res = 600)
p
dev.off()


#mean AUC IFN blood
mean_ifn_blood <- final_merged[final_merged$source == 'Blood', ]
mean_ifn_blood <- mean_ifn_blood[mean_ifn_blood$module == 'IFN response', ]
mean(na.omit(mean_ifn_blood$`Mean AUCs`))
#mean AUC IFN tissue
mean_ifn_tissue <- final_merged[final_merged$source == 'Tissue', ]
mean_ifn_tissue <- mean_ifn_tissue[mean_ifn_tissue$module == 'IFN response', ]
mean(na.omit(mean_ifn_tissue$`Mean AUCs`))



