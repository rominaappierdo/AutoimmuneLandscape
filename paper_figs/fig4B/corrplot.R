setwd("paper_figs/fig4B")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
library(ggplot2)
library(viridis)
library(scales) 


df <- t(data.frame(
  c(0, 0,0,0,0.68,0,0, 0),
  c(0.34, 0, 0, 0, 0, 0.62, 0.2 ,0),
  c(0.12, 0, 0, 0.37, 0.22, 0.27, 0, 0),
  c(0.43, 0, 0.36, 0, 0.25, 0.33, 0, 0),
  c(0.41, 0.3, 0.4, 0.3, 0, 0.42, 0, 0),
  c(0.29, 0.012, 0.29, 0, 0, 0.45, 0, 0.14),
  c(0.39, 0.16, 0.35, 0.12, 0.15, 0.37, 0.12, 0)
))

rownames(df) <- c("Cell types", "Transcription factors", "Cytokines", "KEGG", "Hallmark", "miRNA", "WikiPathways")
colnames(df) <- c("Crohn's disease", "Lupus (Kidney)", "Lupus (Skin)", "Rheumatoid arthritis", 
                  "Scleroderma", "Sjogren's syndrome", "Ulcerative colitis", "Psoriasis")

df <- wide_to_long(df, name_row_column = 'feature', name_column_column = 'disease', 'correlation')

df[df == 0] <- NA
df$feature <- factor(df$feature, levels = c("Cell types", "Transcription factors", "Cytokines", "KEGG", "Hallmark", "miRNA", "WikiPathways"))
df$disease <- factor(df$disease, levels = c("Crohn's disease", "Lupus (Kidney)", "Lupus (Skin)", "Rheumatoid arthritis", 
                                            "Scleroderma", "Sjogren's syndrome", "Ulcerative colitis", "Psoriasis")
)


ggplot(df, aes(x = disease, y = feature, color = correlation))+
  geom_point(size = 8) +
  theme_classic()+
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 1),
        legend.position = 'top') +
  scale_size(range = c(1,15)) +
  # scale_color_gradient2(#midpoint=0, 
  #                       low="blue", 
  #                       #mid="white",
  #                       high="red", 
  #                       limits = c(-1,1),na.value = 'white'
  #                       )
  scale_color_gradient(low="blue", high="red", na.value = 'white', limits = c(-1, 1))+
  #scale_color_gradientn(colors = c("blue", "green", "yellow", "red"), na.value = 'white', limits = c(-1, 1))
  #scale_color_viridis_c(option = "magma", na.value = 'white', limits = c(-1, 1))+
  xlab("")+
  ylab("") +  theme(
    axis.text.y.right = element_text(hjust = 0),  # Adjust the position of y-axis text
    legend.position = "left"  # Move legends to the left
  ) +
  scale_y_discrete(position = "right")  # Move y-axis to the right side


library(RColorBrewer) # For predefined palettes like RdBu
library(ggplot2)

ggplot(df, aes(x = disease, y = feature, color = correlation)) +
  geom_point(size = 8) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "top"
  ) +
  scale_size(range = c(1, 15)) +
  scale_color_gradientn(
    colors = colorRampPalette(c("blue", "lightblue", "pink", "red"))(100), 
    limits = c(-1, 1),
    na.value = "white" # Neutral for NA values
  ) +
  xlab("") +
  ylab("") +
  theme(
    axis.text.y.right = element_text(hjust = 0),
    legend.position = "left"
  ) +
  scale_y_discrete(position = "right")



p <- ggplot(df, aes(x = disease, y = feature, color = correlation)) +
  geom_point(size = 8) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "top"
  ) +
  scale_size(range = c(1, 15)) +
  scale_color_gradientn(
    colors = colorRampPalette(c("darkcyan", "lightgreen", "gold", "orange", "darkorange","red", "darkred"))(100),
    limits = c(-1, 1),
    na.value = "white" # Neutral for NA values
  ) +
  xlab("") +
  ylab("") +
  theme(
    axis.text.y.right = element_text(hjust = 0),
    legend.position = "left"
  ) +
  scale_y_discrete(position = "right")


p
jpeg("fig4B.jpeg", width = 17,height = 20, units = "cm", res = 600)
p
dev.off()


median(na.omit(df$correlation))
