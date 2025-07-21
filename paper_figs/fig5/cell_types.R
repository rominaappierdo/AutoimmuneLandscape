setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig5")
source("/mnt/data3/appierdo/PhD-project/utils_functions/functions_to_process_RNAseq_Data.R")
library(ggplot2)
library(ggtext)

blood_files <- list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/BLOOD/meta_analysis_cell_proportions/", pattern = 'complete', full.names = T, recursive = T)

b <- 0
for (file in blood_files) {
  b <- b+1
  disease <- strsplit(file, "/")[[1]][10]
  df <- read.csv(file, sep ='\t')
  df <- df[df$numStudies >=2, ]
  df <- df[df$effectSizeFDR < 0.05, ]
  df <- df[, c("gene", "effectSize")]
  colnames(df) <- c('gene', disease)
  if (b ==1 ) {
    merged <- df
  } else {
    merged <- merge(merged, df, by = 'gene', all =T)
  }
  
}
rownames(merged) <- merged$gene
merged$gene <- NULL

blood_df <- wide_to_long(merged, name_row_column = 'genes', name_column_column = 'diseases', name_value_column = 'effectSize')


#to get standard errors

b <- 0
for (file in blood_files) {
  b <- b+1
  disease <- strsplit(file, "/")[[1]][10]
  df <- read.csv(file, sep ='\t')
  df <- df[df$numStudies >=2, ]
  df <- df[df$effectSizeFDR < 0.05, ]
  df <- df[, c("gene", "effectSizeStandardError")]
  colnames(df) <- c('gene', disease)
  if (b ==1 ) {
    merged <- df
  } else {
    merged <- merge(merged, df, by = 'gene', all =T)
  }
  
}
rownames(merged) <- merged$gene
merged$gene <- NULL

blood_df_se <- wide_to_long(merged, name_row_column = 'genes', name_column_column = 'diseases', name_value_column = 'SE')

blood_df$SE <- blood_df_se$SE


p1 <- ggplot(blood_df, aes(y = diseases, x = genes)) +
  geom_point(aes(size = effectSize, color = effectSize)) +
  theme_classic()




tissue_files <- list.files("/mnt/data3/appierdo/AutoimmuneLandscape//Analysis/TISSUE/meta_analysis_cell_proportions//", pattern = 'complete', full.names = T, recursive = T)

t <- 0
for (file in tissue_files) {
  t <- t+1
  disease <- strsplit(file, "/")[[1]][12]
  tissue <- strsplit(file, "/|_")[[1]][16]
  disease <- paste0(disease ,"_", tissue)
  df <- read.csv(file, sep ='\t')
  df <- df[df$numStudies >=2, ]
  df <- df[df$effectSizeFDR < 0.05, ]
  df <- df[, c("gene", "effectSize")]
  colnames(df) <- c('gene', disease)
  if (t ==1 ) {
    merged <- df
  } else {
    merged <- merge(merged, df, by = 'gene', all =T)
  }
  
}
rownames(merged) <- merged$gene
merged$gene <- NULL

tissue_df <- wide_to_long(merged, name_row_column = 'genes', name_column_column = 'diseases', name_value_column = 'effectSize')


#to get standard errors
t <- 0
for (file in tissue_files) {
  t <- t+1
  disease <- strsplit(file, "/")[[1]][12]
  tissue <- strsplit(file, "/|_")[[1]][16]
  disease <- paste0(disease ,"_", tissue)
  print(disease)
  df <- read.csv(file, sep ='\t')
  df <- df[df$numStudies >=2, ]
  df <- df[df$effectSizeFDR < 0.05, ]
  df <- df[, c("gene", "effectSizeStandardError")]
  colnames(df) <- c('gene', disease)
  if (t ==1 ) {
    merged <- df
  } else {
    merged <- merge(merged, df, by = 'gene', all =T)
  }
  
}
rownames(merged) <- merged$gene
merged$gene <- NULL

tissue_df2 <- wide_to_long(merged, name_row_column = 'genes', name_column_column = 'diseases', name_value_column = 'SE')

tissue_df$SE <- tissue_df2$SE


p2 <- ggplot(tissue_df, aes(y = diseases, x = genes)) +
  geom_point(aes(size = effectSize, color = effectSize)) +
  theme_classic()



blood_df$source <- rep('Blood', nrow(blood_df))
tissue_df$source <- rep('Tissue', nrow(tissue_df))

m <- rbind(blood_df, tissue_df)
m$diseases[m$diseases == 'crohn'] <- "Crohn's disease"
m$diseases[m$diseases == 'diabetes'] <- "Type I diabetes"
m$diseases[m$diseases == 'ra'] <- "Rheumatoid arthritis"
m$diseases[m$diseases == 'psoriasis'] <- "Psoriasis"
m$diseases[m$diseases == 'sjogren'] <- "Sjogren's syndrome"
m$diseases[m$diseases == 'uc'] <- "Ulcerative colitis"
m$diseases[m$diseases == 'multiple_sclerosis'] <- "Multiple sclerosis"
m$diseases[m$diseases == 'scleroderma'] <- "Scleroderma"
m$diseases[m$diseases == 'jia'] <- "Juvenile idiopathic arthritis"
m$diseases[m$diseases == 'lupus'] <- "Lupus"

m$diseases[m$diseases == 'crohn_GUT'] <- "Crohn's disease (Gut)"
m$diseases[m$diseases == 'ra_SYNOVIUM'] <- "Rheumatoid arthritis (Synovium)"
m$diseases[m$diseases == 'psoriasis_SKIN'] <- "Psoriasis (Skin)"
m$diseases[m$diseases == 'sjogren_SALIVARY.GLAND'] <- "Sjogren's syndrome (Salivary glands)"
m$diseases[m$diseases == 'uc_GUT'] <- "Ulcerative colitis (Gut)"
m$diseases[m$diseases == 'scleroderma_SKIN'] <- "Scleroderma (Skin)"
m$diseases[m$diseases == 'lupus_SKIN'] <- "Lupus (Skin)"
m$diseases[m$diseases == 'lupus_KIDNEY'] <- "Lupus (Kidney)"




m <- m[m$genes == 'plasmacytoid_dendritic_cell',]


estimate <- 10     # The point estimate (e.g., mean)
se <- 2            # Standard error
conf_level <- 0.95


z_score <- qnorm((1 + conf_level) / 2)  # Z critical value for the given confidence level

# Calculate the confidence interval
m$lower_bound <- m$effectSize - z_score * m$SE
m$upper_bound <- m$effectSize + z_score * m$SE

m<- na.omit(m)

p <- ggplot(m, aes(x = effectSize, y = reorder(diseases, ifelse(grepl("Blood", source), 1, 0)),  color = source)) +
  geom_point(size = 4.5) +
  theme_classic() +
  xlim(c(-2, 2)) +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("") +
  theme(
    axis.text.y = element_markdown(),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.text = element_text(size=14),
    legend.title = element_text(size = 16)
  ) +
  ggtitle('Plasmacytoid dendritic cells') +
  scale_y_discrete(labels = c(
    "Ulcerative colitis" = "<span style='color:#ECB602; font-size:20px;'><strong>Ulcerative colitis</strong></strong></span>",
    "Type I diabetes" = "<span style='color:#ECB602; font-size:20px;'><strong>Type I diabetes</strong></span>",
    "Sjogren's syndrome" = "<span style='color:#ECB602; font-size:20px;'><strong>Sjogren's syndrome</strong></span>",
    "Scleroderma" = "<span style='color:#ECB602; font-size:20px;'><strong>Scleroderma</strong></span>",
    "Rheumatoid arthritis" = "<span style='color:#ECB602; font-size:20px;'><strong>Rheumatoid arthritis</strong></span>",
    "Psoriasis" = "<span style='color:#ECB602; font-size:20px;'><strong>Psoriasis</strong></span>",
    "Multiple sclerosis" = "<span style='color:#ECB602; font-size:20px;'><strong>Multiple sclerosis</strong></span>",
    "Lupus" = "<span style='color:#ECB602; font-size:20px;'><strong>Lupus</strong></span>",
    "Juvenile idiopathic arthritis" = "<span style='color:#ECB602; font-size:20px;'><strong>Juvenile idiopathic arthritis</strong></span>",
    "Crohn's disease" = "<span style='color:#ECB602; font-size:20px;'><strong>Crohn's disease</strong></span>",
    "Ulcerative colitis (Gut)" = "<span style='color:#55185D; font-size:20px;'><strong>Ulcerative colitis (Gut)</strong></span>",
    "Sjogren's syndrome (Salivary glands)" = "<span style='color:#55185D; font-size:20px;'><strong>Sjogren's syndrome (Salivary glands)</strong></span>",
    "Scleroderma (Skin)" = "<span style='color:#55185D; font-size:20px;'><strong>Scleroderma (Skin)</strong></span>",
    "Rheumatoid arthritis (Synovium)" = "<span style='color:#55185D; font-size:20px;'><strong>Rheumatoid arthritis (Synovium)</strong></span>",
    "Psoriasis (Skin)" = "<span style='color:#55185D; font-size:20px;'><strong>Psoriasis (Skin)</strong></span>",
    "Lupus (Skin)" = "<span style='color:#55185D; font-size:20px;'><strong>Lupus (Skin)</strong></span>",
    "Lupus (Kidney)" = "<span style='color:#55185D; font-size:20px;'><strong>Lupus (Kidney)</strong></span>",
    "Crohn's disease (Gut)" = "<span style='color:#55185D; font-size:20px;'><span><strong>Crohn's disease (Gut)</strong></span>"
    
  )) +
  geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound), width = 0.2) +
  xlim(-3.1, 3.1) +
  scale_color_manual(values = c('Blood' = '#ECB602', 'Tissue' = '#55185D'))




jpeg("fig5E_plasmacytoids.jpeg", units = 'in', width = 10, height = 8, res = 600)
p
dev.off()


