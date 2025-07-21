setwd("/mnt/data3/appierdo/AutoimmuneLandscape/paper_figs/fig5")
source("../../utils_functions/functions_to_process_RNAseq_Data.R")
options(repr.plot.width=15, repr.plot.height=15)
library(ComplexHeatmap)
library(grid)
library(dendextend)
library(ggtext)

files <- list.files("/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/TISSUE/meta_analysis_pathways_activity/", recursive = T, pattern = '*complete*', full.names = T)
files <- files[grepl("wiki|kegg|hall", files)]

significant_genes<-c() #list of genes that are significant in at least one disease
counter = 0

for (file in files) {
  counter <- counter +1
  t<-read.csv(file, sep='\t')
  t<-t[t["numStudies"]>=2,]
  t <- t[t$heterogeneityPval > 0.05,]
  t<-t[t["effectSizeFDR"]<0.05,]
  t <- t[abs(t$effectSize) > 0.8, ]
  significant_genes<-append(significant_genes, t$gene)
}

tail(sort(table(significant_genes)), 100)

sign_genes_in_50_perc_diseases <- names(table(significant_genes)[table(significant_genes) >=2])


data<-data.frame(matrix(ncol = 3, nrow=0, 
                        dimnames = list(c(), c("effectSize", "gene", "disease"))))

for (file in files) {
  disease<-strsplit(file, "\\/")[[1]][13]
  tissue <- strsplit(file, "\\/|_")[[1]][17]
  disease <- paste0(disease, "_", tissue)
  df<-read.csv(file, sep='\t')
  df<-df[df$gene %in% sign_genes_in_50_perc_diseases, c("effectSize", "gene")]
  df$disease<-rep(disease, nrow(df))
  data<-rbind(data, df)
}

#reshaping
data<-tidyr::spread(data, "disease", "effectSize")
rownames(data)<-data$gene
data$gene<-NULL


n_studies<- c(27, 5, 7, 23, 3, 2, 6, 30)
data<-na.omit(data)
#data <- as.data.frame(t(data))

plot(density(as.matrix(data), na.rm = T))

column_dend = hclust(dist(t(data)))
row_dend = hclust(dist((data)))
#row_dend = color_branches(row_dend, k = 2)

# data[data > 2]<- 2
# data[data < -2]<- -2


column_ha = HeatmapAnnotation(nStudies = anno_barplot(n_studies, gp = gpar(fill='black', col='white', lwd=1), border = F, add_numbers = T, numbers_rot = 0), which = 'column')

colnames(data) <- c("Crohn's disease (Gut)", "Lupus (Kidney)", "Lupus (Skin)", "Psoriasis (Skin)",
                    "Rheumatoid arthritis (Synovium)", "Scleroderma (Skin)", "Sjogren's syndrome (Salivary glands)",
                    "Ulcerative colitis (Gut)")


rownames(data) <- unlist(lapply(
  rownames(data), function(x) {
    ifelse(grepl("ALLOGRAF|CHEMOKINE_SIGNALING_PATHWAY|GLYCOSAMINOGLYCAN_DEGRADATION", x), x, gsub("WP_|HALLMARK_|KEGG_", "", x))
  }
)
)

data


data <- as.data.frame(rowMeans(data))
data$diseases <- rownames(data)
colnames(data) <- c("mean", "diseases")
data$directionality <- ifelse(data$mean > 0, 'down', 'up')


p <- ggplot(data, aes(y = reorder(diseases, mean), x = mean)) +
  geom_bar(stat = 'identity', aes(fill = directionality)) +
  theme_classic() +
  xlab("mean effect size across diseases")+
  theme(axis.title.x = element_text(size = 13),
        axis.text.y = element_markdown()
  ) +
  ylab("") +
  scale_y_discrete(labels = c(
    "INTERFERON_TYPE_I_SIGNALING_PATHWAYS" = "<span style='color:red; font-size:13px;'><strong>INTERFERON_TYPE_I_SIGNALING_PATHWAYS</strong></strong></span>",
    "OVERVIEW_OF_INTERFERONSMEDIATED_SIGNALING_PATHWAY" = "<span style='color:red; font-size:13px;'><strong>OVERVIEW_OF_INTERFERONSMEDIATED_SIGNALING_PATHWAY</strong></strong></span>",
    "TYPE_II_INTERFERON_SIGNALING_IFNG" = "<span style='color:red; font-size:13px;'><strong>TYPE_II_INTERFERON_SIGNALING_IFNG</strong></strong></span>",
    "TYPE_III_INTERFERON_SIGNALING" = "<span style='color:red; font-size:13px;'><strong>TYPE_III_INTERFERON_SIGNALING</strong></strong></span>"
    
    
    
  )) +
  scale_fill_discrete(labels=c('upregulated', 'downregulated'), type=c("#F28E2B","#76B7B2"))


rownames(data)[grepl("ifn|interf", rownames(data), ignore.case = T)]

jpeg("fig5D_ht_path_all.jpeg", units = 'in', width = 15, height = 12, res = 600)
p
dev.off()

